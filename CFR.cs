#region License Information
/*
 * This file is part of HEAL.CFR which is licensed under the MIT license.
 * See the LICENSE file in the project root for more information.
 */
#endregion

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Optimization;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;


namespace HEAL.CFR {
  /// <summary>
  /// Implementation of Continuous Fraction Regression (CFR) as described in 
  /// Pablo Moscato, Haoyuan Sun, Mohammad Nazmul Haque,
  /// Analytic Continued Fractions for Regression: A Memetic Algorithm Approach,
  /// Expert Systems with Applications, Volume 179, 2021, 115018, ISSN 0957-4174,
  /// https://doi.org/10.1016/j.eswa.2021.115018.
  /// </summary>
  public static class CFR {

    /// <summary>
    /// Parameter settings from the paper "Analytic Continued Fractions for Regression: A Memetic Algorithm Approach"
    /// </summary>
    public static Parameters DefaultParameters => new Parameters() {
      depth = 6,
      mutationRate = 0.1,
      numGen = 200,
      stagnatingGens = 5,
      delta = 0.1,
      localSearchParameters = new LocalSearchParameters() {
        iterations = 250,
        restarts = 4, // 5 runs total
        tolerance = 1e-3, // page 12 of the preprint
        minNumRows = 200,
        samplesFrac = 0.2,
      }
    };

    /// <summary>
    /// Parameters for the local search algorithm (Nelder-Mead)
    /// </summary>
    public class LocalSearchParameters {
      /// <summary>
      /// Number of iterations for local search (simplex) (default value 250)
      /// </summary>
      public int iterations;
      /// <summary>
      /// Number of restarts for local search (default value 4)
      /// </summary>
      public int restarts;

      /// <summary>
      /// The minimum number of samples for the local search (default 200)
      /// </summary>
      public int minNumRows;

      /// <summary>
      /// The fraction of samples used for local search. Only used when the number of samples is more than minNumRows (default 20%)
      /// </summary>
      public double samplesFrac;

      /// <summary>
      /// The tolerance value for local search (simplex) (default value: 1e-3) (page 12 of the preprint)
      /// </summary>
      public double tolerance;
    }

    /// <summary>
    /// Parameters for the algorithm. Get default parameter values with CFR.DefaultParameters
    /// </summary>
    public class Parameters {
      /// <summary>
      /// Number of input variables for the model
      /// </summary>
      internal int nVars;

      /// <summary>
      /// Depth of the continued fraction representation (default 6)
      /// </summary>
      public int depth;

      /// <summary>
      /// Mutation rate (default 10%)
      /// </summary>
      public double mutationRate;

      /// <summary>
      /// The maximum number of generations (default 200)
      /// </summary>
      public int numGen;

      /// <summary>
      /// Number of generations after which the population is re-initialized (default value 5)
      /// </summary>
      public int stagnatingGens;

      /// <summary>
      /// The relative weight for the number of variables term in the fitness function (default value: 10%)
      /// </summary>
      public double delta;

      public LocalSearchParameters localSearchParameters;
    }

    /// <summary>
    /// Executes the CFR algorithm.
    /// </summary>
    /// <param name="parameters">The algorithm parameters. Default parameter values can be generated with CFR.DefaultParameters.</param>
    /// <param name="xy">The matrix (n x (k + 1)) of variables. n rows, k input variables, last column is the target variable.</param>
    /// <param name="rand">A random number generator.</param>
    /// <param name="newBestSolutionCallback">A callback to be notified when a new best solution is found. </param>
    /// <param name="iterationCallback">A callback that is called on each iteration.</param>
    /// <returns>The objective value of the best model</returns>
    public static double Run(Parameters parameters, double[,] xy, Random rand,
                             Action<ContinuedFraction, double> newBestSolutionCallback = null,
                             Action<Agent> iterationCallback = null) {
      return Run(parameters, xy, rand, new CancellationTokenSource().Token, newBestSolutionCallback, iterationCallback);
    }

    /// <summary>
    /// Executes the CFR algorithm.
    /// </summary>
    /// <param name="parameters">The algorithm parameters. Default parameter values can be generated with CFR.DefaultParameters.</param>
    /// <param name="xy">The matrix (n x (k + 1)) of variables. n rows, k input variables, last column is the target variable. Scaling all variables to a common range (e.g. [0..1]) may be beneficial.</param>
    /// <param name="rand">A random number generator.</param>
    /// <param name="cancellationToken">Cancellation token to interrupt the algorithm.</param>
    /// <param name="newBestSolutionCallback">A callback to be notified when a new best solution is found. </param>
    /// <param name="iterationCallback">A callback that is called on each iteration.</param>
    /// <returns>The objective value of the best model</returns>
    public static double Run(Parameters parameters, double[,] xy, Random rand,
                           CancellationToken cancellationToken,
                           Action<ContinuedFraction, double> newBestSolutionCallback = null,
                           Action<Agent> iterationCallback = null) {

      parameters.nVars = xy.GetLength(1) - 1;

      /* Algorithm 1 */
      /* Generate initial population by a randomized algorithm */
      var pop = InitialPopulation(parameters, rand, xy);
      var bestObj = pop.pocketObjValue;
      // the best value since the last reset
      var episodeBestObj = pop.pocketObjValue;
      var episodeBestObjGen = 0;
      for (int gen = 1; gen <= parameters.numGen && !cancellationToken.IsCancellationRequested; gen++) {
        /* mutate each current solution in the population */
        var pop_mu = Mutate(parameters, pop, rand, xy);
        /* generate new population by recombination mechanism */
        var pop_r = RecombinePopulation(parameters, pop_mu, rand, xy);

        // Paper:
        // A period of individual search operation is performed every generation on all current solutions.

        // Statement by authors:
        // "We ran the Local Search after Mutation and recombination operations. We executed the local-search only on the Current solutions."
        // "We executed the MaintainInvariant() in the following steps:
        // - After generating the initial population
        // - after resetting the root
        // - after executing the local-search on the whole population.
        // We updated the pocket/ current automatically after mutation and recombination operation."

        /* local search optimization of current solutions */
        foreach (var agent in pop_r.IterateLevels()) {
          LocalSearchSimplex(parameters.localSearchParameters, parameters.delta, agent, xy, rand);
          Debug.Assert(agent.pocketObjValue < agent.currentObjValue);
        }
        foreach (var agent in pop_r.IteratePostOrder()) agent.MaintainInvariant(); // post-order to make sure that the root contains the best model 
        // foreach (var agent in pop_r.IteratePostOrder()) agent.AssertInvariant();

        // for detecting stagnation we track the best objective value since the last reset 
        // and reset if this does not change for stagnatingGens
        if (gen > episodeBestObjGen + parameters.stagnatingGens) {
          Reset(parameters, pop_r, rand, xy);
          episodeBestObj = double.MaxValue;
        }
        if (episodeBestObj > pop_r.pocketObjValue) {
          episodeBestObjGen = gen; // wait at least stagnatingGens until resetting again
          episodeBestObj = pop_r.pocketObjValue;
        }

        /* replace old population with evolved population */
        pop = pop_r;

        /* keep track of the best solution */
        if (bestObj > pop.pocketObjValue) {
          bestObj = pop.pocketObjValue;
          newBestSolutionCallback?.Invoke(pop.pocket, bestObj);
        }

        iterationCallback?.Invoke(pop);
      }
      return bestObj;
    }

    private static Agent InitialPopulation(Parameters parameters, Random rand, double[,] xy) {
      /* instantiate 13 agents in the population */
      var pop = new Agent();
      // see Figure 2
      for (int i = 0; i < 3; i++) {
        pop.children.Add(new Agent());
        for (int j = 0; j < 3; j++) {
          pop.children[i].children.Add(new Agent());
        }
      }

      // Statement by the authors: "Yes, we use post-order traversal here"
      foreach (var agent in pop.IteratePostOrder()) {
        agent.current = new ContinuedFraction(parameters.nVars, parameters.depth, rand);
        agent.pocket = new ContinuedFraction(parameters.nVars, parameters.depth, rand);

        agent.currentObjValue = Evaluate(agent.current, xy, parameters.delta);
        agent.pocketObjValue = Evaluate(agent.pocket, xy, parameters.delta);

        /* within each agent, the pocket solution always holds the better value of guiding
         * function than its current solution
         */
        agent.MaintainInvariant();
      }

      // foreach (var agent in pop.IteratePostOrder()) agent.AssertInvariant();

      return pop;
    }

    // Our criterion for relevance has been fairly strict: if no
    // better model has been produced for five(5) straight generations,
    // then the pocket of the root agent is removed and a new solution is created at random.

    // Statement by the authors: "We only replaced the pocket solution of the root with
    // a randomly generated solution. Then we execute the maintain-invariant process.
    // It does not initialize the solutions in the entire population."
    private static void Reset(Parameters parameters, Agent root, Random rand, double[,] xy) {
      root.pocket = new ContinuedFraction(parameters.nVars, parameters.depth, rand);
      root.pocketObjValue = Evaluate(root.pocket, xy, parameters.delta);

      foreach (var agent in root.IteratePreOrder()) agent.MaintainInvariant(); // Here we use pre-order traversal push the newly created model down the hierarchy. 

      foreach (var agent in root.IteratePostOrder()) agent.AssertInvariant();
    }



    private static Agent RecombinePopulation(Parameters parameters, Agent pop, Random rand, double[,] xy) {
      var l = pop;

      if (pop.children.Count > 0) {
        var s1 = pop.children[0];
        var s2 = pop.children[1];
        var s3 = pop.children[2];

        // Statement by the authors: "we are using recently generated solutions.
        // For an example, in step 1 we got the new current(l), which is being used in
        // Step 2 to generate current(s3). The current(s3) from Step 2 is being used at
        // Step 4. These steps are executed sequentially from 1 to 4. Similarly, in the
        // recombination of lower-level subpopulations, we will have the new current
        // (the supporters generated at the previous level) as the leader of the subpopulation."
        Recombine(parameters, l, s1, SelectRandomOp(rand), rand, xy);
        Recombine(parameters, s3, l, SelectRandomOp(rand), rand, xy);
        Recombine(parameters, s1, s2, SelectRandomOp(rand), rand, xy);
        Recombine(parameters, s2, s3, SelectRandomOp(rand), rand, xy);

        // recombination works from top to bottom
        foreach (var child in pop.children) {
          RecombinePopulation(parameters, child, rand, xy);
        }

      }
      return pop;
    }

    private static ContinuedFraction Recombine(Parameters parameters, Agent a, Agent b, Func<double[], double[], double[]> op, Random rand, double[,] xy) {
      ContinuedFraction p1 = a.pocket;
      ContinuedFraction p2 = b.pocket;
      ContinuedFraction ch = new ContinuedFraction() { h = new Term[p1.h.Length] };
      /* apply a recombination operator chosen uniformly at random on variables of two parents into offspring */
      ch.vars = op(p1.vars, p2.vars);

      /* recombine the coefficients for each term (h) of the continued fraction */
      for (int i = 0; i < p1.h.Length; i++) {
        var coefa = p1.h[i].coef; var varsa = p1.h[i].vars;
        var coefb = p2.h[i].coef; var varsb = p2.h[i].vars;

        /* recombine coefficient values for variables */
        var coefx = new double[parameters.nVars];
        var varsx = new double[parameters.nVars]; // deviates from paper, probably forgotten in the pseudo-code
        for (int vi = 0; vi < parameters.nVars; vi++) {
          if (ch.vars[vi] > 0) {  // CHECK: paper uses featAt()
            if (varsa[vi] * varsb[vi] > 0) {
              coefx[vi] = coefa[vi] + (rand.NextDouble() * 5 - 1) * (coefb[vi] - coefa[vi]) / 3.0;
              varsx[vi] = 1.0;
            } else if (varsa[vi] > 0) {
              coefx[vi] = coefa[vi];
              varsx[vi] = 1.0;
            } else if (varsb[vi] > 0) {
              coefx[vi] = coefb[vi];
              varsx[vi] = 1.0;
            }
          }
        }
        /* update new coefficients of the term in offspring */
        ch.h[i] = new Term() { coef = coefx, vars = varsx };
        /* compute new value of constant (beta) for term hi in the offspring solution ch using 
         * beta of p1.hi and p2.hi */
        ch.h[i].beta = p1.h[i].beta + (rand.NextDouble() * 5 - 1) * (p2.h[i].beta - p1.h[i].beta) / 3.0;
      }

      a.current = ch;
      LocalSearchSimplex(parameters.localSearchParameters, parameters.delta, a, xy, rand);
      return ch;
    }

    private static Agent Mutate(Parameters parameters, Agent pop, Random rand, double[,] xy) {
      foreach (var agent in pop.IterateLevels()) {
        if (rand.NextDouble() < parameters.mutationRate) {
          if (agent.currentObjValue < 1.2 * agent.pocketObjValue ||
              agent.currentObjValue > 2 * agent.pocketObjValue)
            ToggleVariables(agent.current, rand); // major mutation
          else
            ModifyVariable(agent.current, rand); // soft mutation

          // Finally, the local search operation is executed on the mutated solution in order to optimize
          // non-zero coefficients. We do not apply mutation on pocket solutions because we consider them as a "collective memory"
          // of good models visited in the past. They influence the search process via recombination only.
          LocalSearchSimplex(parameters.localSearchParameters, parameters.delta, agent, xy, rand);
        }
      }
      return pop;
    }

    private static void ToggleVariables(ContinuedFraction cfrac, Random rand) {
      double coinToss(double a, double b) {
        return rand.NextDouble() < 0.5 ? a : b;
      }

      /* select a variable index uniformly at random */
      int N = cfrac.vars.Length;
      var vIdx = rand.Next(N);

      /* for each depth of continued fraction, toggle the selection of variables of the term (h) */
      foreach (var h in cfrac.h) {
        /* Case 1: cfrac variable is turned ON: Turn OFF the variable, and either 'Remove' or 
         * 'Remember' the coefficient value at random */
        if (cfrac.vars[vIdx] > 0.0) {
          h.vars[vIdx] = 0.0;
          h.coef[vIdx] = coinToss(0, h.coef[vIdx]);
        } else {
          /* Case 2: term variable is turned OFF: Turn ON the variable, and either 'Remove' 
           * or 'Replace' the coefficient with a random value between -3 and 3 at random */
          if (h.vars[vIdx] == 0.0) {
            h.vars[vIdx] = 1.0;
            h.coef[vIdx] = coinToss(0, rand.NextDouble() * 6 - 3);
          }
        }
      }
      /* toggle the randomly selected variable */
      cfrac.vars[vIdx] = 1.0 - cfrac.vars[vIdx]; // NOT vars[vIdx]
    }

    private static void ModifyVariable(ContinuedFraction cfrac, Random rand) {
      /* randomly select a variable which is turned ON */
      var candVars = new List<int>();
      for (int i = 0; i < cfrac.vars.Length; i++) if (cfrac.vars[i] > 0) candVars.Add(i);
      if (candVars.Count == 0) return; // no variable active
      var vIdx = candVars[rand.Next(candVars.Count)];

      /* randomly select a term (h) of continued fraction */
      var h = cfrac.h[rand.Next(cfrac.h.Length)];

      /* modify the coefficient value */
      if (h.vars[vIdx] > 0) {
        h.coef[vIdx] = 0.0;
      } else {
        h.coef[vIdx] = rand.NextDouble() * 6 - 3;
      }
      /* Toggle the randomly selected variable */
      h.vars[vIdx] = 1.0 - h.vars[vIdx]; // NOT vars[vIdx]
    }

    public static double MeanSquaredError(ContinuedFraction cfrac, double[,] xy) {
      var output = Evaluate(cfrac, xy);
      var yIdx = xy.GetLength(1) - 1;
      var sum = 0.0;
      for (int r = 0; r < output.Length; r++) {
        var resid = xy[r, yIdx] - output[r];
        sum += resid * resid;
      }
      return sum / xy.GetLength(0);
    }

    private static double Evaluate(ContinuedFraction cfrac, double[,] xy, double delta) {
      return MeanSquaredError(cfrac, xy) * (1 + delta * cfrac.vars.Sum()); // count vars_i == 1.0
    }

    // evalutes the cfrac for a single datapoint
    private static double Evaluate(ContinuedFraction cfrac, double[] dataPoint) {
      var res = 0.0;
      for (int i = cfrac.h.Length - 1; i > 1; i -= 2) {
        var hi = cfrac.h[i];
        var hi1 = cfrac.h[i - 1];
        var denom = hi.beta + Dot(hi.vars, hi.coef, dataPoint) + res;
        var numerator = hi1.beta + Dot(hi1.vars, hi1.coef, dataPoint);
        res = numerator / denom;
      }
      var h0 = cfrac.h[0];
      res += h0.beta + Dot(h0.vars, h0.coef, dataPoint);
      return res;
    }

    // vectorized version of Evalute for a set of datapoints
    private static double[] Evaluate(ContinuedFraction cfrac, double[,] xy) {
      var output = new double[xy.GetLength(0)]; // may 
      var numerator = new double[xy.GetLength(0)];
      for (int i = cfrac.h.Length - 1; i > 1; i -= 2) {
        var hi = cfrac.h[i];
        var hi1 = cfrac.h[i - 1];
        Dot(hi.vars, hi.coef, xy, output); // denom = output = vars * coef * xy + output
        Array.Clear(numerator, 0, numerator.Length);
        Dot(hi1.vars, hi1.coef, xy, numerator); // numerator = vars * coeff * xy
        for (int j = 0; j < output.Length; j++)
          output[j] = (numerator[j] + hi1.beta) / (output[j] + hi.beta);
      }
      var h0 = cfrac.h[0];
      Dot(h0.vars, h0.coef, xy, output);
      for (int j = 0; j < output.Length; j++) {
        output[j] += h0.beta;
      }
      return output;
    }

    // we use double instead of bools here because we use the 'vars' as a weight vector 
    private static Func<double[], double[], double[]> SelectRandomOp(Random rand) {
      double[] union(double[] a, double[] b) {
        var res = new double[a.Length];
        for (int i = 0; i < a.Length; i++) res[i] = Math.Sign(a[i] + b[i]); // OR
        return res;
      }
      double[] intersect(double[] a, double[] b) {
        var res = new double[a.Length];
        for (int i = 0; i < a.Length; i++) res[i] = a[i] * b[i]; // AND
        return res;
      }
      double[] symmetricDifference(double[] a, double[] b) {
        var res = new double[a.Length];
        for (int i = 0; i < a.Length; i++) res[i] = a[i] + b[i] == 1.0 ? 1.0 : 0.0; // XOR
        return res;
      }
      return rand.Next(3) switch {
        0 => union,
        1 => intersect,
        2 => symmetricDifference,
        _ => throw new ArgumentException(),
      };
    }



    // The authors used the Nelder Mead solver from https://direct.mit.edu/evco/article/25/3/351/1046/Evolving-a-Nelder-Mead-Algorithm-for-Optimization
    // Using different solvers (e.g. LevMar) is mentioned but not analysed.
    // We use the implementation from MathNET.Numerics
    private static void LocalSearchSimplex(LocalSearchParameters parameters, double delta, Agent a, double[,] xy, Random rand) {
      int maxEvals = parameters.iterations;
      int numSearches = parameters.restarts + 1;
      var numRows = xy.GetLength(0);
      int numSelectedRows = numRows;
      if (numRows > parameters.minNumRows)
        numSelectedRows = (int)(numRows * parameters.samplesFrac);

      var ch = a.current;
      var quality = Evaluate(ch, xy, delta); // get quality with original coefficients

      double[] origCoeff = ExtractCoeff(ch);
      if (origCoeff.Length == 0) return; // no parameters to optimize

      // the best result over all repetitions
      var bestQuality = quality;
      var bestCoeff = CreateVector.Dense(origCoeff);

      // the best result for the repetition
      double currentBestQuality = double.PositiveInfinity;
      Vector<double> currentBestCoeff = null;

      // buffers for selecting random subsets of the training dataset
      double[,] fittingData = null;
      int[] idx = null;

      var objFunc = ObjectiveFunction.Value(curCoeff => {
        SetCoeff(ch, curCoeff);
        var obj = Evaluate(ch, fittingData, delta);
        // we have to track the best result here because MathNET.Numerics does not handle reaching the maximum iterations gracefully.
        if (obj < currentBestQuality) {
          currentBestQuality = obj;
          currentBestCoeff = curCoeff.Clone();
        }
        return obj;
      });

      for (int count = 0; count < numSearches; count++) {
        SelectRandomRows(xy, ref fittingData, ref idx, numSelectedRows, rand);


        var initialGuess = CreateVector.Dense<double>(origCoeff.Length);
        var initialPerturbation = CreateVector.Dense<double>(origCoeff.Length);

        currentBestCoeff = initialGuess;
        currentBestQuality = double.PositiveInfinity;

        for (int i = 0; i < origCoeff.Length; i++) {
          initialGuess[i] = origCoeff[i];
          initialPerturbation[i] = 1.0;
        }
        try {
          // the paper uses early stopping, but this is not possible with the MathNET Numerics implementation
          NelderMeadSimplex.Minimum(objFunc, initialGuess, initialPerturbation, parameters.tolerance, parameters.iterations);
        } catch (MaximumIterationsException e) {
          // MathNET.Numerics does not simply stop when max iterations are reached but throws an exception instead.
          // Therefore we have to track the best result in the objective function instead (see above).
        }
        SetCoeff(ch, currentBestCoeff);

        // the result with the best guiding function value (on the entire dataset) is chosen.
        var newQuality = Evaluate(ch, xy, delta);

        if (newQuality < bestQuality) {
          bestCoeff = currentBestCoeff;
          bestQuality = newQuality;
        }
      }

      SetCoeff(a.current, bestCoeff);
      a.currentObjValue = bestQuality;

      // Unclear what the following means exactly. 
      // 
      // "We remind again that
      // each solution corresponds to a single model, this means that if a current model becomes better than its corresponding
      // pocket model (in terms of the guiding function of the solution), then an individual search optimization step is also
      // performed on the pocket solution/ model before we swap it with the current solution/ model. Individual search can then
      // make a current model better than the pocket model (again, according to the guiding function), and in that case they
      // switch positions within the agent that contains both of them."

      a.MaintainPocketCurrentInvariant();
    }

    // We are reusing selectedDataset and idx to save allocations and GC overhead.
    // selectedDataset is overwritten and only allocated if necessary.
    // idx is shuffled and only allocated if necessary.
    private static void SelectRandomRows(double[,] fullDataset, ref double[,] selectedDataset, ref int[] idx, int numSelectedRows, Random rand) {
      var numRows = fullDataset.GetLength(0);
      var numCols = fullDataset.GetLength(1);
      if (selectedDataset == null || selectedDataset.GetLength(0) != numSelectedRows)
        selectedDataset = new double[numSelectedRows, numCols];

      if (idx == null || idx.Length != numRows) {
        idx = Enumerable.Range(0, numRows).ToArray();
      }
      Shuffle(idx, rand);
      for (int i = 0; i < numSelectedRows; i++) {
        for (int c = 0; c < numCols; c++) {
          selectedDataset[i, c] = fullDataset[idx[i], c];
        }
      }
    }

    private static double[] ExtractCoeff(ContinuedFraction ch) {
      var coeff = new List<double>();
      foreach (var hi in ch.h) {
        coeff.Add(hi.beta);
        for (int vIdx = 0; vIdx < hi.vars.Length; vIdx++) {
          if (hi.vars[vIdx] > 0 && hi.coef[vIdx] != 0) coeff.Add(hi.coef[vIdx]); // paper states twice (for mutation and recombination) that non-zero coefficients are optimized
        }
      }
      return coeff.ToArray();
    }

    private static void SetCoeff(ContinuedFraction ch, Vector<double> curCoeff) {
      int k = 0;
      foreach (var hi in ch.h) {
        hi.beta = curCoeff[k++];
        for (int vIdx = 0; vIdx < hi.vars.Length; vIdx++) {
          if (hi.vars[vIdx] > 0 && hi.coef[vIdx] != 0) hi.coef[vIdx] = curCoeff[k++]; // paper states twice (for mutation and recombination) that non-zero coefficients are optimized
        }
      }
    }

    #region util
    // sum of w_i * x_i * y_i
    private static double Dot(double[] weight, double[] x, double[] y) {
      var s = 0.0;
      for (int i = 0; i < x.Length; i++)
        s += weight[i] * x[i] * y[i];
      return s;
    }

    // vector for each row m_j in m: v_j += w_i * x_i * w_j,i
    private static void Dot(double[] weight, double[] x, double[,] m, double[] v) {
      for (int j = 0; j < v.Length; j++) {
        var s = 0.0;
        for (int i = 0; i < x.Length; i++)
          s += weight[i] * x[i] * m[j, i];
        v[j] += s;
      }
    }

    private static void Shuffle(int[] a, Random rand) {
      // Fisher - Yates shuffle
      for (int i = 0; i < a.Length; i++) {
        var j = rand.Next(i + 1);
        if (i != j) {
          var tmp = a[i];
          a[i] = a[j];
          a[j] = tmp;
        }
      }
    }
    #endregion
  }
}
