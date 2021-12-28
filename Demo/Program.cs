#region License Information
/*
 * This file is part of HEAL.CFR which is licensed under the MIT license.
 * See the LICENSE file in the project root for more information.
 */
#endregion

using HEAL.CFR;
using System;
using System.Diagnostics;

namespace Demo {
  class Program {
    static void Main(string[] args) {
      // Uses CFR to find an approximation to the Gamma function.
      // Example from https://arxiv.org/abs/2001.00624
      // Pablo Moscato, Haoyuan Sun, Mohammad Nazmul Haque.
      // Analytic Continued Fractions for Regression: A Memetic Algorithm Approach.
      // Expert Systems with Applications 179 (2021): 115018, 10.1016/j.eswa.2021.115018

      // generate data for the Gamma function
      // page 13;
      var minX = -2.683;
      var maxX = 4.5;
      var range = maxX - minX;
      int N = 873;
      // generate features: x^1, ... x^6, y
      // we are not using a constant 1 feature because the cfrac representation already contains a constant offset
      var xy = new double[N, 7];
      var y = new double[N];
      for (int step = 0; step < N; step++) {
        var x = minX + step * range / (N - 1);
        xy[step, 0] = x;
        for (int p = 1; p < 6; p++) {
          xy[step, p] = xy[step, p - 1] * x;
        }
        y[step] = MathNet.Numerics.SpecialFunctions.Gamma(x);
        xy[step, 6] = y[step];
      }

      var yVar = MathNet.Numerics.Statistics.ArrayStatistics.Variance(y); // for NMSE = MSE / var(y)

      var param = CFR.DefaultParameters;
      param.depth = 6;
      var randomSeed = 1234;

      var bestObjValue = double.PositiveInfinity;
      var bestMSE = double.PositiveInfinity;
      ContinuedFraction bestModel = null;
      var iter = 0;
      var sw = new Stopwatch();
      void writeCurrentObjectiveValue(Agent obj) {
        iter++;
        if (obj.pocketObjValue < bestObjValue) {
          bestObjValue = obj.pocketObjValue;
          bestModel = obj.pocket.Clone();
        }
        var curMSE = CFR.MeanSquaredError(obj.pocket, xy);
        if (curMSE < bestMSE) bestMSE = curMSE;
        Console.WriteLine($"{iter,5} {sw.Elapsed.TotalSeconds,5:f1} {obj.pocketObjValue,10:e3} {bestObjValue,10:e3} {curMSE,10:e3} {bestMSE,10:e3} {curMSE / yVar,10:e3} {bestMSE / yVar,10:e3}");
      }

      Console.WriteLine($"{"Iter",5} {"Sec",5} {"CurObj",10} {"BestObj",10} {"CurMSE",10} {"BestMSE",10} {"CurNMSE",10} {"BestNMSE",10}");
      sw.Start();
      CFR.Run(param, xy, new Random(randomSeed), null, writeCurrentObjectiveValue);
      Console.WriteLine(bestModel);
      Console.WriteLine($"MSE={CFR.MeanSquaredError(bestModel, xy)}");
    }
  }
}
