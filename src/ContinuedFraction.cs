#region License Information
/*
 * This file is part of HEAL.CFR which is licensed under the MIT license.
 * See the LICENSE file in the project root for more information.
 */
#endregion

using System;
using System.Text;

namespace HEAL.CFR {
  public class ContinuedFraction {

    internal double[] vars; // use a double weight instead of boolean to speed up evaluation
    internal Term[] h;

    public ContinuedFraction() { }
    private ContinuedFraction(ContinuedFraction orig) {
      vars = (double[])orig.vars.Clone();
      h = new Term[orig.h.Length];
      for (int i = 0; i < h.Length; i++) {
        h[i] = new Term() {
          beta = orig.h[i].beta,
          coef = (double[])orig.h[i].coef.Clone(),
          vars = (double[])orig.h[i].vars.Clone()
        };
      }
    }

    public ContinuedFraction(int nVars, int depth, Random rand) {
      vars = new double[nVars];
      for (int i = 0; i < nVars; i++) vars[i] = rand.NextDouble() < 0.3 ? 1.0 : 0.0; // page 12 of the preprint. Each input variable then has a probability p = 1/3 to be present in the whitelist

      h = new Term[depth * 2 + 1];
      for (int i = 0; i < h.Length; i++) {
        h[i] = new Term();
        var hi = h[i];
        hi.vars = (double[])vars.Clone();
        hi.coef = new double[nVars];
        for (int vi = 0; vi < nVars; vi++) {
          if (hi.vars[vi] > 0.0)
            hi.coef[vi] = rand.NextDouble() * 6 - 3;
        }
        hi.beta = rand.NextDouble() * 6 - 3;
      }
    }

    public ContinuedFraction Clone() {
      return new ContinuedFraction(this);
    }

    public override string ToString() {
      var sb = new StringBuilder();
      FormatLinearExpr(sb, h[0]);
      ToInfix(sb, 1);
      return sb.ToString();
    }
    public void ToInfix(StringBuilder sb, int startIdx) {
      if (startIdx > h.Length - 2) return;
      // returns the infix expression
      sb.Append(" + (");
      FormatLinearExpr(sb, h[startIdx]);
      sb.Append(") / (");
      FormatLinearExpr(sb, h[startIdx + 1]);
      ToInfix(sb, startIdx + 2);
      sb.Append(")");
      //         expr1
      // expr0 + -----
      //         expr2

      //           expr1
      // expr0 +  -------------
      //                  expr4
      //          expr3 + -----
      //                  expr5
    }

    private static void FormatLinearExpr(StringBuilder sb, Term term) {
      sb.Append($"{term.beta:g5}");
      for (int i = 0; i < term.coef.Length; i++) {
        if (term.vars[i] * term.coef[i] != 0) { // active and non-zero coeff
          sb.Append($" + {term.coef[i]:g5} * x{i}");
        }
      }
    }

  }
  public class Term {
    public double beta;
    public double[] coef;
    public double[] vars; // use a double weight instead of boolean to speed up evaluation
  }
}