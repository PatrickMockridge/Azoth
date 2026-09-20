// What `SystemSoreideWhitson` does with the six interaction rows that carry a comma.
//
//     javac -proc:none -cp neqsim-3.20.0.jar SoreideWhitsonProbe.java
//     java -cp .:neqsim-3.20.0.jar SoreideWhitsonProbe
//
// `INTER.csv`'s `KIJWhitsonSoriede` column is written with a **decimal comma** on six
// rows - `propane`/`CO2`, `n-butane`/`CO2`, `n-pentane`/`CO2`, `n-hexane`/`CO2`,
// `n-heptane`/`CO2` and `mercury`/`CO2` - against a point on the other 297. Verified
// against NeqSim's own database, `Double.parseDouble` throws on every one of them, and
// `EosMixingRuleHandler:240` is where the column is read:
//
//     intparam[k][l] = Double.parseDouble(dataSet.getString("KIJWhitsonSoriede"));
//
// **And the throw is swallowed.** The read sits inside the per-pair `try` whose
// `catch (Exception ex)` has both of its diagnostics commented out, and none of its
// fallbacks applies to a CO2/alkane pair - so `intparam[k][l]` keeps whatever the earlier
// branch wrote, which is `KIJSRK`. The Whitson-Soreide rule then evaluates those six pairs
// with the ordinary SRK interaction and nothing says so.
//
// Measured, by reading `intparam` off a handler built over a `PhaseSoreideWhitson`: every
// one of the six equals that pair's `KIJSRK` column - `propane`/`CO2` is `0.1018` where the
// Whitson column says `0.1241`, and `mercury`/`CO2` is `0.369` where it says `0.0145`, a
// factor of 25. `CO2`/`methane`, whose row carries a point, is `0.107` - the Whitson value
// - and is the control that says the column is otherwise read correctly.
//
// So this is not dead code and not a crash: it is a fitted parameter silently replaced by a
// different one. Written up to `~/Desktop/neqsim-issue-whitson-comma-rows-silently-fall-back.md`.
//
// The `addSalinity(0, "mole/sec")` call is not incidental: without it every pair throws
// `Index 0 out of bounds for length 0` from the saline parameterisation's empty arrays,
// which is a different matter and was the first reading of this probe's failures.

import neqsim.thermo.system.SystemSoreideWhitson;

public class SoreideWhitsonProbe {
  /** The kij the handler holds for a pair, by index, read off its own array. */
  private static double kij(neqsim.thermo.phase.PhaseInterface phase, int i, int j) {
    try {
      neqsim.thermo.mixingrule.EosMixingRuleHandler handler =
          new neqsim.thermo.mixingrule.EosMixingRuleHandler();
      handler.getMixingRule(11, phase);
      java.lang.reflect.Field field =
          neqsim.thermo.mixingrule.EosMixingRuleHandler.class.getDeclaredField("intparam");
      field.setAccessible(true);
      double[][] intparam = (double[][]) field.get(handler);
      return intparam[i][j];
    } catch (Throwable e) {
      return Double.NaN;
    }
  }

  private static void attempt(String label, String[] names, double[] z) {
    SystemSoreideWhitson s = new SystemSoreideWhitson(298.0, 20.0);
    for (int i = 0; i < names.length; i++) {
      s.addComponent(names[i], z[i], "mole/sec");
    }
    // `addSalinity` is the call that builds the saline parameterisation's arrays; without
    // it every pair throws `Index 0 out of bounds for length 0`, which has nothing to do
    // with the rows this probe is about.
    s.addSalinity(0, "mole/sec");
    s.setTotalFlowRate(15, "mole/sec");
    s.setMixingRule(11);
    System.out.printf("  %-38s ", label);
    try {
      s.init(0);
      s.init(1);
      System.out.println("ok   phase[1] = " + s.getPhase(1).getClass().getSimpleName());
      // **What the phase actually uses for each pair**, against what the table says.
      neqsim.thermo.phase.PhaseInterface p = s.getPhase(1);
      for (int i = 0; i < p.getNumberOfComponents(); i++) {
        for (int j = i + 1; j < p.getNumberOfComponents(); j++) {
          String a = p.getComponent(i).getComponentName();
          String b = p.getComponent(j).getComponentName();
          if (a.equalsIgnoreCase("water") || b.equalsIgnoreCase("water")) {
            continue;
          }
          System.out.printf("        kij(%-10s,%-10s) = %.10g%n", a, b, kij(p, i, j));
        }
      }
    } catch (Throwable e) {
      System.out.println(e.getClass().getSimpleName() + ": " + e.getMessage());
      for (StackTraceElement f : e.getStackTrace()) {
        if (f.getClassName().startsWith("neqsim.")) {
          System.out.println("        at " + f);
        }
      }
    }
  }

  public static void main(String[] args) {
    System.out.println("the shipped test's own mixture:");
    attempt("N2 CO2 C1 C2 water", new String[] {"nitrogen", "CO2", "methane", "ethane", "water"},
        new double[] {0.1, 0.2, 0.3, 0.3, 0.1});
    System.out.println();
    System.out.println("the six pairs whose row carries a comma:");
    for (String alkane : new String[] {"propane", "n-butane", "n-pentane", "n-hexane", "n-heptane", "mercury"}) {
      attempt("CO2 + " + alkane, new String[] {"CO2", alkane, "water"},
          new double[] {0.3, 0.3, 0.4});
    }
    System.out.println();
    System.out.println("and neighbours whose row does not:");
    for (String other : new String[] {"methane", "ethane", "n-octane", "nC10"}) {
      attempt("CO2 + " + other, new String[] {"CO2", other, "water"},
          new double[] {0.3, 0.3, 0.4});
    }
  }
}
