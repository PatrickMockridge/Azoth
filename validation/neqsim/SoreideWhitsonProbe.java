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
// **And this probe does not reach it.** Every mixture below initializes, including the
// six comma pairs and the shipped test's own composition, with the phase confirmed to be
// a `PhaseSoreideWhitson` in each. So the read site is somewhere rule 11 does not go, and
// the column appears to be unread in 3.20.0: the non-aqueous pairs come from
// `WhitsonSoreideMixingRule.getkijWhitsonSoreideNonAqueous`, which falls through to
// `ClassicSRK.getkij` and the `kijsrk`/`kijpr` columns.
//
// That is the state the tranche records rather than a conclusion: the value is
// unparseable and the read site is not reached by the one rule that would use it. Finding
// which dispatch does reach it is what decides whether this is a defect or dead code, and
// until that is answered there is nothing to report upstream.
//
// The `addSalinity(0, "mole/sec")` call is not incidental: without it every pair throws
// `Index 0 out of bounds for length 0` from the saline parameterisation's empty arrays,
// which is a different matter and was the first reading of this probe's failures.

import neqsim.thermo.system.SystemSoreideWhitson;

public class SoreideWhitsonProbe {
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
