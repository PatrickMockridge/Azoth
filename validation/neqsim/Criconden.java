// NeqSim's cricondenbar and cricondentherm, scanned and refined.
//
//     javac -proc:none -cp neqsim-3.20.0.jar Criconden.java
//     java -cp .:neqsim-3.20.0.jar Criconden
//
// Both quantities come two ways, and the driver prints both because they are not the same
// number:
//
//   1. **Scanned** by `calcPTphaseEnvelope(1.0)` and read back with `get("cricondenbar")` /
//      `get("cricondentherm")` - the extreme of the traced locus, accurate to the trace's step.
//   2. **Refined** by `calcCricoP` / `calcCricoT`, which build a `CricondenBarFlash` /
//      `CricondenThermFlash` and Newton the `(ln K, ln T, ln P)` system with a third residual
//      standing in for `dP/dT = 0`. They need a starting `(T, P)` and two starting
//      compositions; the scanned extreme is what a caller has, so that is what this feeds them.
//
// `CricondenbarAnalyser` does the same thing publicly - clone, `calcPTphaseEnvelope(false, 1.)`,
// `get("cricondenbar")[1]` - so the scan is the intended entry and the refinement is the extra
// step.
//
// **What this measures: neither refinement produces anything.** Measured 2026-09-18 on
// methane/n-butane 50/50 and methane/n-butane/n-heptane 40/30/30 under `SystemPrEos`, and
// again on the SRK condensate `CricondenBarThermFlashTest` uses:
//
//   - `calcCricoP` throws `IsNaNException: PhasePrEos:molarVolumeAnalytical - Variable
//     compressibility factor is NaN` under PR, and under SRK returns its seed **unchanged to
//     the last digit** (`274.011887 / 112.349120` in, the same out) - the `run()` revert path
//     at `CricondenBarFlash.java:170`.
//   - `calcCricoT` throws under PR and returns its seed under SRK.
//   - `saturationops.CricondenbarFlash` via `calcCricondenBar()` throws the same
//     `molarVolumeAnalytical` NaN on both fluids.
//
// So NeqSim has three reachable cricondenbar implementations and none of them refines anything.
// The test that exercises them cannot see it: `CricondenBarThermFlashTest` tolerates
// `T=-1, P=-1` ("If not converged, the algorithm falls back to envelope estimate") and asserts
// the result is within **5%** of the envelope estimate, which the reverted value *is* by
// construction. A refinement that fails and a refinement that succeeds both pass.
//
// The third residual is a proxy in any case, and the class says so: "We use the direct
// specification: g_{n+1} = sum_i dgdlnT[i] * dgdlnP[i] as proxy for the cross-sensitivity. The
// proper form is from the bordered matrix approach" - a comment that does not match the line
// below it, in a file that also carries unresolved working notes at `:323-328` and `:358`.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class Criconden {

  static void one(String label, String[] names, double[] moles) {
    SystemInterface fluid = new SystemPrEos(273.15, 1.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);

    ThermodynamicOperations ops = new ThermodynamicOperations(fluid);
    ops.calcPTphaseEnvelope(1.0);

    double[] bar = (double[]) ops.get("cricondenbar");
    double[] therm = (double[]) ops.get("cricondentherm");
    // The trace exposes the phase compositions at both extremes, and `CricondenBarFlash` seeds
    // `ln K` from them. Passing the feed instead - which the first version of this driver did -
    // starts the Newton at `ln K = 0` for a two-phase point, and it does not recover.
    double[] barX0 = (double[]) ops.get("cricondenbarX");
    double[] barY0 = (double[]) ops.get("cricondenbarY");
    double[] thermX0 = (double[]) ops.get("cricondenthermX");
    double[] thermY0 = (double[]) ops.get("cricondenthermY");

    System.out.printf("%n=== %s%n", label);
    System.out.printf("scanned cricondenbar   T=%.10f K  P=%.10f bar%n", bar[0], bar[1]);
    System.out.printf("scanned cricondentherm T=%.10f K  P=%.10f bar%n", therm[0], therm[1]);

    // The refinement, seeded from the scanned extreme *and its compositions*. Seeding with the
    // feed instead starts the Newton at `ln K = 0` for a two-phase point; seeded properly, every
    // refinement here throws or returns its input, which is what this driver exists to show.
    int nc = fluid.getPhase(0).getNumberOfComponents();

    // The scan's own extreme as the seed: what a caller actually has.
    double[] barRefined = new double[3];
    barRefined[0] = bar[0];
    barRefined[1] = bar[1];
    double[] barX = new double[nc];
    double[] barY = new double[nc];
    for (int i = 0; i < nc; i++) {
      barX[i] = barX0[i];
      barY[i] = barY0[i];
    }
    System.out.printf("  seed x=%s%n  seed y=%s%n", java.util.Arrays.toString(barX0),
        java.util.Arrays.toString(barY0));
    try {
      ops.calcCricoP(barRefined, barX, barY);
      System.out.printf("refined cricondenbar   T=%.10f K  P=%.10f bar%n", barRefined[0],
          barRefined[1]);
      for (int i = 0; i < nc; i++) {
        System.out.printf("     x[%s]=%.10g  y[%s]=%.10g%n", names[i], barX[i], names[i], barY[i]);
      }
    } catch (Exception ex) {
      System.out.printf("refined cricondenbar   -> %s: %s%n", ex.getClass().getSimpleName(),
          ex.getMessage());
    }

    double[] thermRefined = new double[3];
    thermRefined[0] = therm[0];
    thermRefined[1] = therm[1];
    double[] thermX = new double[nc];
    double[] thermY = new double[nc];
    for (int i = 0; i < nc; i++) {
      thermX[i] = thermX0[i];
      thermY[i] = thermY0[i];
    }
    try {
      ops.calcCricoT(thermRefined, thermX, thermY);
      System.out.printf("refined cricondentherm T=%.10f K  P=%.10f bar%n", thermRefined[0],
          thermRefined[1]);
      for (int i = 0; i < nc; i++) {
        System.out.printf("     x[%s]=%.10g  y[%s]=%.10g%n", names[i], thermX[i], names[i], thermY[i]);
      }
    } catch (Exception ex) {
      System.out.printf("refined cricondentherm -> %s: %s%n", ex.getClass().getSimpleName(),
          ex.getMessage());
    }
  }

  /// The `saturationops` cricondenbar: a pressure ramp with an inner Newton on `T` that calls a
  /// real `TPflash` at every step. A different algorithm from the Michelsen system, and the only
  /// other reachable cricondenbar NeqSim has.
  static void saturationOps(String label, String[] names, double[] moles) {
    SystemInterface fluid = new SystemPrEos(273.15, 1.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);

    ThermodynamicOperations ops = new ThermodynamicOperations(fluid);
    try {
      ops.calcCricondenBar();
      System.out.printf("saturationops cricondenbar: P=%.10f bar  T=%.10f K%n",
          fluid.getPressure(), fluid.getTemperature());
    } catch (Exception ex) {
      System.out.printf("saturationops cricondenbar -> %s: %s%n", ex.getClass().getSimpleName(),
          ex.getMessage());
    }
  }

  public static void main(String[] args) {
    one("methane/n-butane 50/50", new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5});
    one("methane/n-butane/n-heptane 40/30/30",
        new String[] {"methane", "n-butane", "n-heptane"}, new double[] {0.4, 0.3, 0.3});
    System.out.printf("%n--- the saturationops algorithm, for comparison ---%n");
    saturationOps("methane/n-butane 50/50", new String[] {"methane", "n-butane"},
        new double[] {0.5, 0.5});
    saturationOps("methane/n-butane/n-heptane 40/30/30",
        new String[] {"methane", "n-butane", "n-heptane"}, new double[] {0.4, 0.3, 0.3});
  }
}
