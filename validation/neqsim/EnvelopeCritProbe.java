// Whether a coarse envelope step loses the critical point, and whether anything says so.
//
// `PTPhaseEnvelope` refines the critical point only when the two extreme K-values fall inside
// one window - `Kvallc < 1.05 && Kvalhc > 0.95` - latched once, with no bisection and no
// warning (equinor/neqsim#4078). Nothing bounds K while the step caps bound T and P, so a
// large enough step can walk past the window and `calcCrit()` is never called.
//
// This sweeps the two step caps and reports, for each, the critical point the operation
// publishes and the last point each branch actually reached. The claim is confirmed if a
// coarse cap loses the critical point while a fine cap finds it, and it is *silent* if
// nothing on the returned object distinguishes the two cases.
//
//     javac -proc:none -cp neqsim-f0c7436.jar EnvelopeCritProbe.java
//     java -Xmx2g -cp .:neqsim-f0c7436.jar EnvelopeCritProbe > captures/envelope_crit_probe.tsv

import java.util.Arrays;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.PTPhaseEnvelopeMichelsen;

public class EnvelopeCritProbe {

  static double lastFinite(double[] v) {
    for (int i = v.length - 1; i >= 0; i--) {
      if (Double.isFinite(v[i])) {
        return v[i];
      }
    }
    return Double.NaN;
  }

  static int finiteCount(double[] v) {
    return (int) Arrays.stream(v).filter(Double::isFinite).count();
  }

  /** One run at one pair of step caps. `dt`/`dp` of 0 leaves the class's own default. */
  static void row(String label, String[] names, double[] moles, double dt, double dp) {
    SystemInterface fluid = new SystemPrEos(273.15, 1.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);

    PTPhaseEnvelopeMichelsen op =
        new PTPhaseEnvelopeMichelsen(fluid, "out", 1.0e-10, 1.0, true);
    if (dt > 0.0) {
      op.setDTmax(dt);
    }
    if (dp > 0.0) {
      op.setDPmax(dp);
    }
    op.run();

    double[] bt = op.getBubblePointTemperatures();
    double[] bp = op.getBubblePointPressures();
    double[] dtArr = op.getDewPointTemperatures();
    double[] dpArr = op.getDewPointPressures();

    System.out.println(label);
    System.out.printf("critical_T_K=%.15g%n", op.getCriticalTemperature());
    System.out.printf("critical_P_bar=%.15g%n", op.getCriticalPressure());
    System.out.printf("number_of_critical_points=%d%n", op.getNumberOfCriticalPoints());
    System.out.printf("bubble_points=%d bubble_last_T_K=%.15g bubble_last_P_bar=%.15g%n",
        finiteCount(bt), lastFinite(bt), lastFinite(bp));
    System.out.printf("dew_points=%d dew_last_T_K=%.15g dew_last_P_bar=%.15g%n", finiteCount(dtArr),
        lastFinite(dtArr), lastFinite(dpArr));
    System.out.println();
  }

  /** One fluid, its caps swept from the class's default out to a deliberately absurd one. */
  static void sweep(String fluidName, String[] names, double[] moles) {
    System.out.println();
    System.out.println("################ " + fluidName + " ################");
    System.out.println("=== the class's own default caps (no setter called) ===");
    row("default", names, moles, 0.0, 0.0);

    System.out.println("=== step caps swept ===");
    for (double cap : new double[] { 10.0, 50.0, 100.0, 250.0, 500.0, 1000.0 }) {
      row("dTmax=" + (long) cap + " dPmax=" + (long) cap, names, moles, cap, cap);
    }
  }

  public static void main(String[] args) {
    sweep("methane/n-butane 50/50", new String[] { "methane", "n-butane" }, new double[] { 0.5, 0.5 });
    sweep("methane/n-butane/n-heptane 40/30/30",
        new String[] { "methane", "n-butane", "n-heptane" }, new double[] { 0.4, 0.3, 0.3 });
  }
}
