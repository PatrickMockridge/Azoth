// NeqSim's phase envelope for `eos.pt_phase_envelope`, branch by branch.
//
//     javac -proc:none -cp neqsim-3.20.0.jar PhaseEnvelope.java
//     java -cp .:neqsim-3.20.0.jar PhaseEnvelope
//
// `ThermodynamicOperations.calcPTphaseEnvelope()` builds a `PTPhaseEnvelopeMichelsen` with
// `phasefraction = 1e-10` and `lowPres = 1.0` bar, bubble branch first. That is the entry
// point the spec's validated case names, and the one this drives.
//
// The branches are printed in full because the case is about where they *end*: the item's
// acceptance is the dew branch reaching the critical point, and a summary of the endpoints
// cannot show whether it does. The pressures come back in bar, which is the unit NeqSim
// holds; `NeqSim uses bar throughout` is the hazard this port has been bitten by three
// times, so the unit is in every label below.

import java.util.Arrays;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.PTPhaseEnvelopeMichelsen;

public class PhaseEnvelope {

  static void one(String label, String[] names, double[] moles, double lowPresBar) {
    SystemInterface fluid = new SystemPrEos(273.15, 1.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);

    ThermodynamicOperations ops = new ThermodynamicOperations(fluid);
    ops.calcPTphaseEnvelope(lowPresBar);
    PTPhaseEnvelopeMichelsen op = (PTPhaseEnvelopeMichelsen) ops.getOperation();

    double[] bt = op.getBubblePointTemperatures();
    double[] bp = op.getBubblePointPressures();
    double[] dt = op.getDewPointTemperatures();
    double[] dp = op.getDewPointPressures();

    System.out.printf("%n=== %s   lowPres=%.4f bar%n", label, lowPresBar);
    System.out.printf("critical T=%.10f K   P=%.10f bar%n", op.getCriticalTemperature(),
        op.getCriticalPressure());
    // The arrays arrive with NaN padding at the end of the branch that met the critical
    // point, and the summaries below skip it: a NaN is a slot NeqSim never filled, not a
    // state, and a max over the raw array reports NaN and hides the branch.
    System.out.printf("bubble n=%d  NaN count=%d%n", bt.length, nanCount(bt));
    System.out.printf("dew    n=%d  NaN count=%d%n", dt.length, nanCount(dt));
    System.out.printf("bubble first T=%.10f K P=%.10f bar%n", first(bt), first(bp));
    System.out.printf("dew    first T=%.10f K P=%.10f bar%n", first(dt), first(dp));
    System.out.printf("bubble last  T=%.10f K P=%.10f bar%n", last(bt), last(bp));
    System.out.printf("dew    last  T=%.10f K P=%.10f bar%n", last(dt), last(dp));
    System.out.printf("bubble maxP=%.10f bar at T=%.10f K%n", max(bp), bt[argmax(bp)]);
    System.out.printf("dew    maxP=%.10f bar at T=%.10f K%n", max(dp), dt[argmax(dp)]);
    System.out.printf("bubble maxT=%.10f K at P=%.10f bar%n", max(bt), bp[argmax(bt)]);
    System.out.printf("dew    maxT=%.10f K at P=%.10f bar%n", max(dt), dp[argmax(dt)]);

    System.out.println("--- bubble branch (T K, P bar) ---");
    for (int i = 0; i < bt.length; i++) {
      System.out.printf("  %2d  %.10f  %.10f%n", i, bt[i], bp[i]);
    }
    System.out.println("--- dew branch (T K, P bar) ---");
    for (int i = 0; i < dt.length; i++) {
      System.out.printf("  %2d  %.10f  %.10f%n", i, dt[i], dp[i]);
    }
  }

  /** The largest entry that is a number; NaN is a slot NeqSim never filled. */
  static double max(double[] v) {
    return Arrays.stream(v).filter(Double::isFinite).max().orElse(Double.NaN);
  }

  static int argmax(double[] v) {
    int best = -1;
    for (int i = 0; i < v.length; i++) {
      if (!Double.isFinite(v[i])) {
        continue;
      }
      if (best < 0 || v[i] > v[best]) {
        best = i;
      }
    }
    return best;
  }

  static int nanCount(double[] v) {
    int count = 0;
    for (double value : v) {
      if (!Double.isFinite(value)) {
        count++;
      }
    }
    return count;
  }

  static double first(double[] v) {
    return argmax(v) < 0 ? Double.NaN : v[0];
  }

  /** The last entry before the NaN padding. */
  static double last(double[] v) {
    for (int i = v.length - 1; i >= 0; i--) {
      if (Double.isFinite(v[i])) {
        return v[i];
      }
    }
    return Double.NaN;
  }

  public static void main(String[] args) {
    // The spec's validated case: the pair `eos.pt_phase_envelope`'s own case uses, at the
    // low pressure the default entry point uses.
    one("methane/n-butane 50/50", new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5}, 1.0);
    // The same pair at a lower starting pressure, to see whether the low-pressure end of the
    // trace is the starting pressure or a physical limit.
    one("methane/n-butane 50/50", new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5}, 0.1);
    // A three-component mixture, so the case is not only about a binary.
    one("methane/n-butane/n-heptane 40/30/30",
        new String[] {"methane", "n-butane", "n-heptane"}, new double[] {0.4, 0.3, 0.3}, 1.0);
  }
}
