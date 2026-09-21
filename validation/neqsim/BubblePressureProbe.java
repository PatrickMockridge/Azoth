// The cricondenbar of the envelope case, computed a second way.
//
//     javac -proc:none -cp neqsim-f0c7436.jar BubblePressureProbe.java
//     java -cp .:neqsim-f0c7436.jar BubblePressureProbe
//
// `PhaseEnvelope` reads the cricondenbar off `calcPTphaseEnvelope`'s own bubble branch, which
// is a continuation trace: its maximum is whatever the step sequence happened to sample, and
// the maximum of a smooth curve is the one quantity a trace determines worst. This asks the
// bubble-point flash directly, at a temperature grid, and reports the largest pressure it
// finds - an independent number for the same quantity.
//
// NeqSim holds bar throughout, and every pressure below is labelled.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class BubblePressureProbe {

  public static void main(String[] args) {
    double[] fractions = {0.5, 0.5};
    double bestT = Double.NaN;
    double bestP = 0.0;

    System.out.printf("%8s  %16s%n", "T (K)", "P bubble (bar)");
    for (double t = 296.0; t <= 372.0; t += 2.0) {
      SystemInterface fluid = new SystemPrEos(t, 100.0);
      fluid.addComponent("methane", fractions[0]);
      fluid.addComponent("n-butane", fractions[1]);
      fluid.setMixingRule("classic");
      fluid.setAttractiveTerm(1);
      ThermodynamicOperations ops = new ThermodynamicOperations(fluid);
      try {
        ops.bubblePointPressureFlash(false);
        double p = fluid.getPressure();
        if (p > bestP) {
          bestP = p;
          bestT = t;
        }
        System.out.printf("%8.2f  %16.10f%n", t, p);
      } catch (Exception ex) {
        System.out.printf("%8.2f  %16s%n", t, "failed: " + ex.getMessage());
      }
    }
    System.out.printf("%ncricondenbar (direct bubble flash): P=%.10f bar at T=%.2f K%n", bestP, bestT);
  }
}
