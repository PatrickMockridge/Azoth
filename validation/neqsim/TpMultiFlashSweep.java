// A grid sweep: how often `setMultiPhaseCheck(true)` changes a cubic flash at all.
//
//     javac -proc:none -cp neqsim-3.20.0.jar TpMultiFlashSweep.java
//     java -cp .:neqsim-3.20.0.jar TpMultiFlashSweep
//
// `TpMultiFlashProbe.java` found no state where the multiflash reported three phases and one
// state where it differed from `TPflash` at all - a water/hydrocarbon pair, through the
// aqueous seed. A handful of hand-picked states cannot support that, so this grids the
// temperature-pressure plane for each mixture and counts.
//
// Reports, per mixture: the states run, the states where the phase *count* differs, the
// largest phase-count seen with the flag on, and the largest mole-fraction difference
// between the two runs at any state. A state that neither differs in count nor in
// composition is a state where the flag did nothing.

import java.util.Arrays;
import java.util.Comparator;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class TpMultiFlashSweep {

  static SystemInterface build(double temperatureK, double pressureBar, String[] names, double[] moles) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);
    return fluid;
  }

  /** One state under one setting, or null if the flash threw. */
  static double[][] run(double temperatureK, double pressureBar, String[] names, double[] moles,
      boolean multi, int[] countOut) {
    SystemInterface fluid = build(temperatureK, pressureBar, names, moles);
    fluid.setMultiPhaseCheck(multi);
    try {
      new ThermodynamicOperations(fluid).TPflash();
    } catch (Exception ex) {
      return null;
    }
    int n = fluid.getNumberOfPhases();
    countOut[0] = n;
    // `getNumberOfPhases()` and the phase array disagree after a phase is removed, so a read
    // past the array throws. That is a finding about the state, not a reason to lose the
    // whole sweep: the count is what this driver is counting, so report what can be read.
    double[][] out;
    try {
      out = new double[n][];
      for (int i = 0; i < n; i++) {
        out[i] = new double[fluid.getPhase(i).getNumberOfComponents() + 1];
        out[i][0] = fluid.getPhase(i).getBeta();
        for (int j = 0; j < fluid.getPhase(i).getNumberOfComponents(); j++) {
          out[i][j + 1] = fluid.getPhase(i).getComponent(j).getx();
        }
      }
    } catch (RuntimeException ex) {
      countOut[0] = -n; // readable count, negated: the phases could not all be read
      return null;
    }
    return out;
  }

  static void sweep(String label, String[] names, double[] moles) {
    int states = 0;
    int countDiffers = 0;
    int maxPhasesMulti = 0;
    int maxPhasesPlain = 0;
    double maxXDiff = 0.0;
    double worstT = 0.0;
    double worstP = 0.0;

    double[] pressures = {1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0};
    for (int ti = 0; ti <= 22; ti++) {
      double t = 180.0 + 10.0 * ti;
      for (double p : pressures) {
        int[] nPlain = new int[1];
        int[] nMulti = new int[1];
        double[][] plain = run(t, p, names, moles, false, nPlain);
        double[][] multi = run(t, p, names, moles, true, nMulti);
        if (nPlain[0] == 0 && nMulti[0] == 0) {
          continue; // both threw before reporting a count
        }
        states++;
        maxPhasesPlain = Math.max(maxPhasesPlain, Math.abs(nPlain[0]));
        maxPhasesMulti = Math.max(maxPhasesMulti, Math.abs(nMulti[0]));
        if (nPlain[0] != nMulti[0] || plain == null || multi == null) {
          countDiffers++;
          System.out.printf("      differs: T=%.0f P=%.0f plain=%d multi=%d%s%n", t, p, nPlain[0],
              nMulti[0], (plain == null || multi == null) ? "  (unreadable)" : "");
          continue;
        }
        // Same count: pair the phases by descending beta - the two runs list them in the
        // order the flash happened to build them - and take the largest component gap.
        Comparator<double[]> byBeta = (a, b) -> Double.compare(b[0], a[0]);
        Arrays.sort(plain, byBeta);
        Arrays.sort(multi, byBeta);
        double worst = 0.0;
        for (int i = 0; i < nPlain[0]; i++) {
          for (int k = 1; k < plain[i].length; k++) {
            worst = Math.max(worst, Math.abs(plain[i][k] - multi[i][k]));
          }
        }
        if (worst > maxXDiff) {
          maxXDiff = worst;
          worstT = t;
          worstP = p;
        }
      }
    }
    System.out.printf("%-28s states=%-4d countDiffers=%-3d maxPhases plain=%d multi=%d  maxXDiff=%.3g @ T=%.0f P=%.0f%n",
        label, states, countDiffers, maxPhasesPlain, maxPhasesMulti, maxXDiff, worstT, worstP);
  }

  public static void main(String[] args) {
    sweep("CO2/n-hexane 70/30", new String[] {"CO2", "n-hexane"}, new double[] {0.7, 0.3});
    sweep("CO2/nc10 80/20", new String[] {"CO2", "nc10"}, new double[] {0.8, 0.2});
    sweep("C1/nC4/nc10 60/25/15", new String[] {"methane", "n-butane", "nc10"},
        new double[] {0.6, 0.25, 0.15});
    sweep("C1/nC4 50/50", new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5});
    sweep("N2/C1/nC7 10/60/30", new String[] {"nitrogen", "methane", "n-heptane"},
        new double[] {0.1, 0.6, 0.3});
    sweep("C1/CO2/nC7 40/20/40", new String[] {"methane", "CO2", "n-heptane"},
        new double[] {0.4, 0.2, 0.4});
    sweep("water/n-hexane 50/50", new String[] {"water", "n-hexane"}, new double[] {0.5, 0.5});
    sweep("water/methane 40/60", new String[] {"water", "methane"}, new double[] {0.4, 0.6});
    sweep("C1/C2/C3/nC4/nC5/nC6", new String[] {"methane", "ethane", "propane", "n-butane", "n-pentane",
        "n-hexane"}, new double[] {0.5, 0.15, 0.12, 0.1, 0.08, 0.05});
    // The liquid-liquid-vapour candidates: CO2 with the heaviest ends, cold.
    sweep("CO2/nc14 85/15", new String[] {"CO2", "nc14"}, new double[] {0.85, 0.15});
    sweep("CO2/nc20 90/10", new String[] {"CO2", "nc20"}, new double[] {0.90, 0.10});
    sweep("CO2/C1/nc10 40/30/30", new String[] {"CO2", "methane", "nc10"}, new double[] {0.4, 0.3, 0.3});
    sweep("C1/nC4/nc7/nc16", new String[] {"methane", "n-butane", "n-heptane", "nc16"},
        new double[] {0.45, 0.2, 0.2, 0.15});
    sweep("C2/nC5/nc12 30/30/40", new String[] {"ethane", "n-pentane", "nc12"},
        new double[] {0.3, 0.3, 0.4});
    sweep("H2S/C1/nC6 20/40/40", new String[] {"H2S", "methane", "n-hexane"},
        new double[] {0.2, 0.4, 0.4});
    sweep("N2/CO2/nC8 10/40/50", new String[] {"nitrogen", "CO2", "n-octane"},
        new double[] {0.1, 0.4, 0.5});
  }
}
