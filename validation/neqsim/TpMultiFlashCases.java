// The states where `setMultiPhaseCheck(true)` changes a cubic flash, in full.
//
//     javac -proc:none -cp neqsim-3.20.0.jar TpMultiFlashCases.java
//     java -cp .:neqsim-3.20.0.jar TpMultiFlashCases
//
// `TpMultiFlashSweep.java` grids 16 mixtures over 207 states each and finds the flag changes
// the answer at exactly three kinds of place:
//
//   1. **The trivial-solution escape.** Methane/n-butane 50/50 at 370 K / 40 bar: plain
//      `TPflash` converges to `x = y = z` and reports a single phase, and the multiflash's
//      stability check finds the feed is unstable and adds a liquid at beta ~ 0.01. This is
//      the one state in 207 for that mixture, and it is what `eos.stability_test` is for.
//   2. **The third phase.** CO2/methane/nc10 40/30/30 (12 of 207 states, T 180-220 K) and
//      N2/CO2/n-octane 10/40/50 (9 of 207, T 180-200 K) come back with three phases at states
//      where plain `TPflash` comes back with two. Both mixtures reach three phases under the
//      plain flash somewhere too - the sweep prints `maxPhases plain=3` for exactly these two
//      and no other - so the flag moves *which* states split three ways rather than adding the
//      capability. Both are CO2-rich, cold and at low pressure, which is where a third phase
//      belongs, and neither ever reaches four.
//   3. **The aqueous seed.** Water/n-hexane 50/50 differs in phase *count* at 139 of 207
//      states, and that is the only mixture for which the difference is common. It is
//      `seedHydrocarbonLiquidFromFeed`, gated on a water component.
//
// Every other mixture - 3,300-odd states - is either identical to the last bit or differs by
// less than 3e-9, three mixtures by exactly zero. So the flag is not a second flash on a
// hydrocarbon system; it is the two-phase flash plus a stability check.
//
// This driver prints the states above in the fields a port would report, both ways, so a
// divergence is attributable. `TpMultiFlashProbe.java` is the hand-picked version that came
// first and asked whether there was anything to find at all.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class TpMultiFlashCases {

  static SystemInterface build(String[] names, double[] moles, double temperatureK, double pressureBar,
      boolean multi) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);
    fluid.setMultiPhaseCheck(multi);
    return fluid;
  }

  static void one(String label, String[] names, double[] moles, double temperatureK,
      double pressureBar) {
    System.out.printf("%n=== %s   T=%.4f K   P=%.6f bar%n", label, temperatureK, pressureBar);
    for (int flag = 0; flag < 2; flag++) {
      boolean multi = flag == 1;
      SystemInterface fluid = build(names, moles, temperatureK, pressureBar, multi);
      try {
        new ThermodynamicOperations(fluid).TPflash();
      } catch (Exception ex) {
        System.out.printf("  multiPhaseCheck=%-5s -> %s: %s%n", multi, ex.getClass().getSimpleName(),
            ex.getMessage());
        continue;
      }
      int n = fluid.getNumberOfPhases();
      System.out.printf("  multiPhaseCheck=%-5s  n=%d%n", multi, n);
      for (int i = 0; i < n; i++) {
        try {
          StringBuilder x = new StringBuilder();
          for (int j = 0; j < fluid.getPhase(i).getNumberOfComponents(); j++) {
            x.append(String.format("  %s=%.17g", fluid.getPhase(i).getComponent(j).getName(),
                fluid.getPhase(i).getComponent(j).getx()));
          }
          System.out.printf("     [%d] %-8s beta=%.17g  Z=%.17g  rho=%.10g%s%n", i,
              fluid.getPhase(i).getType(), fluid.getPhase(i).getBeta(), fluid.getPhase(i).getZ(),
              fluid.getPhase(i).getDensity(), x);
        } catch (RuntimeException ex) {
          // `getNumberOfPhases()` and the phase array disagree once a phase is removed, so
          // the count above is not always a readable range. That is a defect upstream, and
          // it is reported rather than swallowed.
          System.out.printf("     [%d] unreadable: %s%n", i, ex.getMessage());
        }
      }
    }
  }

  public static void main(String[] args) {
    // 1. The trivial-solution escape, and the one C1/nC4 state that differs.
    one("C1/nC4 50/50", new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5}, 370.0, 40.0);
    one("C1/nC4 50/50", new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5}, 360.0, 40.0);
    // 2. The third phase, CO2/methane/nc10.
    one("CO2/C1/nc10 40/30/30", new String[] {"CO2", "methane", "nc10"},
        new double[] {0.4, 0.3, 0.3}, 180.0, 2.0);
    one("CO2/C1/nc10 40/30/30", new String[] {"CO2", "methane", "nc10"},
        new double[] {0.4, 0.3, 0.3}, 200.0, 10.0);
    one("CO2/C1/nc10 40/30/30", new String[] {"CO2", "methane", "nc10"},
        new double[] {0.4, 0.3, 0.3}, 220.0, 20.0);
    // 2b. The third phase, N2/CO2/n-octane.
    one("N2/CO2/nC8 10/40/50", new String[] {"nitrogen", "CO2", "n-octane"},
        new double[] {0.1, 0.4, 0.5}, 190.0, 1.0);
    one("N2/CO2/nC8 10/40/50", new String[] {"nitrogen", "CO2", "n-octane"},
        new double[] {0.1, 0.4, 0.5}, 180.0, 10.0);
    // 3. The aqueous seed.
    one("water/nC6 50/50", new String[] {"water", "n-hexane"}, new double[] {0.5, 0.5}, 350.0, 20.0);
  }
}
