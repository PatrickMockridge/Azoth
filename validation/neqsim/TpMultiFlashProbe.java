// What `TPmultiflash` does that `TPflash` does not, on a cubic system.
//
//     javac -proc:none -cp neqsim-f0c7436.jar TpMultiFlashProbe.java
//     java -cp .:neqsim-f0c7436.jar TpMultiFlashProbe
//
// Reached the way a user reaches it: `ThermodynamicOperations.TPflash()` ->
// `TPflash.run()` -> `runInternal()`, which constructs a `TPmultiflash` only under
// `system.doMultiPhaseCheck()`. So the probe is the same flash run twice, once with the flag
// and once without, and the question is whether the flag changes anything.
//
// Written to find out whether `eos.tp_multiflash` has physics to port at all, or whether on
// a cubic it is `eos.pt_flash` with extra branches that a hydrocarbon-only system never
// enters. `TpFlashSweep.java` is the two-phase driver beside it.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class TpMultiFlashProbe {

  static SystemInterface system(String cubic, double temperatureK, double pressureBar, String[] names,
      double[] moles) {
    SystemInterface fluid = cubic.equals("srk") ? new SystemSrkEos(temperatureK, pressureBar)
        : new SystemPrEos(temperatureK, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);
    return fluid;
  }

  static void one(String label, String cubic, double temperatureK, double pressureBar, String[] names,
      double[] moles) {
    System.out.printf("%n=== %s  %s  T=%.2f K  P=%.3f bar%n", label, cubic, temperatureK, pressureBar);
    for (int flag = 0; flag < 2; flag++) {
      boolean multi = flag == 1;
      SystemInterface fluid = system(cubic, temperatureK, pressureBar, names, moles);
      fluid.setMultiPhaseCheck(multi);
      try {
        new ThermodynamicOperations(fluid).TPflash();
      } catch (Exception ex) {
        System.out.printf("  multiPhaseCheck=%-5s -> %s: %s%n", multi, ex.getClass().getSimpleName(),
            ex.getMessage());
        continue;
      }
      StringBuilder line = new StringBuilder();
      for (int i = 0; i < fluid.getNumberOfPhases(); i++) {
        line.append(String.format(" [%s beta=%.10g]", fluid.getPhase(i).getType(),
            fluid.getPhase(i).getBeta()));
      }
      System.out.printf("  multiPhaseCheck=%-5s n=%d %s%n", multi, fluid.getNumberOfPhases(), line);
      if (multi && fluid.getNumberOfPhases() > 2) {
        for (int i = 0; i < fluid.getNumberOfPhases(); i++) {
          StringBuilder x = new StringBuilder();
          for (int j = 0; j < fluid.getPhase(i).getNumberOfComponents(); j++) {
            x.append(String.format(" %s=%.10g", fluid.getPhase(i).getComponent(j).getName(),
                fluid.getPhase(i).getComponent(j).getx()));
          }
          System.out.printf("      %s beta=%.10g x:%s%n", fluid.getPhase(i).getType(),
              fluid.getPhase(i).getBeta(), x);
        }
      }
    }
  }

  public static void main(String[] args) {
    // CO2 + a paraffin: the classic cubic liquid-liquid-vapour candidate.
    one("co2/hexane", "pr", 220.0, 30.0, new String[] {"CO2", "n-hexane"}, new double[] {0.7, 0.3});
    one("co2/hexane", "pr", 250.0, 50.0, new String[] {"CO2", "n-hexane"}, new double[] {0.6, 0.4});
    one("co2/decane", "pr", 230.0, 40.0, new String[] {"CO2", "nc10"}, new double[] {0.8, 0.2});
    // A gas-condensate shape: a light end, an intermediate and a heavy end.
    one("c1/nc4/nc10", "pr", 300.0, 60.0, new String[] {"methane", "n-butane", "nc10"},
        new double[] {0.6, 0.25, 0.15});
    one("c1/nc4/nc10", "pr", 250.0, 40.0, new String[] {"methane", "n-butane", "nc10"},
        new double[] {0.5, 0.3, 0.2});
    // Nitrogen and CO2 with a heavy end, low temperature.
    one("n2/c1/nc7", "pr", 200.0, 60.0, new String[] {"nitrogen", "methane", "n-heptane"},
        new double[] {0.1, 0.6, 0.3});
    // The envelope's own validated pair.
    one("c1/nc4", "pr", 200.0, 20.0, new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5});
    // Water, which is where the aqueous seed lives.
    one("water/nc6", "pr", 350.0, 20.0, new String[] {"water", "n-hexane"}, new double[] {0.5, 0.5});
    one("water/c1", "pr", 320.0, 50.0, new String[] {"water", "methane"}, new double[] {0.4, 0.6});
    // SRK for contrast.
    one("co2/hexane", "srk", 220.0, 30.0, new String[] {"CO2", "n-hexane"}, new double[] {0.7, 0.3});
  }
}
