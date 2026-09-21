// The binary interaction parameters NeqSim's classic mixing rule actually uses for
// methane/n-butane, so a mixture difference between the two libraries can be traced to the
// data rather than to the algebra.
//
//     javac -proc:none -cp neqsim-f0c7436.jar KijProbe.java
//     java -cp .:neqsim-f0c7436.jar KijProbe
//
// NeqSim holds bar throughout.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

public class KijProbe {

  public static void main(String[] args) {
    for (String rule : new String[] {"classic", "no"}) {
      SystemInterface fluid = new SystemPrEos(320.0, 20.0);
      fluid.addComponent("methane", 0.5);
      fluid.addComponent("n-butane", 0.5);
      try {
        fluid.setMixingRule(rule);
      } catch (Exception ex) {
        System.out.println(rule + ": " + ex.getMessage());
        continue;
      }
      fluid.init(0);
      fluid.init(1);
      System.out.printf("%n=== mixing rule %s%n", rule);
      neqsim.thermo.mixingrule.EosMixingRulesInterface mr =
          ((neqsim.thermo.phase.PhaseEos) fluid.getPhase(0)).getMixingRule();
      double[][] kij = mr.getBinaryInteractionParameters();
      for (int i = 0; i < kij.length; i++) {
        System.out.printf("  kij row %d = %s%n", i, java.util.Arrays.toString(kij[i]));
      }
      for (int i = 0; i < 2; i++) {
        System.out.printf("  %-10s Tc=%.6f K  Pc=%.6f bar  omega=%.6f%n",
            fluid.getPhase(0).getComponent(i).getComponentName(),
            fluid.getPhase(0).getComponent(i).getTC(),
            fluid.getPhase(0).getComponent(i).getPC(),
            fluid.getPhase(0).getComponent(i).getAcentricFactor());
      }
    }
  }
}
