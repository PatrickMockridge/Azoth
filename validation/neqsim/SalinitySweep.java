import neqsim.thermo.mixingrule.EosMixingRuleHandler;
import neqsim.thermo.system.SystemSoreideWhitson;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class SalinitySweep {
  private static final double[] CHABAB = {0.43575155, -5.766906744e-2, 8.26464849e-3, 1.29539193e-3,
      -1.6698848e-3, -0.47866096};

  public static void main(String[] args) throws Exception {
    System.out.printf("%14s %22s %22s %22s%n", "NaCl mol/sec", "concentration mol/kg",
        "kij(CO2,water)", "multipK implied");
    for (double salt : new double[] {0.0, 0.02, 0.05, 0.10, 0.14, 0.18, 0.25, 0.40, 1.00}) {
      SystemSoreideWhitson s = new SystemSoreideWhitson(298.0, 20.0);
      s.addComponent("nitrogen", 0.1, "mole/sec");
      s.addComponent("CO2", 0.2, "mole/sec");
      s.addComponent("methane", 0.3, "mole/sec");
      s.addComponent("ethane", 0.3, "mole/sec");
      s.addComponent("water", 0.1, "mole/sec");
      s.addSalinity(salt, "mole/sec");
      s.setTotalFlowRate(15, "mole/sec");
      s.setMixingRule(11);
      new ThermodynamicOperations(s).TPflash();
      s.initProperties();
      double concentration =
          ((neqsim.thermo.phase.PhaseSoreideWhitson) s.getPhase(0)).getSalinityConcentration();

      EosMixingRuleHandler handler = new EosMixingRuleHandler();
      handler.getMixingRule(11, s.getPhase(0));
      java.lang.reflect.Field f = EosMixingRuleHandler.class.getDeclaredField("intparam");
      f.setAccessible(true);
      double[][] intparam = (double[][]) f.get(handler);
      double kij = intparam[1][0];

      // What the correlation gives with multipK = 1.0, so the multiplier is recoverable.
      int water = s.getPhase(1).getComponent("water").getComponentNumber();
      int co2 = s.getPhase(1).getComponent("CO2").getComponentNumber();
      EosMixingRuleHandler h2 = new EosMixingRuleHandler();
      h2.getMixingRule(11, s.getPhase(1));
      double[][] p2 = (double[][]) f.get(h2);
      kij = p2[co2][water];
      double tr = 318.15 / s.getPhase(1).getComponent(co2).getTC();
      double bare = 0.989 * (-0.31092 * (1 + 0.15587 * Math.pow(concentration, 0.75))
          + 0.2358 * (1 + 0.17837 * Math.pow(concentration, 0.98)) * tr
          - 21.2566 * Math.exp(-Math.pow(6.7222, tr) - concentration));
      System.out.printf("%14.4f %22.12g %22.12g %22.12g%n", salt, concentration, kij, kij / bare);
    }
  }
}
