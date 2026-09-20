import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemSoreideWhitson;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class AlphaSalinity {
  private static void report(double salt) throws Exception {
    SystemSoreideWhitson s = new SystemSoreideWhitson(298.0, 20.0);
    s.addComponent("nitrogen", 0.1, "mole/sec");
    s.addComponent("CO2", 0.2, "mole/sec");
    s.addComponent("methane", 0.3, "mole/sec");
    s.addComponent("ethane", 0.3, "mole/sec");
    s.addComponent("water", 0.1, "mole/sec");
    s.addSalinity(salt, "mole/sec");
    s.setTotalFlowRate(15, "mole/sec");
    s.setMixingRule(11);
    s.setTemperature(45.0, "C");
    s.setPressure(40.0, "bara");
    new ThermodynamicOperations(s).TPflash();
    s.initProperties();

    PhaseInterface aqueous = s.getPhase(1);
    double concentration =
        ((neqsim.thermo.phase.PhaseSoreideWhitson) aqueous).getSalinityConcentration();
    int w = aqueous.getComponent("water").getComponentNumber();
    neqsim.thermo.component.attractiveeosterm.AttractiveTermInterface term = aqueous.getComponent(w).getAttractiveTerm();
    double alpha = term.alpha(s.getTemperature());
    // The same component at zero salinity, for the comparison that matters.
    double alphaAtZero = new neqsim.thermo.component.attractiveeosterm.AttractiveTermSoreideWhitson(
        (neqsim.thermo.component.ComponentEosInterface) aqueous.getComponent(w))
        .alpha(s.getTemperature());
    // The CO2/water kij, which is the other thing the salinity reaches.
    neqsim.thermo.mixingrule.EosMixingRuleHandler handler =
        new neqsim.thermo.mixingrule.EosMixingRuleHandler();
    handler.getMixingRule(11, aqueous);
    java.lang.reflect.Field f = neqsim.thermo.mixingrule.EosMixingRuleHandler.class
        .getDeclaredField("intparam");
    f.setAccessible(true);
    double[][] intparam = (double[][]) f.get(handler);
    int co2 = aqueous.getComponent("CO2").getComponentNumber();
    double kij = intparam[co2][w];
    double tr = s.getTemperature() / aqueous.getComponent(co2).getTC();
    double bare = 0.989 * (-0.31092 * (1 + 0.15587 * Math.pow(concentration, 0.75))
        + 0.2358 * (1 + 0.17837 * Math.pow(concentration, 0.98)) * tr
        - 21.2566 * Math.exp(-Math.pow(6.7222, tr) - concentration));
    System.out.printf("  salt = %6.3f   c = %18.12g   alpha(water) = %.12f   kij(CO2,water)"
            + " = %20.12g   multipK = %s%n",
        salt, concentration, alpha, kij, (Math.abs(bare) > 1e-12 && kij != 0.1896)
            ? String.format("%.6f", kij / bare) : "not this branch");
  }

  public static void main(String[] args) throws Exception {
    System.out.println("water's alpha in a Soreide-Whitson aqueous phase, against salinity:");
    for (double salt : new double[] {0.0, 0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20}) {
      report(salt);
    }
  }
}
