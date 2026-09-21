import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemSoreideWhitson;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

// Which of a Soreide-Whitson phase's salinity-dependent quantities actually respond.
//
//     javac -proc:none -cp neqsim-f0c7436.jar AlphaSalinity.java
//     java -cp .:neqsim-f0c7436.jar AlphaSalinity
//
// The tranche's third F instrument, and the one that separates two things that look alike.
// Over eight brines from zero to `7.59 mol/kg`:
//
//   * **`alpha(water)` is `1.576491803452` at every one of them** - the value a freshly
//     constructed `AttractiveTermSoreideWhitson` gives with no salinity set. The push in
//     `SystemSoreideWhitson.calcSalinity` is guarded by `comp.getClass().getName().equals(
//     "...ComponentEosInterface")`, which compares a concrete class name to an interface
//     name, so it is false for every component and `setSalinityFromPhase` never runs.
//   * **`ln phi(CO2)` in the aqueous phase moves from `4.517210` to `6.005369`.** The
//     salinity-dependent CO2/water correlation *is* wired up.
//
// **And the second of those is why a first reading of this was wrong.** The correlation
// lives in `EosMixingRuleHandler.WhitsonSoreideMixingRule.calcA`, where it multiplies `aij`;
// it never writes the `intparam` matrix, which is where `getkij` reads and where an earlier
// probe looked. Reading that matrix showed a constant `0.1896` at every salinity - the
// interaction table's value - and looked like a dead correlation. It is not one, and this
// probe measures the quantity that carries it instead.
//
// `calcA` also reads the salinity **only when water's mole fraction exceeds 0.8**, so the
// brine's concentration reaches the aqueous phase and not the gas. That is coherent - the
// concentration is set on the aqueous phase by `calcSalinity` and means a brine property.

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
    // **What the salinity actually reaches**: `calcA`'s aqueous branch multiplies `aij` by
    // `(1 - kijWhitsonSoreideAqueous(...))`, which never touches `intparam` - so the CO2
    // fugacity coefficient in the brine is the observable, not that matrix.
    double lnPhiCO2 = Math.log(aqueous.getComponent(co2).getFugacityCoefficient());
    System.out.printf("  salt = %6.3f   c = %18.12g   alpha(water) = %.12f   "
            + "ln phi(CO2, aqueous) = %.12f%n", salt, concentration, alpha, lnPhiCO2);
  }

  public static void main(String[] args) throws Exception {
    System.out.println("water's alpha in a Soreide-Whitson aqueous phase, against salinity:");
    for (double salt : new double[] {0.0, 0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20}) {
      report(salt);
    }
  }
}
