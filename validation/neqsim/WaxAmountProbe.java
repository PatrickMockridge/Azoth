// The wax amount on a fluid azoth can **build**, which the characterised one is not.
//
// `WaxProbe` is NeqSim's own test fluid, and it is made from TBP and plus fractions - so the
// fifteen components it contains exist only because `characterisePlusFraction` invented them,
// and a port has to be handed them as data. That is fine for capturing the *model*, and it is
// not a state a caller of azoth's API can express.
//
// **The 16 n-alkanes from n-hexane to nc24 carry `waxformer = 1` in `COMP.csv`** with a heat
// of fusion and a triple point beside it, so a fluid of *named* components can precipitate
// wax too - and the same `TPmultiflashWAX` runs on it. That makes this the family's
// representative: the amount-solve on a fluid both libraries can name.
//
//   javac -proc:none -cp neqsim-3.20.0.jar WaxAmountProbe.java
//   java -cp .:neqsim-3.20.0.jar WaxAmountProbe > captures/wax_amount_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class WaxAmountProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    double[] temperatures = {290.0, 275.0, 265.0, 255.0, 245.0, 235.0};
    for (double temperatureK : temperatures) {
      report(temperatureK, 5.0);
    }
  }

  private static SystemInterface build() {
    SystemInterface system = new SystemSrkEos(320.0, 5.0);
    system.addComponent("methane", 70.0);
    system.addComponent("n-heptane", 10.0);
    system.addComponent("nc14", 10.0);
    system.addComponent("nc20", 10.0);
    system.setMixingRule(2);
    system.addSolidComplexPhase("wax");
    system.setMultiphaseWaxCheck(true);
    system.setMultiPhaseCheck(true);
    system.init(0);
    system.init(1);
    return system;
  }

  private static void report(double temperatureK, double pressureBara) {
    System.out.printf("# named wax fluid at T = %.15g K, P = %.15g bara%n", temperatureK, pressureBara);
    try {
      SystemInterface fluid = build();
      for (int i = 0; i < fluid.getPhase(0).getNumberOfComponents(); i++) {
        ComponentInterface component = fluid.getPhase(0).getComponent(i);
        System.out.printf("component[%d] = %s%n", i, component.getName());
        row("component[" + i + "].mw", component.getMolarMass());
        row("component[" + i + "].tc", component.getTC());
        row("component[" + i + "].pc", component.getPC());
        row("component[" + i + "].acentric", component.getAcentricFactor());
        row("component[" + i + "].moles", component.getNumberOfmoles());
        row("component[" + i + "].wax_former", component.isWaxFormer() ? 1.0 : 0.0);
        row("component[" + i + "].heat_of_fusion", component.getHeatOfFusion());
        row("component[" + i + "].triple_point_temperature", component.getTriplePointTemperature());
      }

      fluid.setTemperature(temperatureK);
      fluid.setPressure(pressureBara);
      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      operations.TPflash();

      row("phases", fluid.getNumberOfPhases());
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        System.out.printf("phase[%d] = %s%n", p, phase.getType());
        row("phase[" + p + "].beta", phase.getBeta());
        for (int i = 0; i < phase.getNumberOfComponents(); i++) {
          String name = phase.getComponent(i).getName();
          row("phase[" + p + "].x[" + name + "]", phase.getComponent(i).getx());
          row(
              "fugcoef[" + p + "][" + name + "]",
              phase.getComponent(i).getFugacityCoefficient());
        }
      }

      // The material balance, as the other two probes print it: a wax fraction is one.
      int n = fluid.getPhase(0).getNumberOfComponents();
      for (int i = 0; i < n; i++) {
        double inPhases = 0.0;
        for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
          inPhases += fluid.getPhase(p).getBeta() * fluid.getPhase(p).getComponent(i).getx();
        }
        String name = fluid.getPhase(0).getComponent(i).getName();
        row("balance_error[" + name + "]", inPhases - fluid.getPhase(0).getComponent(i).getz());
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
