// NeqSim's GE activity coefficients, read the only way that works.
//
//     javac -proc:none -cp neqsim-3.20.0.jar GeGamma.java
//     java -cp .:neqsim-3.20.0.jar GeGamma
//
// `PhaseGE.getActivityCoefficient` returns the component's **cached** `gamma` field, and
// only a GE evaluation sets it - `ComponentGENRTLmodifiedHV` assigns `gamma =
// exp(lngamma)` inside its nine-argument `getGamma`, which the phase calls from
// `getExcessGibbsEnergy`. A probe that skips that call reads whatever was in the field:
// measured, it returns `gamma = 1.000002` for water in a CO2/water mixture under the
// Wong-Sandler rule, which is a plausible-looking number and not the model's answer. So
// the `getExcessGibbsEnergy` line below is load-bearing, and it is why this driver exists
// beside `FlashTp.java` rather than inside it: every `getActivityCoefficient` in that
// file's `nrtl()` and `wilson()` reads a phase a flash has already initialised.
//
// `SystemSrkEos.setMixingRule` numbers: 3 = CLASSIC_HV, 4 = Huron-Vidal with DijT,
// 5 = Wong-Sandler.

import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermo.mixingrule.EosMixingRuleHandler;
import neqsim.thermo.component.ComponentGEInterface;
import neqsim.thermo.phase.PhaseGENRTLmodifiedHV;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;

public class GeGamma {
  /** One GE rule's activity coefficients at a fixed composition, read after evaluating. */
  static void one(String first, String second, double temperatureK, double pressureBar, int rule,
      String label) {
    SystemSrkEos fluid = new SystemSrkEos(temperatureK, pressureBar);
    fluid.addComponent(first, 0.5);
    fluid.addComponent(second, 0.5);
    fluid.init(0);
    fluid.setMixingRule(rule);

    PhaseInterface phase = fluid.getPhase(0);
    EosMixingRuleHandler.SRKHuronVidal2 mixRule =
        (EosMixingRuleHandler.SRKHuronVidal2) phase.getMixingRule();
    PhaseGENRTLmodifiedHV gePhase = (PhaseGENRTLmodifiedHV) mixRule.getGEPhase();
    int n = phase.getNumberOfComponents();
    gePhase.getExcessGibbsEnergy(phase, n, temperatureK, pressureBar, PhaseType.LIQUID);

    StringBuilder lnGamma = new StringBuilder();
    for (int i = 0; i < n; i++) {
      if (i > 0) {
        lnGamma.append(", ");
      }
      lnGamma.append(Math.log(((ComponentGEInterface) gePhase.getComponent(i)).getGamma()));
    }
    System.out.println(label);
    System.out.println("  " + first + "/" + second + "  T=" + temperatureK + " P=" + pressureBar
        + "  rule=" + rule);
    System.out.println("  ln_gamma   [" + lnGamma + "]");
    System.out.println("  HValpha[0][1]=" + mixRule.getHValphaParameter(0, 1)
        + "  HVgij[0][1]=" + mixRule.getHVDijParameter(0, 1)
        + "  HVgijT[0][1]=" + mixRule.getHVDijTParameter(0, 1)
        + "  HVgijT[1][0]=" + mixRule.getHVDijTParameter(1, 0));
  }

  public static void main(String[] args) {
    one("water", "ethanol", 300.0, 1.0, 3, "water/ethanol, CLASSIC_HV");
    one("water", "ethanol", 300.0, 1.0, 4, "water/ethanol, Huron-Vidal with DijT");
    one("CO2", "water", 350.0, 5.0, 3, "CO2/water, CLASSIC_HV");
    one("CO2", "water", 350.0, 5.0, 4, "CO2/water, Huron-Vidal with DijT");
    one("CO2", "water", 350.0, 5.0, 5, "CO2/water, Wong-Sandler");
  }
}
