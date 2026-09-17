// What `PhaseGENRTL` reports for a liquid, for `eos.ge_nrtl_phase`.
//
//     javac -proc:none -cp neqsim-3.20.0.jar GeNrtlPhi.java
//     java -cp .:neqsim-3.20.0.jar GeNrtlPhi
//
// The phase's fugacity coefficient is `gamma_i * P0_i / P` - see `ComponentGE.fugcoef` -
// so it is the composition of the NRTL activity coefficient, already ported, with the
// Antoine vapour pressure, already ported. This prints all three so a port can be
// checked against each separately rather than only against their product.
//
// **Read from the GE phase, not from a flash.** At these states `SystemNRTL`'s TPflash
// converges to a single phase and leaves the SRK one, so `getPhase(0)` is a cubic with a
// cubic's fugacity coefficient - measured, `phi = 0.9746` where `gamma P0 / P` is
// `0.1693`. The GE phase is built by the system's constructor and is reachable before
// any flash, which is also how `GeGamma.java` reaches the mixing rules' inner phase.
//
// The evaluation before the read is the same requirement `GeGamma.java` documents: a
// component's `gamma` is a *cached* field, and `fugcoef` reads it.

import neqsim.thermo.system.SystemNRTL;
import neqsim.thermo.component.ComponentGEInterface;
import neqsim.thermo.phase.PhaseGEInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;

public class GeNrtlPhi {
  static void one(String first, String second, double temperatureK, double pressureBar) {
    SystemNRTL system = new SystemNRTL(temperatureK, pressureBar);
    system.addComponent(first, 0.5);
    system.addComponent(second, 0.5);
    system.createDatabase(true);
    system.setMixingRule("classic");
    system.init(0);

    PhaseInterface ge = null;
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      if (system.getPhase(p) instanceof PhaseGEInterface) {
        ge = system.getPhase(p);
      }
    }
    if (ge == null) {
      System.out.println("no GE phase was built");
      return;
    }

    ((PhaseGEInterface) ge).getExcessGibbsEnergy(ge, ge.getNumberOfComponents(), temperatureK,
        pressureBar, PhaseType.LIQUID);

    System.out.println(first + "/" + second + "  T=" + temperatureK + " P=" + pressureBar
        + "  (" + ge.getClass().getSimpleName() + ")");
    for (int i = 0; i < ge.getNumberOfComponents(); i++) {
      // `fugcoef` is what reads the cached gamma and builds the coefficient; without it
      // the coefficient is whatever the last evaluation left behind.
      ge.getComponent(i).fugcoef(ge);
      double gamma = ((ComponentGEInterface) ge.getComponent(i)).getGamma();
      double p0 = ge.getComponent(i).getAntoineVaporPressure(temperatureK);
      double phi = ge.getComponent(i).getFugacityCoefficient();
      System.out.println("  " + ge.getComponent(i).getName()
          + "  x=" + ge.getComponent(i).getx()
          + "  gamma=" + gamma
          + "  P0_bar=" + p0
          + "  phi=" + phi
          + "  gamma*P0/P=" + (gamma * p0 / pressureBar));
    }
  }

  public static void main(String[] args) {
    one("methanol", "water", 298.15, 1.0);
    one("ethanol", "water", 350.0, 1.0);
  }
}
