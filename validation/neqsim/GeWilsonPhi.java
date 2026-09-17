// What a `PhaseGEWilson` liquid *should* report, for `eos.ge_wilson_phase`.
//
//     javac -proc:none -cp neqsim-3.20.0.jar GeWilsonPhi.java
//     java -cp .:neqsim-3.20.0.jar GeWilsonPhi
//
// **Why the driver composes the coefficient rather than reading it.**
// `PhaseGEWilson.getExcessGibbsEnergy` sums `x_i ln(gamma_i)` and never writes `gamma`
// back to the component, so `ComponentGE.getGamma()` returns its initial zero and
// `ComponentGE.fugcoef` sets `phi_i = 0 * P0_i / P = 0`. Measured: every component of
// every state below reports `gamma = 0.0` and `phi = 0.0`. `ComponentGEWilson`'s
// eight-argument `getGamma` override - the one every other GE component implements - is
// stubbed to `return 0.0` as well.
//
// So `SystemGEWilson` has no fugacity coefficient to compare against, and the two factors
// it *does* compute are both public: `getWilsonActivityCoefficient` and
// `getAntoineVaporPressure`. This prints those, and their product over `P`, which is what
// `ComponentGE.fugcoef` would set if the field were written. The activity coefficient is
// the same number `validation/eos/n_butane_nc12_wilson_against_neqsim.json` validates
// `eos.wilson_activity_coefficients` against, read through the same direct method and for
// the same reason.
//
// What is NeqSim's here is both factors and the phase's group of components. What is the
// driver's is the one multiplication, and it is stated rather than hidden because it is
// the line upstream leaves out.

import neqsim.thermo.system.SystemGEWilson;
import neqsim.thermo.component.ComponentGEWilson;
import neqsim.thermo.phase.PhaseGEInterface;
import neqsim.thermo.phase.PhaseInterface;

public class GeWilsonPhi {

  static void one(String first, String second, double temperatureK, double pressureBar, double xFirst) {
    SystemGEWilson system = new SystemGEWilson(temperatureK, pressureBar);
    system.addComponent(first, xFirst);
    system.addComponent(second, 1.0 - xFirst);
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

    System.out.println(first + "/" + second + "  T=" + temperatureK + " P=" + pressureBar
        + "  x=" + xFirst);
    for (int i = 0; i < ge.getNumberOfComponents(); i++) {
      ComponentGEWilson component = (ComponentGEWilson) ge.getComponent(i);
      double gamma = component.getWilsonActivityCoefficient(ge);
      double p0 = component.getAntoineVaporPressure(temperatureK);
      component.fugcoef(ge);
      System.out.println("  " + component.getName()
          + "  x=" + component.getx()
          + "  gamma=" + gamma
          + "  P0_bar=" + p0
          + "  gamma*P0/P=" + (gamma * p0 / pressureBar)
          + "  [the phase reports phi=" + component.getFugacityCoefficient() + "]");
    }
  }

  public static void main(String[] args) {
    // Solvent-tagged pairs only: for a component NeqSim tags a Henry's-law solute the
    // phase's own `fugcoef` takes the other branch, so a Raoult comparison would not be
    // a comparison. `n-butane` is one of those, which is why the pair the Wilson
    // *activity* model is validated on (`validation/eos/n_butane_nc12_...`) is not the
    // pair this phase can be stated against.
    // And with both members carrying a real vapour-pressure row. `nc12` and every heavier
    // pseudo-component share one default row in NeqSim's `COMP.csv` - `(-7.76451, 1.45838,
    // -2.7758, -1.23303, 0.0)`, identical for `nc12`, `nc14`, `nc16`, `nc20` and even
    // `nc6-benzene` - so a case stated over one of those is stated over filler.
    one("nc10", "nc12", 298.15, 1.0, 0.5);
    one("n-octane", "nc10", 350.0, 1.0, 0.4);
    one("n-heptane", "n-nonane", 320.0, 1.0, 0.6);
  }
}
