// NeqSim's own fugacity coefficients for `eos.ge_wilson_phase`, after #3774.
//
//     javac -proc:none -cp neqsim-3.20.0.jar GeWilsonPhi.java
//     java -cp .:neqsim-3.20.0.jar GeWilsonPhi
//
// **This driver used to compose the coefficient. It no longer has to.**
// `PhaseGEWilson.getExcessGibbsEnergy` summed `x_i ln(gamma_i)` and never wrote `gamma` back
// to the component, so `ComponentGE.getGamma()` returned its initial zero and
// `ComponentGE.fugcoef` set `phi_i = 0 * P0_i / P = 0`; `ComponentGEWilson`'s eight-argument
// `getGamma` override was stubbed to `return 0.0` as well. With nothing published, the
// driver multiplied `getWilsonActivityCoefficient` by `getAntoineVaporPressure` and divided
// by `P`, which is what `fugcoef` would have set had the field been written.
//
// NeqSim commit `c5ec5fb` (PR #3774, closing #3770) stores `lngamma` and `gamma` in
// `getWilsonActivityCoefficient` and delegates the eight-argument `getGamma` to it. The
// sequence below is the library's own, not the driver's: `system.init(1)` reaches
// `PhaseGE.init(totalNumberOfMoles, numberOfComponents, initType, pt, beta)`, which calls the
// five-argument `getExcessGibbsEnergy` - the only overload that evaluates
// `getWilsonActivityCoefficient` and therefore the only one that publishes anything. Then
// `fugcoef` reads the field it published. Run `fugcoef` without the phase having initialised
// and the zero is back, which is exactly the defect #3774 fixed and the reason this driver
// exists.
//
// The pair choices are the same as before and for the same reason: solvent-tagged members
// only, because `ComponentGE.fugcoef` takes a Henry's-law branch for a solute, and members
// with a real vapour-pressure row. `nC12`'s row used to carry coefficients shared with
// `nc14`, `nc16`, `nc20` and `nc6-benzene` - filler - and commit `83b64e5` (PR #3775) has
// since replaced them with `none`, so that state's second component now reports `NaN` and is
// not an oracle for anything.

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
    // The library's sequence: this is what publishes `gamma` on every component.
    system.init(1);

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

    System.out.printf("%s/%s  T=%.4f P=%.4f  x=%.4f%n", first, second, temperatureK, pressureBar,
        xFirst);
    for (int i = 0; i < ge.getNumberOfComponents(); i++) {
      ComponentGEWilson component = (ComponentGEWilson) ge.getComponent(i);
      component.fugcoef(ge);
      double phi = component.getFugacityCoefficient();
      System.out.printf("  %-10s x=%.17g gamma=%.17g ln_gamma=%.17g P0_bar=%.17g phi=%.17g ln_phi=%.17g%n",
          component.getName(), component.getx(), component.getGamma(),
          component.getLnGamma(), component.getAntoineVaporPressure(temperatureK), phi,
          Math.log(phi));
    }
  }

  public static void main(String[] args) {
    one("nc10", "nc12", 298.15, 1.0, 0.5);
    one("n-octane", "nc10", 350.0, 1.0, 0.4);
    one("n-heptane", "n-nonane", 320.0, 1.0, 0.6);
  }
}
