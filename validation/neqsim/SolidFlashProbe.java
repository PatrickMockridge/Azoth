// The solid-phase flash's layers, printed so each can be ported and checked.
//
// `ThermodynamicOperations.TPSolidflash()` builds a `SolidFlash1`: a flash of the fluid phases
// plus **one pure solid**, whose component is the one `setSolidPhaseCheck(name)` selected. The
// solid's fugacity is `ComponentSolid.fugcoef2` - a high-pressure *liquid* reference phase of
// the host's own class, multiplied by an exponential of a fusion term, a heat-capacity term
// and a volume term, all three from the component's **tabulated** properties rather than from
// a correlation.
//
// **Water is the representative**, and NeqSim says so itself: `fugcoef2` overrides
// `deltaCpSL = 37.12` for the component named `water` and computes it from the table for every
// other, so ice is the case the class was written around.
//
//   javac -proc:none -cp neqsim-3.20.0.jar SolidFlashProbe.java
//   java -cp .:neqsim-3.20.0.jar SolidFlashProbe > captures/solid_flash_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.component.ComponentSolid;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class SolidFlashProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    // Below water's triple point at every pressure here, so the ice forms at all three.
    for (double temperatureC : new double[] {0.0, -10.0, -20.0}) {
      report(temperatureC, 10.0);
    }
  }

  private static void report(double temperatureC, double pressureBara) {
    System.out.printf("# water/methane at %.15g C, %.15g bara%n", temperatureC, pressureBara);
    try {
      SystemInterface fluid = new SystemSrkEos(273.15 + temperatureC, pressureBara);
      fluid.addComponent("water", 0.5);
      fluid.addComponent("methane", 0.5);
      fluid.setMixingRule(2);
      fluid.setMultiPhaseCheck(true);
      // The one component allowed to precipitate: `TPSolidflash` checks all of them when this
      // is not set, and refuses where more than one would.
      fluid.setSolidPhaseCheck("water");
      fluid.init(0);
      fluid.init(1);

      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      operations.TPSolidflash();

      row("temperature_K", fluid.getTemperature());
      row("phases", fluid.getNumberOfPhases());
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        System.out.printf("phase[%d] = %s%n", p, phase.getType());
        row("phase[" + p + "].beta", phase.getBeta());
        for (int i = 0; i < phase.getNumberOfComponents(); i++) {
          String name = phase.getComponent(i).getName();
          row("phase[" + p + "].x[" + name + "]", phase.getComponent(i).getx());
          // **The solid model's own coefficient and the numbers it is built from**, read off
          // whichever phase's component is a `ComponentSolid` - the tabulated route's whole
          // arithmetic, so a port can be checked against it without the flash in between.
          ComponentInterface component = phase.getComponent(i);
          if (component instanceof ComponentSolid solid) {
            row("solid_fugacity_coefficient[" + name + "]", solid.getFugacityCoefficient());
            row("solid_heat_of_fusion[" + name + "]", solid.getHeatOfFusion());
            row("solid_triple_point[" + name + "]", solid.getTriplePointTemperature());
            row(
                "solid_cp_liquid[" + name + "]",
                solid.getPureComponentCpLiquid(solid.getTriplePointTemperature()));
            row(
                "solid_cp_solid[" + name + "]",
                solid.getPureComponentCpSolid(solid.getTriplePointTemperature()));
            row(
                "solid_density_liquid[" + name + "]",
                solid.getPureComponentLiquidDensity(solid.getTriplePointTemperature()));
            row(
                "solid_density_solid[" + name + "]",
                solid.getPureComponentSolidDensity(solid.getTriplePointTemperature()));
            row("solid_do_check[" + name + "]", solid.doSolidCheck() ? 1.0 : 0.0);
          }
        }
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
