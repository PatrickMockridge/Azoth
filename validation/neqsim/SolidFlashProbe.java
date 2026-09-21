// The solid-phase flash's two routes, run side by side so each can be ported and checked.
//
// `ThermodynamicOperations.TPSolidflash()` builds a `SolidFlash1`. A plain
// `ThermodynamicOperations.TPflash()` on a system whose solid check is on reaches
// `Flash.solidPhaseFlash()`, which builds a `SolidFlash` instead - and that is the route
// `AsphalteneOnsetPressureFlash` takes, because it calls `TPflash` and not `TPSolidflash`.
//
// The two classes share the equilibrium equations and differ in the arithmetic around them:
// `SolidFlash` normalises each phase (`setXY` ends in `normalize()`) and damps its Newton
// step by `(iter+1)/(10+iter)`; `SolidFlash1` normalises nothing and damps by a line search
// on `Q`. The states below are where that shows.
//
// **Water is the representative**, and NeqSim says so itself: `ComponentSolid.fugcoef2`
// overrides `deltaCpSL = 37.12` for the component named `water` and computes it from the
// table for every other, so ice is the case the class was written around.
//
//   javac -proc:none -cp neqsim-f0c7436.jar SolidFlashProbe.java
//   java -cp .:neqsim-f0c7436.jar SolidFlashProbe > captures/solid_flash_probe.tsv

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
    for (double temperatureC : new double[] {0.0, -10.0, -20.0, 5.0}) {
      report(new String[] {"water", "methane"}, new double[] {0.5, 0.5}, temperatureC, 10.0);
    }
    // A solid that is not water: the heat-capacity difference comes off the tables instead of
    // the class's override, and the liquid density polynomial **is** zero on CO2's row, so the
    // reference-phase fallback is the branch this state measures.
    for (double temperatureC : new double[] {-80.0, -100.0}) {
      report(new String[] {"CO2", "methane"}, new double[] {0.5, 0.5}, temperatureC, 10.0);
    }
  }

  private static void report(String[] names, double[] feed, double temperatureC, double pressureBara) {
    for (String route : new String[] {"solid_flash_1", "solid_flash"}) {
      System.out.printf(
          "# %s: %s at %.15g C, %.15g bara%n",
          route, String.join("/", names), temperatureC, pressureBara);
      try {
        SystemInterface fluid = new SystemSrkEos(273.15 + temperatureC, pressureBara);
        for (int i = 0; i < names.length; i++) {
          fluid.addComponent(names[i], feed[i]);
        }
        fluid.setMixingRule(2);
        fluid.setMultiPhaseCheck(true);
        // The one component allowed to precipitate: `TPSolidflash` checks all of them when
        // this is not set, and refuses where more than one would.
        fluid.setSolidPhaseCheck(names[0]);
        fluid.init(0);
        fluid.init(1);

        ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
        if (route.equals("solid_flash_1")) {
          operations.TPSolidflash();
        } else {
          operations.TPflash();
        }

        row("temperature_K", fluid.getTemperature());
        row("phases", fluid.getNumberOfPhases());
        // **The material balance the reported state actually satisfies**, per component:
        // `sum_k beta_k x_ik` beside the feed's `z_i`. A state whose phases do not sum to one
        // is a state whose recovered feed is not the one it was given.
        double betaSum = 0.0;
        double[] recovered = new double[fluid.getPhase(0).getNumberOfComponents()];
        for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
          betaSum += fluid.getPhase(p).getBeta();
          for (int i = 0; i < recovered.length; i++) {
            recovered[i] += fluid.getPhase(p).getBeta() * fluid.getPhase(p).getComponent(i).getx();
          }
        }
        row("beta_sum", betaSum);
        for (int i = 0; i < recovered.length; i++) {
          row("recovered[" + fluid.getPhase(0).getComponent(i).getName() + "]", recovered[i]);
          row("feed[" + fluid.getPhase(0).getComponent(i).getName() + "]",
              fluid.getPhase(0).getComponent(i).getz());
        }
        for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
          PhaseInterface phase = fluid.getPhase(p);
          System.out.printf("phase[%d] = %s%n", p, phase.getType());
          row("phase[" + p + "].beta", phase.getBeta());
          double sum = 0.0;
          for (int i = 0; i < phase.getNumberOfComponents(); i++) {
            String name = phase.getComponent(i).getName();
            row("phase[" + p + "].x[" + name + "]", phase.getComponent(i).getx());
            sum += phase.getComponent(i).getx();
            // **The solid model's own coefficient and the numbers it is built from**, read
            // off whichever phase's component is a `ComponentSolid` - the tabulated route's
            // whole arithmetic, so a port can be checked against it without the flash between.
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
          // A sum that is not one is a phase whose composition the class never normalised.
          row("phase[" + p + "].x_sum", sum);
        }
      } catch (Exception error) {
        System.out.printf(
            "# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
      }
      System.out.println();
    }
  }
}
