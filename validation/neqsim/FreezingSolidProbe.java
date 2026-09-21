// The freezing point on the **tabulated** solid route, printed a layer at a time.
//
// `FreezingPointTemperatureFlash` has two solid routes and an `instanceof` chooses: a
// `PhaseSolidHelmholtzEos` takes the molar-Gibbs difference that `FreezingProbe` already
// captures, and every other solid takes the tabulated one, whose residual is the multiphase
// appearance condition over the fluid's own phases:
//
//   residual = ln z_k - ln( sum_p beta_p phi_solid_k / phi_kp )
//
// NeqSim loops the components its caller enabled a solid check for and keeps the **highest**
// freezing point, because a fluid at a temperature where any of its substances freezes has
// frozen. So the controlling component is part of the answer and is printed.
//
// **NeqSim's own test for this state only bounds the answer** -
// `FreezingPointTemperatureFlashTest.testLNGFreezingPointFlashAfterFluidOnlyTPFlash` asserts
// `90 < T < 220` - so the number below is what this probe is for.
//
// **`setSolidPhaseCheck("CO2")` also turns the multiphase check on** as a side effect, which is
// why the phase count is printed: the same call is the defect in
// `neqsim-issue-solid-check-collapses-the-flash.md`, and whether it bites here is a
// measurement rather than an assumption.
//
// A capture, not a test: run it against the pinned jar and commit what it prints.
//   javac -proc:none -cp neqsim-f0c7436.jar FreezingSolidProbe.java
//   java -cp .:neqsim-f0c7436.jar FreezingSolidProbe > captures/freezing_solid_probe.tsv

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.flashops.saturationops.FreezingPointResult;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class FreezingSolidProbe {

  /** The feed NeqSim's own test uses, with CO2 the only component enabled for a solid. */
  private static final String[] NAMES = {"CO2", "nitrogen", "methane", "ethane", "propane"};

  private static final double[] MOLES = {0.17, 1.1011731548, 0.324, 0.274, 0.0306};

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    for (double pressureBara : new double[] {5.0, 20.0, 50.0}) {
      report(pressureBara);
    }
  }

  private static void report(double pressureBara) {
    System.out.printf("# LNG feed at 120.35 K, %.15g bara, CO2 the solid candidate%n",
        pressureBara);
    try {
      SystemInterface fluid = new SystemSrkEos(120.35, pressureBara);
      for (int i = 0; i < NAMES.length; i++) {
        fluid.addComponent(NAMES[i], MOLES[i]);
      }
      fluid.setMixingRule(2);
      fluid.setSolidPhaseCheck("CO2");

      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      FreezingPointResult result = operations.freezingPointTemperatureFlashResult();

      row("converged", result.isConverged() ? 1.0 : 0.0);
      row("temperature_K", result.getTemperature("K"));
      row("iterations", result.getIterations());
      row("residual", result.getResidual());
      System.out.printf("component = %s%n", result.getComponentName());
      if (result.getFailureReason() != null) {
        System.out.printf("# failure = %s%n", result.getFailureReason());
      }

      // **The state the residual is built from.** The flash re-flashes at its answer, so the
      // phase set and the coefficients below are the ones the reported residual was formed
      // from - which is what a port has to reproduce, not the ones the feed arrived with.
      row("phases", fluid.getNumberOfPhases());
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        System.out.printf("phase[%d] = %s%n", p, phase.getType());
        row("phase[" + p + "].beta", phase.getBeta());
        if (phase.hasComponent("CO2")) {
          row("phase[" + p + "].ln_phi_CO2",
              phase.getComponent("CO2").getLogFugacityCoefficient());
        }
      }
      for (int i = 0; i < MOLES.length; i++) {
        row("z[" + NAMES[i] + "]", fluid.getPhase(0).getComponent(NAMES[i]).getz());
      }
      row("total_moles", fluid.getTotalNumberOfMoles());
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
