// The freezing-point operation's layers, printed one at a time so each can be ported and
// checked.
//
// `FreezingPointTemperatureFlash` has two solid routes and the `instanceof` decides which: a
// `PhaseSolidHelmholtzEos` (para-hydrogen, argon) takes the reference-EOS route, and anything
// else - in practice the `PhasePureComponentSolid` that `addSolidPhase()` builds - takes a
// fluid-only TP flash first and then a residual over the component's tabulated solid
// properties. para-hydrogen is the first route and the one this prints, because both halves
// of it are already ported here: `eos.hydrogen_phase` is the Leachman equation the test uses
// and `eos.parahydrogen_solid_phase` is the solid.
//
// The states are NeqSim's own `FreezingPointTemperatureFlashTest`'s: the calibrated triple
// point at 13.6 K and 0.07042 bara, then the two pressures whose invalid bracket trials must
// not abort a valid melting state.
//
// A capture, not a test: run it against the pinned jar and commit what it prints.
//   javac -proc:none -cp neqsim-3.20.0.jar FreezingProbe.java
//   java -cp .:neqsim-3.20.0.jar FreezingProbe > captures/freezing_probe.tsv

import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.flashops.saturationops.FreezingPointResult;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemLeachmanEos;

public class FreezingProbe {

  /** One `key = value` row, which is the shape `tools/neqsim_layer_diff.py` reads. */
  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    double[][] states = {{13.6, 0.07042}, {13.8, 3.512706909625152}, {13.8, 18.76432785899884}};
    for (double[] state : states) {
      report(state[0], state[1]);
    }
  }

  private static void report(double temperatureK, double pressureBara) {
    System.out.printf("# para-hydrogen at T = %.15g K, P = %.15g bara%n", temperatureK, pressureBara);
    try {
      SystemInterface system = new SystemLeachmanEos(temperatureK, pressureBara, "para-hydrogen", true);
      system.setSolidPhaseCheck("para-hydrogen");
      ThermodynamicOperations operations = new ThermodynamicOperations(system);

      FreezingPointResult result = operations.freezingPointTemperatureFlashResult();

      row("converged", result.isConverged() ? 1.0 : 0.0);
      row("temperature_K", result.getTemperature("K"));
      row("iterations", result.getIterations());
      row("residual", result.getResidual());
      System.out.printf("component = %s%n", result.getComponentName());
      // What the walk left the system at, so the operation can be checked on the fluid side
      // and not only at its answer.
      row("system_temperature_K", system.getTemperature());
      row("system_pressure_bara", system.getPressure());
      row("phases", system.getNumberOfPhases());
      row("z", phaseZ(system));
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }

  /** The first phase's compressibility factor, or NaN where there is no phase to read it from. */
  private static double phaseZ(SystemInterface system) {
    try {
      return system.getPhase(0).getZ();
    } catch (Throwable error) {
      return Double.NaN;
    }
  }
}
