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
      // **The two Gibbs energies the residual is the difference of.** At the answer they are
      // equal by construction, so printing them says which of the two a port has to get right
      // rather than only that their difference is zero: a port whose solid is 3 J/mol off and
      // whose fluid compensates would reproduce the temperature and not the physics.
      row("g_fluid_J_per_mol", molarGibbs(system.getPhase(0)));
      // The fluid's whole property set, for the same reason as the solid's: the calibration
      // the freezing operation applies needs the liquid's Gibbs energy *and* its entropy at
      // the triple point, so a port that has one right and the other wrong would still land
      // on a temperature.
      row("fluid_v", system.getPhase(0).getMolarVolume());
      row("fluid_u_J_per_mol", molarProperty(system.getPhase(0), "u"));
      row("fluid_h_J_per_mol", molarProperty(system.getPhase(0), "h"));
      row("fluid_s_J_per_molK", molarProperty(system.getPhase(0), "s"));
      row("fluid_cv_J_per_molK", molarProperty(system.getPhase(0), "cv"));
      row("fluid_cp_J_per_molK", molarProperty(system.getPhase(0), "cp"));
      row("solid_phase_index", solidPhaseIndex(system));
      row("g_solid_J_per_mol", molarGibbs(solidPhase(system)));
      // The solid's whole property set, because the Gibbs energy is a rearrangement of the
      // same Helmholtz evaluation its `u`, `h`, `s` and `cp` come from: a port that has one
      // of them wrong usually has the others wrong too, and printing all of them says which.
      neqsim.thermo.phase.PhaseInterface solid = solidPhase(system);
      row("solid_v", solid == null ? Double.NaN : solid.getMolarVolume());
      row("solid_u_J_per_mol", molarProperty(solid, "u"));
      row("solid_h_J_per_mol", molarProperty(solid, "h"));
      row("solid_s_J_per_molK", molarProperty(solid, "s"));
      row("solid_cv_J_per_molK", molarProperty(solid, "cv"));
      row("solid_cp_J_per_molK", molarProperty(solid, "cp"));
      row("solid_z", solid == null ? Double.NaN : solid.getZ());
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

  /** The phase of the configured solid type, or null. */
  private static neqsim.thermo.phase.PhaseInterface solidPhase(SystemInterface system) {
    for (neqsim.thermo.phase.PhaseInterface phase : system.getPhases()) {
      if (phase != null && phase.getType() == neqsim.thermo.phase.PhaseType.SOLID) {
        return phase;
      }
    }
    return null;
  }

  private static double solidPhaseIndex(SystemInterface system) {
    for (int i = 0; i < system.getPhases().length; i++) {
      if (system.getPhases()[i] != null
          && system.getPhases()[i].getType() == neqsim.thermo.phase.PhaseType.SOLID) {
        return i;
      }
    }
    return Double.NaN;
  }

  /** A phase's molar Gibbs energy, or NaN where it has no moles to divide by. */
  private static double molarGibbs(neqsim.thermo.phase.PhaseInterface phase) {
    return molarProperty(phase, "g");
  }

  /** One of a phase's extensive properties, per mole. */
  private static double molarProperty(neqsim.thermo.phase.PhaseInterface phase, String which) {
    if (phase == null) {
      return Double.NaN;
    }
    try {
      double moles = phase.getNumberOfMolesInPhase();
      if (!(moles > 0.0)) {
        return Double.NaN;
      }
      double total =
          switch (which) {
            case "g" -> phase.getGibbsEnergy();
            case "u" -> phase.getInternalEnergy();
            case "h" -> phase.getEnthalpy();
            case "s" -> phase.getEntropy();
            case "cv" -> phase.getCv();
            case "cp" -> phase.getCp();
            default -> Double.NaN;
          };
      return total / moles;
    } catch (Throwable error) {
      return Double.NaN;
    }
  }
}
