// The asphaltene onset's layers, printed so each can be ported and checked.
//
// `AsphalteneOnsetPressureFlash` is a *scan* rather than a solve: it starts at the reservoir
// pressure, steps down by `pressureStep` (5 bara), runs a `TPflash` with
// `setSolidPhaseCheck("asphaltene")` at each step, and bisects between the last pressure
// where the fluid was whole and the first where an asphaltene phase appeared. So the layer
// underneath it is `SolidFlash` over the `asphaltene` databank row - the same `ComponentSolid`
// arithmetic `eos.tp_solid_flash` ports - and what this probe prints is both the scan's own
// steps and the composition of the fluid the solid appears in.
//
// **NeqSim's own test for this class never reaches the onset.** Every assertion in
// `AsphalteneOnsetFlashTest` is `Double.isNaN(x) || x > 0`, which is true of every answer the
// method can return, and its fluid carries no asphaltene component at all - so the numbers
// below are the first external check the class has had.
//
//   javac -proc:none -cp neqsim-3.20.0.jar AsphalteneOnsetProbe.java
//   java -cp .:neqsim-3.20.0.jar AsphalteneOnsetProbe > captures/asphaltene_onset_probe.tsv

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class AsphalteneOnsetProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    // A live oil: light ends and a C7+ tail, with the generic asphaltene row at 2 mol% - the
    // component `AsphalteneCharacterization` adds and then overrides the molar mass of.
    String[] names = {"methane", "ethane", "propane", "n-butane", "n-heptane", "nC10", "asphaltene"};
    double[] feed = {0.30, 0.05, 0.05, 0.03, 0.05, 0.50, 0.02};
    for (double temperatureC : new double[] {60.0, 100.0, 140.0}) {
      report(names, feed, 273.15 + temperatureC);
    }
  }

  private static void report(String[] names, double[] feed, double temperatureK) {
    System.out.printf("# heavy oil at %.15g K, scanned from 150 bara%n", temperatureK);
    // **The solid check changes the fluid flash**, which is the finding this family turns on.
    // The same state, flashed three ways: plainly; with a solid candidate that *cannot*
    // precipitate; and with the asphaltene candidate. Where the plain flash splits, the two
    // checked ones return the feed as a single phase.
    for (double pressure : new double[] {150.0, 100.0, 60.0, 20.0, 5.0}) {
      for (String candidate : new String[] {null, "methane", "asphaltene"}) {
        SystemInterface fluid = fresh(names, feed, temperatureK, pressure, candidate);
        new ThermodynamicOperations(fluid).TPflash();
        StringBuilder line = new StringBuilder(String.format("  seed[%.12g]", pressure));
        line.append(" candidate=").append(candidate == null ? "none" : candidate).append(" ->");
        for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
          line.append(' ')
              .append(fluid.getPhase(p).getType())
              .append(':')
              .append(String.format("%.12g", fluid.getPhase(p).getBeta()))
              .append(" x_me=")
              .append(String.format("%.6g", fluid.getPhase(p).getComponent(0).getx()));
        }
        System.out.println(line);
      }
    }
    // **A fresh system at every pressure**, because a `TPflash` mutates the system it runs on
    // and NeqSim's own scan carries that mutation from step to step: a sweep that reuses one
    // system measures the search's path rather than the states it passes through. Both are
    // printed, because both are things a port has to reproduce.
    for (double pressure : new double[] {150.0, 100.0, 60.0, 40.0, 30.0, 25.0, 20.0, 15.0, 10.0, 5.0}) {
      try {
        SystemInterface fluid = fresh(names, feed, temperatureK, pressure);
        ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
        operations.TPflash();
        StringBuilder line = new StringBuilder(String.format("  state[%.12g] =", pressure));
        double solid = 0.0;
        for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
          PhaseInterface phase = fluid.getPhase(p);
          line.append(' ').append(phase.getType()).append(':').append(phase.getBeta());
          if (phase.getPhaseTypeName().toLowerCase().contains("asphaltene")
              || phase.getPhaseTypeName().toLowerCase().contains("solid")) {
            solid += phase.getBeta();
          }
        }
        System.out.println(line);
        row("  solid_at[" + pressure + "]", solid);
        if (solid > 0.0) {
          for (int i = 0; i < names.length; i++) {
            row("  x_oil[" + pressure + "][" + names[i] + "]",
                fluid.getPhase(0).getComponent(i).getx());
          }
        }
      } catch (Exception error) {
        System.out.printf("  # %.15g bara failed: %s%n", pressure, error.getMessage());
      }
    }

    // **The search from a pressure the fluid is whole at**, which is the only way to get a
    // non-degenerate answer: `startPressure` defaults to the system's own, so a scan begun
    // where the solid is already there reports the starting pressure and bisects nothing.
    for (double start : new double[] {60.0, 80.0}) {
      try {
        SystemInterface fluid = fresh(names, feed, temperatureK, start);
        ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
        double onset = operations.asphalteneOnsetPressure(start, 1.0);
        row("onset_pressure_bara_from[" + start + "]", onset);
        row("onset_found_from[" + start + "]", Double.isNaN(onset) ? 0.0 : 1.0);
        if (!Double.isNaN(onset)) {
          // The precipitated share the class reports, taken one bar below the onset.
          SystemInterface below = fresh(names, feed, temperatureK, onset - 1.0);
          ThermodynamicOperations ops2 = new ThermodynamicOperations(below);
          ops2.TPflash();
          double solid = 0.0;
          for (int p = 0; p < below.getNumberOfPhases(); p++) {
            String type = below.getPhase(p).getPhaseTypeName().toLowerCase();
            if (type.contains("asphaltene") || type.contains("solid")) {
              solid += below.getPhase(p).getBeta();
            }
          }
          row("solid_below_onset_from[" + start + "]", solid);
        }
      } catch (Exception error) {
        System.out.printf("# search from %.15g failed: %s: %s%n",
            start, error.getClass().getSimpleName(), error.getMessage());
      }
    }

    try {
      // And the temperature axis, at the pressure the system was built at.
      SystemInterface fluid = fresh(names, feed, temperatureK, 60.0);
      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      double onset = operations.asphalteneOnsetTemperature();
      row("onset_temperature_K", onset);
      row("onset_temperature_found", Double.isNaN(onset) ? 0.0 : 1.0);
    } catch (Exception error) {
      System.out.printf("# temperature search failed: %s: %s%n",
          error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }

  private static SystemInterface fresh(
      String[] names, double[] feed, double temperatureK, double pressureBara) {
    return fresh(names, feed, temperatureK, pressureBara, "asphaltene");
  }

  private static SystemInterface fresh(
      String[] names, double[] feed, double temperatureK, double pressureBara, String candidate) {
    SystemInterface fluid = new SystemSrkEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], feed[i]);
    }
    fluid.setMixingRule(2);
    if (candidate != null) {
      fluid.setSolidPhaseCheck(candidate);
    }
    fluid.init(0);
    fluid.init(1);
    return fluid;
  }
}
