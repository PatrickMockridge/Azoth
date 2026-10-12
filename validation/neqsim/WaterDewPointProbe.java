// Water dew point and wax appearance, for `eos.water_dew_point_temperature`,
// `eos.water_dew_point_multiphase`, `eos.water_dew_point_line` and
// `eos.wax_appearance_temperature`.
//
//     javac -proc:none -cp neqsim-f0c7436.jar WaterDewPointProbe.java
//     java -cp .:neqsim-f0c7436.jar WaterDewPointProbe > captures/water_dew_point_probe.tsv
//
// # What it drives
//
// `ThermodynamicOperations.waterDewPointTemperatureFlash()`,
// `...waterDewPointTemperatureMultiphaseFlash()`, `...waterDewPointLine(min, max)` and
// `...calcWAT()` - the four public entry points, so nothing here reaches a class the
// caller of the library could not. All four are plain cubics: this library can flash
// neither an associating water nor a MEG, the divergence `eos.hydrate_inhibitor_concentration`
// already records.
//
// # The gas
//
// NeqSim's own `WaterDewPointTest` fixture: CO2 0.02, nitrogen 0.01, methane 0.9 - x,
// ethane 0.05, propane 0.01, i-butane 0.005, n-butane 0.005 and water x, on a
// `SystemSrkEos` at 70 bara, with `x` in mole fraction. That test's two assertions are
// the only numbers NeqSim publishes for these operations - `-12.0953140 C` at 35 ppm on
// the prescribed route - and the rest of this probe's rows exist because the multiphase
// route has **no plain-cubic assertion upstream at all**: its own test runs on
// `SystemSrkCPAstatoil`, which this library cannot build.
//
// # The two routes disagree, and that is the point
//
// The prescribed route states its second phase: it seeds a nearly pure water phase and
// steps the temperature to the fugacity equality. The multiphase route asks a multiphase
// flash to *find* an aqueous phase and brackets on where it first appears. Upstream keeps
// them as separate entry points; the port does too, and this capture is what says whether
// they answer the same question.
//
// # The wax row
//
// `WATcalc` is a *wax* appearance temperature - its own javadoc says so - and its fixture
// in `WaxFlashTest` needs an `addTBPWax()` plus fraction. This probe uses a **nameable**
// feed instead, because the port's wax formers are the component databank's own
// `waxformer` rows (`n-hexane`, `n-heptane`, `nc11`-`nc24`); the plus-fraction route is a
// feed construction this library does not build and the case records that.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.util.database.NeqSimDataBase;

public class WaterDewPointProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  /** The `WaterDewPointTest` gas, on a plain cubic. `ppmWater` is a mole fraction in ppm. */
  static SystemInterface gas(double ppmWater, double pressureBara) {
    double x = ppmWater * 1e-6;
    SystemInterface system = new SystemSrkEos(270.15, pressureBara);
    system.addComponent("CO2", 0.02);
    system.addComponent("nitrogen", 0.01);
    system.addComponent("methane", 0.9 - x);
    system.addComponent("ethane", 0.05);
    system.addComponent("propane", 0.01);
    system.addComponent("i-butane", 0.005);
    system.addComponent("n-butane", 0.005);
    system.addComponent("water", x);
    system.createDatabase(true);
    system.setMixingRule(2);
    return system;
  }

  /**
   * The wax feed. `WaxAmountProbe`'s own - methane 70, n-heptane 10, nc14 10 and nc20 10, the
   * three heavy ends being databank `waxformer` rows. That probe's capture shows this feed
   * **is** three-phase at 265 K and two-phase at 275 K, so the onset is bracketable here.
   *
   * A *lean* gas is not: the same probe's richer cousins show the wax comes from the heavy
   * ends' mole fraction, and a few hundred ppm of one former has no onset at all.
   */
  static SystemInterface waxFeed() {
    SystemInterface system = new SystemSrkEos(320.0, 5.0);
    system.addComponent("methane", 70.0);
    system.addComponent("n-heptane", 10.0);
    system.addComponent("nc14", 10.0);
    system.addComponent("nc20", 10.0);
    system.setMixingRule(2);
    system.addSolidComplexPhase("wax");
    system.setMultiphaseWaxCheck(true);
    system.setMultiPhaseCheck(true);
    system.init(0);
    system.init(1);
    return system;
  }

  /**
   * A heavier variant of the same feed, so the onset lands well above the region where the
   * trial flash goes unstable: the tree's own feed onsets at about 270 K, which is exactly
   * where `WaxAmountProbe`'s capture shows the flash becoming marginal.
   */
  static SystemInterface richWaxFeed() {
    SystemInterface system = new SystemSrkEos(320.0, 1.0);
    system.addComponent("methane", 40.0);
    system.addComponent("n-heptane", 10.0);
    system.addComponent("nc14", 10.0);
    system.addComponent("nc20", 20.0);
    system.addComponent("nc24", 20.0);
    system.setMixingRule(2);
    system.addSolidComplexPhase("wax");
    system.setMultiphaseWaxCheck(true);
    system.setMultiPhaseCheck(true);
    system.init(0);
    system.init(1);
    return system;
  }

  /** A lean gas with one former in it: the state `calcWAT` refuses rather than answers. */
  static SystemInterface leanWaxFeed() {
    SystemInterface system = new SystemSrkEos(273.0 + 30, 50.0);
    system.addComponent("CO2", 0.018);
    system.addComponent("nitrogen", 0.333);
    system.addComponent("methane", 96.702);
    system.addComponent("ethane", 1.773);
    system.addComponent("propane", 0.496);
    system.addComponent("i-butane", 0.099);
    system.addComponent("n-butane", 0.115);
    system.addComponent("i-pentane", 0.004);
    system.addComponent("n-pentane", 0.024);
    system.addComponent("n-heptane", 0.324);
    system.createDatabase(true);
    system.setMixingRule(2);
    system.addSolidComplexPhase("wax");
    system.setMultiphaseWaxCheck(true);
    return system;
  }

  /** `WaxFlashTest`'s fixture verbatim - the plus fraction, the characterisation, the TBP wax. */
  static SystemInterface waxFixture() {
    NeqSimDataBase.setCreateTemporaryTables(true);
    SystemInterface system = new SystemSrkEos(273.0 + 30, 50.0);
    system.addComponent("CO2", 0.018);
    system.addComponent("nitrogen", 0.333);
    system.addComponent("methane", 96.702);
    system.addComponent("ethane", 1.773);
    system.addComponent("propane", 0.496);
    system.addComponent("i-butane", 0.099);
    system.addComponent("n-butane", 0.115);
    system.addComponent("i-pentane", 0.004);
    system.addComponent("n-pentane", 0.024);
    system.addComponent("n-heptane", 0.324);
    system.addPlusFraction("C9", 0.095, 207.0 / 1000.0, 0.8331);
    system.getCharacterization().characterisePlusFraction();
    system.getWaxModel().addTBPWax();
    system.createDatabase(true);
    system.setMixingRule(2);
    system.addSolidComplexPhase("wax");
    system.setMultiphaseWaxCheck(true);
    NeqSimDataBase.setCreateTemporaryTables(false);
    return system;
  }

  static void prescribed(String label, double ppmWater, double pressureBara) {
    SystemInterface system = gas(ppmWater, pressureBara);
    System.out.printf("# %s%n", label);
    row("water_mole_fraction", ppmWater * 1e-6);
    row("pressure_bara", system.getPressure());
    try {
      new ThermodynamicOperations(system).waterDewPointTemperatureFlash();
      row("temperature_K", system.getTemperature("K"));
      row("temperature_C", system.getTemperature("C"));
    } catch (Exception failure) {
      System.out.printf("error = %s%n", failure.toString());
    }
    System.out.println();
  }

  static void multiphase(String label, double ppmWater, double pressureBara) {
    SystemInterface system = gas(ppmWater, pressureBara);
    System.out.printf("# %s%n", label);
    row("water_mole_fraction", ppmWater * 1e-6);
    row("pressure_bara", system.getPressure());
    try {
      new ThermodynamicOperations(system).waterDewPointTemperatureMultiphaseFlash();
      row("temperature_K", system.getTemperature("K"));
      row("temperature_C", system.getTemperature("C"));
      row("has_aqueous_at_answer", system.hasPhaseType("aqueous") ? 1.0 : 0.0);
      row("number_of_phases_at_answer", system.getNumberOfPhases());
    } catch (Exception failure) {
      System.out.printf("error = %s%n", failure.toString());
    }
    System.out.println();
  }

  static void line(double ppmWater, double minimumPressure, double maximumPressure) {
    SystemInterface system = gas(ppmWater, maximumPressure);
    System.out.printf("# water dew point line, %g to %g bara%n", minimumPressure, maximumPressure);
    row("water_mole_fraction", ppmWater * 1e-6);
    row("minimum_pressure_bara", minimumPressure);
    row("maximum_pressure_bara", maximumPressure);
    try {
      ThermodynamicOperations operations = new ThermodynamicOperations(system);
      operations.waterDewPointLine(minimumPressure, maximumPressure);
      double[][] points = operations.getOperation().getPoints(0);
      for (int i = 0; i < points[0].length; i++) {
        row(String.format("point_%d_temperature_K", i), points[0][i]);
        row(String.format("point_%d_pressure_bara", i), points[1][i]);
      }
    } catch (Exception failure) {
      System.out.printf("error = %s%n", failure.toString());
    }
    System.out.println();
  }

  /**
   * A plain TP flash down the same temperatures `WATcalc`'s bracket walks, printing the wax
   * phase fraction beside the phases. `WATcalc` refuses these states; this is what it refuses.
   */
  static void waxSweep(double pressureBara) {
    System.out.printf("# plain TP flash sweep of the named wax feed at %g bara%n", pressureBara);
    row("pressure_bara", pressureBara);
    for (double temperatureK : new double[] {300.0, 290.0, 280.0, 275.0, 270.0, 265.0, 260.0,
        255.0, 245.0, 235.0}) {
      SystemInterface system = waxFeed();
      system.setPressure(pressureBara);
      system.setTemperature(temperatureK);
      System.out.printf("temperature_K = %.15g%n", temperatureK);
      try {
        new ThermodynamicOperations(system).TPflash();
        row("phases", system.getNumberOfPhases());
        for (int p = 0; p < system.getNumberOfPhases(); p++) {
          System.out.printf("phase_%d_type = %s%n", p, system.getPhase(p).getType());
          row(String.format("phase_%d_beta", p), system.getBeta(p));
        }
        row("has_wax_phase", system.hasPhaseType("wax") ? 1.0 : 0.0);
        if (system.hasPhaseType("wax")) {
          row("wax_mass_fraction", system.getPhaseFraction("wax", "mass"));
          row("wax_mole_fraction", system.getPhaseFraction("wax", "mole"));
          // The mechanism: a mass fraction of ~0 for a phase with real moles is a degenerate
          // molar mass, and this is which side of the ratio is degenerate.
          int wax = system.getPhaseIndex("wax");
          row("wax_phase_molar_mass", system.getPhase(wax).getMolarMass());
          row("wax_phase_number_of_moles", system.getPhase(wax).getNumberOfMolesInPhase());
          row("wax_phase_beta", system.getBeta(wax));
          row("fluid_molar_mass", system.getMolarMass());
          row("fluid_total_moles", system.getTotalNumberOfMoles());
        }
      } catch (Exception failure) {
        System.out.printf("error = %s%n", failure.toString());
      }
    }
    System.out.println();
  }

  static void wax(String label, SystemInterface system, double pressureBara) {
    system.setPressure(pressureBara);
    System.out.printf("# %s%n", label);
    row("pressure_bara", pressureBara);
    try {
      new ThermodynamicOperations(system).calcWAT();
      row("temperature_K", system.getTemperature("K"));
      row("temperature_C", system.getTemperature("C"));
    } catch (Exception failure) {
      System.out.printf("error = %s%n", failure.toString());
    }
    System.out.println();
  }

  public static void main(String[] args) {
    neqsim.util.database.NeqSimDataBase.createTemporaryTables();

    prescribed("prescribed route (WaterDewPointTemperatureFlash), SRK, 35 ppm at 70 bara",
        35.0, 70.0);
    prescribed("prescribed route, SRK, 22 ppm at 70 bara", 22.0, 70.0);
    prescribed("prescribed route, SRK, 35 ppm at 20 bara", 35.0, 20.0);

    multiphase("multiphase route (WaterDewPointTemperatureMultiphaseFlash), SRK, 22 ppm at 70 bara",
        22.0, 70.0);
    multiphase("multiphase route, SRK, 35 ppm at 70 bara", 35.0, 70.0);
    multiphase("multiphase route, SRK, 200 ppm at 70 bara", 200.0, 70.0);

    line(35.0, 1.0, 200.0);

    waxSweep(5.0);

    wax("wax appearance temperature, the named feed methane/n-heptane/nc14/nc20 at 5 bara",
        waxFeed(), 5.0);
    wax("wax appearance temperature, a heavier named feed at 1 bara", richWaxFeed(), 1.0);
    wax("wax appearance temperature, a lean gas with one former", leanWaxFeed(), 50.0);
    wax("wax appearance temperature, WaxFlashTest's own fixture (plus fraction + addTBPWax)",
        waxFixture(), 50.0);
  }
}
