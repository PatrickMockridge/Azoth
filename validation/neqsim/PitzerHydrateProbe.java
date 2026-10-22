// PitzerHydrateProbe - the Pitzer-brine hydrate saturation at the pinned NeqSim commit.
//
// The class is `PitzerHydrateFlash`, reached both directly and through the public
// `ThermodynamicOperations` entry points. Each state prints the solved temperature or
// pressure, the stable structure, the residual, the water activity and the convergence flag.
// The fluid is built the way `PitzerHydrateFlashTest` builds it: a one-kilogram water basis,
// excess guest, ions in mol per kilogram, the classic mixing rule, multiphase on, and for a
// CO2 + chloride brine the PHREEQC CO2-chloride parameters with the test's zero zeta map.
//
//   javac -proc:none -cp neqsim-f0c7436.jar PitzerHydrateProbe.java
//   java -cp .:neqsim-f0c7436.jar PitzerHydrateProbe > captures/pitzer_hydrate_probe.tsv

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Supplier;
import neqsim.thermo.phase.PhaseHydrate;
import neqsim.thermo.system.SystemPitzer;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.flashops.saturationops.PitzerHydrateFlash;

public class PitzerHydrateProbe {

  private static final double WATER_MOLES = 1.0 / 0.01801528;
  private static final String[] CATIONS = {"Na+", "K+", "Ca++", "Mg++"};
  private static final double[] SALT_MOLAR_MASS = {0.05844, 0.0745513, 0.11098, 0.095211};

  static Map<String, Double> zeroZeta() {
    Map<String, Double> values = new LinkedHashMap<String, Double>();
    for (String cation : CATIONS) {
      values.put(cation, 0.0);
    }
    return values;
  }

  /**
   * A charge-balanced brine on a one-kilogram water basis: `saltWtPct[i]` is the mass of cation
   * `i`'s chloride salt per hundred grams of solution.
   */
  static SystemPitzer brine(String[] guests, double[] guestMoles, double[] saltWtPct, double pressure) {
    SystemPitzer fluid = new SystemPitzer(280.0, pressure);
    for (int i = 0; i < guests.length; i++) {
      fluid.addComponent(guests[i], guestMoles[i]);
    }
    fluid.addComponent("water", WATER_MOLES);
    double solute = 0.0;
    for (double value : saltWtPct) {
      solute += value;
    }
    double waterPart = 100.0 - solute;
    double chloride = 0.0;
    for (int i = 0; i < CATIONS.length; i++) {
      double moles = saltWtPct[i] / waterPart / SALT_MOLAR_MASS[i];
      fluid.addComponent(CATIONS[i], moles);
      chloride += moles * (i < 2 ? 1.0 : 2.0);
    }
    fluid.addComponent("Cl-", chloride);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    boolean hasCo2 = false;
    for (String guest : guests) {
      hasCo2 |= "CO2".equals(guest);
    }
    if (hasCo2 && solute > 0.0) {
      fluid.applyPhreeqcCo2ChlorideParameters(zeroZeta(), 0.0);
    }
    return fluid;
  }

  static SystemPitzer waterAnd(String guest, double guestMoles, double pressure) {
    return brine(new String[] {guest}, new double[] {guestMoles}, new double[4], pressure);
  }

  static SystemPitzer waterOnly(double pressure) {
    SystemPitzer fluid = new SystemPitzer(280.0, pressure);
    fluid.addComponent("water", WATER_MOLES);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    return fluid;
  }

  static SystemPitzer guestOnly(String guest, double pressure) {
    SystemPitzer fluid = new SystemPitzer(280.0, pressure);
    fluid.addComponent(guest, 10.0);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    return fluid;
  }

  static SystemPitzer singleSalt(String guest, int salt, double weightPercent, double pressure) {
    double[] parts = new double[4];
    parts[salt] = weightPercent;
    return brine(new String[] {guest}, new double[] {10.0}, parts, pressure);
  }

  static String structure(SystemPitzer fluid) {
    return String.valueOf(((PhaseHydrate) fluid.getPhases()[4]).getStableHydrateStructure());
  }

  static void temperatureRow(String label, SystemPitzer fluid, int requested) {
    PitzerHydrateFlash flash = new PitzerHydrateFlash(fluid, false);
    if (requested != 0) {
      flash.setStructure(requested);
    }
    try {
      flash.run();
    } catch (RuntimeException failure) {
      System.out.println("  " + label + " status=refused reason=" + failure.getClass().getSimpleName()
          + ":" + failure.getMessage());
      return;
    }
    System.out.println("  " + label + " status=solved temperature_k=" + fluid.getTemperature()
        + " structure=" + structure(fluid) + " residual=" + flash.getResidual()
        + " water_activity=" + flash.getWaterActivity() + " converged=" + flash.isConverged());
  }

  static void pressureRow(String label, SystemPitzer fluid, int requested) {
    PitzerHydrateFlash flash = new PitzerHydrateFlash(fluid, true);
    if (requested != 0) {
      flash.setStructure(requested);
    }
    try {
      flash.run();
    } catch (RuntimeException failure) {
      System.out.println("  " + label + " status=refused reason=" + failure.getClass().getSimpleName()
          + ":" + failure.getMessage());
      return;
    }
    System.out.println("  " + label + " status=solved pressure_bara=" + fluid.getPressure()
        + " structure=" + structure(fluid) + " residual=" + flash.getResidual()
        + " water_activity=" + flash.getWaterActivity() + " converged=" + flash.isConverged());
  }

  static void refusalRow(String label, Supplier<SystemPitzer> build) {
    try {
      SystemPitzer fluid = build.get();
      PitzerHydrateFlash flash = new PitzerHydrateFlash(fluid, false);
      flash.run();
      System.out.println("  " + label + " status=solved temperature_k=" + fluid.getTemperature());
    } catch (RuntimeException failure) {
      System.out.println("  " + label + " status=refused reason=" + failure.getClass().getSimpleName()
          + ":" + failure.getMessage());
    }
  }

  public static void main(String[] args) throws Exception {
    System.out.println("pitzer_hydrate_flash");

    System.out.println(" temperature_solve");
    temperatureRow("co2_water_30bara", waterAnd("CO2", 10.0, 30.0), 0);
    temperatureRow("co2_water_60bara", waterAnd("CO2", 10.0, 60.0), 0);
    temperatureRow("co2_water_100bara", waterAnd("CO2", 10.0, 100.0), 0);
    temperatureRow("co2_nacl_1molal_30bara", singleSalt("CO2", 0, 5.5213, 30.0), 0);
    temperatureRow("co2_nacl_5wt_30bara", singleSalt("CO2", 0, 5.0, 30.0), 0);
    temperatureRow("co2_nacl_10wt_30bara", singleSalt("CO2", 0, 10.0, 30.0), 0);
    temperatureRow("co2_nacl_20wt_30bara", singleSalt("CO2", 0, 20.0, 30.0), 0);
    temperatureRow("co2_mgcl2_5wt_30bara", singleSalt("CO2", 3, 5.0, 30.0), 0);
    temperatureRow("co2_mgcl2_15wt_30bara", singleSalt("CO2", 3, 15.0, 30.0), 0);
    temperatureRow("methane_water_50bara", waterAnd("methane", 10.0, 50.0), 0);
    temperatureRow("methane_water_100bara", waterAnd("methane", 10.0, 100.0), 0);
    temperatureRow("co2_methane_water_50bara",
        brine(new String[] {"CO2", "methane"}, new double[] {5.0, 5.0}, new double[4], 50.0), 0);

    System.out.println(" pressure_solve");
    for (double temperature : new double[] {274.5, 278.0, 280.0, 283.0}) {
      SystemPitzer fluid = waterAnd("CO2", 10.0, 30.0);
      fluid.setTemperature(temperature);
      pressureRow("co2_water_" + temperature + "k", fluid, 0);
    }
    SystemPitzer salty = singleSalt("CO2", 0, 5.0, 30.0);
    salty.setTemperature(278.0);
    pressureRow("co2_nacl_5wt_278k", salty, 0);

    System.out.println(" structure_request");
    temperatureRow("co2_water_30bara_auto", waterAnd("CO2", 10.0, 30.0), 0);
    temperatureRow("co2_water_30bara_si", waterAnd("CO2", 10.0, 30.0), 1);
    temperatureRow("co2_water_30bara_sii", waterAnd("CO2", 10.0, 30.0), 2);
    temperatureRow("co2_water_60bara_auto", waterAnd("CO2", 10.0, 60.0), 0);
    temperatureRow("co2_water_60bara_si", waterAnd("CO2", 10.0, 60.0), 1);
    temperatureRow("co2_water_60bara_sii", waterAnd("CO2", 10.0, 60.0), 2);

    System.out.println(" public_operations_api");
    SystemPitzer api = singleSalt("CO2", 0, 5.0, 30.0);
    ThermodynamicOperations operations = new ThermodynamicOperations(api);
    try {
      operations.hydrateFormationTemperature();
      PitzerHydrateFlash flash = (PitzerHydrateFlash) operations.getOperation();
      System.out.println("  hydrate_formation_temperature status=solved temperature_k=" + api.getTemperature()
          + " structure=" + structure(api) + " residual=" + flash.getResidual()
          + " water_activity=" + flash.getWaterActivity());
      double solvedTemperature = api.getTemperature();
      api.setPressure(20.0);
      operations.hydrateFormationPressure();
      System.out.println("  hydrate_formation_pressure status=solved pressure_bara=" + api.getPressure()
          + " temperature_k=" + api.getTemperature()
          + " temperature_change_k=" + (api.getTemperature() - solvedTemperature));
    } catch (RuntimeException failure) {
      System.out.println("  public_operations_api status=refused reason=" + failure.getClass().getSimpleName()
          + ":" + failure.getMessage());
    }

    System.out.println(" refusals");
    refusalRow("water_only_30bara", () -> waterOnly(30.0));
    refusalRow("co2_only_30bara", () -> guestOnly("CO2", 30.0));
  }
}
