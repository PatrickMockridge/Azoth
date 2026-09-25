// The gas-phase binary diffusivity of `FullerSchettlerGiddingsDiffusivity`, on the states
// `DiffusivityExperimentalValidationTest` validates against published data.
//
// Compile and run from this directory:
//
//     javac -cp .:neqsim-f0c7436.jar -d . FullerProbe.java
//     java -cp .:neqsim-f0c7436.jar FullerProbe > captures/fuller_probe.tsv
//
// **The two experimental rows are the oracle and the probe cannot state them itself**: the class
// is a correlation, so what a port can be held to is the number the class answers for a pair whose
// measured value is published - methane/nitrogen at 2.2e-5 m2/s and CO2/nitrogen at 1.67e-5, both
// Marrero & Mason (1972) at 298 K and one atmosphere, which the class's own test asserts to 25 per
// cent. The rows below print the class's answer beside the inputs a port needs, so the port is held
// to the class's arithmetic rather than to the band: the volume each component's diffusion volume
// is resolved to is the port's own ladder (the special-molecule table, then `0.285 * Vc`), and a
// port whose table disagrees answers a different number here.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class FullerProbe {

  public static void main(String[] args) {
    // The class's own experimental rows: Marrero & Mason (1972) at 298 K, one atmosphere.
    row("methane_nitrogen_298K_1atm", "methane", "nitrogen", 298.15, 1.01325);
    row("co2_nitrogen_298K_1atm", "CO2", "nitrogen", 298.15, 1.01325);
    // The two scaling rows its test asserts: D ~ 1/P and D ~ T^1.75.
    row("methane_nitrogen_298K_5bar", "methane", "nitrogen", 298.15, 5.0);
    row("methane_nitrogen_300K_1atm", "methane", "nitrogen", 300.0, 1.01325);
    row("methane_nitrogen_400K_1atm", "methane", "nitrogen", 400.0, 1.01325);
    // **A pair where one name is in the special table and the other is not**, so the row exercises
    // both rungs of the ladder at once: methane is listed, and `ammonia` is resolved from its
    // critical volume because the table's key for it is the formula `NH3` and the lookup is by
    // name. `n-dodecane` is the second such name.
    row("methane_ammonia_298K_1atm", "methane", "ammonia", 298.15, 1.01325);
    row("methane_n_dodecane_298K_1atm", "methane", "n-dodecane", 298.15, 1.01325);
    // The states the rate-based column runs at, so the model's own inputs are captured too.
    row("methane_co2_313K_50bar", "methane", "CO2", 313.15, 50.0);
  }

  static void row(String label, String first, String second, double temperatureK,
      double pressureBara) {
    try {
      rowOf(label, first, second, temperatureK, pressureBara);
    } catch (RuntimeException error) {
      // A name the database does not carry would otherwise take the whole probe down with it.
      System.out.println(label);
      System.out.println("error=" + error.getMessage());
      System.out.println();
    }
  }

  static void rowOf(String label, String first, String second, double temperatureK,
      double pressureBara) {
    SystemInterface system = new SystemSrkEos(temperatureK, pressureBara);
    system.addComponent(first, 0.5);
    system.addComponent(second, 0.5);
    system.createDatabase(true);
    system.setMixingRule(2);
    new ThermodynamicOperations(system).TPflash();
    system.initPhysicalProperties();

    if (!system.hasPhaseType("gas")) {
      System.out.println(label);
      System.out.println("has_gas_phase=false");
      System.out.println();
      return;
    }
    system.getPhase("gas").getPhysicalProperties()
        .setDiffusionCoefficientModel("Fuller-Schettler-Giddings");
    double d = system.getPhase("gas").getPhysicalProperties().diffusivityCalc
        .calcBinaryDiffusionCoefficient(0, 1, 0);

    System.out.println(label);
    System.out.println("cubic=srk");
    System.out.println("first=" + first);
    System.out.println("second=" + second);
    System.out.println("temperature_K=" + system.getPhase("gas").getTemperature());
    System.out.println("pressure_bara=" + system.getPhase("gas").getPressure());
    System.out.println("first_molar_mass_g_per_mol="
        + system.getPhase("gas").getComponent(0).getMolarMass() * 1000.0);
    System.out.println("second_molar_mass_g_per_mol="
        + system.getPhase("gas").getComponent(1).getMolarMass() * 1000.0);
    System.out.println("first_critical_volume_cm3_per_mol="
        + system.getPhase("gas").getComponent(0).getCriticalVolume());
    System.out.println("second_critical_volume_cm3_per_mol="
        + system.getPhase("gas").getComponent(1).getCriticalVolume());
    System.out.println("d_m2_per_s=" + d);
    System.out.println();
  }
}
