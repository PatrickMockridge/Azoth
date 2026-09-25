// The transport properties a *phase* answers by default, and the models behind them.
//
// Compile and run from this directory:
//
//     javac -cp .:neqsim-f0c7436.jar -d . PhaseTransportProbe.java
//     java -cp .:neqsim-f0c7436.jar PhaseTransportProbe > captures/phase_transport_probe.tsv
//
// **The defaults are not the models the class's own names suggest.** `GasPhysicalProperties`
// comments the Chung methods out and assigns `PFCTViscosityMethodHeavyOil` and
// `PFCTConductivityMethodMod86` for a gas as well as for a liquid, and its diffusivity is the
// base `Diffusivity` class - Chapman-Enskog on Lennard-Jones parameters - not
// Fuller-Schettler-Giddings, which is only reached by `setDiffusionCoefficientModel`. The
// rate-based column sets no model at all, so these defaults are what its segment fluxes use.
//
// What the probe prints is the whole input set of each correlation: the phase's T and P, its
// components' molar masses and Lennard-Jones parameters, and - per pair - the binary matrix the
// diffusivity answers. A port that reads any of those differently lands somewhere else here.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class PhaseTransportProbe {

  public static void main(String[] args) {
    row("methane_nitrogen_298K_1atm", new String[] { "methane", "nitrogen" },
        new double[] { 0.5, 0.5 }, 298.15, 1.01325);
    row("methane_co2_313K_50bar", new String[] { "methane", "CO2" }, new double[] { 0.5, 0.5 },
        313.15, 50.0);
    row("methane_nbutane_300K_20bar", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0);
    row("co2_water_313K_50bar", new String[] { "CO2", "water" }, new double[] { 0.02, 0.98 },
        313.15, 50.0);
  }

  static void row(String label, String[] names, double[] z, double temperatureK,
      double pressureBara) {
    SystemInterface system = new SystemSrkEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.createDatabase(true);
    system.setMixingRule(2);
    new ThermodynamicOperations(system).TPflash();
    system.initPhysicalProperties();

    System.out.println(label);
    System.out.println("components=" + String.join(" ", names));
    System.out.println("phase_count=" + system.getNumberOfPhases());
    StringBuilder molarMasses = new StringBuilder("molar_mass_kg_per_mol=");
    StringBuilder diameters = new StringBuilder("lennard_jones_diameter_angstrom=");
    StringBuilder energies = new StringBuilder("lennard_jones_energy_k=");
    for (int c = 0; c < names.length; c++) {
      molarMasses.append(system.getPhase(0).getComponent(c).getMolarMass()).append(" ");
      diameters.append(system.getPhase(0).getComponent(c).getLennardJonesMolecularDiameter())
          .append(" ");
      energies.append(system.getPhase(0).getComponent(c).getLennardJonesEnergyParameter())
          .append(" ");
    }
    System.out.println(molarMasses.toString().trim());
    System.out.println(diameters.toString().trim());
    System.out.println(energies.toString().trim());

    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      System.out.println("phase" + p + "_type=" + system.getPhase(p).getType());
      System.out.println("phase" + p + "_temperature_K=" + system.getPhase(p).getTemperature());
      System.out.println("phase" + p + "_pressure_bara=" + system.getPhase(p).getPressure());
      System.out.println(
          "phase" + p + "_density_kg_per_m3=" + system.getPhase(p).getPhysicalProperties().getDensity());
      System.out.println("phase" + p + "_molar_mass_kg_per_mol=" + system.getPhase(p).getMolarMass());
      System.out.println("phase" + p + "_viscosity_Pa_s="
          + system.getPhase(p).getPhysicalProperties().getViscosity());
      System.out.println("phase" + p + "_conductivity_W_per_mK="
          + system.getPhase(p).getPhysicalProperties().getConductivity());
      StringBuilder composition = new StringBuilder("phase" + p + "_x=");
      for (int c = 0; c < system.getPhase(p).getNumberOfComponents(); c++) {
        composition.append(names[c]).append(":").append(system.getPhase(p).getComponent(c).getx())
            .append(" ");
      }
      System.out.println(composition.toString().trim());
      // The matrix and the effective vector are filled by `calcDiffusionCoefficients`, which
      // reads the model's own binary coefficients first.
      system.getPhase(p).getPhysicalProperties().diffusivityCalc.calcDiffusionCoefficients(0, 1);
      // The default binary diffusivity matrix, at the phase's own state.
      for (int i = 0; i < names.length; i++) {
        for (int j = 0; j < names.length; j++) {
          if (i == j) {
            continue;
          }
          System.out.println("phase" + p + "_d_" + i + j + "_m2_per_s=" + system.getPhase(p)
              .getPhysicalProperties().diffusivityCalc.calcBinaryDiffusionCoefficient(i, j, 0));
        }
      }
      for (int i = 0; i < names.length; i++) {
        System.out.println("phase" + p + "_d_effective_" + i + "_m2_per_s="
            + system.getPhase(p).getPhysicalProperties().getEffectiveDiffusionCoefficient(i));
      }
    }
    // **The selected `Chapman-Enskog` model is not the default one.** The string sets
    // `useDiffusionLJOverride(true)`, which swaps the database's Lennard-Jones parameters for
    // Poling/BSL's textbook values - so the same arithmetic answers two different numbers, and
    // the rate-based column, which selects nothing, takes the first.
    if (system.getNumberOfPhases() > 0) {
      system.getPhase(0).getPhysicalProperties().setDiffusionCoefficientModel("Chapman-Enskog");
      System.out.println("selected_model_d_01_m2_per_s=" + system.getPhase(0).getPhysicalProperties()
          .diffusivityCalc.calcBinaryDiffusionCoefficient(0, 1, 0));
    }
    System.out.println();
  }
}
