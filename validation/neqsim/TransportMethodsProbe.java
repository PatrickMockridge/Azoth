// The two selectable transport methods P1 still owes: Filippov conductivity and friction-theory
// viscosity. Both are direct calculations over a phase; neither needs a flash.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar TransportMethodsProbe.java
//     java -cp .:neqsim-f0c7436.jar TransportMethodsProbe > captures/transport_methods_probe.tsv
//
// The model strings are `PhysicalProperties.setConductivityModel("Filippov")` and
// `setViscosityModel("friction theory")` (`PhysicalProperties.java:185-240`). Both are per phase,
// so the probe picks the densest phase of each state and records which one it drove.
//
// The units are the reason this capture prints every intermediate. `FrictionTheoryViscosityMethod`
// takes `Pc` in **bar**, `Vc` in cm3/mol, `criticalViscosity` in Pa*s scaled by 1e7 to micropoise,
// and returns Pa*s after a final 1e-7. `Filippov`'s Missenard correction uses `Tc`/`Pc` in the
// component's own units (K, bar) against a pressure in bara, and its javadoc says `P - Psat`
// while the code uses `P`.

import neqsim.physicalproperties.methods.commonphasephysicalproperties.viscosity.FrictionTheoryViscosityMethod;
import neqsim.physicalproperties.methods.liquidphysicalproperties.conductivity.FilippovConductivityMethod;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class TransportMethodsProbe {

  public static void main(String[] args) {
    frictionRow("srk_methane_nbutane_300K_20bara", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, false);
    frictionRow("pr_methane_nbutane_300K_20bara", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, true);
    frictionRow("pr_methane_300K_10bara", new String[] { "methane" }, new double[] { 1.0 }, 300.0,
        10.0, true);

    filippovRow("aqueous_co2_water_313K_50bara", new String[] { "CO2", "water" },
        new double[] { 0.02, 0.98 }, 313.15, 50.0, false);
    filippovRow("oil_methane_nbutane_300K_20bara", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, false);
    filippovRow("water_300K_1000bara", new String[] { "water" }, new double[] { 1.0 }, 300.0, 1000.0,
        false);
    filippovRow("c6_c10_313K_5bara", new String[] { "n-hexane", "n-decane" },
        new double[] { 0.5, 0.5 }, 313.15, 5.0, false);
  }

  static void frictionRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, boolean pengRobinson) {
    SystemInterface system = pengRobinson ? new SystemPrEos(temperatureK, pressureBara)
        : new SystemSrkEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.createDatabase(true);
    system.setMixingRule(2);
    new ThermodynamicOperations(system).TPflash();
    system.initPhysicalProperties();
    int index = densestPhase(system);

    System.out.println(label);
    System.out.println("eos=" + system.getClass().getSimpleName());
    System.out.println("phase_count=" + system.getNumberOfPhases());
    printState(system, names);
    system.getPhase(index).getPhysicalProperties().setViscosityModel("friction theory");
    system.initPhysicalProperties();

    double mu = system.getPhase(index).getPhysicalProperties().getViscosity();
    Object calc = system.getPhase(index).getPhysicalProperties().viscosityCalc;
    System.out.println("viscosity_model_class=" + calc.getClass().getSimpleName());
    System.out.println("viscosity_Pa_s=" + mu);
    System.out.println("phase_index=" + index);
    System.out.println("phase_type=" + system.getPhase(index).getType());
    System.out.println("phase_density_kg_per_m3="
        + system.getPhase(index).getPhysicalProperties().getDensity());
    System.out.println("phase_molar_volume_m3_per_mol=" + system.getPhase(index).getMolarVolume());
    try {
      PhaseEosInterface eos = (PhaseEosInterface) system.getPhase(index);
      System.out.println("pressure_repulsive_Pa=" + eos.getPressureRepulsive());
      System.out.println("pressure_attractive_Pa=" + eos.getPressureAttractive());
    } catch (RuntimeException ex) {
      System.out.println("pressure_split_error=" + ex.getClass().getSimpleName());
    }
    if (calc instanceof FrictionTheoryViscosityMethod) {
      FrictionTheoryViscosityMethod friction = (FrictionTheoryViscosityMethod) calc;
      for (int c = 0; c < system.getPhase(index).getNumberOfComponents(); c++) {
        System.out.println("component_" + c + "_x=" + system.getPhase(index).getComponent(c).getx());
        System.out.println("component_" + c + "_tc_K="
            + system.getPhase(index).getComponent(c).getTC());
        System.out.println("component_" + c + "_pc_bar="
            + system.getPhase(index).getComponent(c).getPC());
        System.out.println("component_" + c + "_critical_volume_cm3_per_mol="
            + system.getPhase(index).getComponent(c).getCriticalVolume());
        System.out.println("component_" + c + "_omega="
            + system.getPhase(index).getComponent(c).getAcentricFactor());
        System.out.println("component_" + c + "_molar_mass_kg_per_mol="
            + system.getPhase(index).getComponent(c).getMolarMass());
        System.out.println("component_" + c + "_critical_viscosity_Pa_s="
            + system.getPhase(index).getComponent(c).getCriticalViscosity());
        System.out.println("component_" + c + "_viscosity_correction_factor="
            + system.getPhase(index).getComponent(c).getViscosityCorrectionFactor());
        System.out.println("component_" + c + "_dipole_debye="
            + system.getPhase(index).getComponent(c).getDebyeDipoleMoment());
        System.out.println("component_" + c + "_friction_k="
            + system.getPhase(index).getComponent(c).getViscosityFrictionK());
        System.out.println("component_" + c + "_pure_viscosity_micropoise="
            + friction.pureComponentViscosity[c]);
        System.out.println("component_" + c + "_fc=" + friction.Fc[c]);
        System.out.println("component_" + c + "_omega_visc=" + friction.omegaVisc[c]);
      }
    }
    System.out.println();
  }

  static void filippovRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, boolean pengRobinson) {
    SystemInterface system = pengRobinson ? new SystemPrEos(temperatureK, pressureBara)
        : new SystemSrkEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.createDatabase(true);
    system.setMixingRule(2);
    new ThermodynamicOperations(system).TPflash();
    system.initPhysicalProperties();
    int index = densestPhase(system);

    System.out.println(label);
    System.out.println("eos=" + system.getClass().getSimpleName());
    System.out.println("phase_count=" + system.getNumberOfPhases());
    printState(system, names);
    system.getPhase(index).getPhysicalProperties().setConductivityModel("Filippov");
    system.initPhysicalProperties();

    Object calc = system.getPhase(index).getPhysicalProperties().conductivityCalc;
    System.out.println("conductivity_model_class=" + calc.getClass().getSimpleName());
    System.out.println("phase_index=" + index);
    System.out.println("phase_type=" + system.getPhase(index).getType());
    System.out.println("phase_density_kg_per_m3="
        + system.getPhase(index).getPhysicalProperties().getDensity());
    System.out.println("temperature_K=" + system.getPhase(index).getTemperature());
    System.out.println("pressure_bar=" + system.getPhase(index).getPressure());
    double tcMix = 0.0;
    double pcMix = 0.0;
    double wtSum = 0.0;
    for (int c = 0; c < system.getPhase(index).getNumberOfComponents(); c++) {
      tcMix += system.getPhase(index).getComponent(c).getx()
          * system.getPhase(index).getComponent(c).getTC();
      pcMix += system.getPhase(index).getComponent(c).getx()
          * system.getPhase(index).getComponent(c).getPC();
      wtSum += system.getPhase(index).getWtFrac(c);
    }
    System.out.println("mixture_tc_K=" + tcMix);
    System.out.println("mixture_pc_bar=" + pcMix);
    System.out.println("tr=" + system.getPhase(index).getTemperature() / tcMix);
    System.out.println("pr=" + system.getPhase(index).getPressure() / pcMix);
    System.out.println("weight_fraction_sum=" + wtSum);
    for (int c = 0; c < system.getPhase(index).getNumberOfComponents(); c++) {
      System.out.println("component_" + c + "_x=" + system.getPhase(index).getComponent(c).getx());
      System.out.println("component_" + c + "_w=" + system.getPhase(index).getWtFrac(c));
      System.out.println("component_" + c + "_tc_K="
          + system.getPhase(index).getComponent(c).getTC());
      System.out.println("component_" + c + "_pc_bar="
          + system.getPhase(index).getComponent(c).getPC());
      System.out.println("component_" + c + "_liquid_conductivity_0="
          + system.getPhase(index).getComponent(c).getLiquidConductivityParameter(0));
      System.out.println("component_" + c + "_liquid_conductivity_1="
          + system.getPhase(index).getComponent(c).getLiquidConductivityParameter(1));
      System.out.println("component_" + c + "_liquid_conductivity_2="
          + system.getPhase(index).getComponent(c).getLiquidConductivityParameter(2));
      if (calc instanceof FilippovConductivityMethod) {
        System.out.println("component_" + c + "_pure_conductivity_W_per_mK="
            + ((FilippovConductivityMethod) calc).pureComponentConductivity[c]);
      }
    }
    double withCorrection = system.getPhase(index).getPhysicalProperties().getConductivity();
    System.out.println("conductivity_with_pressure_correction_W_per_mK=" + withCorrection);
    if (calc instanceof FilippovConductivityMethod) {
      FilippovConductivityMethod filippov = (FilippovConductivityMethod) calc;
      System.out.println("filippov_coefficient=" + filippov.getFilippovCoefficient());
      System.out.println("use_pressure_correction=" + filippov.isUsePressureCorrection());
      filippov.setUsePressureCorrection(false);
      system.getPhase(index).getPhysicalProperties().init(system.getPhase(index));
      System.out.println("conductivity_without_pressure_correction_W_per_mK="
          + system.getPhase(index).getPhysicalProperties().getConductivity());
      filippov.setUsePressureCorrection(true);
    }
    System.out.println();
  }

  static int densestPhase(SystemInterface system) {
    int index = 0;
    double density = -1.0;
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      double candidate = system.getPhase(p).getPhysicalProperties().getDensity();
      if (candidate > density) {
        density = candidate;
        index = p;
      }
    }
    return index;
  }

  static void printState(SystemInterface system, String[] names) {
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      System.out.println("phase" + p + "_type=" + system.getPhase(p).getType());
      System.out.println(
          "phase" + p + "_density_kg_per_m3=" + system.getPhase(p).getPhysicalProperties().getDensity());
      System.out.println("phase" + p + "_molar_mass_kg_per_mol=" + system.getPhase(p).getMolarMass());
      StringBuilder composition = new StringBuilder("phase" + p + "_x=");
      for (int c = 0; c < system.getPhase(p).getNumberOfComponents(); c++) {
        composition.append(names[c]).append(":").append(system.getPhase(p).getComponent(c).getx())
            .append(" ");
      }
      System.out.println(composition.toString().trim());
    }
  }
}
