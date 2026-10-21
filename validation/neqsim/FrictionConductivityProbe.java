// `FrictionTheoryConductivityMethod`, reached the only way the tree reaches it - an
// explicit `setConductivityModel("friction theory")` on a phase. The class reads the
// phase's *existing* viscosity, so that is printed beside the answer, with each
// component's Cv0 and the class's own pure-component arrays.
//
//     javac -proc:none -cp neqsim-f0c7436.jar FrictionConductivityProbe.java
//     java -cp .:neqsim-f0c7436.jar FrictionConductivityProbe > captures/friction_conductivity_probe.tsv

import java.lang.reflect.Field;
import neqsim.physicalproperties.methods.commonphasephysicalproperties.conductivity.FrictionTheoryConductivityMethod;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class FrictionConductivityProbe {

  public static void main(String[] args) {
    row("srk_methane_nbutane_300K_20bara", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, false);
    row("pr_methane_nbutane_300K_20bara", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, true);
    row("pr_methane_300K_10bara", new String[] { "methane" }, new double[] { 1.0 }, 300.0, 10.0,
        true);
    row("aqueous_co2_water_313K_50bara", new String[] { "CO2", "water" },
        new double[] { 0.02, 0.98 }, 313.15, 50.0, false);
    row("npentane_300K_1bara", new String[] { "n-pentane" }, new double[] { 1.0 }, 300.0, 1.0,
        false);
  }

  static void row(String label, String[] names, double[] z, double temperatureK,
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
    System.out.println("phase_index=" + index);
    System.out.println("phase_type=" + system.getPhase(index).getType());
    System.out.println("temperature_K=" + system.getPhase(index).getTemperature());
    System.out.println("molar_mass_kg_per_mol=" + system.getPhase(index).getMolarMass());

    // The class reads `phase.getViscosity()`, so the phase's own model answers that.
    double viscosity = system.getPhase(index).getPhysicalProperties().getViscosity();
    System.out.println("viscosity_Pa_s=" + viscosity);
    system.getPhase(index).getPhysicalProperties().setConductivityModel("friction theory");
    system.initPhysicalProperties();
    double conductivity = system.getPhase(index).getPhysicalProperties().getConductivity();
    Object calc = system.getPhase(index).getPhysicalProperties().conductivityCalc;
    System.out.println("conductivity_model_class=" + calc.getClass().getSimpleName());
    System.out.println("conductivity_W_per_m_K=" + conductivity);
    System.out.println("viscosity_after_Pa_s="
        + system.getPhase(index).getPhysicalProperties().getViscosity());

    for (int c = 0; c < system.getPhase(index).getNumberOfComponents(); c++) {
      System.out.println("component_" + c + "_x=" + system.getPhase(index).getComponent(c).getx());
      System.out.println("component_" + c + "_tc_K="
          + system.getPhase(index).getComponent(c).getTC());
      System.out.println("component_" + c + "_critical_volume_cm3_per_mol="
          + system.getPhase(index).getComponent(c).getCriticalVolume());
      System.out.println("component_" + c + "_omega="
          + system.getPhase(index).getComponent(c).getAcentricFactor());
      System.out.println("component_" + c + "_dipole_debye="
          + system.getPhase(index).getComponent(c).getDebyeDipoleMoment());
      System.out.println("component_" + c + "_molar_mass_kg_per_mol="
          + system.getPhase(index).getComponent(c).getMolarMass());
      System.out.println("component_" + c + "_cv0_J_per_mol_K="
          + system.getPhase(index).getComponent(c).getCv0(system.getPhase(index).getTemperature()));
    }
    if (calc instanceof FrictionTheoryConductivityMethod) {
      FrictionTheoryConductivityMethod method = (FrictionTheoryConductivityMethod) calc;
      for (int c = 0; c < system.getPhase(index).getNumberOfComponents(); c++) {
        System.out.println("component_" + c + "_pure_conductivity_W_per_m_K="
            + method.pureComponentConductivity[c]);
        System.out.println("component_" + c + "_fc=" + method.fc[c]);
        System.out.println("component_" + c + "_omegaCond=" + method.omegaCond[c]);
        try {
          Field field =
              FrictionTheoryConductivityMethod.class.getDeclaredField("pureComponentViscosity");
          field.setAccessible(true);
          double[] pureViscosity = (double[]) field.get(method);
          System.out.println(
              "component_" + c + "_pure_viscosity_Pa_s=" + pureViscosity[c]);
        } catch (ReflectiveOperationException ex) {
          System.out.println("pure_viscosity_error=" + ex.getClass().getSimpleName());
        }
      }
    }
    System.out.println();
  }

  /** NeqSim orders phases gas-first; the densest one is the liquid this probe wants. */
  static int densestPhase(SystemInterface system) {
    int index = 0;
    double densest = -1.0;
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      double density = system.getPhase(p).getPhysicalProperties().getDensity();
      if (density > densest) {
        densest = density;
        index = p;
      }
    }
    return index;
  }
}
