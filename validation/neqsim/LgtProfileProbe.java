// The gradient-theory profile `LGTSurfaceTension` integrates, and the influence parameters
// it integrates. The final sigma is in `surface_tension_solvers_probe.tsv`; this probe
// records the profile itself, because a port that only matches the last number can hide a
// compensating pair of unit errors.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar LgtProfileProbe.java
//     java -cp .:neqsim-f0c7436.jar LgtProfileProbe > captures/lgt_profile_probe.tsv
//
// The densities are the class's own, in the unit its `1e5` factor produces; the influence
// parameter's unit is NeqSim's, and the two are what the port has to reproduce.

import neqsim.physicalproperties.interfaceproperties.surfacetension.LGTSurfaceTension;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class LgtProfileProbe {

  public static void main(String[] args) {
    row("methane_nbutane_300K_20bara", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0);
    row("methane_150K_10bara", new String[] { "methane" }, new double[] { 1.0 }, 150.0, 10.0);
    row("co2_water_313K_50bara", new String[] { "CO2", "water" }, new double[] { 0.02, 0.98 },
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
    System.out.println("phase_count=" + system.getNumberOfPhases());
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      System.out.println("phase" + p + "_type=" + system.getPhase(p).getType());
      System.out.println(
          "phase" + p + "_density_kg_per_m3=" + system.getPhase(p).getPhysicalProperties().getDensity());
      System.out.println("phase" + p + "_molar_volume_m3_per_mol="
          + system.getPhase(p).getMolarVolume());
      System.out.println("phase" + p + "_molar_mass_kg_per_mol=" + system.getPhase(p).getMolarMass());
      for (int c = 0; c < system.getPhase(p).getNumberOfComponents(); c++) {
        System.out.println("phase" + p + "_component_" + c + "_x="
            + system.getPhase(p).getComponent(c).getx());
        System.out.println("phase" + p + "_component_" + c + "_nbp_K="
            + system.getPhase(p).getComponent(c).getNormalBoilingPoint());
        System.out.println("phase" + p + "_component_" + c
            + "_influence_parameter=" + system.getPhase(p).getComponent(c)
                .getSurfaceTenisionInfluenceParameter(system.getPhase(p).getTemperature()));
        neqsim.thermo.component.ComponentEosInterface componentEos =
            (neqsim.thermo.component.ComponentEosInterface) system.getPhase(p).getComponent(c);
        System.out.println("phase" + p + "_component_" + c + "_aT=" + componentEos.getaT());
        System.out.println("phase" + p + "_component_" + c + "_b=" + componentEos.getb());
        System.out.println("phase" + p + "_component_" + c + "_tc_K="
            + system.getPhase(p).getComponent(c).getTC());
        System.out.println("phase" + p + "_component_" + c + "_pc_bar="
            + system.getPhase(p).getComponent(c).getPC());
        System.out.println("phase" + p + "_component_" + c + "_omega="
            + system.getPhase(p).getComponent(c).getAcentricFactor());
      }
    }
    if (system.getNumberOfPhases() < 2) {
      System.out.println("no_interface=true");
      System.out.println();
      return;
    }

    LGTSurfaceTension lgt = new LGTSurfaceTension(system);
    double sigma = lgt.calcSurfaceTension(0, 1);
    System.out.println("sigma_N_per_m=" + sigma);
    int count = system.getPhase(0).getNumberOfComponents();
    for (int step = 0; step < 200; step++) {
      StringBuilder densities = new StringBuilder("step" + step + "_molar_density=");
      for (int c = 0; c < count; c++) {
        densities.append(lgt.getMolarDensity(c)[step]).append(" ");
      }
      System.out.println(densities.toString().trim());
      System.out.println("step" + step + "_pressure_Pa=" + lgt.getPressure()[step]);
      System.out.println("step" + step + "_z=" + lgt.getz()[step]);
    }
    System.out.println();
  }
}
