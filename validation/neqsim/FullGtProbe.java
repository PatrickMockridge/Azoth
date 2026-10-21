// `GTSurfaceTensionFullGT.runcase` on several states. The solver-probe capture has one
// FullGT number (model 1's gas-oil slot on methane/n-butane); this gives the port more than
// one state to be wrong about.
//
//     javac -proc:none -cp neqsim-f0c7436.jar FullGtProbe.java
//     java -cp .:neqsim-f0c7436.jar FullGtProbe > captures/full_gt_probe.tsv
//
// The class exposes no intermediate state - runcase is the whole surface - so this records
// the bulk phases, the influence parameters and the answer.

import neqsim.physicalproperties.interfaceproperties.surfacetension.GTSurfaceTensionFullGT;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class FullGtProbe {

  public static void main(String[] args) {
    row("methane_nbutane_300K_20bara", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0);
    row("methane_nheptane_313K_30bara", new String[] { "methane", "n-heptane" },
        new double[] { 0.5, 0.5 }, 313.15, 30.0);
    row("co2_water_313K_50bara", new String[] { "CO2", "water" }, new double[] { 0.02, 0.98 },
        313.15, 50.0);
    row("methane_150K_10bara", new String[] { "methane" }, new double[] { 1.0 }, 150.0, 10.0);
    row("nheptane_water_300K_1bara", new String[] { "n-heptane", "water" },
        new double[] { 0.5, 0.5 }, 300.0, 1.0);
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
      System.out.println("phase" + p + "_molar_volume=" + system.getPhase(p).getMolarVolume());
      for (int c = 0; c < system.getPhase(p).getNumberOfComponents(); c++) {
        System.out.println("phase" + p + "_component_" + c + "_x="
            + system.getPhase(p).getComponent(c).getx());
        System.out.println("phase" + p + "_component_" + c + "_influence_parameter="
            + system.getPhase(p).getComponent(c)
                .getSurfaceTenisionInfluenceParameter(system.getPhase(p).getTemperature()));
      }
    }
    if (system.getNumberOfPhases() >= 2) {
      try {
        GTSurfaceTensionFullGT full = new GTSurfaceTensionFullGT(system, 0, 1);
        System.out.println("sigma_N_per_m=" + full.runcase());
      } catch (RuntimeException ex) {
        System.out.println("error=" + ex.getClass().getSimpleName() + ":" + ex.getMessage());
      }
    } else {
      System.out.println("no_interface=true");
    }
    System.out.println();
  }
}
