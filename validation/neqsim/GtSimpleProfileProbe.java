// The profile `GTSurfaceTensionSimple` walks, and the `dmu/dn` matrix its Newton solves
// against. LGT's capture gave the linear path; this one adds the inner solve, which is the
// only new machinery between the two ports.
//
//     javac -proc:none -cp neqsim-f0c7436.jar GtSimpleProfileProbe.java
//     java -cp .:neqsim-f0c7436.jar GtSimpleProfileProbe > captures/gt_simple_profile_probe.tsv
//
// The class's stored densities are mol/m³, its pressures Pa, and its `getDmudn2` holds the
// last step's matrix. The interaction rule in this class assigns 0.0 in both arms of its
// water condition, so every pair here is zero; the probe records the class's own answer
// for a water pair so that dead branch is visible rather than assumed.

import neqsim.physicalproperties.interfaceproperties.surfacetension.GTSurfaceTensionSimple;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class GtSimpleProfileProbe {

  public static void main(String[] args) {
    row("methane_nbutane_300K_20bara", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0);
    row("co2_water_313K_50bara", new String[] { "CO2", "water" }, new double[] { 0.02, 0.98 },
        313.15, 50.0);
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
    if (system.getNumberOfPhases() < 2) {
      System.out.println("no_interface=true");
      System.out.println();
      return;
    }

    GTSurfaceTensionSimple gt = new GTSurfaceTensionSimple(system);
    double sigma = gt.calcSurfaceTension(0, 1);
    System.out.println("sigma_N_per_m=" + sigma);
    int count = system.getPhase(0).getNumberOfComponents();
    for (int step = 0; step < 200; step++) {
      StringBuilder densities = new StringBuilder("step" + step + "_molar_density=");
      for (int c = 0; c < count; c++) {
        densities.append(gt.getMolarDensity(c)[step]).append(" ");
      }
      System.out.println(densities.toString().trim());
      System.out.println("step" + step + "_pressure_Pa=" + gt.getPressure()[step]);
      System.out.println("step" + step + "_z=" + gt.getz()[step]);
    }
    double[][][] dmu = gt.getDmudn2();
    if (dmu != null) {
      for (int i = 0; i < count; i++) {
        for (int k = 0; k < count; k++) {
          System.out.println("last_dmudn_" + i + k + "=" + dmu[199][i][k]);
        }
      }
    }
    System.out.println();
  }
}
