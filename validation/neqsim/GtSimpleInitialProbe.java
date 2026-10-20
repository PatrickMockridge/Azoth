// The local system `GTSurfaceTensionSimple` starts from: the chemical potentials and the
// `dmu/dn` matrix at the first phase's density, before its initial gradient solve. The
// profile capture gives the step that solve produces; this gives the numbers it solves on,
// so a port can see which side of the equation differs.
//
//     javac -proc:none -cp neqsim-f0c7436.jar GtSimpleInitialProbe.java
//     java -cp .:neqsim-f0c7436.jar GtSimpleInitialProbe > captures/gt_simple_initial_probe.tsv

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class GtSimpleInitialProbe {

  public static void main(String[] args) {
    SystemInterface system = new SystemSrkEos(300.0, 20.0);
    system.addComponent("methane", 0.5);
    system.addComponent("n-butane", 0.5);
    system.createDatabase(true);
    system.setMixingRule(2);
    new ThermodynamicOperations(system).TPflash();
    system.initPhysicalProperties();

    int count = system.getPhase(0).getNumberOfComponents();
    SystemInterface local = system.clone();
    for (int i = 0; i < count; i++) {
      local.addComponent(local.getPhase(0).getComponent(i).getName(),
          -local.getPhase(0).getComponent(i).getNumberOfmoles());
      local.addComponent(local.getPhase(0).getComponent(i).getName(),
          system.getPhase(0).getComponent(i).getx() / system.getPhase(0).getMolarVolume());
    }
    local.init(0);
    local.setUseTVasIndependentVariables(true);
    local.setNumberOfPhases(1);
    local.getPhase(0).setTotalVolume(1.0);
    local.useVolumeCorrection(false);
    local.init_x_y();
    local.init(3);

    System.out.println("pressure_bar=" + local.getPhase(0).getPressure());
    System.out.println("molar_volume=" + local.getPhase(0).getMolarVolume());
    for (int i = 0; i < count; i++) {
      System.out.println("component_" + i + "_x=" + local.getPhase(0).getComponent(i).getx());
      System.out.println("component_" + i + "_mu="
          + local.getPhase(0).getComponent(i).getChemicalPotential(local.getPhase(0)));
      for (int k = 0; k < count; k++) {
        System.out.println("dmudn_" + i + k + "="
            + local.getPhase(0).getComponent(i).getChemicalPotentialdNTV(k, local.getPhase(0)));
      }
    }
  }
}
