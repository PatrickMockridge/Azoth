// What `LGTSurfaceTension`'s local one-phase system answers at the profile steps where the
// cubic's pressure is negative. The class records a negative `pressure_interface` and keeps
// going, so its `getChemicalPotential` is finite there - this probe rebuilds that local
// system outside the class and prints the potentials, so a port knows what to reproduce
// rather than guessing at NeqSim's entropy convention.
//
//     javac -proc:none -cp neqsim-f0c7436.jar LgtLocalProbe.java
//     java -cp .:neqsim-f0c7436.jar LgtLocalProbe > captures/lgt_local_probe.tsv

import neqsim.physicalproperties.interfaceproperties.surfacetension.LGTSurfaceTension;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class LgtLocalProbe {

  public static void main(String[] args) {
    SystemInterface system = new SystemSrkEos(300.0, 20.0);
    system.addComponent("methane", 0.5);
    system.addComponent("n-butane", 0.5);
    system.createDatabase(true);
    system.setMixingRule(2);
    new ThermodynamicOperations(system).TPflash();
    system.initPhysicalProperties();

    int count = system.getPhase(0).getNumberOfComponents();
    double[] den1 = new double[count];
    double[] den2 = new double[count];
    double[] delta = new double[count];
    for (int i = 0; i < count; i++) {
      den1[i] = 1e5 * system.getPhase(0).getComponent(i).getx()
          / system.getPhase(0).getMolarVolume();
      den2[i] = 1e5 * system.getPhase(1).getComponent(i).getx()
          / system.getPhase(1).getMolarVolume();
      delta[i] = (den2[i] - den1[i]) / 200.0;
    }

    // The class's own setup, reproduced: a clone turned into a single TV phase.
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

    double[] current = den1.clone();
    for (int j = 1; j < 200; j++) {
      for (int i = 0; i < count; i++) {
        current[i] += delta[i];
        local.addComponent(local.getPhase(0).getComponent(i).getName(), delta[i] / 1.0e5);
      }
      local.init_x_y();
      local.init(3);
      if (j == 1 || j == 52 || j == 53 || j == 60 || j == 100 || j == 150 || j == 199) {
        System.out.println("step=" + j);
        System.out.println("pressure_bar=" + local.getPhase(0).getPressure());
        System.out.println("pressure_Pa=" + 1e5 * local.getPhase(0).getPressure());
        System.out.println("molar_volume=" + local.getPhase(0).getMolarVolume());
        System.out.println("moles=" + local.getPhase(0).getNumberOfMolesInPhase());
        for (int i = 0; i < count; i++) {
          System.out.println("component_" + i + "_x=" + local.getPhase(0).getComponent(i).getx());
          System.out.println("component_" + i + "_mu_J_per_mol="
              + local.getPhase(0).getComponent(i).getChemicalPotential(local.getPhase(0)));
          System.out.println("component_" + i + "_gibbs="
              + local.getPhase(0).getComponent(i).getGibbsEnergy(local.getPhase(0).getTemperature(),
                  local.getPhase(0).getPressure()));
        }
      }
    }

    // The class's own view, for the same steps.
    LGTSurfaceTension lgt = new LGTSurfaceTension(system);
    lgt.calcSurfaceTension(0, 1);
    System.out.println("class_sigma=" + lgt.calcSurfaceTension(0, 1));
  }
}
