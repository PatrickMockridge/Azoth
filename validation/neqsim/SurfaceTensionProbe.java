// The interfacial-tension models `InterfaceProperties` can select, and the default answer.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar SurfaceTensionProbe.java
//     java -cp .:neqsim-f0c7436.jar SurfaceTensionProbe > captures/surface_tension_probe.tsv
//     java -cp .:neqsim-f0c7436.jar SurfaceTensionProbe solvers > captures/surface_tension_solvers_probe.tsv
//
// With no argument the probe records model 0 - the field's initial value and the only default -
// and the two calculators it installs, Parachor and Firozabadi-Ramley, called directly. Model 0
// is parchor for gas/oil and gas/aqueous and Firozabadi for liquid/liquid
// (`InterfaceProperties.java:237-240`).
//
// `solvers` additionally runs models 1, 2 and 4, which build GT, LGT and GT-simple density-profile
// solvers. Those are the next tranche's expected values, and they are kept out of the no-argument
// run so the weekly oracle sweep does not pay for a solver it does not test.
//
// The pair is selected by phase type, not by the arguments' order: GAS/OIL goes to the gas-liquid
// slot, GAS/AQUEOUS to the gas-aqueous one, and everything else - including a reversed oil/gas
// pair - to the liquid-liquid slot (`InterfaceProperties.java:168-176`). The probe prints the
// phase inventory so which slot ran is visible in the capture.

import neqsim.physicalproperties.interfaceproperties.surfacetension.FirozabadiRamleyInterfaceTension;
import neqsim.physicalproperties.interfaceproperties.surfacetension.ParachorSurfaceTension;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class SurfaceTensionProbe {

  public static void main(String[] args) {
    boolean solvers = args.length > 0 && args[0].equals("solvers");
    row("methane_nbutane_300K_20bara", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, false, solvers);
    row("co2_water_313K_50bara", new String[] { "CO2", "water" }, new double[] { 0.02, 0.98 },
        313.15, 50.0, false, solvers);
    row("nheptane_water_300K_1bara", new String[] { "n-heptane", "water" },
        new double[] { 0.5, 0.5 }, 300.0, 1.0, true, solvers);
  }

  static void row(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, boolean multiPhase, boolean solvers) {
    SystemInterface system = new SystemSrkEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.createDatabase(true);
    system.setMixingRule(2);
    if (multiPhase) {
      system.setMultiPhaseCheck(true);
    }
    new ThermodynamicOperations(system).TPflash();
    system.initPhysicalProperties();

    System.out.println(label);
    System.out.println("components=" + String.join(" ", names));
    System.out.println("temperature_K=" + temperatureK);
    System.out.println("pressure_bara=" + pressureBara);
    System.out.println("phase_count=" + system.getNumberOfPhases());
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
    if (system.getNumberOfPhases() < 2) {
      System.out.println("no_interface=true");
      System.out.println();
      return;
    }

    // The components' critical temperatures and parachors, and phase 1's pseudo-critical
    // temperature. Firozabadi-Ramley's `Tr` is that phase's own sum(x_i Tc_i), so a port that
    // takes a pre-computed value has to say which phase it came from.
    StringBuilder tcs = new StringBuilder("component_tc_K=");
    StringBuilder parachors = new StringBuilder("component_parachor=");
    double pseudoCritical = 0.0;
    for (int c = 0; c < system.getPhase(0).getNumberOfComponents(); c++) {
      double tc = system.getPhase(0).getComponent(c).getTC();
      tcs.append(tc).append(" ");
      parachors.append(system.getPhase(0).getComponent(c).getParachorParameter()).append(" ");
      pseudoCritical += system.getPhase(0).getComponent(c).getx() * tc;
    }
    System.out.println(tcs.toString().trim());
    System.out.println(parachors.toString().trim());
    System.out.println("phase0_pseudo_critical_K=" + pseudoCritical);
    double deltaDensity = Math.abs(system.getPhase(1).getPhysicalProperties().getDensity()
        - system.getPhase(0).getPhysicalProperties().getDensity());
    System.out.println("delta_density_kg_per_m3=" + deltaDensity);
    System.out.println("delta_density_over_1000=" + deltaDensity / 1000.0);

    // Model 0 is the default: `interfacialTensionModel` starts at 0 and `init()` re-applies it
    // whenever a slot is null.
    system.getInterphaseProperties().setInterfacialTensionModel(0);
    System.out.println("model_0_sigma_N_per_m="
        + system.getInterphaseProperties().getSurfaceTension(0, 1));

    // The two calculators model 0 installs, called directly on the same system.
    System.out.println("parachor_direct_sigma_N_per_m="
        + new ParachorSurfaceTension(system).calcSurfaceTension(0, 1));
    System.out.println("firozabadi_direct_sigma_N_per_m="
        + new FirozabadiRamleyInterfaceTension(system).calcSurfaceTension(0, 1));

    if (solvers) {
      for (int model : new int[] { 1, 2, 4 }) {
        long started = System.currentTimeMillis();
        try {
          system.getInterphaseProperties().setInterfacialTensionModel(model);
          System.out.println("model_" + model + "_sigma_N_per_m="
              + system.getInterphaseProperties().getSurfaceTension(0, 1));
        } catch (RuntimeException ex) {
          System.out.println("model_" + model + "_error=" + ex.getClass().getSimpleName() + ":"
              + ex.getMessage());
        }
        System.out.println("model_" + model + "_milliseconds="
            + (System.currentTimeMillis() - started));
      }
    }
    System.out.println();
  }
}
