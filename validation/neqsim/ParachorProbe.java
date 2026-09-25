// The interface (mixture) surface tension of `ParachorSurfaceTension.calcSurfaceTension`, which
// `InterfaceProperties.getSurfaceTension(phase1, phase2)` dispatches to for a gas/oil and for a
// gas/aqueous pair alike.
//
// Compile and run from this directory:
//
//     javac -cp .:neqf0c7436.jar -d . ParachorProbe.java
//     java -cp .:neqsim-f0c7436.jar ParachorProbe > captures/parachor_probe.tsv
//
// **The pair form is not the pure component form.** `eos.parachor_surface_tension` ports
// `calcPureComponentSurfaceTension`, where the mole fractions divide out; the interface form sums
// `P_i * (rho_2/M_2 * x_2i - rho_1/M_1 * x_1i)` over the components - a *molar density* difference
// per component - and raises the total to the fourth. The probe prints every input that sum needs,
// so a port is held to the arithmetic and not to a band: the per-phase density and molar mass, both
// compositions, and each component's own parachor.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class ParachorProbe {

  public static void main(String[] args) {
    // A gas/oil pair, which is `gasLiquidSurfaceTensionCalc`'s branch: the binary column's own feed
    // state, two-phase at 300 K and 20 bara.
    row("methane_nbutane_gas_oil_300K_20bar", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0);
    // A three-component oil, so the sum is over more than two terms.
    row("methane_nbutane_nheptane_gas_oil_300K_20bar",
        new String[] { "methane", "n-butane", "n-heptane" }, new double[] { 0.4, 0.4, 0.2 }, 300.0,
        20.0);
    // A gas/aqueous pair, which is `gasAqueousSurfaceTensionCalc`'s branch - the same
    // `ParachorSurfaceTension` calculator by default, and the branch the rate-based column's own
    // CO2 absorber state takes.
    row("co2_water_gas_aqueous_313K_50bar", new String[] { "CO2", "water" },
        new double[] { 0.02, 0.98 }, 313.15, 50.0);
    // A methane-rich gas over water, which is what a dehydration or hydrate state looks like.
    row("methane_water_gas_aqueous_300K_50bar", new String[] { "methane", "water" },
        new double[] { 0.05, 0.95 }, 300.0, 50.0);
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
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      System.out.println("phase" + i + "_type=" + system.getPhase(i).getType());
      System.out.println("phase" + i + "_temperature_K=" + system.getPhase(i).getTemperature());
      System.out.println("phase" + i + "_pressure_bara=" + system.getPhase(i).getPressure());
      System.out.println(
          "phase" + i + "_density_kg_per_m3=" + system.getPhase(i).getPhysicalProperties().getDensity());
      System.out.println("phase" + i + "_molar_mass_kg_per_mol=" + system.getPhase(i).getMolarMass());
      StringBuilder composition = new StringBuilder("phase" + i + "_x=");
      for (int c = 0; c < system.getPhase(i).getNumberOfComponents(); c++) {
        composition.append(names[c]).append(":").append(system.getPhase(i).getComponent(c).getx())
            .append(" ");
      }
      System.out.println(composition.toString().trim());
    }
    StringBuilder parachors = new StringBuilder("parachors=");
    for (int c = 0; c < names.length; c++) {
      parachors.append(names[c]).append(":")
          .append(system.getPhase(0).getComponent(c).getParachorParameter()).append(" ");
    }
    System.out.println(parachors.toString().trim());

    if (system.getNumberOfPhases() < 2) {
      System.out.println("surface_tension_N_per_m=none");
      System.out.println();
      return;
    }
    int gas = -1;
    int liquid = -1;
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      String type = system.getPhase(i).getType().toString();
      if ("GAS".equals(type) && gas < 0) {
        gas = i;
      } else if (liquid < 0) {
        liquid = i;
      }
    }
    System.out.println("gas_phase_index=" + gas);
    System.out.println("liquid_phase_index=" + liquid);
    double sigma = system.getInterphaseProperties().getSurfaceTension(gas, liquid);
    System.out.println("surface_tension_N_per_m=" + sigma);
    System.out.println();
  }
}
