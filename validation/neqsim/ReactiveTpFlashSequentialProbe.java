// The sequential chemical dispatch inside `TPflash`, for `reactions.reactive_flash_sequential`.
//
// A chemical system gets a different flash from a physical one: `TPflash`'s successive
// substitution converges on K-values, and then - while the loop is still running - it hands
// every phase **but the first** to `ChemicalReactionOperations.solveChemEq`, measures how far
// that moved the composition, and repeats until the movement stops changing. The result is a
// phase split that is at *both* phase and chemical equilibrium, reached by alternating two
// operations rather than by solving them together.
//
// What the probe prints is the outcome: the phases, their compositions and their fractions.
// The interior - `chemdev` and `diffChem` per round - is local to the class's loop and is not
// exposed, so the port is held to the answer rather than to the path.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ReactiveTpFlashSequentialProbe.java
//     java -cp .:neqsim-f0c7436.jar ReactiveTpFlashSequentialProbe > captures/reactive_tp_flash_sequential_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class ReactiveTpFlashSequentialProbe {
  public static void main(String[] args) {
    // **Which fluids the dispatch is even reachable for.** The loop below fires only when
    // `system.isChemicalSystem()`, so what `chemicalReactionInit` makes of a fluid decides
    // whether the sequential path exists for it - and every fluid it reacts is a fluid it has
    // added ions to.
    neutrality();
    one("co2-water-298K", 298.15, 1.01325, new String[] { "CO2", "water" },
        new double[] { 0.01, 10.0 });
    one("co2-h2s-water-298K", 298.15, 1.01325, new String[] { "CO2", "H2S", "water" },
        new double[] { 0.01, 0.01, 10.0 });
    one("co2-water-313K-10bar", 313.15, 10.0, new String[] { "CO2", "water" },
        new double[] { 0.1, 10.0 });
  }

  /// One line per fluid: how many components the reaction tables added, whether the system is
  /// a chemical one, and whether any of its components is charged.
  ///
  /// Five of these are reactions a chemist would expect - ammonia synthesis, the reverse
  /// water-gas shift, cracking, methanol synthesis, combustion - and **the tables carry none of
  /// them**: a fluid becomes chemical only where water is, and water brings `OH-`, `H3O+` and
  /// the carbonate ions with it.
  static void neutrality() {
    String[][] fluids = {
      { "CO", "water", "CO2", "hydrogen" },
      { "CO2", "water" },
      { "CO2", "hydrogen", "CO", "water" },
      { "nitrogen", "hydrogen" },
      { "nitrogen", "hydrogen", "ammonia" },
      { "CO2", "hydrogen" },
      { "methane", "ethane", "propane" },
      { "methanol", "CO", "hydrogen" },
      { "oxygen", "nitrogen", "CO2" },
    };
    for (String[] names : fluids) {
      SystemInterface system = new SystemSrkEos(298.15, 1.01325);
      for (String name : names) {
        system.addComponent(name, 1.0);
      }
      system.chemicalReactionInit();
      StringBuilder out = new StringBuilder();
      boolean charged = false;
      int count = system.getPhase(0).getNumberOfComponents();
      for (int i = 0; i < count; i++) {
        out.append(system.getPhase(0).getComponent(i).getName()).append(" ");
        if (system.getPhase(0).getComponent(i).getIonicCharge() != 0) {
          charged = true;
        }
      }
      System.out.println("neutrality " + String.join("+", names) + " -> components=" + count
          + " chemical=" + system.isChemicalSystem() + " charged=" + charged + " | "
          + out.toString().trim());
    }
    System.out.println();
  }

  static void one(String label, double temperature, double pressure, String[] names,
      double[] moles) {
    System.out.println("fluid=" + label);
    System.out.println("temperature_K=" + temperature);
    System.out.println("pressure_bara=" + pressure);

    SystemInterface system = new SystemSrkEos(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      System.out.println("  feed[" + names[i] + "]=" + moles[i]);
      system.addComponent(names[i], moles[i]);
    }
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule(2);
    system.init(0);
    System.out.println("is_chemical_system=" + system.isChemicalSystem());

    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    operations.TPflash();

    System.out.println("phases=" + system.getNumberOfPhases());
    System.out.println("total_moles=" + system.getNumberOfMoles());
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      PhaseInterface held = system.getPhase(phase);
      System.out.println("  phase[" + phase + "]=" + held.getPhaseTypeName()
          + " beta=" + held.getBeta() + " moles=" + held.getNumberOfMolesInPhase());
      StringBuilder composition = new StringBuilder("    x=");
      for (int i = 0; i < held.getNumberOfComponents(); i++) {
        ComponentInterface component = held.getComponent(i);
        composition.append(component.getName()).append(":").append(component.getx()).append(" ");
      }
      System.out.println(composition.toString().trim());
    }
    System.out.println();
  }
}
