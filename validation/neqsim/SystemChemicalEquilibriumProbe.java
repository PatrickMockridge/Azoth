// The 100-line operation, for `reactions.system_chemical_equilibrium`.
//
// `thermodynamicoperations/chemicalequilibrium/ChemicalEquilibrium` - **a different class from the
// solver of the same name in `chemicalreactions/`** - loops the phases calling `solveChemEq`
// until the composition stops moving, and reports the reaction heat of `HCO3-` before and after.
//
// Three things are measured here, and the third is the one the plan calls the quirk:
//
// 1. **On a neutral fluid the operation does nothing at all.** Its whole body is inside
//    `if (system.isChemicalSystem())`, and a fluid only becomes chemical where water is - which
//    is where the ions are.
// 2. What it leaves behind on a fluid it does act on: the composition and `getDeltaReactionHeat`.
// 3. **`solveChemEq` discards the phase index the loop hands it**
//    (`if (phaseNum != reactivePhase) phaseNum = reactivePhase`), so the loop over phases
//    re-solves the *same* phase every pass. Driving `solveChemEq(1)` and `solveChemEq(0)` on the
//    same two-phase system is what shows it: both move phase 1, and phase 0 never moves.
//
//     javac -proc:none -cp neqsim-f0c7436.jar SystemChemicalEquilibriumProbe.java
//     java -cp .:neqsim-f0c7436.jar SystemChemicalEquilibriumProbe > captures/system_chemical_equilibrium_probe.tsv

import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.chemicalequilibrium.ChemicalEquilibrium;

public class SystemChemicalEquilibriumProbe {

  public static void main(String[] args) {
    neutral();
    one("co2-water-298K", 298.15, 1.01325, new String[] { "CO2", "water" },
        new double[] { 0.01, 10.0 });
    retarget("co2-water-298K-retarget", 298.15, 1.01325, new String[] { "CO2", "water" },
        new double[] { 0.01, 10.0 });
  }

  /// The operation on a fluid that is not chemical: its body never runs.
  static void neutral() {
    SystemInterface system = new SystemSrkEos(298.15, 1.01325);
    system.addComponent("methane", 0.5);
    system.addComponent("CO2", 0.3);
    system.addComponent("hydrogen", 0.2);
    system.chemicalReactionInit();
    system.setMixingRule(2);
    system.init(0);
    System.out.println("neutral_fluid=methane+CO2+hydrogen");
    System.out.println("neutral_is_chemical_system=" + system.isChemicalSystem());
    ChemicalEquilibrium operation = new ChemicalEquilibrium(system);
    operation.run();
    System.out.println("neutral_delta_reaction_heat=" + system.getChemicalReactionOperations()
        .getDeltaReactionHeat());
    System.out.println();
  }

  /// The operation on a fluid it acts on, after a flash has made the phases.
  static void one(String label, double temperature, double pressure, String[] names,
      double[] moles) {
    System.out.println("fluid=" + label);
    SystemInterface system = build(temperature, pressure, names, moles);

    ThermodynamicOperations flash = new ThermodynamicOperations(system);
    flash.TPflash();
    System.out.println("phases_before=" + system.getNumberOfPhases());
    printPhase("before", system);

    double molesInReactivePhase =
        system.getPhase(1).getNumberOfMolesInPhase();
    System.out.println("reactive_phase_moles=" + molesInReactivePhase);

    ChemicalEquilibrium operation = new ChemicalEquilibrium(system);
    operation.run();

    printPhase("after", system);
    System.out.println("delta_reaction_heat=" + system.getChemicalReactionOperations()
        .getDeltaReactionHeat());
    System.out.println();
  }

  /// **Does the loop's phase index reach the solve at all?** `solveChemEq(i)` is called for
  /// `i = 0` and `i = 1` on the same system, and each call's effect on each phase is printed:
  /// the class re-targets to the reactive phase, so the index is a fiction.
  static void retarget(String label, double temperature, double pressure, String[] names,
      double[] moles) {
    System.out.println("fluid=" + label);
    SystemInterface system = build(temperature, pressure, names, moles);
    ThermodynamicOperations flash = new ThermodynamicOperations(system);
    flash.TPflash();
    System.out.println("phases=" + system.getNumberOfPhases());

    // **The claim is read from the source and this only corroborates it.** `solveChemEq`
    // re-targets (`if (phaseNum != reactivePhase) phaseNum = reactivePhase`), so neither call
    // should move phase 0. The flash has already solved the chemistry, so what moves is small -
    // so this is a consistency check rather than a demonstration, and it is labelled as one.
    for (int requested = 0; requested < system.getNumberOfPhases(); requested++) {
      system.init(1);
      double[] before = snapshot(system);
      ChemicalReactionOperations operations = system.getChemicalReactionOperations();
      operations.solveChemEq(requested, 1);
      double[] after = snapshot(system);
      StringBuilder moved = new StringBuilder("requested_phase=" + requested + " moved=");
      for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
        double change = 0.0;
        for (int i = 0; i < names.length; i++) {
          change += Math.abs(after[phase * names.length + i] - before[phase * names.length + i]);
        }
        moved.append(" phase").append(phase).append(":=").append(change).append(" ");
      }
      System.out.println(moved.toString().trim());
    }
    System.out.println();
  }

  static SystemInterface build(double temperature, double pressure, String[] names,
      double[] moles) {
    SystemInterface system = new SystemSrkEos(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule(2);
    system.init(0);
    return system;
  }

  /// Every phase's every mole fraction, flattened, so a call's effect is measurable.
  static double[] snapshot(SystemInterface system) {
    int nc = system.getPhase(0).getNumberOfComponents();
    double[] out = new double[system.getNumberOfPhases() * nc];
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      for (int i = 0; i < nc; i++) {
        out[phase * nc + i] = system.getPhase(phase).getComponent(i).getx();
      }
    }
    return out;
  }

  static void printPhase(String label, SystemInterface system) {
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      StringBuilder composition = new StringBuilder("  " + label + "_phase[" + phase + "]="
          + system.getPhase(phase).getPhaseTypeName() + " beta="
          + system.getPhase(phase).getBeta() + " x=");
      for (int i = 0; i < names(system); i++) {
        ComponentInterface component = system.getPhase(phase).getComponent(i);
        composition.append(component.getName()).append(":").append(component.getx()).append(" ");
      }
      System.out.println(composition.toString().trim());
    }
  }

  static int names(SystemInterface system) {
    return system.getPhase(0).getNumberOfComponents();
  }
}
