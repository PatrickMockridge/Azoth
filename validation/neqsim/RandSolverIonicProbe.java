// The RAND solver's **ionic branch**, for the port of `reactions.reactive_tp_flash`'s refusal.
//
// `ModifiedRANDSolver.computeG0` has two branches for an ion, and which one runs is a property
// of the *system class*: `isElectrolyteEOS = (system instanceof SystemFurstElectrolyteEos)`,
// and only there is `lnPhiRef[i] = ph.getLogInfiniteDiluteFugacity(i, solvent)` computed. On
// every other system the correction is zero and `g0[ion] = dGfAq / (R T)` stands alone - which
// is what the class's own comment says, and what this probe measures rather than quotes.
//
// Three things are printed per fluid, because the port needs all three and only the first is
// visible in a flash's answer:
//
//   * `is_electrolyte_eos`, which selects the branch;
//   * `ln_phi_ref[i]` and `g0[i]` for every component, which are the branch's own numbers;
//   * the answer, so that a branch which runs and changes nothing is distinguishable from one
//     that does not run.
//
// **The ions are the ones the reaction tables carry**, and the fluid is the class's own: water
// with sodium chloride, which is what `SystemFurstElectrolyteEos` is for. A neutral control on
// the same solver shows the other branch.
//
//     javac -proc:none -cp neqsim-f0c7436.jar RandSolverIonicProbe.java
//     java -cp .:neqsim-f0c7436.jar RandSolverIonicProbe > captures/rand_solver_ionic_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemFurstElectrolyteEos;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.flashops.reactiveflash.FormulaMatrix;
import neqsim.thermodynamicoperations.flashops.reactiveflash.ModifiedRANDSolver;

public class RandSolverIonicProbe {

  public static void main(String[] args) {
    one("nacl-brine-furst", brine());
    one("wgs-srk-control", waterGasShift());
  }

  /// Water with sodium chloride: the fluid `SystemFurstElectrolyteEos` exists for.
  static SystemInterface brine() {
    SystemFurstElectrolyteEos system = new SystemFurstElectrolyteEos(298.15, 1.01325);
    system.addComponent("water", 55.5);
    system.addComponent("Na+", 1.0);
    system.addComponent("Cl-", 1.0);
    system.chemicalReactionInit();
    system.setMixingRule(10);
    system.setMultiPhaseCheck(false);
    return system;
  }

  /// The neutral control: the water-gas shift the port's own capture already holds.
  static SystemInterface waterGasShift() {
    SystemSrkEos system = new SystemSrkEos(600.0, 1.0);
    system.addComponent("CO", 0.25);
    system.addComponent("water", 0.25);
    system.addComponent("CO2", 0.25);
    system.addComponent("hydrogen", 0.25);
    system.setMixingRule("classic");
    return system;
  }

  static void one(String label, SystemInterface system) {
    System.out.println("fluid=" + label);
    system.init(0);
    system.init(1);
    System.out.println("system_class=" + system.getClass().getSimpleName());
    System.out.println("is_chemical_system=" + system.isChemicalSystem());
    System.out.println("has_ions=" + system.hasIons());

    FormulaMatrix matrix = new FormulaMatrix(system);
    String[] names = matrix.getComponentNames();
    System.out.println("reactive_components=" + String.join(" ", names));
    System.out.println("has_ionic_species=" + matrix.hasIonicSpecies());
    StringBuilder charges = new StringBuilder();
    for (int i = 0; i < names.length; i++) {
      charges.append(matrix.getIonicCharges()[i]).append(i + 1 < names.length ? " " : "");
    }
    System.out.println("ionic_charges=" + charges);

    ModifiedRANDSolver solver = new ModifiedRANDSolver(system, matrix);
    System.out.println("is_electrolyte_eos=" + solver.isElectrolyteEOS());

    // The pure aqueous Gibbs energy each ionic `g0` is corrected from, so the branch's two
    // operands are both on the record rather than only their difference.
    StringBuilder aqueous = new StringBuilder();
    for (int i = 0; i < names.length; i++) {
      aqueous.append(system.getPhase(0).getComponent(i).getGibbsEnergyOfFormation())
          .append(i + 1 < names.length ? " " : "");
    }
    System.out.println("gibbs_formation=" + aqueous);

    // `computeG0` runs inside the solve - the arrays are null until it does - so the branch's
    // own numbers are read afterwards and labelled as the state they are a property of.
    solver.solve();
    System.out.println("ln_phi_ref=" + join(solver.getLnPhiRef()));
    System.out.println("g0=" + join(solver.getG0()));
    System.out.println("ran_to=" + solver.getIterationsUsed() + " iterations");
    System.out.println("final_residual=" + solver.getFinalResidual());
    System.out.println("final_element_residual=" + solver.getFinalElementResidual());
    System.out.println("total_moles=" + solver.getTotalMoles());
    for (int phase = 0; phase < solver.getMoleFractions().length; phase++) {
      System.out.println("phase_x[" + phase + "]=" + join(solver.getMoleFractions()[phase]));
    }
    for (int phase = 0; phase < solver.getMoles().length; phase++) {
      System.out.println("phase_moles[" + phase + "]=" + join(solver.getMoles()[phase]));
    }
    System.out.println();
  }

  static String join(double[] values) {
    StringBuilder line = new StringBuilder();
    for (int i = 0; i < values.length; i++) {
      line.append(values[i]).append(i + 1 < values.length ? " " : "");
    }
    return line.toString();
  }
}
