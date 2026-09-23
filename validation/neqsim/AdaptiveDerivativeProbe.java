// The adaptive-derivative refinement, for `reactions.reactive_phase_equilibrium`'s second pass.
//
// `ChemicalReactionOperations.solveChemEq` runs **two** refinements and switches
// `useAdaptiveDerivatives` on from the second, so the second is only reached where the first
// failed the certificate - which is every captured fluid. What the second pass changes is the
// **M matrix**: with the switch on, `M[i][k]` is `dfugdn(i, k) + 1/n_t` rather than the ideal
// `delta_ik / n_i`.
//
// So this probe prints three things per fluid:
//
// 1. `getdfugdn(i, k)` at the phase's own composition and moles - **the derivative surface the
//    port has to reproduce**, and the quantity whose *scaling* is the open question (azoth's
//    `PhaseDerivatives::d_ln_phi_dn` is at a total mole number of one, NeqSim's is at the
//    phase's actual moles);
// 2. the ideal solve (`setUseAdaptiveDerivatives(false)`), which is what the library ports;
// 3. the adaptive solve (`setUseAdaptiveDerivatives(true)`), which is the second refinement.
//
// The last two are printed beside each other because the question the plan asks is whether the
// second refinement changes the answer or only the path.
//
//     javac -proc:none -cp neqsim-f0c7436.jar AdaptiveDerivativeProbe.java
//     java -cp .:neqsim-f0c7436.jar AdaptiveDerivativeProbe > captures/adaptive_derivative_probe.tsv

import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.chemicalreactions.chemicalequilibrium.ChemicalEquilibrium;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class AdaptiveDerivativeProbe {
  public static void main(String[] args) {
    one("co2-water", new String[] { "CO2", "water" }, new double[] { 0.01, 10.0 });
    one("co2-h2s-water", new String[] { "CO2", "H2S", "water" },
        new double[] { 0.01, 0.01, 10.0 });
    one("water-meg", new String[] { "water", "MEG" }, new double[] { 1.0, 1.0 });
  }

  static SystemSrkEos build(String[] names, double[] moles) {
    SystemSrkEos system = new SystemSrkEos(298.15, 1.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    // The same three calls the chemical-equilibrium probe makes, in the same order:
    // `chemicalReactionInit` is what populates the reaction operations, and `createDatabase`
    // re-sizes the interaction matrices for the components it added.
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule(2);
    system.init(0);
    return system;
  }

  static void one(String label, String[] names, double[] moles) {
    System.out.println("fluid=" + label);

    // The reaction operation builds the matrix and the amounts, and its readers are public -
    // so the setup is the class's own rather than a re-derivation, exactly as the chemical
    // equilibrium probe does it.
    SystemSrkEos setup = build(names, moles);
    ChemicalReactionOperations operations = setup.getChemicalReactionOperations();
    operations.setReactiveComponents(0);
    ComponentInterface[] reactive = operations.getComponents();
    double[][] aMatrix = operations.calcAmatrix();
    double[] bVector = operations.calcBVector();
    System.out.println("reactive_components=" + reactive.length);

    // The components the matrix is indexed by are the operation's, in its own order, and the
    // solve is handed the same list.
    String[] solved = new String[reactive.length];
    for (int i = 0; i < reactive.length; i++) {
      solved[i] = reactive[i].getName();
    }
    System.out.println("solved_components=" + String.join(" ", solved));

    SystemInterface system = build(names, moles);
    system.init(1);
    system.init(2);
    int phase = 0;

    double molesInPhase = system.getPhase(phase).getNumberOfMolesInPhase();
    System.out.println("phase_moles=" + molesInPhase);

    // **The derivative surface**, at the phase's own moles and composition. **`init(3)` comes
    // first**, because the class's own `chemSolve` calls it before reading a derivative and a
    // reader taken before it is zero for every pair.
    system.init(3, phase);
    StringBuilder line = new StringBuilder("dfugdn=");
    for (int i = 0; i < reactive.length; i++) {
      for (int k = 0; k < reactive.length; k++) {
        line.append(system.getPhase(phase).getComponent(reactive[i].getName()).getdfugdn(k));
        line.append(" ");
      }
    }
    System.out.println(line.toString().trim());

    // The ideal solve: this library's path, and the first refinement.
    double[] ideal = solveOnce(aMatrix, bVector, system, reactive, phase, false);
    printMoles("ideal_moles", ideal);

    // **The adaptive solve runs on the same system**, because that is what the class's second
    // refinement does: `solveChemEq` builds a new solver over a system the first refinement has
    // already left at its answer. An adaptive solve from the *feed* is a different run, and one
    // this probe's first draft took by mistake.
    double[] adaptive = solveOnce(aMatrix, bVector, system, reactive, phase, true);
    printMoles("adaptive_moles", adaptive);

    // And the production path, which runs both refinements: `solveChemEq(phaseNum, type)`, the
    // same call the operations probe makes.
    SystemSrkEos production = build(names, moles);
    production.init(1);
    ChemicalReactionOperations ops = production.getChemicalReactionOperations();
    boolean converged = ops.solveChemEq(phase, 1);
    System.out.println("solveChemEq_converged=" + converged);
    double[] fromProduction = new double[reactive.length];
    for (int i = 0; i < reactive.length; i++) {
      fromProduction[i] = production.getPhase(phase)
          .getComponent(reactive[i].getName()).getNumberOfMolesInPhase();
    }
    printMoles("solveChemEq_moles", fromProduction);
    System.out.println();
  }

  /// One `ChemicalEquilibrium` solve, ideal or adaptive, with its own numbers printed.
  static double[] solveOnce(double[][] aMatrix, double[] bVector, SystemInterface system,
      ComponentInterface[] components, int phase, boolean adaptive) {
    ChemicalEquilibrium solver =
        new ChemicalEquilibrium(aMatrix, bVector, system, components, phase);
    solver.setUseAdaptiveDerivatives(adaptive);
    boolean converged = solver.solve();
    String tag = adaptive ? "adaptive" : "ideal";
    System.out.println(tag + "_converged=" + converged);
    System.out.println(tag + "_iterations=" + solver.getLastIterationCount());
    System.out.println(tag + "_error=" + solver.getLastError());
    return solver.getMoles();
  }

  static void printMoles(String label, double[] values) {
    StringBuilder out = new StringBuilder(label + "=");
    for (int i = 0; i < values.length; i++) {
      out.append(values[i]).append(i + 1 < values.length ? " " : "");
    }
    System.out.println(out);
  }
}
