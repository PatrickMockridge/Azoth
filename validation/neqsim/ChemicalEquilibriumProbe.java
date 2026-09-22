// The Smith-Missen equilibrium solve, for `reactions.chemical_equilibrium`.
//
// Two things are printed per fluid, because they are two different measurements.
//
// **The direct solve** constructs `ChemicalEquilibrium` and calls `solve()` once. That is
// the *ideal* path - `M[i][k] = delta_ik / n_i`, no fugacity derivatives - which is what
// `solveChemEq` gets from its first refinement, and it is what this library ports. Its
// inputs are printed as well as its answer: the element matrix, the element amounts, the
// starting moles, the reduced reference potentials and the activity coefficients, so a
// case can be built from the capture rather than from a re-derivation.
//
// `chem_ref` and `logactivityVec` are private on the class, so they are recomputed here
// the way `calcRefPot` computes them - `getReferencePotential() / (R T)`, and the phase's
// own `getLogActivityCoefficient` where `calcActivity()` says there is one. That is the
// same call and not a second implementation.
//
// **`solveChemEq`** is the production entry point and is printed beside it: it runs an LP
// initial estimate and then a refinement loop that switches fugacity derivatives on from
// the second pass. The two can disagree, and one of these fluids is where they do.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ChemicalEquilibriumProbe.java
//     java -cp .:neqsim-f0c7436.jar ChemicalEquilibriumProbe > captures/chemical_equilibrium_probe.tsv

import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.chemicalreactions.chemicalequilibrium.ChemicalEquilibrium;
import neqsim.thermo.ThermodynamicConstantsInterface;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemSrkEos;

public class ChemicalEquilibriumProbe {
  public static void main(String[] args) {
    one("co2-water", new String[] { "CO2", "water" }, new double[] { 0.01, 10.0 });
    one("co2-h2s-water", new String[] { "CO2", "H2S", "water" },
        new double[] { 0.01, 0.01, 10.0 });
  }

  static SystemSrkEos build(String[] names, double[] moles) {
    SystemSrkEos system = new SystemSrkEos(298.15, 1.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule(2);
    system.init(0);
    return system;
  }

  static void one(String label, String[] names, double[] moles) {
    System.out.println("fluid=" + label);

    // --- the direct ideal-path solve -------------------------------------
    SystemSrkEos direct = build(names, moles);
    ChemicalReactionOperations operations = direct.getChemicalReactionOperations();
    operations.setReactiveComponents(0);
    ComponentInterface[] components = operations.getComponents();
    int phaseNum = 0;

    System.out.println("components=" + components.length);
    for (int i = 0; i < components.length; i++) {
      System.out.println("  component[" + i + "]=" + components[i].getName()
          + " reference_state=" + components[i].getReferenceStateType());
    }

    double[][] a = operations.calcAmatrix();
    operations.calcBVector();
    double[] b = operations.calcBVector();
    double[] n = new double[components.length];
    for (int i = 0; i < components.length; i++) {
      n[i] = direct.getPhase(phaseNum).getComponent(components[i].getComponentNumber())
          .getNumberOfMolesInPhase();
    }

    printMatrix("A", a);
    printVector("b", b);
    printVector("n0", n);

    // Recomputed the way `calcRefPot` does it, because the fields are private. **The
    // `init(1)` comes first**, because the constructor does it first: the activity
    // coefficients are a level-1 quantity, and reading them before this gives `-Infinity`
    // for every species.
    direct.init(1, phaseNum);
    PhaseInterface phase = direct.getPhase(phaseNum);
    double temperature = phase.getTemperature();
    double[] chemRef = new double[components.length];
    double[] logActivity = new double[components.length];
    int waterNumb = 0;
    for (int i = 0; i < components.length; i++) {
      if (components[i].getComponentName().equals("water")) {
        waterNumb = i;
        break;
      }
    }
    for (int i = 0; i < components.length; i++) {
      chemRef[i] = components[i].getReferencePotential()
          / (ThermodynamicConstantsInterface.R * temperature);
      logActivity[i] = 0.0;
      if (components[i].calcActivity()) {
        logActivity[i] = phase.getLogActivityCoefficient(components[i].getComponentNumber(),
            components[waterNumb].getComponentNumber());
      }
    }
    System.out.println("temperature=" + temperature);
    System.out.println("moles_in_phase=" + phase.getNumberOfMolesInPhase());
    printVector("chem_ref_reduced", chemRef);
    printVector("log_activity", logActivity);

    ChemicalEquilibrium solver =
        new ChemicalEquilibrium(a, b, direct, components, phaseNum);
    boolean converged = solver.solve();
    System.out.println("direct_converged=" + converged);
    System.out.println("direct_iterations=" + solver.getLastIterationCount());
    System.out.println("direct_error=" + solver.getLastError());
    System.out.println("direct_tolerance=" + solver.getConvergenceTolerance());
    printVector("direct_moles", solver.getMoles());

    // --- the production entry point --------------------------------------
    SystemSrkEos production = build(names, moles);
    ChemicalReactionOperations ops2 = production.getChemicalReactionOperations();
    ops2.setReactiveComponents(0);
    boolean solveChemEq = ops2.solveChemEq(0, 0);
    System.out.println("solveChemEq_converged=" + solveChemEq);
    System.out.println("solveChemEq_charge_mols=" + ops2.getReactivePhaseChargeMoles());
    System.out.println("solveChemEq_max_element_residual="
        + ops2.getMaximumAbsoluteElementBalanceResidual());
    System.out.println("solveChemEq_max_reaction_log_residual="
        + ops2.getMaximumAbsoluteReactionLogResidual());
  }

  static void printMatrix(String name, double[][] m) {
    System.out.println(name + "=" + m.length + "x" + (m.length == 0 ? 0 : m[0].length));
    for (int r = 0; r < m.length; r++) {
      StringBuilder row = new StringBuilder("  " + name + "[" + r + "]=");
      for (int c = 0; c < m[r].length; c++) {
        row.append(m[r][c]).append(c + 1 < m[r].length ? " " : "");
      }
      System.out.println(row);
    }
  }

  static void printVector(String name, double[] v) {
    StringBuilder out = new StringBuilder(name + "=");
    for (int i = 0; i < v.length; i++) {
      out.append(v[i]).append(i + 1 < v.length ? " " : "");
    }
    System.out.println(out);
  }
}
