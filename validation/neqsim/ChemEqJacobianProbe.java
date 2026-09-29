// The Newton matrix `ChemicalEquilibrium` installs, in each of its five modes - and whether
// the solve converges in them.
//
// equinor/neqsim#4081 claims the default installs `delta_ik / n_i` and nothing else, and that
// the adaptive path reaches `dfugdN + 1/n_t` for *every* `i, k`: the diagonal cancelled, the
// sign of `1/n_t` flipped, leaving a rank-one singular matrix wherever `dfugdN` is near zero.
// It also proposes using `getFugacitydN` for the entry rather than hand-rolling it.
//
// `M_matrix` is package-private, so it is read reflectively. Each mode is built at the same
// state from a fresh system, converged, and then asked for the matrix it installs - with the
// matrix's rank and singular values printed beside it, because "rank one" is a claim about the
// singular values and the second one is what separates a collapsed matrix from a dominant but
// full-rank one.
//
// **The convergence column is the load-bearing one and the matrices are not.** A mode that
// diverges leaves its matrix describing a diverged state, so the rank numbers for those modes
// are read at whichever state the iteration abandoned - printed for completeness, not as a
// property of the mode at a physical composition.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ChemEqJacobianProbe.java
//     java -Xmx2g -cp .:neqsim-f0c7436.jar ChemEqJacobianProbe > captures/chem_eq_jacobian_probe.tsv

import java.lang.reflect.Field;
import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.chemicalreactions.chemicalequilibrium.ChemicalEquilibrium;
import neqsim.thermo.ThermodynamicConstantsInterface;
import neqsim.thermo.component.Component;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemSrkEos;

public class ChemEqJacobianProbe {

  static double readDouble(Object o, String name) throws Exception {
    Field f = o.getClass().getDeclaredField(name);
    f.setAccessible(true);
    return f.getDouble(o);
  }

  static double[] readVector(Object o, String name) throws Exception {
    Field f = o.getClass().getDeclaredField(name);
    f.setAccessible(true);
    return (double[]) f.get(o);
  }

  static double[][] readMatrix(Object o, String name) throws Exception {
    Field f = o.getClass().getDeclaredField(name);
    f.setAccessible(true);
    return (double[][]) f.get(o);
  }

  static void printMatrix(String name, double[][] m) {
    System.out.println(name + "=" + m.length + "x" + (m.length == 0 ? 0 : m[0].length));
    for (int r = 0; r < m.length; r++) {
      StringBuilder row = new StringBuilder("  " + name + "[" + r + "]=");
      for (int c = 0; c < m[r].length; c++) {
        row.append(String.format("%.15g", m[r][c])).append(c + 1 < m[r].length ? " " : "");
      }
      System.out.println(row);
    }
  }

  static void printVector(String name, double[] v) {
    StringBuilder row = new StringBuilder(name + "=");
    for (int i = 0; i < v.length; i++) {
      row.append(String.format("%.15g", v[i])).append(i + 1 < v.length ? " " : "");
    }
    System.out.println(row);
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

  /** One mode: build the solver at the same state, install the mode, expose its matrix. */
  static void mode(String label, String[] names, double[] moles, boolean full, boolean fug,
      boolean adaptive) throws Exception {
    SystemSrkEos system = build(names, moles);
    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    operations.setReactiveComponents(0);
    ComponentInterface[] components = operations.getComponents();
    int phaseNum = 0;
    double[][] a = operations.calcAmatrix();
    double[] b = operations.calcBVector();

    system.init(1, phaseNum);
    ChemicalEquilibrium solver = new ChemicalEquilibrium(a, b, system, components, phaseNum);
    solver.setUseFullMMatrix(full);
    solver.setUseFugacityDerivatives(fug);
    solver.setUseAdaptiveDerivatives(adaptive);
    solver.calcRefPot();

    // `step()` installs the matrix but needs the solver initialised first - it reads
    // `x_solve`, which the flash-side setup inside `solve()` is what creates. So converge
    // once and then read the matrix the mode installs at the converged state.
    boolean converged = solver.solve();

    double[] nMol = readVector(solver, "n_mol");
    double nTot = readDouble(solver, "n_t");
    System.out.println();
    System.out.println("--- " + label + " ---");
    System.out.printf("useFullMMatrix=%b useFugacityDerivatives=%b useAdaptiveDerivatives=%b%n", full,
        fug, adaptive);
    System.out.printf("n_t=%.15g%n", nTot);
    printVector("n_mol", nMol);
    printVector("chem_ref", readVector(solver, "chem_ref"));
    printVector("logactivity", readVector(solver, "logactivityVec"));

    System.out.printf("solve_converged=%b iterations=%d error=%.15g%n", converged,
        solver.getLastIterationCount(), solver.getLastError());

    // `step()` is what installs a matrix; the residual it is the Jacobian of is printed
    // beside it so the two can be compared rather than assumed.
    double error = solver.step();
    System.out.printf("step_error=%.15g%n", error);
    double[][] m = readMatrix(solver, "M_matrix");
    printMatrix("M_matrix", m);

    // **The rank, because the claim is that the matrix collapses to rank one.** Singular
    // values are printed beside it: a matrix that is *nearly* rank one has a second value
    // far below its first rather than exactly zero, and the ratio is what separates the
    // two cases. Jama is on the classpath because NeqSim uses it.
    Jama.Matrix jm = new Jama.Matrix(m);
    System.out.println("matrix_rank=" + jm.rank());
    double[][] svdS = jm.svd().getS().getArrayCopy();
    StringBuilder svs = new StringBuilder("singular_values=");
    for (int i = 0; i < svdS.length; i++) {
      svs.append(String.format("%.6g", svdS[i][i])).append(i + 1 < svdS.length ? " " : "");
    }
    System.out.println(svs);
    printVector("next_moles", solver.getMoles());

    // The two candidate formulas the claim is written in, from the same state.
    PhaseInterface phase = system.getPhase(phaseNum);
    int nc = components.length;
    StringBuilder dfug = new StringBuilder("dfugdN=");
    for (int i = 0; i < nc; i++) {
      dfug.append(String.format("%.15g", ((Component) components[i]).getFugacitydN(0, phase)))
          .append(i + 1 < nc ? " " : "");
    }
    System.out.println(dfug);
  }

  public static void main(String[] args) throws Exception {
    String[] names = { "CO2", "water" };
    double[] moles = { 0.01, 10.0 };

    System.out.println("fluid=co2-water");
    System.out.printf("R=%.15g%n", ThermodynamicConstantsInterface.R);
    mode("default (diagonal only)", names, moles, false, false, false);
    mode("full M matrix", names, moles, true, false, false);
    mode("fugacity derivatives on", names, moles, false, true, false);
    mode("adaptive derivatives on", names, moles, false, false, true);
    mode("adaptive + full", names, moles, true, false, true);
  }
}
