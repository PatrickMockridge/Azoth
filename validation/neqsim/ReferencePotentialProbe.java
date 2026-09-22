// The independent reaction basis and the reference potentials, for
// `reactions.reference_potentials`.
//
// `ChemicalReactionList.calcReferencePotentials` decides which components' standard-state
// potentials to solve for directly and which to propagate from the rest. Both halves are
// greedy rank tests on a `Jama.Matrix`, so the answer depends on JAMA's `rank()` and on
// nothing else - and JAMA's rank is a singular-value count against
// `max(m, n) * s[0] * 2**-52`.
//
// **The singular values are printed beside that tolerance, because the size of the gap
// decides how much of JAMA has to be reproduced.** If every singular value the test
// rejects is orders of magnitude below the cutoff, then any correct rank routine selects
// the same basis and the tolerance is not load-bearing; if one sits near it, the
// Golub-Reinsch SVD is part of the port and there is no shortcut.
//
// The matrices are the real ones: three fluids on a plain `SystemSrkEos`, whose
// `chemicalReactionInit` loads the standard source's active rows and whose components are
// the ones the element table carries. Three, because the size of the gap is the whole
// question and one fluid answers it only for that fluid's matrix. Nothing here
// reimplements the correlation - the matrix, the rank calls and the potentials are
// NeqSim's own.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ReferencePotentialProbe.java
//     java -cp .:neqsim-f0c7436.jar ReferencePotentialProbe > captures/reference_potential_probe.tsv

import java.util.ArrayList;
import Jama.Matrix;
import Jama.SingularValueDecomposition;
import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.chemicalreactions.chemicalreaction.ChemicalReaction;
import neqsim.chemicalreactions.chemicalreaction.ChemicalReactionList;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemSrkEos;

public class ReferencePotentialProbe {
  static final double EPS = Math.pow(2.0, -52.0);

  public static void main(String[] args) {
    // Several fluids, because the size of the gap between the smallest rejected singular
    // value and the cutoff is the whole question. More reactions means a larger matrix.
    one("co2-water", new String[] { "CO2", "water" }, new double[] { 0.01, 10.0 });
    one("co2-h2s-water", new String[] { "CO2", "H2S", "water" }, new double[] { 0.01, 0.01, 10.0 });
    one("co2-h2s-mdea-water", new String[] { "CO2", "H2S", "MDEA", "water" },
        new double[] { 0.01, 0.01, 1.0, 10.0 });
  }

  static void one(String label, String[] names, double[] moles) {
    System.out.println("fluid=" + label);
    SystemSrkEos system = new SystemSrkEos(298.15, 1.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule(2);
    system.init(0);
    system.setNumberOfPhases(1);
    system.setMaxNumberOfPhases(1);

    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    ChemicalReactionList list = operations.getReactionList();
    operations.setReactiveComponents(0);
    ComponentInterface[] components = operations.getComponents();
    PhaseInterface phase = system.getPhase(0);

    System.out.println("components=" + components.length);
    for (int i = 0; i < components.length; i++) {
      System.out.println("  component[" + i + "]=" + components[i].getName());
    }

    ArrayList<ChemicalReaction> reactions = list.getChemicalReactionList();
    System.out.println("reactions=" + reactions.size());
    for (int i = 0; i < reactions.size(); i++) {
      ChemicalReaction reaction = reactions.get(i);
      StringBuilder nameList = new StringBuilder();
      for (String name : reaction.getNames()) {
        nameList.append(name).append(" ");
      }
      System.out.println("  reaction[" + i + "]=" + reaction.getName() + " names=" + nameList.toString().trim());
    }

    list.createReactionMatrix(phase, components);
    double[][] g = list.getReactionGMatrix();
    int nRows = g.length;
    int nCols = nRows == 0 ? 0 : g[0].length - 1;
    System.out.println("gmatrix=" + nRows + "x" + (nCols + 1));

    for (int r = 0; r < nRows; r++) {
      StringBuilder row = new StringBuilder("  G[" + r + "]=");
      for (int c = 0; c < g[r].length; c++) {
        row.append(g[r][c]).append(c + 1 < g[r].length ? " " : "");
      }
      System.out.println(row);
    }

    // The greedy column selection, replayed so the rank at each step and the singular
    // values behind it are in the capture. Same loop as `calcReferencePotentials`.
    ArrayList<Integer> independent = new ArrayList<Integer>();
    ArrayList<Integer> dependent = new ArrayList<Integer>();
    Matrix current = null;

    for (int j = 0; j < nCols; j++) {
      Matrix next;
      if (current == null) {
        next = new Matrix(nRows, 1);
        for (int i = 0; i < nRows; i++) {
          next.set(i, 0, g[i][j]);
        }
      } else {
        next = new Matrix(nRows, current.getColumnDimension() + 1);
        next.setMatrix(0, nRows - 1, 0, current.getColumnDimension() - 1, current);
        for (int i = 0; i < nRows; i++) {
          next.set(i, current.getColumnDimension(), g[i][j]);
        }
      }

      int currentRank = current == null ? 0 : current.rank();
      int nextRank = next.rank();

      double[] singular = new SingularValueDecomposition(next).getSingularValues();
      double tolerance = Math.max(next.getRowDimension(), next.getColumnDimension())
          * singular[0] * EPS;

      System.out.println("  step[col=" + j + "] current_rank=" + currentRank + " next_rank="
          + nextRank + " tolerance=" + tolerance + " smallest_sv=" + singular[singular.length - 1]
          + " largest_sv=" + singular[0]);

      if (nextRank > currentRank) {
        current = next;
        independent.add(j);
        if (independent.size() == nRows) {
          for (int k = j + 1; k < nCols; k++) {
            dependent.add(k);
          }
          break;
        }
      } else {
        dependent.add(j);
      }
    }

    StringBuilder indep = new StringBuilder();
    for (int c : independent) {
      indep.append(c).append(" ");
    }
    StringBuilder dep = new StringBuilder();
    for (int c : dependent) {
      dep.append(c).append(" ");
    }
    System.out.println("independent_columns=" + indep.toString().trim());
    System.out.println("dependent_columns=" + dep.toString().trim());
    System.out.println("basis_rank_too_low=" + (independent.size() < nRows));

    double[] potentials = list.calcReferencePotentials();
    System.out.println("potentials=" + potentials.length);
    for (int i = 0; i < potentials.length; i++) {
      System.out.println("  potential[" + i + "]=" + potentials[i]);
    }
  }
}
