// `ChemicalReactionList.calcReferencePotentials`' fallbacks, for
// `reactions.reference_potentials`.
//
// The class answers from a component's Gibbs energy of formation in two places: when the
// component rank falls below the reaction count it returns null and its caller
// `ChemicalReactionOperations.calcChemRefPot` reads that as every component's Gibbs energy
// of formation; and when the propagation deadlocks it seeds the first uncomputed dependent
// component with that same number and carries on from it.
//
// **The deadlock needs a component order, and the order is the caller's.**
// `createReactionMatrix(phase, components)` stores the array it is handed and
// `calcReferencePotentials` indexes it, so this probe passes the order under test
// explicitly. `2 H2O = H3O+ + OH-` is one reaction over three components: water is the
// only column that raises the rank to the reaction count, so it is solved for, and the
// two ions have no reaction in which every other component is known - which is the
// deadlock, exactly as `reactions.reference_potentials`' own sweep finds it.
//
// The seed is printed beside each component's `getGibbsEnergyOfFormation()`, because the
// claim is that they are the same number taken raw and not negated: azoth reads that
// column out of the same `COMP.csv` row.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ReferenceFallbackProbe.java
//     java -cp .:neqsim-f0c7436.jar ReferenceFallbackProbe > captures/reference_fallback_probe.tsv

import java.util.ArrayList;
import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.chemicalreactions.chemicalreaction.ChemicalReaction;
import neqsim.chemicalreactions.chemicalreaction.ChemicalReactionList;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemSrkEos;

public class ReferenceFallbackProbe {

  public static void main(String[] args) {
    // The one order that deadlocks on the standard source: the two ions are dependent and
    // neither has a reaction whose other components are all known.
    one("water-h3o-oh", new String[] { "water", "H3O+", "OH-" }, new double[] { 10.0, 1.0e-7, 1.0e-7 });
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
    ComponentInterface[] loaded = operations.getComponents();
    PhaseInterface phase = system.getPhase(0);

    System.out.println("loaded_components=" + loaded.length);
    for (int i = 0; i < loaded.length; i++) {
      System.out.println("  loaded[" + i + "]=" + loaded[i].getName());
    }

    // The order under test, which is the caller's and is what the deadlock depends on.
    ComponentInterface[] components = new ComponentInterface[names.length];
    for (int i = 0; i < names.length; i++) {
      for (ComponentInterface candidate : loaded) {
        if (candidate.getName().equals(names[i])) {
          components[i] = candidate;
        }
      }
      if (components[i] == null) {
        System.out.println("missing_component=" + names[i]);
        return;
      }
      System.out.println("  component[" + i + "]=" + components[i].getName());
    }

    for (int i = 0; i < components.length; i++) {
      System.out.println("  gibbs_formation[" + components[i].getName() + "]="
          + components[i].getGibbsEnergyOfFormation());
    }

    ArrayList<ChemicalReaction> reactions = list.getChemicalReactionList();
    System.out.println("reactions=" + reactions.size());
    for (int i = 0; i < reactions.size(); i++) {
      ChemicalReaction reaction = reactions.get(i);
      StringBuilder nameList = new StringBuilder();
      for (String name : reaction.getNames()) {
        nameList.append(name).append(" ");
      }
      System.out.println("  reaction[" + i + "]=" + reaction.getName() + " names="
          + nameList.toString().trim());
    }

    // The matrix this order produces, so the capture shows there is one independent
    // column and two dependents rather than asserting it.
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

    double[] potentials = list.calcReferencePotentials();
    System.out.println("potentials_null=" + (potentials == null));
    if (potentials != null) {
      for (int i = 0; i < potentials.length; i++) {
        System.out.println("  potential[" + components[i].getName() + "]=" + potentials[i]);
      }
    }
  }
}
