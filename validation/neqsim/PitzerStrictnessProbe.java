// The Pitzer source's validation gate, for `reactions.reference_potentials`.
//
// `ChemicalReactionList.requireValidatedEvidenceForActiveReactions` throws for any
// **active** row whose `validationstatus` is not `VALIDATED`, and only the pitzer source
// requires it (`ChemicalReactionDataSource.PITZER` is the one enum constant with the flag
// set). The gate runs at `ChemicalReactionOperations.java:179`, **after**
// `removeJunkReactions` and `removeDependentReactions`, so what it judges is the survivor
// set - not the table.
//
// Because it sits inside the operations constructor, a fluid that trips it cannot be
// initialised at all: this probe builds each fluid as a `SystemPitzer`, calls
// `chemicalReactionInit()` inside a try/catch, and prints whether the gate threw and what
// it named. For the fluids that survive it also prints **every loaded row with its
// status**, so the survivor set itself is oracled rather than inferred from the refusal.
//
//     43 of the table's 47 rows are not `VALIDATED`, and the four that are - `CO2water`,
//     `waterreac`, `carbonate` and `water-H2S` - are the whole of what a plain CO2/water
//     fluid can run. So CO2/water is the fluid where the gate fires on nothing, and the
//     other fluids here are the ones where it fires.
//
//     javac -proc:none -cp neqsim-f0c7436.jar PitzerStrictnessProbe.java
//     java -cp .:neqsim-f0c7436.jar PitzerStrictnessProbe > captures/pitzer_strictness_probe.tsv

import java.util.ArrayList;

import neqsim.chemicalreactions.chemicalreaction.ChemicalReaction;
import neqsim.chemicalreactions.chemicalreaction.ChemicalReactionList;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPitzer;

public class PitzerStrictnessProbe {

  public static void main(String[] args) {
    fluid("co2-water", 298.15, 10.0, new String[] { "CO2", "water" }, new double[] { 0.5, 55.5 });
    fluid("h2s-water", 298.15, 10.0, new String[] { "H2S", "water" }, new double[] { 0.5, 55.5 });
    fluid("co2-h2s-water", 298.15, 10.0, new String[] { "CO2", "H2S", "water" },
        new double[] { 0.5, 0.5, 55.5 });
    fluid("methane-oxygen-water", 298.15, 10.0, new String[] { "methane", "oxygen", "water" },
        new double[] { 0.5, 0.5, 55.5 });
    fluid("mdea-water-co2", 313.15, 5.0, new String[] { "MDEA", "water", "CO2" },
        new double[] { 1.0, 10.0, 0.1 });
    fluid("methane-water-nacl", 313.15, 10.0, new String[] { "methane", "water", "Na+", "Cl-" },
        new double[] { 5.0, 55.5, 1.0, 1.0 });
  }

  /// One fluid: does `chemicalReactionInit` throw, and what does the survivor list hold?
  static void fluid(String label, double temperature, double pressure, String[] names,
      double[] moles) {
    System.out.println("fluid=" + label);
    SystemInterface system = new SystemPitzer(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }

    boolean threw = false;
    String message = "";
    try {
      system.chemicalReactionInit();
    } catch (RuntimeException e) {
      threw = true;
      message = e.getMessage() == null ? e.getClass().getName() : e.getMessage();
    }
    System.out.println("gate_threw=" + threw);
    if (threw) {
      System.out.println("gate_message=" + message.replace('\n', ' '));
    }

    // **`system.getComponentNames()` after the attempt is the list the gate judged**, and
    // that is why it is printed for every fluid: the constructor's loop adds the ions its
    // chemistry needs and *re-reads* them, so the refusal above is raised on the augmented
    // pass and not on the feed. A caller who wants to reproduce the verdict must supply
    // these names, not the ones the fluid was built from.
    System.out.println("names=" + String.join(" ", system.getComponentNames()));

    if (!threw) {
      system.createDatabase(true);
      system.setMixingRule("classic");
      System.out.println("survivors="
          + system.getChemicalReactionOperations().getReactionList()
              .getChemicalReactionList().size());
      for (ChemicalReaction reaction : system.getChemicalReactionOperations().getReactionList()
          .getChemicalReactionList()) {
        System.out.println("  survivor[" + reaction.getName() + "]="
            + reaction.getValidationStatus());
      }
    }

    // **The survivor set the gate judged, replayed on the augmented names.** On the fluids
    // the gate passes this is a cross-check of the list above; on the one it refuses it is
    // the only evidence of what was active when it refused.
    ChemicalReactionList replayed = new ChemicalReactionList();
    replayed.readReactions(system);
    replayed.removeJunkReactions(system.getComponentNames());
    replayed.removeDependentReactions();
    ArrayList<ChemicalReaction> rows = replayed.getChemicalReactionList();
    System.out.println("replayed_survivors=" + rows.size());
    for (ChemicalReaction reaction : rows) {
      System.out.println("  replayed[" + reaction.getName() + "]="
          + reaction.getValidationStatus());
    }
    System.out.println("replayed_source=" + replayed.getReactionDataSource().getIdentifier());
    System.out.println();
  }
}
