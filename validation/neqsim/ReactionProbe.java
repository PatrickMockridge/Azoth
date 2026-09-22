// The reaction tables' equilibrium constant, its temperature derivative and the rate
// factor, for `reactions.equilibrium_constant`.
//
// The reaction is read through `ChemicalReactionFactory.getChemicalReaction`, which is
// NeqSim's own public lookup: it queries `reactiondata` by name, builds the
// `ChemicalReaction` from the row and its `stoccoefdata` coefficients, and hands it back.
// So the arithmetic below is the class's and not a re-implementation of the correlation.
//
// `getK` and `getReactionHeat` take a `PhaseInterface` and read only its temperature, so
// the phase here is an empty `PhaseSrkEos` whose temperature is set per state. Attaching a
// fluid would check a flash as well and this has no flash in it.
//
// **`getKineticRateLaw()` is printed because the answer turns on it.** A reaction built
// from the table is put on `LEGACY_TEMPERATURE_CORRELATION` by
// `ChemicalReactionList.readReactions`, whose comment is "preserve database kinetic
// behavior until each rate's units and provenance are qualified" - so the stored
// `ACTENERGY` is not evaluated, and the legacy law ignores it. `methanecombustion` is
// included to show the dormant-flag path as well: its `usereaction` is 0.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ReactionProbe.java
//     java -cp .:neqsim-f0c7436.jar ReactionProbe > captures/reaction_probe.tsv

import neqsim.chemicalreactions.chemicalreaction.ChemicalReaction;
import neqsim.chemicalreactions.chemicalreaction.ChemicalReactionFactory;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseSrkEos;

public class ReactionProbe {
  static final String[] REACTIONS = {
      "CO2water", "waterreac", "carbonate", "water-H2S", "MDEAprot", "DEAprot", "methanecombustion",
  };

  static final double[] TEMPERATURES = { 273.15, 298.15, 323.15, 373.15, 423.15 };

  public static void main(String[] args) {
    PhaseInterface phase = new PhaseSrkEos();

    for (String name : REACTIONS) {
      ChemicalReaction reaction = ChemicalReactionFactory.getChemicalReaction(name);
      if (reaction == null) {
        System.out.println("reaction=" + name + " absent");
        continue;
      }

      double[] k = reaction.getK();
      System.out.println("reaction=" + name);
      System.out.println("  k1=" + k[0] + " k2=" + k[1] + " k3=" + k[2] + " k4=" + k[3]);
      System.out.println("  rate_law=" + reaction.getKineticRateLaw());
      System.out.println("  stored_rate=" + reaction.getRateFactor() + " stored_activation_energy="
          + reaction.getActivationEnergy() + " reference_temperature="
          + reaction.getReferenceTemperature());

      for (double t : TEMPERATURES) {
        phase.setTemperature(t);
        double value = reaction.getK(phase);
        System.out.println("  T=" + t
            + " lnK=" + Math.log(value)
            + " K=" + value
            + " reaction_heat=" + reaction.getReactionHeat(phase)
            + " rate_factor=" + reaction.getRateFactor(phase));
      }
    }
  }
}
