// The multi-salt precipitation's layers, printed so each can be ported and checked.
//
// `MultiSaltPrecipitation.solve()` is a **complementarity loop**: it holds one solid amount per
// mineral and repeatedly takes the mineral whose `|SR - 1|` is largest, precipitating where the
// ratio is over one and dissolving where it is under, until every mineral is at saturation.
// `CalcSaltSatauration` is the per-mineral half - `precipitate()` and `dissolve()` - and the
// loop refuses rather than returns when it stalls or when its ledger does not close.
//
//   javac -proc:none -cp neqsim-f0c7436.jar MultiSaltProbe.java
//   java -cp .:neqsim-f0c7436.jar MultiSaltProbe > captures/multi_salt_pitzer_probe.tsv
//   java -cp .:neqsim-f0c7436.jar MultiSaltProbe three-salts
//
// **`three-salts` refuses on every revision pinned so far**, so it is run to see the refusal
// and not to reproduce a capture. `multi_salt_probe.tsv` is the same brine from before
// `PitzerParameterCoverage` existed and is the loop's only oracle; the header above says what
// it is and why it cannot be rebuilt.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPitzer;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.flashops.saturationops.MultiSaltPrecipitationResult;
import neqsim.thermodynamicoperations.flashops.saturationops.SaltPrecipitationResult;

public class MultiSaltProbe {

  /** The halite brine, whose one mineral the loop settles in a single round. */
  private static final String[] HALITE = {"water", "Na+", "Cl-", "CO3--", "HCO3-"};

  /** The brine three minerals can form out of, so the complementarity loop takes three. */
  private static final String[] THREE_SALTS = {
      "water", "Na+", "Cl-", "Ca++", "SO4--", "CO3--", "HCO3-"
  };

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    // `precipitateScales`, and not the older `calcMultiSaltPrecipitation`: the operation was
    // renamed, and the probe's first compile against the pinned jar is what said so.
    if (args.length > 0 && args[0].equals("three-salts")) {
      // **The complementarity loop's own state, and the only one that exercises it.** Three
      // minerals take three rounds (`updates = 3`), so this is the oracle for a loop azoth
      // does not carry yet; the halite brine below settles in one and would not show it.
      //
      // **NeqSim refuses this brine and so does azoth, for the same reason**: `Ca++|Na+` and
      // `Cl-|SO4--` carry no theta in the legacy Pitzer dataset and both phases refuse a set
      // whose parameters are missing rather than running it with zeros. So the run ends in
      // the refusal, which is itself the measurement - and the reason `multi_salt_probe.tsv`
      // cannot be rebuilt from a pinned revision.
      report(25.0, 10.0, THREE_SALTS, new String[] {"CaCO3", "CaSO4_A", "NaCl"},
          new double[] {1.0, 0.05, 0.05, 0.01, 0.005, 0.004, 0.002});
      return;
    }
    // **A brine at halite's saturation**, because a mineral under it precipitates nothing
    // and exercises none of the solve: `NaCl`'s product is about `38` and its molality at
    // saturation about `6`, so the run that matters is the one that crosses that.
    for (double temperatureC : new double[] {25.0, 60.0, 90.0}) {
      report(temperatureC, 10.0, HALITE, new String[] {"NaCl"},
          new double[] {1.0, 0.12, 0.12, 0.004, 0.002});
    }
  }

  /** The brine the last two arguments describe, ending at the minerals it can form. */
  private static void report(double temperatureC, double pressureBara, String[] components,
      String[] minerals, double[] moles) {
    System.out.printf("# brine at %.15g C, %.15g bara%n", temperatureC, pressureBara);
    StringBuilder names = new StringBuilder();
    for (String mineral : minerals) {
      names.append(" ").append(mineral);
    }
    System.out.printf("# minerals:%s%n", names);
    try {
      // **A Pitzer brine, and not a cubic one.** The operation reads
      // `getActivityCoefficient` off the phase it is pointed at, so the numbers a port is
      // compared on are the *phase model's* - and `SystemSrkEos` answers with the cubic
      // phase's, which is a different model from the electrolyte one azoth carries.
      SystemInterface fluid = new SystemPitzer(273.15 + temperatureC, pressureBara);
      for (int i = 0; i < components.length; i++) {
        fluid.addComponent(components[i], moles[i]);
      }
      fluid.setMixingRule(2);
      fluid.setMultiPhaseCheck(true);
      fluid.init(0);
      fluid.init(1);

      // **The activity coefficients the ratios are built from**, and the phase's own mole
      // fractions: a port that read different gammas would produce a plausible answer for a
      // different brine, so they are the model's inputs and are captured as such.
      for (int phase = 0; phase < fluid.getNumberOfPhases(); phase++) {
        if (!fluid.getPhase(phase).hasComponent("water")) {
          continue;
        }
        System.out.printf("phase[%d] = %s%n", phase, fluid.getPhase(phase).getType());
        for (int i = 0; i < fluid.getPhase(phase).getNumberOfComponents(); i++) {
          String name = fluid.getPhase(phase).getComponent(i).getName();
          row("x[" + name + "]", fluid.getPhase(phase).getComponent(i).getx());
          row("gamma[" + name + "]", fluid.getPhase(phase).getActivityCoefficient(i, 0));
        }
      }

      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      MultiSaltPrecipitationResult result = operations.precipitateScales(minerals);

      row("updates", result.getEquilibriumUpdates());
      row("final_violation", result.getMaximumComplementarityViolation());
      row("maximum_absolute_residual_moles", result.getMaximumComponentBalanceResidualMoles());
      row("maximum_normalized_residual", result.getMaximumNormalizedBalanceResidual());
      for (String mineral : minerals) {
        SaltPrecipitationResult mineralResult = result.getMineralResult(mineral);
        row("initial_sr[" + mineral + "]", mineralResult.getInitialSaturationRatio());
        row("final_sr[" + mineral + "]", mineralResult.getFinalSaturationRatio());
        row("precipitated_moles[" + mineral + "]", mineralResult.getPrecipitatedMoles());
        row("precipitated_grams[" + mineral + "]", mineralResult.getPrecipitatedMassGrams());
        row("complementarity[" + mineral + "]", mineralResult.getComplementarityViolation());
        row("balance_residual[" + mineral + "]", mineralResult.getMaximumIonBalanceResidualMoles());
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
