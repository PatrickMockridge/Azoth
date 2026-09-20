// What `PhasePitzer.getPitzerParameterCoverage` reports, and what it refuses.
//
//     javac -proc:none -cp neqsim-3.20.0.jar PitzerCoverageProbe.java
//     java -cp .:neqsim-3.20.0.jar PitzerCoverageProbe
//
// The audit answers the question the *selection* rule does not: the selection rule asks
// whether the PHREEQC catalogue covers the topology and falls back if not, while this asks
// whether whichever dataset loaded covers it. The two can disagree - the legacy CSV is
// chosen precisely because the catalogue was incomplete, and the audit then reports what
// the CSV is missing in turn.
//
// Four things the port needs that are not in the class's javadoc:
//
//   1. **What the defined sets are**, per dataset. `definedBinaryPairs`, `definedThetaPairs`
//      and `definedPsiTuples` are filled as parameters are added, so what counts as
//      "defined" is the loaded dataset's row set, not the audit's opinion.
//   2. **Which of the three families goes missing for which topology** - the theta/psi
//      pairs are same-sign, so a single-salt brine has neither and a mixed brine has both.
//   3. **That an ion below `ACTIVE_ION_MOLALITY` leaves the topology**, so a trace
//      component does not make the audit refuse.
//   4. **What the refusal says**, verbatim, since the port's own refusal must not read
//      differently from NeqSim's when the same topology is audited.
//
// The reaction-species variant (`getPitzerReactionParameterCoverage`) is a separate
// observable: it does not exclude `H3O+`, `OH-`, `HCO3-` and `CO3--`, so the same brine can
// be complete under one and not the other.

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhasePitzer;
import neqsim.thermo.phase.PitzerParameterCoverage;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPitzer;

public class PitzerCoverageProbe {

  private static SystemInterface build(String[] names, double[] z) {
    SystemInterface s = new SystemPitzer(298.15, 1.0);
    for (int i = 0; i < names.length; i++) {
      s.addComponent(names[i], z[i]);
    }
    s.setMixingRule("classic");
    s.init(0);
    s.init(1);
    for (int i = 0; i < s.getNumberOfPhases(); i++) {
      if (s.getPhase(i) instanceof PhasePitzer) {
        return s;
      }
    }
    throw new IllegalStateException("no PhasePitzer in " + java.util.Arrays.toString(names));
  }

  private static PhasePitzer pitz(SystemInterface s) {
    for (int i = 0; i < s.getNumberOfPhases(); i++) {
      if (s.getPhase(i) instanceof PhasePitzer) {
        return (PhasePitzer) s.getPhase(i);
      }
    }
    throw new IllegalStateException("no PhasePitzer");
  }

  private static void audit(String label, String[] names, double[] z) {
    SystemInterface s;
    try {
      s = build(names, z);
    } catch (IllegalStateException e) {
      System.out.printf("%n=== %s ===%n", label);
      System.out.printf("  init()         = REFUSED: %s%n", e.getMessage());
      return;
    }
    PhasePitzer p = pitz(s);
    PitzerParameterCoverage c = p.getPitzerParameterCoverage();
    System.out.printf("%n=== %s ===%n", label);
    System.out.printf("  dataset        = %s%n", c.getDatasetId());
    System.out.printf("  activeCations  = %s%n", c.getActiveCations());
    System.out.printf("  activeAnions   = %s%n", c.getActiveAnions());
    System.out.printf("  missingBinary  = %s%n", c.getMissingBinaryPairs());
    System.out.printf("  missingTheta   = %s%n", c.getMissingThetaPairs());
    System.out.printf("  missingPsi     = %s%n", c.getMissingPsiTuples());
    System.out.printf("  isComplete     = %s%n", c.isComplete());
    try {
      p.requireCompletePitzerParameterCoverage();
      System.out.printf("  require()      = accepted%n");
    } catch (IllegalStateException e) {
      System.out.printf("  require()      = REFUSED: %s%n", e.getMessage());
    }
    PitzerParameterCoverage reaction = p.getPitzerReactionParameterCoverage();
    System.out.printf("  reaction isComplete = %s%n", reaction.isComplete());
    if (!reaction.isComplete()) {
      System.out.printf("  reaction missingBinary = %s%n", reaction.getMissingBinaryPairs());
      System.out.printf("  reaction missingTheta  = %s%n", reaction.getMissingThetaPairs());
      System.out.printf("  reaction missingPsi    = %s%n", reaction.getMissingPsiTuples());
    }
  }

  public static void main(String[] args) {
    // A single salt: one cation, one anion, so no same-sign pair exists and only the
    // binary family can go missing.
    audit("water + Na+ + Cl-", new String[] {"water", "Na+", "Cl-"},
        new double[] {0.98, 0.01, 0.01});

    // Two cations and one anion: theta(Na+,K+) and psi(Na+,K+,Cl-) are required.
    audit("water + Na+ + K+ + Cl-", new String[] {"water", "Na+", "K+", "Cl-"},
        new double[] {0.97, 0.01, 0.01, 0.01});

    // Two anions and one cation: the psi tuple is (anion, anion, cation).
    audit("water + Na+ + Cl- + SO4--", new String[] {"water", "Na+", "Cl-", "SO4--"},
        new double[] {0.97, 0.01, 0.01, 0.01});

    // Both directions at once.
    audit("water + Na+ + K+ + Cl- + SO4--",
        new String[] {"water", "Na+", "K+", "Cl-", "SO4--"},
        new double[] {0.95, 0.01, 0.01, 0.01, 0.01});

    // Divalent against monovalent, which the catalogue covers from its own tests.
    audit("water + Na+ + Ca++ + Cl-", new String[] {"water", "Na+", "Ca++", "Cl-"},
        new double[] {0.97, 0.01, 0.01, 0.01});

    // A topology that falls back to the CSV, so the audit runs against a dataset with
    // only cation-anion rows: `HCO3-` is in no catalogue family.
    audit("water + Na+ + HCO3-", new String[] {"water", "Na+", "HCO3-"},
        new double[] {0.98, 0.01, 0.01});

    // The same ions under the reaction-species variant, which does not exclude them.
    audit("water + Na+ + HCO3- + CO3--",
        new String[] {"water", "Na+", "HCO3-", "CO3--"},
        new double[] {0.97, 0.01, 0.01, 0.01});

    // A mixed brine forced onto the legacy CSV, which carries no same-sign rows at all.
    audit("water + Na+ + K+ + Cl- + HCO3-",
        new String[] {"water", "Na+", "K+", "Cl-", "HCO3-"},
        new double[] {0.95, 0.01, 0.01, 0.01, 0.01});

    // An ion with no row in either dataset: `li+` is a component and appears in neither.
    audit("water + Na+ + Li+ + Cl-", new String[] {"water", "Na+", "Li+", "Cl-"},
        new double[] {0.97, 0.01, 0.01, 0.01});

    // And an anion the same way.
    audit("water + Na+ + Cl- + Br-", new String[] {"water", "Na+", "Cl-", "Br-"},
        new double[] {0.97, 0.01, 0.01, 0.01});

    // A single salt with no row anywhere.
    audit("water + NH4+ + Cl-", new String[] {"water", "NH4+", "Cl-"},
        new double[] {0.98, 0.01, 0.01});

    // A trace ion: the audit's own threshold is 1e-8 mol/kg, and this is below it.
    // x = 1e-11 with this much water is a molality of 5.7e-10.
    audit("water + Na+ + Cl- + trace K+ (molality 5.7e-10)",
        new String[] {"water", "Na+", "Cl-", "K+"},
        new double[] {0.979_999_999_99, 0.01, 0.01, 1.0e-11});
  }
}
