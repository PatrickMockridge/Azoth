// What a built `SystemPitzer` actually reads, before anything is ported.
//
//     javac -proc:none -cp neqsim-f0c7436.jar PitzerProbe.java
//     java -cp .:neqsim-f0c7436.jar PitzerProbe
//
// The port needs to know three things this prints, and none of them is in the class's
// javadoc:
//
//   1. **Which parameter dataset is in force.** `PhasePitzer.loadParametersFromDatabase`
//      first tries `PitzerParameterDatasets.tryApplyCompletePhreeqcPitzerCatalog`, and
//      returns if it succeeds; only then does it read the `pitzerparameters` table, which
//      is `PitzerParameters.csv` - the 30 rows azoth vendors. `getParameterDatasetId()`
//      names the one that won.
//   2. **What the parameters are**, through the same getters the model's own arithmetic
//      uses, at a state rather than from a table.
//   3. **Which phase slot is Pitzer's.** `SystemPitzer` is a hybrid system, so the phase
//      array is an EoS gas, the GE aqueous phase, and an EoS oil.
//
// What it measured, 2026-09-20:
//
//   * **Both datasets are live, and the topology chooses.** `water + Na+ + Cl-` and
//     `water + Na+ + Ca++ + Cl-` apply the PHREEQC catalogue; `water + CO2` and
//     `water + Na+ + HCO3-` fall back to `neqsim-legacy-pitzer-parameters-v1`, which is
//     `PitzerParameters.csv` - the fallback happens when a required row is absent from
//     the catalogue, so it is per-mixture rather than per-build.
//   * **The two disagree on a pair they share.** Na+/Cl- is `0.07534 / 0.2769 / 0.00148`
//     under PHREEQC against the CSV's `0.0765 / 0.2664 / 0.00127` - Cphi by 16%.
//   * **The catalogue carries a `B2` the CSV does not**: `Ca++/Cl-` is `-1.13` there and
//     zero in the CSV, whose four `beta2` rows are the divalent sulphates.
//   * **The `LAMBDA`/`ZETA` families reach a neutral solute** - `water + CO2` is the
//     topology where the catalogue's neutral interactions matter.
//   * The molality is the textbook `n_i / m_water`: `Na+` at `x = 0.1` reports `6.1677`
//     mol/kg, which is `0.1 / (0.9 * 0.018015)`.
//
// The PHREEQC catalogue is
// `src/main/resources/neqsim/thermo/phase/phreeqc/pitzer-<commit>.dat` - 280 lines,
// eleven parameter families, up to six temperature coefficients per row, and PHREEQC
// species names (`Ba+2`, `B(OH)4-`) rather than NeqSim's (`Ba++`).

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhasePitzer;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPitzer;

public class PitzerProbe {

  private static SystemInterface build(double t, double pBar, String[] names, double[] z) {
    SystemInterface s = new SystemPitzer(t, pBar);
    for (int i = 0; i < names.length; i++) {
      s.addComponent(names[i], z[i]);
    }
    s.setMixingRule("classic");
    return s;
  }

  private static void report(String label, SystemInterface s) {
    System.out.printf("%n=== %s ===%n", label);
    s.init(0);
    s.init(1);
    System.out.printf("phases = %d%n", s.getNumberOfPhases());
    for (int i = 0; i < s.getNumberOfPhases(); i++) {
      PhaseInterface phase = s.getPhase(i);
      System.out.printf("  [%d] %-12s %s%n", i, phase.getType(),
          phase.getClass().getSimpleName());
      if (phase instanceof PhasePitzer) {
        PhasePitzer pitz = (PhasePitzer) phase;
        System.out.printf("       dataset id     = %s%n", pitz.getParameterDatasetId());
        System.out.printf("       common-ion     = %s%n", pitz.isPhreeqcCommonIonTermsActive());
        System.out.printf("       ionic strength = %.10g mol/kg%n", pitz.getIonicStrength());
        System.out.printf("       solvent weight = %.10g kg%n", pitz.getSolventWeight());
        System.out.printf("       osmotic (water) = %.10g%n", pitz.getOsmoticCoefficientOfWater());
        System.out.printf("       osmotic (molal) = %.10g%n",
            pitz.getOsmoticCoefficientOfWaterMolality());
        System.out.printf("       density         = %.10g kg/m3%n", pitz.getDensity());
        for (int a = 0; a < phase.getNumberOfComponents(); a++) {
          for (int b = a + 1; b < phase.getNumberOfComponents(); b++) {
            System.out.printf("       b0(%s,%s)  = %.10g   b1 = %.10g   cphi = %.10g   b2 = %.10g%n",
                phase.getComponent(a).getComponentName(),
                phase.getComponent(b).getComponentName(),
                pitz.getBeta0ij(a, b, s.getTemperature()),
                pitz.getBeta1ij(a, b, s.getTemperature()),
                pitz.getCphiij(a, b, s.getTemperature()),
                pitz.getBeta2ij(a, b, s.getTemperature()));
          }
        }
      }
    }
    for (int i = 0; i < s.getNumberOfComponents(); i++) {
      System.out.printf("  x[%d] %-10s = %.10g  charge = %+.0f  molality = %.10g%n", i,
          s.getPhase(1).getComponent(i).getComponentName(),
          s.getPhase(1).getComponent(i).getx(),
          s.getPhase(1).getComponent(i).getIonicCharge(),
          s.getPhase(1).getComponent(i).getMolality(s.getPhase(1)));
    }
  }

  public static void main(String[] args) {
    // The simplest brine the model is fitted for.
    report("water + NaCl",
        build(298.15, 1.0, new String[] {"water", "Na+", "Cl-"}, new double[] {0.9, 0.1, 0.1}));

    // A pair the 30-row CSV fits and the PHREEQC catalogue also carries, so a divergence
    // between the two datasets would show.
    report("water + NaCl + CaCl2",
        build(298.15, 1.0, new String[] {"water", "Na+", "Ca++", "Cl-", "Cl-"},
            new double[] {0.85, 0.05, 0.02, 0.08, 0.0}));

    // A neutral solute, which is where `LAMBDA`/`ZETA` and the Henry reference live.
    report("water + CO2",
        build(298.15, 1.0, new String[] {"water", "CO2"}, new double[] {0.99, 0.01}));

    // And one the CSV fits but PHREEQC may not: a carbonate system.
    report("water + Na+ + HCO3-",
        build(298.15, 1.0, new String[] {"water", "Na+", "HCO3-"}, new double[] {0.9, 0.05, 0.05}));
  }
}
