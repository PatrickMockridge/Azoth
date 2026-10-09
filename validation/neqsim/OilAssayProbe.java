// The oil assay's mass-fraction resolution, as `OilAssayCharacterisation` runs it.
//
// `getResolvedMassFractions` is the one method of the class that answers without touching the
// thermodynamic system: it validates that every cut carries exactly one of a mass or a volume
// fraction, that they share one basis, that they sum to one within `1e-3`, normalises them, and -
// for a volume basis - converts with each cut's specific gravity and renormalises.
//
// The rest of `apply` is system mutation: it resolves a molar mass per cut, adds standard
// components out of the databank, and requires the reconstructed mass to close on the configured
// total within `1e-10`. That is a different id's scope, and the reason is measured rather than
// assumed - see the capture's last rows, which record the two refusals it raises.
//
//   javac -proc:none -cp /path/to/neqsim-f0c7436.jar OilAssayProbe.java
//   java -cp .:/path/to/neqsim-f0c7436.jar OilAssayProbe > captures/oil_assay_probe.tsv

import java.util.Arrays;
import neqsim.thermo.characterization.OilAssayCharacterisation;
import neqsim.thermo.characterization.OilAssayCharacterisation.AssayCut;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class OilAssayProbe {

  private static SystemInterface empty() {
    SystemInterface system = new SystemSrkEos(298.15, 1.0);
    system.setMixingRule(2);
    return system;
  }

  private static String vector(double[] values) {
    StringBuilder line = new StringBuilder();
    for (int i = 0; i < values.length; i++) {
      if (i > 0) {
        line.append(", ");
      }
      line.append(String.format("%.15g", values[i]));
    }
    return line.toString();
  }

  /** One assay, its cuts' declared fractions and gravities, and the mass fractions it resolves. */
  private static void sweep(String label, String basis, OilAssayCharacterisation assay,
      double[] declared) {
    System.out.printf("# assay = %s%n", label);
    System.out.printf("# basis = %s%n", basis);
    System.out.printf("declared_fraction = [%s]%n", vector(declared));
    System.out.printf("cut_count = %d%n", declared.length);
    double[] densities = new double[declared.length];
    for (int i = 0; i < declared.length; i++) {
      densities[i] = assay.getCuts().get(i).resolveDensity() * 1000.0;
    }
    System.out.printf("density = [%s]%n", vector(densities));
    try {
      System.out.printf("mass_fraction = [%s]%n", vector(assay.getResolvedMassFractions()));
      System.out.printf("total_declared_fraction = %.15g%n", Arrays.stream(declared).sum());
      System.out.printf("bulk_specific_gravity = %.15g%n", assay.getBulkSpecificGravity());
    } catch (RuntimeException e) {
      System.out.printf("# refused: %s%n",
          e.getMessage() == null ? e.toString() : e.getMessage().split("\n")[0]);
    }
    System.out.println();
  }

  public static void main(String[] args) {
    // A four-cut TBP assay on a volume basis, which is the shape an atmospheric assay arrives in.
    OilAssayCharacterisation tbp = new OilAssayCharacterisation(empty());
    tbp.addTBPCutBoundariesCelsius("TBP",
        new double[] {0.0, 25.0, 55.0, 80.0, 100.0},
        new double[] {60.0, 120.0, 200.0, 320.0, 420.0},
        new double[] {0.68, 0.78, 0.85, 0.92});
    sweep("tbp_volume_basis", "volume", tbp, new double[] {0.25, 0.30, 0.25, 0.20});

    // The same cuts stated on a mass basis, which skips the gravity conversion entirely.
    OilAssayCharacterisation mass = new OilAssayCharacterisation(empty());
    mass.addCut(new AssayCut("M1").withMassFraction(0.4).withSpecificGravity(0.70)
        .withAverageBoilingPointKelvin(350.0));
    mass.addCut(new AssayCut("M2").withMassFraction(0.35).withSpecificGravity(0.82)
        .withAverageBoilingPointKelvin(450.0));
    mass.addCut(new AssayCut("M3").withMassFraction(0.25).withSpecificGravity(0.90)
        .withAverageBoilingPointKelvin(560.0));
    sweep("mass_basis", "mass", mass, new double[] {0.4, 0.35, 0.25});

    // Fractions that miss one by less than the class's `1e-3` tolerance, and by more: the
    // boundary is inclusive of nothing - `1.001` is outside as much as `1.003` is.
    OilAssayCharacterisation loose = new OilAssayCharacterisation(empty());
    loose.addCut(new AssayCut("L1").withMassFraction(0.5).withSpecificGravity(0.75));
    loose.addCut(new AssayCut("L2").withMassFraction(0.5005).withSpecificGravity(0.85));
    sweep("inside_the_closure_tolerance", "mass", loose, new double[] {0.5, 0.5005});

    OilAssayCharacterisation broken = new OilAssayCharacterisation(empty());
    broken.addCut(new AssayCut("B1").withMassFraction(0.5).withSpecificGravity(0.75));
    broken.addCut(new AssayCut("B2").withMassFraction(0.51).withSpecificGravity(0.85));
    sweep("outside_the_closure_tolerance", "mass", broken, new double[] {0.5, 0.51});

    // One cut with both bases stated, and one pair with neither: the two ways the class refuses a
    // cut that does not say what its fraction is a fraction *of*.
    OilAssayCharacterisation both = new OilAssayCharacterisation(empty());
    both.addCut(new AssayCut("X1").withMassFraction(0.5).withVolumeFraction(0.5).withSpecificGravity(0.8));
    both.addCut(new AssayCut("X2").withMassFraction(0.5).withSpecificGravity(0.8));
    sweep("a_cut_with_both_bases", "mass", both, new double[] {0.5, 0.5});

    OilAssayCharacterisation neither = new OilAssayCharacterisation(empty());
    neither.addCut(new AssayCut("N1").withSpecificGravity(0.8));
    neither.addCut(new AssayCut("N2").withMassFraction(1.0).withSpecificGravity(0.8));
    sweep("a_cut_with_neither_basis", "mass", neither, new double[] {0.0, 1.0});
  }
}
