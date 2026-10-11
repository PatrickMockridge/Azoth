// `OilAssayCharacterisation.apply` - each cut's molar mass and moles, and the fluid they build.
//
// Three routes reach one cut's molar mass, and the class takes them in this order:
//
//   * a **standard component**, whose molar mass comes out of the component database by name;
//   * a **stated molar mass**, used as it is;
//   * otherwise the closed form `5.805e-5 * Tb^2.3776 / SG^0.9371`, on a boiling point the cut
//     states or one resolved from its Watson factor as `(W * SG)^3 / 1.8`, and a gravity it
//     states or one resolved from its API gravity as `141.5 / (API + 131.5)`.
//
// The moles follow as `totalAssayMass * massFraction / molarMass`, a cut with no positive mass
// fraction is skipped, and the class then checks that what it rebuilt is the total it started
// from to `1e-10` relative.
//
//   javac -proc:none -cp /path/to/neqsim-f0c7436.jar OilAssayDerivationProbe.java
//   java -cp .:/path/to/neqsim-f0c7436.jar OilAssayDerivationProbe > captures/oil_assay_derivation_probe.tsv

import neqsim.thermo.characterization.OilAssayCharacterisation;
import neqsim.thermo.characterization.OilAssayCharacterisation.AssayCut;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class OilAssayDerivationProbe {

  /** One assay, and the fluid `apply` builds from it. */
  private static void sweep(String label, double totalAssayMass, AssayCut[] cuts) {
    SystemInterface system = new SystemSrkEos(298.15, 1.0);
    OilAssayCharacterisation assay = new OilAssayCharacterisation(system);
    assay.setTotalAssayMass(totalAssayMass);
    for (AssayCut cut : cuts) {
      assay.addCut(cut);
    }

    System.out.printf("# assay = %s%n", label);
    System.out.printf("cut_count = %d%n", cuts.length);
    System.out.printf("cut_names = [%s]%n", cutNames(cuts));
    System.out.printf("mass_fraction = [%s]%n", vector(assay.getResolvedMassFractions()));

    try {
      assay.apply();
    } catch (Exception ex) {
      System.out.printf("# refused: %s%n", ex);
      return;
    }

    system.init(0);
    System.out.printf("component_names = [%s]%n", componentNames(system));
    System.out.printf("component_molar_mass = [%s]%n", componentMolarMasses(system));
    System.out.printf("component_moles = [%s]%n", componentMoles(system));
    System.out.printf("total_moles = %.15g%n", system.getTotalNumberOfMoles());
    // The class's own gate, read back off the fluid it built rather than off a number it
    // publishes - it keeps the check to itself, so this is the same arithmetic from outside.
    double rebuilt = 0.0;
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      rebuilt += system.getPhase(0).getComponent(i).getNumberOfmoles()
          * system.getPhase(0).getComponent(i).getMolarMass();
    }
    System.out.printf("mass_balance_error = %.15g%n", Math.abs(rebuilt - totalAssayMass) / totalAssayMass);
  }

  private static String cutNames(AssayCut[] cuts) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < cuts.length; i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(cuts[i].getName());
    }
    return out.toString();
  }

  private static String componentNames(SystemInterface system) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(system.getPhase(0).getComponent(i).getComponentName());
    }
    return out.toString();
  }

  private static String componentMolarMasses(SystemInterface system) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(String.format("%.15g", system.getPhase(0).getComponent(i).getMolarMass()));
    }
    return out.toString();
  }

  private static String componentMoles(SystemInterface system) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(String.format("%.15g", system.getPhase(0).getComponent(i).getNumberOfmoles()));
    }
    return out.toString();
  }

  private static String vector(double[] values) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < values.length; i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(String.format("%.15g", values[i]));
    }
    return out.toString();
  }

  public static void main(String[] args) {
    // **All three routes in one mass-basis assay**: a gravity and a boiling point the cut
    // states, a standard component read out of the database by name, and a stated molar mass.
    sweep(
        "a_mass_basis_slate",
        1000.0,
        new AssayCut[] {
          new AssayCut("cut_a").withMassFraction(0.4).withDensity(0.78)
              .withAverageBoilingPointKelvin(400.0),
          new AssayCut("cut_b").withMassFraction(0.3).withStandardComponent("n-heptane"),
          new AssayCut("cut_c").withMassFraction(0.3).withDensity(0.85).withMolarMass(0.25),
        });

    // A volume basis, so each cut's density comes from its own API gravity rather than from a
    // stated gravity, and the volume fractions are converted before anything else happens.
    sweep(
        "a_volume_basis_slate",
        500.0,
        new AssayCut[] {
          new AssayCut("light").withVolumeFraction(0.6).withApiGravity(55.0)
              .withAverageBoilingPointKelvin(340.0),
          new AssayCut("heavy").withVolumeFraction(0.4).withApiGravity(20.0)
              .withAverageBoilingPointKelvin(560.0),
        });

    // A Watson factor and no boiling point, which is the resolution `(W * SG)^3 / 1.8`.
    sweep(
        "a_watson_factor_cut",
        250.0,
        new AssayCut[] {
          new AssayCut("w1").withMassFraction(0.5).withDensity(0.80)
              .withWatsonCharacterizationFactor(11.5),
          new AssayCut("w2").withMassFraction(0.5).withDensity(0.88).withMolarMass(0.30),
        });
  }
}
