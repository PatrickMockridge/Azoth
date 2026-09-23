// The kinetic rate law and the mass-transfer rate matrix, for `reactions.chemical_kinetics`.
//
// `ChemicalReaction.getRateFactor(phase)` is **two laws behind one selector**, and the plan's
// note on it recorded only the first:
//
// * `LEGACY_TEMPERATURE_CORRELATION`, which is what an object with no selector gets:
//   `2.576e9 * exp(-6024 / T) / 1000` - a literal, the same for every reaction at a given
//   temperature, which is why `ACTENERGY` is never read in this branch;
// * `REFERENCE_ARRHENIUS`, which needs a reference rate, an activation energy and a reference
//   temperature and is the branch the class documents as the migration target:
//   `rateFactor * exp(-Ea/R * (1/T - 1/refT))`.
//
// Both are printed, at two temperatures, beside the reaction's own names and stoichiometric
// coefficients - because the law is a property of the *reaction* and the numbers have to be
// read with it.
//
// The mass-transfer matrix that consumes the rate factor is `Kinetics.calcReacMatrix`, whose
// inputs are the phase pair and an effective-diffusion vector; those are printed for the
// aqueous phase of a chemical system.
//
//     javac -proc:none -cp neqsim-f0c7436.jar KineticsProbe.java
//     java -cp .:neqsim-f0c7436.jar KineticsProbe > captures/kinetics_probe.tsv

import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.chemicalreactions.chemicalreaction.ChemicalReaction;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class KineticsProbe {

  public static void main(String[] args) {
    rateLaw("co2-water-298K", 298.15, new String[] { "CO2", "water" }, new double[] { 0.01, 10.0 });
    rateMatrix("co2-water-298K", 298.15, new String[] { "CO2", "water" },
        new double[] { 0.01, 10.0 });
    freshKinetics("co2-water-298K", 298.15, new String[] { "CO2", "water" },
        new double[] { 0.01, 10.0 });
  }

  /// **`phiInfinite` and `isIrreversible` are fields of the object, and neither is reset.**
  ///
  /// `Kinetics.calcReacMatrix` assigns `phiInfinite` per sibling and sets `isIrreversible`
  /// once any reaction's scaled `1/K` comes in under `1e-3`; both are read back off the
  /// object afterwards, so on a shared instance the value a caller sees is the *last* one
  /// any call produced. The port answers per call, so the two are only comparable against
  /// a fresh object - which is what this block is for.
  ///
  /// The shared object's values are printed first, beside them, because a reader who takes
  /// `operations.getKinetics()` gets those and the difference is the point.
  static void freshKinetics(String label, double temperature, String[] names, double[] moles) {
    System.out.println("fluid=" + label + "-fresh-kinetics");
    SystemInterface system = build(temperature, names, moles);
    ThermodynamicOperations flash = new ThermodynamicOperations(system);
    flash.TPflash();
    system.init(1);
    PhaseInterface aqueous = system.getPhase(1);
    ChemicalReactionOperations operations = system.getChemicalReactionOperations();

    // The shared instance, after every row has been taken from it.
    for (int comp = 0; comp < aqueous.getNumberOfComponents(); comp++) {
      operations.solveKinetics(1, aqueous, comp);
    }
    for (int comp = 0; comp < aqueous.getNumberOfComponents(); comp++) {
      String name = aqueous.getComponent(comp).getName();
      System.out.println("shared_phi_infinite[" + name + "]="
          + operations.getKinetics().getPhiInfinite());
      System.out.println("shared_irreversible[" + name + "]="
          + operations.getKinetics().isIrreversible());
    }

    // A fresh object per component: the per-call answer.
    for (int comp = 0; comp < aqueous.getNumberOfComponents(); comp++) {
      String name = aqueous.getComponent(comp).getName();
      neqsim.chemicalreactions.kinetics.Kinetics fresh =
          new neqsim.chemicalreactions.kinetics.Kinetics(operations);
      fresh.calcReacMatrix(aqueous, aqueous, comp);
      System.out.println("fresh_phi_infinite[" + name + "]=" + fresh.getPhiInfinite());
      System.out.println("fresh_irreversible[" + name + "]=" + fresh.isIrreversible());
    }
    System.out.println();
  }

  /// **The two rate laws**, per reaction, at two temperatures.
  static void rateLaw(String label, double temperature, String[] names, double[] moles) {
    System.out.println("fluid=" + label);
    SystemInterface system = build(temperature, names, moles);
    system.init(1);
    PhaseInterface phase = system.getPhase(1);

    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    System.out.println("reactions=" + operations.getReactionList().getChemicalReactionList().size());

    int index = 0;
    for (ChemicalReaction reaction : operations.getReactionList().getChemicalReactionList()) {
      System.out.println("  reaction[" + index + "]=" + String.join(" ",
          reaction.getNames()));
      StringBuilder coefficients = new StringBuilder("    stoc_coefs=");
      for (double coefficient : reaction.getStocCoefs()) {
        coefficients.append(coefficient).append(" ");
      }
      System.out.println(coefficients.toString().trim());
      System.out.println("    rate_law=" + reaction.getKineticRateLaw());
      System.out.println("    stored_rate_factor=" + reaction.getRateFactor());
      for (double t : new double[] { 298.15, 373.15 }) {
        phase.setTemperature(t);
        System.out.println("    legacy_rate_factor_at_" + t + "=" + reaction.getRateFactor(phase));
      }
      // The Arrhenius branch, on the same reaction, for the comparison the selector exists for.
      reaction.setReferenceKinetics(1.0e-3, 50_000.0, 298.15);
      for (double t : new double[] { 298.15, 373.15 }) {
        phase.setTemperature(t);
        System.out.println("    arrhenius_rate_factor_at_" + t + "=" + reaction.getRateFactor(phase));
      }
      phase.setTemperature(temperature);
      index++;
    }
    System.out.println();
  }

  /// The rate matrix and the diffusion vector it consumes.
  static void rateMatrix(String label, double temperature, String[] names, double[] moles) {
    System.out.println("fluid=" + label + "-kinetics");
    SystemInterface system = build(temperature, names, moles);
    ThermodynamicOperations flash = new ThermodynamicOperations(system);
    flash.TPflash();
    system.init(1);

    System.out.println("phases=" + system.getNumberOfPhases());
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      System.out.println("  phase[" + phase + "]=" + system.getPhase(phase).getPhaseTypeName());
    }

    PhaseInterface aqueous = system.getPhase(1);
    aqueous.getPhysicalProperties().calcEffectiveDiffusionCoefficients();
    StringBuilder effective = new StringBuilder("effective_diffusion=");
    for (int i = 0; i < aqueous.getNumberOfComponents(); i++) {
      effective.append(aqueous.getPhysicalProperties().getEffectiveDiffusionCoefficient(i)).append(" ");
    }
    System.out.println(effective.toString().trim());

    // **The binary matrix the effective vector is assembled from**, by index: `D_eff_i` is
    // `(1 - x_i) / sum_{j != i} x_j / D_ij`, so the vector is not reproducible without it.
    // **The binary matrix, and the reason the effective vector cannot be reproduced from it.**
    // `getFickDiffusionCoefficient` is declared on `PhysicalProperties` and not on the interface
    // the phase hands back - so a caller needs the concrete cast - and what it returns is
    // **diagonal**: only the `i == i` entries are non-zero, so the assembly's
    // `sum_{j != i} x_j / D_ij` divides by zero and the recomputation below comes out empty.
    // The array `calcEffectiveDiffusionCoefficients` consumed is the model's own private one, and
    // no public getter returns it.
    neqsim.physicalproperties.system.PhysicalProperties properties =
        (neqsim.physicalproperties.system.PhysicalProperties) aqueous.getPhysicalProperties();
    StringBuilder binaryLine = new StringBuilder("binary_diffusion=");
    for (int i = 0; i < aqueous.getNumberOfComponents(); i++) {
      for (int j = 0; j < aqueous.getNumberOfComponents(); j++) {
        binaryLine.append(properties.getFickDiffusionCoefficient(i, j)).append(" ");
      }
    }
    System.out.println(binaryLine.toString().trim());

    // **The assembly, recomputed from that matrix**, beside the vector the class holds: if the
    // two agree then the matrix printed above is the one `calcEffectiveDiffusionCoefficients`
    // consumed, and the formula is `(1 - x_i) / sum_{j != i} x_j / D_ij`.
    StringBuilder recomputed = new StringBuilder("effective_recomputed=");
    for (int i = 0; i < aqueous.getNumberOfComponents(); i++) {
      double sum = 0.0;
      for (int j = 0; j < aqueous.getNumberOfComponents(); j++) {
        if (i != j) {
          sum += aqueous.getComponent(j).getx() / properties.getFickDiffusionCoefficient(i, j);
        }
      }
      recomputed.append((1.0 - aqueous.getComponent(i).getx()) / sum).append(" ");
    }
    System.out.println(recomputed.toString().trim());

    // **The context the matrix is built from**: the phases' densities, each component's mole
    // fraction and molar mass, and each reaction's equilibrium constant at the state. Without
    // them the matrix is a number with nothing behind it.
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      PhaseInterface held = system.getPhase(phase);
      System.out.println("density[" + phase + "]=" + held.getPhysicalProperties().getDensity());
      for (int i = 0; i < held.getNumberOfComponents(); i++) {
        System.out.println("  component[" + phase + "][" + held.getComponent(i).getName()
            + "]_x=" + held.getComponent(i).getx()
            + " molar_mass=" + held.getComponent(i).getMolarMass());
      }
    }
    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    for (ChemicalReaction reaction : operations.getReactionList().getChemicalReactionList()) {
      System.out.println("reaction_k[" + String.join("_", reaction.getNames()) + "]="
          + reaction.getK(aqueous));
      System.out.println("reaction_rate_factor[" + String.join("_", reaction.getNames()) + "]="
          + reaction.getRateFactor(aqueous));
    }

    for (int comp = 0; comp < aqueous.getNumberOfComponents(); comp++) {
      String name = aqueous.getComponent(comp).getName();
      double matrix = operations.solveKinetics(1, aqueous, comp);
      System.out.println("reac_matrix[" + name + "]=" + matrix);
      System.out.println("  phi_infinite=" + operations.getKinetics().getPhiInfinite());
    }
    System.out.println();
  }

  static SystemInterface build(double temperature, String[] names, double[] moles) {
    SystemInterface system = new SystemSrkEos(temperature, 1.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule(2);
    system.init(0);
    return system;
  }
}
