// Every TBP model NeqSim offers, printed over the same cuts, so they can be ported.
//
// `TbpFractionProbe` captures `PedersenSRK` alone - the default `Characterise` selects. This
// probe takes `TBPfractionModel.getAvailableModels()` and runs the same twelve cuts through
// each of the ten, so the port has every model's equations *and* the values they give.
//
// The cuts cross every branch a model has: `calcTB` switches at MW 540, and the Pedersen
// coefficient set at MW 1120, `RiaziDaubert` switches to its parent above MW 300, and
// `CavettModel`'s API correction applies only below API 30 (SG above ~0.876). A port that
// picked the wrong side of any of them would still be right at a single cut.
//
//   javac -proc:none -cp neqsim-f0c7436.jar CharacterizationProbe.java
//   java -cp .:neqsim-f0c7436.jar CharacterizationProbe > captures/characterization_probe.tsv

import neqsim.thermo.characterization.TBPfractionModel;
import neqsim.thermo.characterization.TBPModelInterface;

public class CharacterizationProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  /** One model's properties over every cut. */
  private static void sweep(String name, double[][] cuts) {
    TBPModelInterface model = new TBPfractionModel().getModel(name);
    System.out.printf("# model = %s%n", name);
    for (int i = 0; i < cuts.length; i++) {
      double molarMass = cuts[i][0];
      double density = cuts[i][1];
      System.out.printf("# cut[%d] MW = %.15g g/mol, density = %.15g g/cm3%n", i, molarMass, density);
      row("cut[" + i + "].tc", model.calcTC(molarMass, density));
      row("cut[" + i + "].pc", model.calcPC(molarMass, density));
      row("cut[" + i + "].tb", model.calcTB(molarMass, density));
      row("cut[" + i + "].acentric", model.calcAcentricFactor(molarMass, density));
      // `calcm` throws for the five models that do not use the cubic's alpha exponent, which is
      // a statement the port owes rather than a gap in this probe: `isCalcm()` is what the caller
      // reads, and it is false for exactly those.
      if (model.isCalcm()) {
        row("cut[" + i + "].calcm", model.calcm(molarMass, density));
      } else {
        System.out.printf("cut[%d].calcm = # not applicable (isCalcm false)%n", i);
      }
      row("cut[" + i + "].critical_volume", model.calcCriticalVolume(molarMass, density));
      // **The models do not share one unit for `molarMass`, and `addTBPfraction` compensates per
      // call site.** There the local is kg/mol, and it passes `molarMass * 1000.0` to
      // `calcCriticalViscosity` (g/mol) but the bare `molarMass` to `calcParachorParameter`, which
      // multiplies by 1000 inside. This probe works in g/mol throughout, so parachor takes the
      // thousandth and the others do not.
      row("cut[" + i + "].parachor", model.calcParachorParameter(molarMass / 1000.0, density));
      row("cut[" + i + "].critical_viscosity", model.calcCriticalViscosity(molarMass, density));
      row("cut[" + i + "].watson_k", model.calcWatsonCharacterizationFactor(molarMass, density));
      System.out.println();
    }
  }

  public static void main(String[] args) {
    // MW in g/mol and density in g/cm3, which are the units `addTBPfraction` hands over.
    double[][] cuts = {
      {100.0, 0.70},
      {150.0, 0.75},
      {200.0, 0.80},
      {300.0, 0.85},
      {500.0, 0.88},
      {539.0, 0.89},
      {541.0, 0.90},
      {600.0, 0.90},
      {900.0, 0.95},
      {1119.0, 0.99},
      {1121.0, 1.00},
      {2000.0, 1.05},
    };

    for (String name : TBPfractionModel.getAvailableModels()) {
      sweep(name, cuts);
    }

    // **`calcRacketZ` is not here**, as in `TbpFractionProbe`: it takes a *flashed* reference
    // system - the Peneloux shift is the difference between a cut's measured molar volume and the
    // cubic's - so it is not a function of the two numbers a cut is described by.
    System.out.println("# racketZ needs a flashed reference system and is not a function of (MW, density)");
  }
}
