// The TBP pseudo-component's correlations, printed one pair at a time so they can be ported.
//
// A wax fluid is built from TBP and plus fractions, and `addTBPfraction(moles, MW, density)`
// turns each into a component through `Characterise.TBPFractionModelName`, whose default is
// **PedersenSRK**. Everything below is what that class computes from the two numbers a cut is
// described by, so a port has the equations *and* the values they give.
//
// The pairs are chosen to cross every branch: `calcTB` switches at MW 540 and `calcTC` and
// `calcPC` switch their coefficient set at MW 1120, and a port that picked the wrong side of
// either would still be right at a single pair.
//
//   javac -proc:none -cp neqsim-f0c7436.jar TbpFractionProbe.java
//   java -cp .:neqsim-f0c7436.jar TbpFractionProbe > captures/tbp_fraction_probe.tsv

import neqsim.thermo.characterization.TBPfractionModel;
import neqsim.thermo.characterization.TBPModelInterface;

public class TbpFractionProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    TBPModelInterface model = new TBPfractionModel().getModel("PedersenSRK");
    System.out.printf("# model = %s%n", model.getName());

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

    for (int i = 0; i < cuts.length; i++) {
      double molarMass = cuts[i][0];
      double density = cuts[i][1];
      System.out.printf("# cut[%d] MW = %.15g g/mol, density = %.15g g/cm3%n", i, molarMass, density);
      row("cut[" + i + "].mw", molarMass);
      row("cut[" + i + "].density", density);
      row("cut[" + i + "].tc", model.calcTC(molarMass, density));
      row("cut[" + i + "].pc", model.calcPC(molarMass, density));
      row("cut[" + i + "].tb", model.calcTB(molarMass, density));
      row("cut[" + i + "].acentric", model.calcAcentricFactor(molarMass, density));
      row("cut[" + i + "].calcm", model.isCalcm() ? model.calcm(molarMass, density) : Double.NaN);
      row("cut[" + i + "].critical_volume", model.calcCriticalVolume(molarMass, density));
      row("cut[" + i + "].parachor", model.calcParachorParameter(molarMass, density));
      row("cut[" + i + "].critical_viscosity", model.calcCriticalViscosity(molarMass, density));
      row(
          "cut[" + i + "].watson_k",
          model.calcWatsonCharacterizationFactor(molarMass, density));
      System.out.println();
    }

    // **`calcRacketZ` is not here**, and that is a statement rather than an omission: it takes
    // a *flashed* reference system - the Peneloux shift is the difference between a cut's
    // measured molar volume and the cubic's - so it is not a function of the two numbers a cut
    // is described by and it cannot be captured from this probe's inputs.
    System.out.println("# racketZ needs a flashed reference system and is not a function of (MW, density)");
  }
}
