// The two temperature derivatives of the Soreide-Whitson water alpha.
//
//     javac -proc:none -cp neqsim-3.20.0.jar SoreideAlphaDerivatives.java
//     java -cp .:neqsim-3.20.0.jar SoreideAlphaDerivatives
//
// `AttractiveTermSoreideWhitson` provides `alpha`, `diffalphaT` and `diffdiffalphaT`, which
// are `d(alpha)/dT` and `d2(alpha)/dT2` **in kelvin, not in reduced temperature** - both
// carry a `1/Tc` per derivative. `eos.soreide_whitson_alpha` reports only the first, so the
// other two are captured here before being ported.
//
// The salinity is set by reflection, because `setSalinityFromPhase` is the method
// `SystemSoreideWhitson.calcSalinity` never reaches - the same defect reported separately -
// so a phase can never drive these at a real brine.

import neqsim.thermo.component.ComponentSoreideWhitson;
import neqsim.thermo.component.attractiveeosterm.AttractiveTermSoreideWhitson;

public class SoreideAlphaDerivatives {
  public static void main(String[] args) throws Exception {
    // **A bare component, not a phase.** The term reads only its `Tc` and its name, and a
    // one-component system has no pairs for `calcA` to fill, so building one throws before
    // the term is ever exercised.
    ComponentSoreideWhitson water = new ComponentSoreideWhitson("water", 1.0, 1.0, 0);
    double tc = water.getTC();
    System.out.printf("water Tc = %.10g K%n", tc);
    System.out.printf("%10s %10s %24s %24s %24s%n", "salinity", "T / K", "alpha", "dalpha/dT",
        "d2alpha/dT2");
    for (double salinity : new double[] {0.0, 1.0, 4.0, 8.0}) {
      AttractiveTermSoreideWhitson term =
          (AttractiveTermSoreideWhitson) water.getAttractiveTerm();
      java.lang.reflect.Method set =
          AttractiveTermSoreideWhitson.class.getDeclaredMethod("setSalinityFromPhase", double.class);
      set.setAccessible(true);
      set.invoke(term, salinity);
      for (double t : new double[] {273.15, 298.15, 323.15, 373.15, 473.15, 573.15}) {
        System.out.printf("%10.4f %10.2f %24.15g %24.15g %24.15g%n", salinity, t,
            term.alpha(t), term.diffalphaT(t), term.diffdiffalphaT(t));
      }
    }
  }
}
