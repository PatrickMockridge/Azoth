// The three pure-component vapour pressures of the water-nitric-sulfuric acid system, for
// `eos.nitric_sulfuric_acid_vapor_pressure`.
//
//     javac -proc:none -cp neqsim-f0c7436.jar AcidVaporPressure.java
//     java -cp .:neqsim-f0c7436.jar AcidVaporPressure
//
// `NitricSulfuricAcidVaporPressure` is a static utility: no system, no flash, no tuning
// between the caller and the number. So the oracle for this calculation is the most direct
// in the port - each of the three functions is called at a state and its pascal value
// printed, and there is nothing else in the path.
//
// The three are different equations, not one form with three coefficient sets:
//
//   water          log10(P0/mbar) = 8.42926609 - 1827.17843/T - 71208.271/T^2,  x 100 Pa
//   nitric acid    log10(P0/torr) = 7.57628 - 1470.385/(T - 43.0),              x 133.322368421 Pa
//   sulfuric acid  ln(P0/atm)     = -10156.0/T + 16.259,                        x 101325 Pa
//
// The nitric-acid pair is an engineering adjustment rather than the paper's: NeqSim's
// javadoc records A/B refitted from Pennington (1951)'s 7.61628/1486.238 to reproduce the
// 83 C boiling point and to raise P0 at 273.15 K by 6.9 % against Vandoni (1944)'s ternary
// salting-out data. C = 43.0 is unchanged. The port carries NeqSim's values.
//
// The validity is per species rather than per class: water and sulfuric acid are stated
// for roughly 190-298 K, nitric acid's Antoine for roughly 190-400 K, and the class's own
// MIN/MAX_RECOMMENDED_TEMPERATURE_K of 190/298 is the conservative envelope.

import neqsim.thermo.util.empiric.NitricSulfuricAcidVaporPressure;

public class AcidVaporPressure {

  static void one(double temperatureK) {
    System.out.println("T=" + temperatureK
        + "  water=" + NitricSulfuricAcidVaporPressure.pureVaporPressureWater(temperatureK)
        + "  nitric=" + NitricSulfuricAcidVaporPressure.pureVaporPressureNitricAcid(temperatureK)
        + "  sulfuric=" + NitricSulfuricAcidVaporPressure.pureVaporPressureSulfuricAcid(temperatureK));
  }

  public static void main(String[] args) {
    one(190.0);
    one(220.0);
    one(250.0);
    one(273.15);
    one(298.0);
    one(340.0);
  }
}
