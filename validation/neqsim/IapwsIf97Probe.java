// NeqSim's `Iapws_if97` steam tables, so azoth's IF97 kernel has an oracle this repo re-runs.
//
//     javac -proc:none -cp neqsim-f0c7436.jar IapwsIf97Probe.java
//     java -cp .:neqsim-f0c7436.jar IapwsIf97Probe > captures/iapws_if97_probe.tsv
//
// `crates/azoth-eos/src/iapws_if97.rs` and `eos.water_phase` were landed with the literals of
// NeqSim's own shipped `Iapws_if97Test` copied by hand - four states and two saturation points.
// That is an oracle of a sort, but it is not this repository's: a pinned-jar move could not say
// whether it moved, because there is nothing here to re-run. This probe prints the class's
// public surface over a spread wide enough to cross both regions and both saturation
// directions, and the capture becomes the source the Rust test and the model's cases read.
//
// **Both saturation directions are printed, deliberately.** `p4_T` and `T4_p` are separate fits
// in the class and are not inverses of one another, so the pair can be checked against each
// other as well as against azoth.
//
// **No near-critical state.** The class implements Regions 1 and 2 and the Region 4 curve;
// Region 3 is not in it, so a state above Region 1's 623.15 K ceiling at a pressure inside the
// dome would be answered by Region 1's coefficients outside their own domain. The states below
// stay inside the two regions the formulation defines. The most interesting of them is
// `0.101325 MPa / 373.15 K`, which is **vapour**: the saturation temperature at that pressure
// is 373.1243 K, so the boiling point of water at one atmosphere is twelve hundredths of a
// kelvin below the 373.15 K a reader would reach for.
//
// Units are the class's own and not SI: pressure in MPa, temperature in K, specific volume in
// m^3/kg, enthalpy and entropy per kg, heat capacity per kg, and the speed of sound in m/s.

import neqsim.thermo.util.steam.Iapws_if97;

public class IapwsIf97Probe {

  /** The Region 4 curve is only defined between the triple and the critical point. */
  private static final double[] SATURATION_TEMPERATURES = {273.16, 300.0, 373.15, 500.0, 640.0};

  private static final double[] SATURATION_PRESSURES = {0.000611, 0.101325, 1.0, 5.0, 10.0, 22.064};

  /** `{p [MPa], T [K]}`, the four `eos.water_phase` states first. */
  private static final double[][] STATES = {
      {1.0, 300.0},
      {5.0, 500.0},
      {0.1, 400.0},
      {10.0, 600.0},
      {0.101325, 373.15},
      {3.0, 500.0},
      {22.0, 373.15},
      {1.0, 773.15},
  };

  public static void main(String[] args) {
    System.out.println("# azoth IapwsIf97Probe - NeqSim's thermo.util.steam.Iapws_if97.");
    System.out.println("# Pressure in MPa, temperature in K, and the class's own units for");
    System.out.println("# everything else: m^3/kg, kJ/kg, kJ/(kg.K), m/s.");

    System.out.println();
    System.out.println("# Region 4, forward: the saturation pressure at a temperature");
    for (double t : SATURATION_TEMPERATURES) {
      System.out.printf("psat_t(%s) = %.15g%n", trim(t), Iapws_if97.psat_t(t));
    }

    System.out.println();
    System.out.println("# Region 4, backward: the saturation temperature at a pressure");
    for (double p : SATURATION_PRESSURES) {
      System.out.printf("tsat_p(%s) = %.15g%n", trim(p), Iapws_if97.T4_p(p));
    }

    System.out.println();
    System.out.println("# the properties at (p, T), one state per line so the region each");
    System.out.println("# falls in can be read beside the number it produces");
    for (double[] state : STATES) {
      double p = state[0];
      double t = state[1];
      System.out.printf(
          "p_MPa=%s T_K=%s region=%d v=%.15g h=%.15g s=%.15g cp=%.15g w=%.15g%n",
          trim(p),
          trim(t),
          t <= Iapws_if97.T4_p(p) ? 1 : 2,
          Iapws_if97.v_pt(p, t),
          Iapws_if97.h_pt(p, t),
          Iapws_if97.s_pt(p, t),
          Iapws_if97.cp_pt(p, t),
          Iapws_if97.w_pt(p, t));
    }
  }

  /** A state's own number with no trailing zeros, so a key reads as the argument it was given. */
  private static String trim(double value) {
    return java.math.BigDecimal.valueOf(value).stripTrailingZeros().toPlainString();
  }
}
