import neqsim.thermodynamicoperations.flashops.RachfordRice;

/**
 * Drives NeqSim's `RachfordRice` on the states `eos.rachford_rice` carries cases for.
 *
 * <p>Both solvers are called, so the two can be compared at the same input, and the
 * feeds that have no root are included so the clamp is measured rather than assumed.
 */
public class RachfordRiceProbe {

  private static void probe(String label, double[] z, double[] K) {
    RachfordRice rr = new RachfordRice();
    String nielsen;
    String michelsen;
    try {
      nielsen = String.format("%.17g", rr.calcBetaNielsen2023(K, z));
    } catch (Exception e) {
      nielsen = e.getClass().getSimpleName();
    }
    try {
      michelsen = String.format("%.17g", rr.calcBetaMichelsen2001(K, z));
    } catch (Exception e) {
      michelsen = e.getClass().getSimpleName();
    }
    System.out.printf("%-34s nielsen=%-24s michelsen=%s%n", label, nielsen, michelsen);
  }

  public static void main(String[] args) {
    // A real two-phase state, and a second one so the answer is not one number.
    probe("CH4/nC4 z=0.6/0.4", new double[] {0.6, 0.4},
        new double[] {7.304244305324782, 0.33749596785762953});
    probe("CH4/nC4 z=0.6/0.4 (300 K)", new double[] {0.6, 0.4},
        new double[] {5.799172708809654, 0.14913889826410118});

    // No root: every K on one side of one.
    probe("all liquid z=0.5/0.5", new double[] {0.5, 0.5}, new double[] {0.2, 0.3});
    probe("all vapour z=0.5/0.5", new double[] {0.5, 0.5}, new double[] {3.0, 5.0});

    // The reformulation's own case: K spanning twelve orders of magnitude.
    probe("wide K z=0.3/0.7", new double[] {0.3, 0.7}, new double[] {1.0e6, 1.0e-6});

    // A root below zero - this library's negative flash - reached both ways round.
    probe("subcooled z=0.1/0.9", new double[] {0.1, 0.9},
        new double[] {5.799172708809655, 0.14913889826410245});
    probe("superheated z=0.5/0.5", new double[] {0.5, 0.5}, new double[] {1.5, 0.9});

    // An ion alongside a splitting pair.
    probe("ion z=0.1/0.45/0.45", new double[] {0.1, 0.45, 0.45},
        new double[] {1.0e-40, 7.304244305324782, 0.33749596785762953});
  }
}
