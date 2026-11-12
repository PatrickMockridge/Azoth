// `LMTDcorrectionFactor` driven directly, on states of its own.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-118d9af.jar LmtdCorrectionProbe.java
//     java -cp .:neqsim-118d9af.jar LmtdCorrectionProbe > captures/lmtd_correction_probe.tsv
//
// **The four temperature inputs are the class's own**, and `R`, `P` and `F_t` are its
// intermediates. The multi-pass path is the Bowman-Mueller-Nagle conversion: the overall
// P becomes a per-shell P1 before the single-shell formula runs.

import neqsim.process.mechanicaldesign.heatexchanger.LMTDcorrectionFactor;

public class LmtdCorrectionProbe {

  public static void main(String[] args) {
    // The class's own R/P test state: R = (150-90)/(80-30) = 1.2, P = 50/120.
    row("class_rp_state", 150.0, 90.0, 30.0, 80.0);
    // R = 1: equal temperature changes on both sides.
    row("r_equals_one", 150.0, 90.0, 30.0, 90.0);
    // R = 1 at a small effectiveness: the R = 1 branch is not clamped, so F_t can pass one.
    row("r_equals_one_small_p", 150.0, 140.0, 30.0, 40.0);
    // No heat transfer: the cold side does not move, so P is zero and F_t is one.
    row("no_heat_transfer", 150.0, 150.0, 30.0, 30.0);
    // Equal hot and cold inlets: the denominator of P is zero, which the class answers P = 0.
    row("same_inlets", 150.0, 90.0, 150.0, 80.0);
    // A temperature cross: R*P >= 1, so the single-shell formula answers zero.
    row("temperature_cross", 150.0, 20.0, 30.0, 80.0);
    // A large effectiveness that needs three shell passes to clear 0.75.
    row("multi_pass_need", 150.0, 54.0, 30.0, 126.0);
    // A cross that no number of shell passes up to six clears.
    row("unachievable", 150.0, 31.2, 30.0, 148.8);
    System.out.flush();
  }

  static void row(String label, double hotIn, double hotOut, double coldIn, double coldOut) {
    System.out.println(label);
    System.out.println("hot_inlet_k=" + hotIn);
    System.out.println("hot_outlet_k=" + hotOut);
    System.out.println("cold_inlet_k=" + coldIn);
    System.out.println("cold_outlet_k=" + coldOut);
    System.out.println("r=" + LMTDcorrectionFactor.calcR(hotIn, hotOut, coldIn, coldOut));
    System.out.println("p=" + LMTDcorrectionFactor.calcP(hotIn, hotOut, coldIn, coldOut));
    System.out.println(
        "ft_1_shell_pass=" + LMTDcorrectionFactor.calcFt1ShellPass(hotIn, hotOut, coldIn, coldOut));
    System.out.println(
        "ft_2_shell_passes=" + LMTDcorrectionFactor.calcFt2ShellPass(hotIn, hotOut, coldIn, coldOut));
    for (int passes = 1; passes <= 6; passes++) {
      System.out.println(
          "ft_" + passes + "_passes="
              + LMTDcorrectionFactor.calcFt(hotIn, hotOut, coldIn, coldOut, passes));
    }
    System.out.println(
        "required_shell_passes="
            + LMTDcorrectionFactor.requiredShellPasses(hotIn, hotOut, coldIn, coldOut));
    System.out.println();
  }
}
