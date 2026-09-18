import neqsim.thermo.component.ComponentEos;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

/**
 * Whether NeqSim's cubic `dFdN` is the derivative of its own `F`.
 *
 * <p>
 * `CpaSweep` measures this for an *associating* mixture and finds that it is not: at 356 K
 * and 1 bara water/methanol gives `dFdN = [-0.0901941, -0.1211108]` against a finite
 * difference of `[-0.0419683, -0.0743505]`, while `ComponentSrkCPA.dFCPAdN` is a true
 * derivative. This asks the same question of a **plain** cubic, which decides how far the
 * finding reaches.
 *
 * <p>
 * It matters because NeqSim's ordinary cubic flashes *do* agree with azoth's, and azoth's
 * `ln phi` is the textbook Soave-Redlich-Kwong expression. If a plain cubic's `dFdN` also
 * fails to be its own `F`'s derivative, then NeqSim's cubic `ln phi` is not `integral
 * dFdN` either - and the error would have to cancel in `ln phi_L - ln phi_V`, which is
 * what a flash compares.
 *
 * <p>
 * `F` is a function of `(T, V, n)`, so along any path
 * `dF = sum_i (dF/dn_i)|_{T,V} dn_i + (dF/dV)|_{T,n} dV`. Perturbing one mole number at
 * constant `T` and `P` therefore isolates `(dF/dn_i)|_{T,V}` once the volume's share is
 * subtracted, and no prescribed-volume state is needed.
 *
 * <p>
 * Usage:
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar CubicFdProbe.java
 * java -cp .:neqsim-3.20.0.jar CubicFdProbe
 * </pre>
 */
public final class CubicFdProbe {

  private CubicFdProbe() {}

  private static final double T_K = 330.0;
  private static final double P_BARA = 25.0;
  private static final double D = 1.0e-7;

  /** One PR methane/n-butane state at explicit mole numbers. */
  private static PhaseInterface at(double methane, double butane) {
    SystemInterface system = new SystemPrEos(T_K, P_BARA);
    system.addComponent("methane", methane);
    system.addComponent("n-butane", butane);
    system.setMixingRule("classic");
    system.setAttractiveTerm(1);
    system.init(0);
    system.init(1);
    return system.getPhase(0);
  }

  private static void row(String label, double methane, double butane) {
    PhaseInterface p = at(methane, butane);
    double f = ((PhaseEos) p).getF();
    double fv = p.FV();
    double v = p.getTotalVolume();
    System.out.printf("%-18s F=%.17g FV=%.17g V=%.17g%n", label, f, fv, v);
    for (int i = 0; i < p.getNumberOfComponents(); i++) {
      ComponentEos c = (ComponentEos) p.getComponent(i);
      System.out.printf("  dFdN[%d]=%.17g lnphi=%.17g%n", i,
          c.dFdN(p, p.getNumberOfComponents(), T_K, P_BARA),
          Math.log(c.getFugacityCoefficient()));
    }
  }

  /** The isolated `dF/dn_i` at constant `T` and `V`, from two constant-`T,P` solves. */
  private static void difference(String label, double methane, double butane, int which) {
    double dM = which == 0 ? D : 0.0;
    double dB = which == 1 ? D : 0.0;
    PhaseInterface up = at(methane + dM, butane + dB);
    PhaseInterface down = at(methane - dM, butane - dB);
    PhaseInterface base = at(methane, butane);
    double fv = base.FV();
    double df = ((PhaseEos) up).getF() - ((PhaseEos) down).getF();
    double dv = up.getTotalVolume() - down.getTotalVolume();
    System.out.printf("%-18s dFdN_fd[%d]=%.17g   (dF=%.17g FV*dV=%.17g)%n",
        label, which, (df - fv * dv) / (2.0 * D), df, fv * dv);
  }

  public static void main(String[] args) {
    double methane = 0.6;
    double butane = 0.4;
    row("base", methane, butane);
    difference("base", methane, butane, 0);
    difference("base", methane, butane, 1);
    // The extensivity test: doubling every mole number at the same T and P doubles an
    // extensive `F` and leaves a molar one alone.
    System.out.printf("F_scale_two=%.17g%n", ((PhaseEos) at(2.0 * methane, 2.0 * butane)).getF());
  }
}
