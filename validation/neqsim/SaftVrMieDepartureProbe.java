import java.lang.reflect.Field;
import neqsim.thermo.component.ComponentSAFTVRMie;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseSAFTVRMie;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSAFTVRMie;

/**
 * NeqSim's SAFT-VR-Mie temperature derivatives at one state, in the order the departure is
 * built from them.
 *
 * <p>
 * {@code eos.saft_vr_mie_phase} reports a state but not an enthalpy. This is the oracle for
 * that port: each component's diameter and its derivative, the packing fraction's own
 * derivative, then each Helmholtz term's derivative, then the departure NeqSim assembles
 * from them.
 *
 * <p>
 * <b>Every key is a {@code d/dT}, not a {@code T d/dT}.</b> NeqSim carries the plain
 * derivative and the departure multiplies by the temperature where it needs to; azoth's
 * reduced form carries {@code T d/dT} so that no absolute temperature is needed. The header
 * says which is which on every line, because the two differ by a factor of 350 at an
 * ordinary state and nothing about the number says which it is.
 *
 * <p>
 * <b>Unlike PC-SAFT, this class's {@code dFdT} is correct</b> - {@code PhaseSAFTVRMie}'s
 * diameter derivative is the plain chain rule, where {@code PhasePCSAFT.getdDSAFTdT}
 * multiplies it by an extra {@code 3 d^2}. So the numbers below are the oracle, and the
 * finite difference at the bottom is a second opinion rather than the only one: it walks the
 * constant-pressure path and subtracts the volume's own contribution, which uses none of the
 * derivative code above it.
 *
 * <p>
 * The dispersion's temperature derivatives are <b>central differences in {@code T}</b>
 * ({@code PhaseSAFTVRMie.java:1225}, step {@code T 1e-5}), and its {@code eta} derivatives
 * are central differences in {@code eta}. So an analytic port should land on these to about
 * {@code 1e-11} relative, and a disagreement at {@code 1e-6} is one of the two steps rather
 * than a term. The individual {@code da1/2/3DispDT} fields have no getter and are read by
 * reflection.
 *
 * <p>
 * <b>{@code dF_HC_SAFTdT} does not survive that check for any fluid with a chain</b>, and
 * the {@code lnGCE} lines are why: they evaluate {@code calcLnGChainEffective} by hand at
 * {@code eta} and {@code eta +- |eta| 1e-4} and reproduce {@code dgHSSAFTdN} exactly, while
 * the same block's {@code gHS} line gives the derivative of the contact value the chain
 * term actually uses - twenty-five times larger. So the difference quotient is faithful and
 * the value it is differencing is not the one kept. See
 * {@code ~/Desktop/neqsim-saft-vr-mie-chain-contact-value-temperature.md}.
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar SaftVrMieDepartureProbe.java
 * java -cp .:neqsim-3.20.0.jar SaftVrMieDepartureProbe [T_K P_bara name:z ...]
 * </pre>
 */
public final class SaftVrMieDepartureProbe {

  private SaftVrMieDepartureProbe() {}

  private static double field(PhaseSAFTVRMie phase, String name) {
    try {
      Field f = PhaseSAFTVRMie.class.getDeclaredField(name);
      f.setAccessible(true);
      return f.getDouble(phase);
    } catch (ReflectiveOperationException e) {
      throw new IllegalStateException(name, e);
    }
  }

  /**
   * {@code {F, F_hc, F_disp}} at a temperature with the molar volume pinned, so that the
   * difference is taken at constant volume.
   */
  private static double[] fAtFixedVolume(String[] names, double[] z, double t, double pBar,
      double vm, int n) {
    SystemInterface s = new SystemSAFTVRMie(t, pBar);
    for (int i = 0; i < n; i++) {
      s.addComponent(names[i], z[i]);
    }
    s.setMixingRule("classic");
    s.init(0);
    s.init(1);
    s.init(2);
    PhaseSAFTVRMie ph = (PhaseSAFTVRMie) s.getPhase(0);
    ph.setMolarVolume(vm);
    ph.volInit();
    return new double[] {ph.getF(), ph.F_HC_SAFT(), ph.F_DISP_SAFT(), ph.getGhsSAFT()};
  }

  private static double lnGChainEffective(PhaseSAFTVRMie phase, double eta, double alpha) {
    try {
      java.lang.reflect.Method m =
          PhaseSAFTVRMie.class.getDeclaredMethod("calcLnGChainEffective", double.class, double.class);
      m.setAccessible(true);
      return (double) m.invoke(phase, eta, alpha);
    } catch (ReflectiveOperationException e) {
      throw new IllegalStateException(e);
    }
  }

  private static void probe(double t, double pBar, String[] names, double[] z) {
    SystemInterface system = new SystemSAFTVRMie(t, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    system.init(2);

    PhaseInterface phase = system.getPhase(0);
    PhaseSAFTVRMie mie = (PhaseSAFTVRMie) phase;
    int n = phase.getNumberOfComponents();

    System.out.printf("%n# %s at T=%g K, P=%g bara%n", String.join("/", names), t, pBar);
    for (int i = 0; i < n; i++) {
      ComponentSAFTVRMie c = (ComponentSAFTVRMie) phase.getComponent(i);
      System.out.printf("component[%d] %s d=%.15g ddSAFTidT=%.15g%n", i, names[i], c.getdSAFTi(),
          c.getDdSAFTidT());
    }
    double nmoles = phase.getNumberOfMolesInPhase();
    System.out.printf("mSAFT=%.15g mmin1SAFT=%.15g%n", mie.getmSAFT(), mie.getMmin1SAFT());
    System.out.printf("dSAFT=%.15g dDSAFTdT=%.15g%n", mie.getDSAFT(), mie.getdDSAFTdT());
    System.out.printf("nSAFT=%.15g dNSAFTdT=%.15g  (both d/dT)%n", mie.getNSAFT(),
        field(mie, "dNSAFTdT"));
    System.out.printf("aHS=%.15g gHS=%.15g%n", mie.getAHSSAFT(), mie.getGhsSAFT());

    System.out.printf("A1Disp=%.15g da1DispDT=%.15g da1DispDeta=%.15g%n", mie.getA1Disp(),
        field(mie, "da1DispDT"), mie.getDa1DispDeta());
    System.out.printf("A2Disp=%.15g da2DispDT=%.15g da2DispDeta=%.15g%n", mie.getA2Disp(),
        field(mie, "da2DispDT"), field(mie, "da2DispDeta"));
    System.out.printf("A3Disp=%.15g da3DispDT=%.15g da3DispDeta=%.15g%n", mie.getA3Disp(),
        field(mie, "da3DispDT"), field(mie, "da3DispDeta"));

    System.out.printf("F_hc=%.15g dF_HC_SAFTdT=%.15g%n", mie.F_HC_SAFT(), mie.dF_HC_SAFTdT());
    System.out.printf("F_disp=%.15g dF_DISP_SAFTdT=%.15g%n", mie.F_DISP_SAFT(),
        mie.dF_DISP_SAFTdT());
    System.out.printf("F=%.15g dFdT=%.15g%n", mie.getF(), mie.dFdT());

    // The departure NeqSim assembles: `AresTV + T SresTV + P V - n R T`, with
    // `SresTV = -(dAres/dT)_V`. So `Hres = -T d(Ares)/dT + P V - n R T` and, per mole,
    // `Hres/(R T) = -T dFdT + Z - 1` when `Ares = R T F`.
    System.out.printf("Z=%.15g v=%.15g%n", phase.getZ(),
        phase.getVolume() / phase.getNumberOfMolesInPhase());
    System.out.printf("HresTP/n=%.15g SresTP/n=%.15g%n", phase.getHresTP() / nmoles,
        phase.getSresTP() / nmoles);

    // **The (dF/dT)_V the port is checked against.** NeqSim's own `dFdT` is usable here, so
    // this is a second opinion: it perturbs the temperature, pins the molar volume back to
    // the converged one and re-runs `volInit` - the same call `dFdVdVdV` makes - so it reads
    // `F` at `(T +- h, v)` and uses none of the `d/dT` code above.
    double vm = phase.getMolarVolume();
    double[] hs = {0.01, 0.05, 0.2};
    for (int k = 0; k < 3; k++) {
      double[] p = fAtFixedVolume(names, z, t + hs[k], pBar, vm, n);
      double[] m = fAtFixedVolume(names, z, t - hs[k], pBar, vm, n);
      System.out.printf(
          "dFdT_at_V_fd[h=%g] F=%.15g F_hc=%.15g F_disp=%.15g gHS=%.15g%n", hs[k],
          (p[0] - m[0]) / (2.0 * hs[k]), (p[1] - m[1]) / (2.0 * hs[k]),
          (p[2] - m[2]) / (2.0 * hs[k]), (p[3] - m[3]) / (2.0 * hs[k]));
    }
    System.out.printf("daHSSAFTdN=%.15g dgHSSAFTdN=%.15g dg/dEta_from_gHS_fd=%.15g%n",
        field(mie, "daHSSAFTdN"), field(mie, "dgHSSAFTdN"),
        (fAtFixedVolume(names, z, t + 0.01, pBar, vm, n)[3]
            - fAtFixedVolume(names, z, t - 0.01, pBar, vm, n)[3]) / (2.0 * 0.01)
            / field(mie, "dNSAFTdT"));
    // **The function `dgHSSAFTdN` claims to be the derivative of.** `ghsSAFT` is
    // `exp(calcLnGChainEffective(eta, alpha))`, and the class differences that same call at
    // `eta +- |eta| 1e-4`. Evaluating it by hand at the same three points settles whether the
    // value or the difference is wrong.
    double eta = mie.getNSAFT();
    double dEta = Math.max(Math.abs(eta) * 1.0e-4, 1.0e-12);
    for (double alpha : new double[] {1.0, 0.0}) {
      double ln0 = lnGChainEffective(mie, eta, alpha);
      double lnP = lnGChainEffective(mie, eta + dEta, alpha);
      double lnM = lnGChainEffective(mie, Math.max(eta - dEta, 1.0e-15), alpha);
      System.out.printf("lnGCE[alpha=%g] eta=%.15g etaP=%.15g etaM=%.15g fd=%.15g%n", alpha, ln0,
          lnP, lnM, (Math.exp(lnP) - Math.exp(lnM)) / (2.0 * dEta));
    }
    System.out.printf("dFdT_at_V_neqsim=%.15g T_dFdT_at_V_neqsim=%.15g%n", mie.dFdT(),
        t * mie.dFdT());

    // `dF/dv` two ways, because this class scales its volume derivatives by `1e-5` where
    // PC-SAFT's are in m^3 and the convention is not visible in the number.
    double dv = vm * 1.0e-5;
    double dfdvFd = (fAtFixedVolume(names, z, t, pBar, vm + dv, n)[0]
        - fAtFixedVolume(names, z, t, pBar, vm - dv, n)[0]) / (2.0 * dv);
    System.out.printf("dFdV=%.15g dFdV_x1e5=%.15g dFdv_fd=%.15g%n", mie.dFdV(),
        mie.dFdV() * 1.0e5, dfdvFd);
  }

  public static void main(String[] args) {
    System.out.println("# azoth SaftVrMieDepartureProbe - NeqSim 3.20.0's SystemSAFTVRMie.");
    System.out.println("# Every dT key is d/dT; azoth carries T*d/dT, so multiply by T.");
    if (args.length >= 4) {
      double t = Double.parseDouble(args[0]);
      double p = Double.parseDouble(args[1]);
      String[] names = new String[(args.length - 2) / 2];
      double[] z = new double[names.length];
      for (int i = 0; i < names.length; i++) {
        names[i] = args[2 + 2 * i];
        z[i] = Double.parseDouble(args[3 + 2 * i]);
      }
      probe(t, p, names, z);
      return;
    }
    probe(350.0, 30.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4});
    probe(300.0, 50.0, new String[] {"methane"}, new double[] {1.0});
  }
}
