import neqsim.thermo.component.ComponentPCSAFT;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhasePCSAFT;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPCSAFT;

/**
 * NeqSim's PC-SAFT temperature derivatives at one state, in the order the departure is
 * built from them.
 *
 * <p>
 * {@code eos.pcsaft_rahmat_phase} reports a state but not an enthalpy, because the
 * temperature derivative of its Helmholtz energy was never ported. This is the oracle for
 * that port: the segment diameter and its derivative, then each Helmholtz term's
 * derivative, then the departure NeqSim assembles from them.
 *
 * <p>
 * <b>Every key is a {@code d/dT}, not a {@code T d/dT}.</b> NeqSim carries the plain
 * derivative and the departure multiplies by the temperature where it needs to; azoth's
 * reduced form carries {@code T d/dT} so that no absolute temperature is needed. The
 * header says which is which on every line, because the two differ by a factor of 298 at
 * an ordinary state and nothing about the number says which it is.
 *
 * <pre>
 * javac -proc:none -cp neqsim-f0c7436.jar PcsaftDepartureProbe.java
 * java -cp .:neqsim-f0c7436.jar PcsaftDepartureProbe [T_K P_bara name:z ...]
 * </pre>
 */
public final class PcsaftDepartureProbe {

  private PcsaftDepartureProbe() {}

  private static void probe(double t, double pBar, String[] names, double[] z) {
    SystemInterface system = new SystemPCSAFT(t, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    system.init(2);

    PhaseInterface phase = system.getPhase(0);
    PhasePCSAFT safT = (PhasePCSAFT) phase;
    int n = phase.getNumberOfComponents();

    System.out.printf("%n# %s at T=%g K, P=%g bara%n", String.join("/", names), t, pBar);
    for (int i = 0; i < n; i++) {
      ComponentPCSAFT c = (ComponentPCSAFT) phase.getComponent(i);
      System.out.printf("d_i[%d]=%.15g%n", i, c.getdSAFTi());
    }
    System.out.printf("dSAFT=%.15g dmeanSAFT=%.15g mSAFT=%.15g%n", safT.getDSAFT(),
        safT.getDmeanSAFT(), safT.getmSAFT());
    System.out.printf("dDSAFTdT=%.15g  (d/dT of dSAFT)%n", safT.getdDSAFTdT());
    System.out.printf("nSAFT=%.15g dnSAFTdV=%.15g%n", safT.getNSAFT(), safT.getDnSAFTdV());

    System.out.printf("F_hc=%.15g dF_hc_dT=%.15g%n", safT.F_HC_SAFT(), safT.dF_HC_SAFTdT());
    System.out.printf("F_disp1=%.15g dF_disp1_dT=%.15g%n", safT.F_DISP1_SAFT(),
        safT.dF_DISP1_SAFTdT());
    System.out.printf("F_disp2=%.15g dF_disp2_dT=%.15g%n", safT.F_DISP2_SAFT(),
        safT.dF_DISP2_SAFTdT());
    System.out.printf("F=%.15g dFdT=%.15g%n", safT.getF(), safT.dFdT());

    // The departure NeqSim assembles: `AresTV + T SresTV + P V - n R T`, with
    // `SresTV = -(dAres/dT)_V`. So `Hres = -T d(Ares)/dT + P V - n R T` and, per mole,
    // `Hres/(R T) = -T dFdT + Z - 1` when `Ares = R T F`.
    PhaseEos eos = (PhaseEos) phase;
    double moles = phase.getNumberOfMolesInPhase();
    double pv = phase.getPressure() * phase.getMolarVolume();
    System.out.printf("Z=%.15g%n", phase.getZ());
    System.out.printf("AresTV/(nRT)=%.15g T*SresTV/(nRT)=%.15g PV/(nRT)=%.15g%n",
        eos.getAresTV() / (moles * 8.3144621 * t), t * eos.getSresTV() / (moles * 8.3144621 * t),
        pv / (8.3144621 * t));
    System.out.printf("HresTP/n=%.15g SresTP/n=%.15g%n", phase.getHresTP() / moles,
        phase.getSresTP() / moles);

    // **The (dF/dT)_V the port is checked against.** NeqSim's own `dFdT` cannot be used for
    // PC-SAFT: `getdDSAFTdT` multiplies the correct chain rule by `3 d**2`, so every
    // eta-dependent part of it is ~1e-18 too small. This instead walks the constant-pressure
    // path and subtracts the volume's own contribution, which uses no derivative code of
    // NeqSim's that this defect reaches.
    double h = 0.05;
    double[] f = new double[2];
    double[] vol = new double[2];
    for (int k = 0; k < 2; k++) {
      double tk = t + (k == 0 ? h : -h);
      SystemInterface s2 = new SystemPCSAFT(tk, pBar);
      for (int i = 0; i < n; i++) {
        s2.addComponent(names[i], z[i]);
      }
      s2.setMixingRule("classic");
      s2.init(0);
      s2.init(1);
      s2.init(2);
      PhasePCSAFT ph = (PhasePCSAFT) s2.getPhase(0);
      f[k] = ph.getF();
      vol[k] = ph.getMolarVolume();
    }
    double dfdT_p = (f[0] - f[1]) / (2.0 * h);
    double dvdT_p = (vol[0] - vol[1]) / (2.0 * h);
    double dfdV = safT.dFdV();
    double dfdT_v = dfdT_p - dfdV * dvdT_p;
    System.out.printf("F_fd_plus=%.15g F_fd_minus=%.15g v_plus=%.15g v_minus=%.15g%n",
        f[0], f[1], vol[0], vol[1]);
    System.out.printf("dFdV=%.15g dFdT_at_P=%.15g dVdT_at_P=%.15g%n", dfdV, dfdT_p, dvdT_p);
    System.out.printf("dFdT_at_V=%.15g  T*dFdT_at_V=%.15g%n", dfdT_v, t * dfdT_v);
  }

  public static void main(String[] args) {
    System.out.println("# azoth PcsaftDepartureProbe - NeqSim master's SystemPCSAFT.");
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
    probe(300.0, 50.0, new String[] {"methane"}, new double[] {1.0});
    probe(350.0, 30.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4});
  }
}
