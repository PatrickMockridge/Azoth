import neqsim.thermo.component.ComponentPCSAFT;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhasePCSAFT;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPCSAFT;

/**
 * NeqSim's PC-SAFT layers at one state, in the order the model is built.
 *
 * <p>
 * The lesson this tranche paid for twice: a port that disagrees with its oracle needs the
 * oracle's <em>intermediates</em>, not its answer, or the difference gets attributed to
 * whichever quantity was compared last. {@code CpaSweep} exists because of that for the
 * associating model; this is the same instrument for PC-SAFT, one state per invocation.
 *
 * <p>
 * Every key names the NeqSim method it comes from. The layers are the component's own
 * parameters, then the mixture's, then the packing fraction, then the hard-sphere terms,
 * then each dispersion sub-quantity, then the three Helmholtz terms, then the state.
 *
 * <p>
 * <b>The default mixing rule is the interesting one.</b> NeqSim reads no interaction
 * column for PC-SAFT unless a caller asks for the classic rule - {@code PhaseEos}'s
 * constructor sets {@code ClassicVdW} with every {@code k_ij} zero - so {@code KIJPCSAFT}
 * is carried and unused by default. The {@code kij} keys print what the phase actually
 * uses, which is the number a port has to match rather than the one in the table.
 *
 * <pre>
 * javac -proc:none -cp neqsim-f0c7436.jar PcsaftProbe.java
 * java -cp .:neqsim-f0c7436.jar PcsaftProbe [T_K P_bara name:z ...]
 * </pre>
 */
public final class PcsaftProbe {

  private PcsaftProbe() {}

  private static void probe(double t, double pBar, String[] names, double[] z) {
    SystemInterface system = new SystemPCSAFT(t, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    // **The classic rule, because NeqSim's own default cannot run a mixture.** Leaving
    // the mixing rule alone leaves `PhaseEos`'s `ClassicVdW` with a null `intparam`, and
    // `calcF1dispSumTerm` dereferences it for every pair - so a binary throws a
    // NullPointerException rather than using a zero `k_ij`. A pure fluid never asks for a
    // pair and runs either way, which is why the crash is a mixture's. And the rule that
    // works is the one whose branch, keyed on the phase class, reads `KIJPCSAFT`: the
    // zero-default story is what the source reads like, not what a caller can do.
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
      System.out.printf("component[%d] %s m=%.15g sigma=%.15g epsik=%.15g d=%.15g%n", i, names[i],
          c.getmSAFTi(), c.getSigmaSAFTi(), c.getEpsikSAFT(), c.getdSAFTi());
    }
    System.out.printf("mixture m=%.15g mmin1=%.15g d=%.15g dmean=%.15g md=%.15g%n",
        safT.getmSAFT(), safT.getMmin1SAFT(), safT.getDSAFT(), safT.getDmeanSAFT(), safT.getmdSAFT());
    System.out.printf("volumeSAFT=%.15g nSAFT=%.15g%n", safT.getVolumeSAFT(), safT.getNSAFT());
    System.out.printf("aHS=%.15g gHS=%.15g%n", safT.getAHSSAFT(), safT.getGhsSAFT());
    System.out.printf("f1vol=%.15g f1sum=%.15g f2sum=%.15g I1=%.15g I2=%.15g C1=%.15g%n",
        safT.getF1dispVolTerm(), safT.getF1dispSumTerm(), safT.getF2dispSumTerm(),
        safT.getF1dispI1(), safT.getF2dispI2(), safT.getF2dispZHC());
    System.out.printf("F_hc=%.15g F_disp1=%.15g F_disp2=%.15g F=%.15g%n", safT.F_HC_SAFT(),
        safT.F_DISP1_SAFT(), safT.F_DISP2_SAFT(), safT.getF());
    // `getMolarVolume(String)` divides a molar mass by a density that a PC-SAFT phase does
    // not populate, and returns `Infinity`; the field it does hold is in the unit the
    // phase's own `volumeSAFT` converts by `1.0e-5`.
    System.out.printf("Z=%.15g v=%.15g%n", phase.getZ(),
        phase.getMolarVolume() * 1.0e-5 / phase.getNumberOfMolesInPhase());
    for (int i = 0; i < n; i++) {
      System.out.printf("lnPhi[%d]=%.15g%n", i,
          Math.log(phase.getComponent(i).getFugacityCoefficient()));
    }
  }

  public static void main(String[] args) {
    System.out.println("# azoth PcsaftProbe - NeqSim master's SystemPCSAFT (PhasePCSAFTRahmat).");
    System.out.println("# The default mixing rule is NeqSim's own, so every k_ij this phase");
    System.out.println("# uses is zero unless a caller sets the classic rule.");
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
