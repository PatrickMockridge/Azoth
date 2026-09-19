import neqsim.thermo.component.ComponentSAFTVRMie;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseSAFTVRMie;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSAFTVRMie;

/**
 * NeqSim's SAFT-VR-Mie layers at one state, in the order the model is built.
 *
 * <p>
 * The third instrument of its kind, after {@code CpaSweep} and {@code PcsaftProbe}, and for
 * the same reason: a port that disagrees with its oracle needs the oracle's
 * <em>intermediates</em>, not its answer, or the difference gets attributed to whichever
 * quantity was compared last.
 *
 * <p>
 * Every key names the NeqSim method it comes from. The layers are each component's
 * parameters, then the mixture's, then the packing fraction, then the hard-sphere terms,
 * then the three dispersion terms of the Mie expansion with their derivatives, then the
 * Helmholtz energy, then the state.
 *
 * <p>
 * <b>This phase reads no interaction parameter.</b> There is no
 * {@code getBinaryInteractionParameter} call anywhere in {@code PhaseSAFTVRMie}, so unlike
 * PC-SAFT it needs no mixing rule to run a mixture - its dispersion is built from the
 * component sets alone. The probe still sets the classic rule, so that what it prints is
 * what a caller running the model the same way would get.
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar SaftVrMieProbe.java
 * java -cp .:neqsim-3.20.0.jar SaftVrMieProbe [T_K P_bara name:z ...]
 * </pre>
 */
public final class SaftVrMieProbe {

  private SaftVrMieProbe() {}

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
      System.out.printf(
          "component[%d] %s m=%.15g lambdaR=%.15g lambdaA=%.15g sigma=%.15g epsik=%.15g d=%.15g%n",
          i, names[i], c.getmSAFTi(), c.getLambdaRSAFTVRMie(), c.getLambdaASAFTVRMie(),
          c.getSigmaSAFTi(), c.getEpsikSAFT(), c.getdSAFTi());
    }
    System.out.printf("mixture m=%.15g mmin1=%.15g d=%.15g%n", mie.getmSAFT(),
        mie.getMmin1SAFT(), mie.getDSAFT());
    System.out.printf("volumeSAFT=%.15g nSAFT=%.15g%n", mie.getVolumeSAFT(), mie.getNSAFT());
    System.out.printf("aHS=%.15g gHS=%.15g dDSAFTdT=%.15g%n", mie.getAHSSAFT(), mie.getGhsSAFT(),
        mie.getdDSAFTdT());

    // **The chain RDF's own layers**, because the mixture's `gHS` is not the
    // Carnahan-Starling contact value: with `mmin1SAFT = sum_i x_i (m_i - 1) > 0` the class
    // replaces it with `exp(sum_i w_i ln g_Mie_ii / W)` over `w_i = x_i (m_i - 1)`, each
    // component's `g_Mie_ii` being `gHS0_i exp(alpha tau (g1_i + tau g2_i)/gHS0_i)`. So a
    // port that stops at the contact value reproduces a pure fluid - where `mmin1` is zero
    // and the block is skipped, which is why methane's `gHS` above is exactly the CS
    // value - and is wrong for every mixture.
    double blend = mie.getGhsSAFT();
    if (mie.getMmin1SAFT() <= 1.0e-10) {
      System.out.println("chain: skipped, mmin1SAFT = " + mie.getMmin1SAFT()
          + ", so gHS above is the Carnahan-Starling contact value");
    } else {
      // The blend fraction has no getter, so it is recovered from the printed `gHS` and
      // the weight-averaged form of the line above - printing it beside the g1 and g2 it
      // came from is what makes the recovery checkable rather than asserted.
      double logBlend = Math.log(blend);
      double weightedLog0 = 0.0;
      double weightedSlope = 0.0;
      for (int i = 0; i < n; i++) {
        ComponentSAFTVRMie c = (ComponentSAFTVRMie) phase.getComponent(i);
        double w = phase.getComponent(i).getNumberOfMolesInPhase() / phase.getNumberOfMolesInPhase()
            * (c.getmSAFTi() - 1.0);
        double sigmaI = c.getSigmaSAFTi();
        double dI = c.getdSAFTi();
        double x0 = (dI > 0) ? sigmaI / dI : 1.0;
        double epsOvKT = c.getEpsikSAFT() / t;
        double lr = c.getLambdaRSAFTVRMie();
        double la = c.getLambdaASAFTVRMie();
        double cMie = ComponentSAFTVRMie.calcMiePrefactor(lr, la);
        double zetaSt = mie.getNSAFT() * x0 * x0 * x0;
        double gHS0 = PhaseSAFTVRMie.calcGHS_x0(mie.getNSAFT(), x0);
        double g1 = PhaseSAFTVRMie.calcG1Chain(mie.getNSAFT(), lr, la, cMie, x0);
        double g2 = PhaseSAFTVRMie.calcG2Chain(mie.getNSAFT(), zetaSt, lr, la, epsOvKT, cMie, x0);
        weightedLog0 += w * Math.log(gHS0);
        weightedSlope += w * epsOvKT * (g1 + epsOvKT * g2) / gHS0;
        System.out.printf("chain[%d] w=%.15g x0=%.15g zetaSt=%.15g cMie=%.15g alphaLafitte=%.15g%n",
            i, w, x0, zetaSt, cMie, PhaseSAFTVRMie.calcMieAlpha(lr, la));
        System.out.printf(
            "chain[%d] aS1Bare(la)=%.15g B(la)=%.15g aS1Bare(lr)=%.15g B(lr)=%.15g%n", i,
            PhaseSAFTVRMie.calcAS1Bare(mie.getNSAFT(), la),
            PhaseSAFTVRMie.calcBBare(mie.getNSAFT(), la, x0),
            PhaseSAFTVRMie.calcAS1Bare(mie.getNSAFT(), lr),
            PhaseSAFTVRMie.calcBBare(mie.getNSAFT(), lr, x0));
        System.out.printf("chain[%d] khs=%.15g gHS0=%.15g g1=%.15g g2=%.15g%n", i,
            PhaseSAFTVRMie.calcKHS(mie.getNSAFT()), gHS0, g1, g2);
      }
      double mmin1 = mie.getMmin1SAFT();
      System.out.printf("chain alpha: weightedLog0=%.15g weightedSlope=%.15g ln(gHS)=%.15g%n",
          weightedLog0 / mmin1, weightedSlope / mmin1, logBlend);
      double alpha = weightedSlope == 0.0 ? 0.0
          : (logBlend - weightedLog0 / mmin1) / (weightedSlope / mmin1);
      System.out.printf("chain alphaRecovered=%.15g from gHS=%.15g%n", alpha, blend);
    }
    System.out.printf("A1Disp=%.15g A2Disp=%.15g A3Disp=%.15g%n", mie.getA1Disp(), mie.getA2Disp(),
        mie.getA3Disp());
    System.out.printf("dA1dEta=%.15g dA2dEta=%.15g dA3dEta=%.15g%n", mie.getDa1DispDeta(),
        mie.getDa2DispDeta(), mie.getDa3DispDeta());
    System.out.printf("F_hc=%.15g F_disp=%.15g F=%.15g%n", mie.F_HC_SAFT(), mie.F_DISP_SAFT(),
        mie.getF());
    System.out.printf("Z=%.15g v=%.15g%n", phase.getZ(), phase.getVolume()
        / phase.getNumberOfMolesInPhase());
    for (int i = 0; i < n; i++) {
      System.out.printf("lnPhi[%d]=%.15g%n", i,
          Math.log(phase.getComponent(i).getFugacityCoefficient()));
    }
  }

  public static void main(String[] args) {
    System.out.println("# azoth SaftVrMieProbe - NeqSim 3.20.0's SystemSAFTVRMie (PhaseSAFTVRMie).");
    System.out.println("# The phase reads no interaction parameter, so a mixture runs with");
    System.out.println("# the component sets alone; no k_ij appears in these layers.");
    if (args.length >= 5) {
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
