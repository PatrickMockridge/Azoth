import neqsim.thermo.phase.PhasePCSAFT;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPCSAFT;
import neqsim.thermo.system.SystemPCSAFTa;

/**
 * The second volume derivative of the hard-sphere and chain term, from the two PC-SAFT
 * phase classes that disagree about it.
 *
 * <p>
 * `PhasePCSAFTa` extends `PhasePCSAFT` and inherits its `dF_HC_SAFTdVdV`, whose chain term
 * carries one `dnSAFT/dV` where the second derivative needs two. `SystemPCSAFT` builds
 * `PhasePCSAFTRahmat`, which overrides it. This prints both at one state, the chain term
 * as the base class writes it, and the factor it is short of - so the difference is a
 * number rather than an expression.
 *
 * <pre>
 * javac -proc:none -cp neqsim-f0c7436.jar PcsaftDVdVProbe.java
 * java -cp .:neqsim-f0c7436.jar PcsaftDVdVProbe
 * </pre>
 */
public final class PcsaftDVdVProbe {

  private PcsaftDVdVProbe() {}

  private static PhasePCSAFT solve(SystemInterface system) {
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    return (PhasePCSAFT) system.getPhase(0);
  }

  private static void show(String label, PhasePCSAFT phase) {
    // The chain term as the base class writes it: one `dnSAFTdV` where the second
    // derivative needs two, so the correct value is this times `dnSAFTdV`.
    double written = phase.getMmin1SAFT() * Math.pow(phase.getDgHSSAFTdN(), 2.0)
        / Math.pow(phase.getGhsSAFT(), 2.0) * phase.getDnSAFTdV();

    System.out.printf("%s%n", label);
    System.out.printf("   class                    = %s%n", phase.getClass().getName());
    System.out.printf("   dFdVdV                   = %.15g%n", phase.dFdVdV());
    System.out.printf("   mmin1, ghS, dghSdN       = %.15g, %.15g, %.15g%n",
        phase.getMmin1SAFT(), phase.getGhsSAFT(), phase.getDgHSSAFTdN());
    System.out.printf("   dnSAFTdV                 = %.15g%n", phase.getDnSAFTdV());
    System.out.printf("   chain term as written    = %.15g%n", written);
    System.out.printf("   the correct chain term   = %.15g (written * dnSAFTdV)%n",
        written * phase.getDnSAFTdV());
  }

  public static void main(String[] args) {
    System.out.println("# azoth PcsaftDVdVProbe - PhasePCSAFT's second volume derivative.");

    SystemInterface rahmat = new SystemPCSAFT(350.0, 30.0);
    rahmat.addComponent("methane", 0.6);
    rahmat.addComponent("n-butane", 0.4);
    PhasePCSAFT a = solve(rahmat);
    show("SystemPCSAFT (builds PhasePCSAFTRahmat)", a);

    SystemInterface associating = new SystemPCSAFTa(350.0, 30.0);
    associating.addComponent("methane", 0.6);
    associating.addComponent("n-butane", 0.4);
    PhasePCSAFT b = solve(associating);
    show("SystemPCSAFTa (builds PhasePCSAFTa, which inherits the base)", b);

    System.out.printf("%nthe base class's dFdVdV is %.6g%% from the Rahmat one%n",
        100.0 * (b.dFdVdV() / a.dFdVdV() - 1.0));
  }
}
