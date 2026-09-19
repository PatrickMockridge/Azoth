import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

/**
 * `(dH/dT)_P` at a two-phase state, against the `Cp` NeqSim's PH-flash Jacobian puts there.
 *
 * <p>
 * {@code SysNewtonRhapsonPHflash} carries two schemes and both of them use the same
 * derivative for the temperature entry. Its single-phase loop takes {@code dValdT =
 * system.getCp()} outright, and its two-phase Jacobian takes {@code Ett = -system.getCp() /
 * R} for the enthalpy row at constant pressure and composition. **In one phase that
 * derivative is right and in two it is not**: across a phase boundary the enthalpy carries
 * the latent heat of the phase that is appearing, so
 *
 * <pre>
 * (dH/dT)_P = Cp + (H_vap - H_liq) (dbeta/dT)_P
 * </pre>
 *
 * and the second term is not small where the phase fractions are moving.
 *
 * <p>
 * This prints both, at a two-phase state and at a single-phase one, because the claim is
 * not "the derivative is wrong" - it is "wrong <em>across a phase boundary</em>", and a
 * probe that only ran the two-phase state could not tell that from a library whose Cp is
 * wrong everywhere. The state is the one {@code eos.ph_flash}'s own case uses.
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar PhFlashJacobianProbe.java
 * java -cp .:neqsim-3.20.0.jar PhFlashJacobianProbe
 * </pre>
 */
public final class PhFlashJacobianProbe {

  private PhFlashJacobianProbe() {}

  /** The quality and the total enthalpy at a temperature, from a fresh TP flash. */
  private static double[] stateAt(SystemInterface system, double t) {
    system.setTemperature(t);
    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    operations.TPflash();
    system.init(3);
    return new double[] {system.getBeta(), system.getEnthalpy(), system.getCp(),
        system.getEntropy()};
  }

  private static void probe(double t, double pBar, String[] names, double[] z, double delta) {
    SystemInterface system = new SystemSrkEos(t, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("classic");

    // A fresh flash fixes the number of phases at this state; the perturbations are taken
    // on that phase structure, so a state that split does not silently become single.
    ThermodynamicOperations initial = new ThermodynamicOperations(system);
    initial.TPflash();
    system.init(3);
    int phases = system.getNumberOfPhases();

    double[] here = stateAt(system, t);
    double[] up = stateAt(system, t + delta);
    double[] down = stateAt(system, t - delta);

    double dHdT = (up[1] - down[1]) / (2.0 * delta);
    double dSdT = (up[3] - down[3]) / (2.0 * delta);
    double dbeta = (up[0] - down[0]) / (2.0 * delta);

    System.out.printf("%n# %s at T=%g K, P=%g bara: %d phase(s)%n", String.join("/", names), t,
        pBar, phases);
    System.out.printf("   beta                     = %.15g%n", here[0]);
    System.out.printf("   dbeta/dT                 = %.15g 1/K%n", dbeta);
    System.out.printf("   (dH/dT)_P  measured      = %.15g J/(mol K)%n", dHdT);
    System.out.printf("   Cp        the Jacobian   = %.15g J/(mol K)%n", here[2]);
    System.out.printf("   latent contribution      = %.15g J/(mol K)%n", dHdT - here[2]);
    System.out.printf("   the Jacobian is out by   = %.6g%%%n",
        100.0 * (here[2] / dHdT - 1.0));
    // The S-flash's Jacobian is the same statement with `Cp/T`, and the entropy carries
    // the same phase-change term, so it is measured here rather than assumed.
    System.out.printf("   (dS/dT)_P  measured      = %.15g J/(mol K^2)%n", dSdT);
    System.out.printf("   Cp/T      the Jacobian   = %.15g J/(mol K^2)%n", here[2] / t);
    System.out.printf("   the Jacobian is out by   = %.6g%%%n",
        100.0 * ((here[2] / t) / dSdT - 1.0));
  }

  public static void main(String[] args) {
    System.out.println("# azoth PhFlashJacobianProbe - NeqSim 3.20.0's PH-flash Jacobian.");
    System.out.println("# Two-phase state first, then a single-phase control: the claim is");
    System.out.println("# that dH/dT = Cp is wrong ACROSS a phase boundary, not everywhere.");

    // `eos.ph_flash`'s own case: methane/n-butane at 2 MPa, 300 K, where it splits.
    probe(300.0, 20.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, 0.5);
    // The control: the same fluid and pressure above the dew point, one phase.
    probe(400.0, 20.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, 0.5);
  }
}
