import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSAFTVRMie;

/**
 * NeqSim's SAFT-VR-Mie flash, reached the way a caller reaches it.
 *
 * <p>
 * <b>The point of the probe is the route, not the answer.</b> {@code ThermodynamicOperations}
 * dispatches on {@code instanceof SystemSAFTVRMie} to {@code TPflashSAFT}, which is the one
 * place in the library where the SAFT-VR-Mie family is <em>not</em> flashed by the ordinary
 * {@code TPflash} - the generic flash reaches its K-values through a Newton solve on the
 * derivative surface, and this model's is taken by successive substitution instead. So a port
 * that runs the generic flash here answers a different algorithm, and this prints what the
 * library's own algorithm gives.
 *
 * <p>
 * The scheme, from {@code TPflashSAFT}: Wilson K-values to seed, then successive substitution
 * to a relative K change of {@code 1e-6} or fifty iterations, with Rachford-Rice inside for
 * the vapour fraction, and **each phase's fugacity coefficients from a separate single-phase
 * system built at that trial composition**. The vapour fraction is bounded at
 * {@code 1e-10}; outside it the flash reports a single phase, chosen between gas and liquid by
 * whichever of their Gibbs energies is lower - so the single-phase branch is a comparison of
 * two explicit solves rather than a stability test.
 *
 * <p>
 * <b>And the class's answer depends on whether the caller initialised the system first.</b>
 * {@code init(0)} before the flash finds the split; without it the flash reports a single gas
 * phase - which is what {@code ThermodynamicOperations.TPflash()} itself does, because it does
 * not initialise. So this runs all three: the class without a prior {@code init(0)}, the class
 * with one, and the dispatched call. The first and third are the same failure and the second is
 * the algorithm's answer, so the numbers a port is checked against are the second column's.
 *
 * <pre>
 * javac -proc:none -cp neqsim-f0c7436.jar SaftVrMieFlashProbe.java
 * java -cp .:neqsim-f0c7436.jar SaftVrMieFlashProbe [T_K P_bara name:z ...]
 * </pre>
 */
public final class SaftVrMieFlashProbe {

  private SaftVrMieFlashProbe() {}

  private static void probe(double t, double pBar, String[] names, double[] z) {
    SystemInterface system = new SystemSAFTVRMie(t, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("classic");

    System.out.printf("%n# %s at T=%g K, P=%g bara%n", String.join("/", names), t, pBar);
    System.out.printf("doMultiPhaseCheck=%s%n", system.doMultiPhaseCheck());

    try {
      new ThermodynamicOperations(system).TPflash();
    } catch (Exception e) {
      System.out.printf("TPflash THREW %s: %s%n", e.getClass().getSimpleName(), e.getMessage());
      return;
    }

    int phases = system.getNumberOfPhases();
    System.out.printf("numberOfPhases=%d beta=%.15g%n", phases, system.getBeta());
    System.out.printf("phaseType[0]=%s%n", system.getPhase(0).getType());
    for (int i = 0; i < phases; i++) {
      PhaseInterface phase = system.getPhase(i);
      System.out.printf("phase[%d] type=%s beta=%.15g Z=%.15g%n", i, phase.getType(),
          phase.getBeta(), phase.getZ());
    }
    for (int i = 0; i < phases; i++) {
      PhaseInterface phase = system.getPhase(i);
      System.out.printf("phase[%d] x = [", i);
      for (int c = 0; c < names.length; c++) {
        System.out.printf("%.15g%s", phase.getComponent(c).getx(),
            c + 1 < names.length ? ", " : "");
      }
      System.out.println("]");
    }
    for (int i = 0; i < phases; i++) {
      PhaseInterface phase = system.getPhase(i);
      System.out.printf("phase[%d] lnPhi = [", i);
      for (int c = 0; c < names.length; c++) {
        System.out.printf("%.15g%s",
            Math.log(phase.getComponent(c).getFugacityCoefficient()),
            c + 1 < names.length ? ", " : "");
      }
      System.out.println("]");
    }

    // **The same class, two ways.** The count is the whole difference.
    for (boolean initialise : new boolean[] {false, true}) {
      SystemInterface other = build(t, pBar, names, z);
      if (initialise) {
        other.init(0);
      }
      new neqsim.thermodynamicoperations.flashops.TPflashSAFT(other).run();
      System.out.printf("TPflashSAFT, %-18s phases=%d beta=%.15g%n",
          initialise ? "init(0) first" : "no init(0) first", other.getNumberOfPhases(),
          other.getBeta());
      // The converged compositions and K-values: what a port is checked against when
      // the class is given a system that was initialised.
      for (int i = 0; i < other.getNumberOfPhases(); i++) {
        PhaseInterface phase = other.getPhase(i);
        System.out.printf("   %s type=%s x = [", initialise ? "init(0)" : "no-init",
            phase.getType());
        for (int c = 0; c < names.length; c++) {
          System.out.printf("%.15g%s", phase.getComponent(c).getx(),
              c + 1 < names.length ? ", " : "");
        }
        System.out.println("]");
      }
      if (other.getNumberOfPhases() == 2) {
        System.out.printf("   K = [");
        for (int c = 0; c < names.length; c++) {
          System.out.printf("%.15g%s", other.getPhase(0).getComponent(c).getK(),
              c + 1 < names.length ? ", " : "");
        }
        System.out.println("]");
      }
    }
  }

  private static SystemInterface build(double t, double pBar, String[] names, double[] z) {
    SystemInterface system = new SystemSAFTVRMie(t, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("classic");
    return system;
  }

  public static void main(String[] args) {
    System.out.println("# azoth SaftVrMieFlashProbe - NeqSim master's TPflashSAFT, by dispatch.");
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
    probe(350.0, 30.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4});
    probe(250.0, 30.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4});
    probe(200.0, 30.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4});
  }
}
