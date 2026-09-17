import neqsim.thermo.system.SystemPrEos;

/**
 * Prints NeqSim's fugacity-coefficient derivatives at the state
 * {@code crates/azoth-eos/tests/mixture.rs} differentiates: methane/n-butane
 * 0.6/0.4 at 330 K and 25 bar, under Peng-Robinson with the classic mixing rule.
 *
 * <p>The three families are {@code getdfugdx(j)} - the composition derivative, in the
 * mole-fraction variable NeqSim's second-order TP flash is written in -
 * {@code getdfugdt()} and {@code getdfugdp()}. They are printed raw, because the
 * comparison they exist for is a comparison of <em>conventions</em> as much as of
 * values: NeqSim builds its composition derivative from {@code dFdNdN} and
 * {@code dFdNdV}, which is the constant-volume frame, while the surface azoth
 * differentiates is the constant-pressure one.
 *
 * <p>{@code getdfugdx(i)} is {@code getdfugdn(i) * numberOfMolesInPhase}, so the two
 * are one derivative in two variables.
 */
public class FugacityDerivativeProbe {

  private static void report(String label, SystemPrEos system) throws Exception {
    // `init(3)` is the flag that fills `logfugcoefdT`, `logfugcoefdP` and
    // `logfugcoefdN`; `init(1)` leaves every one of them at its placeholder.
    new neqsim.thermodynamicoperations.ThermodynamicOperations(system).TPflash();
    system.init(3);
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      var p = system.getPhase(phase);
      System.out.printf("%s phase %d (%s) T=%.4f P=%.4f%n", label, phase, p.getType(),
          p.getTemperature(), p.getPressure());
      for (int i = 0; i < p.getNumberOfComponents(); i++) {
        var c = p.getComponent(i);
        System.out.printf("  %-10s x=%.17g z=%.17g dfugdt=%.17g dfugdp=%.17g%n",
            c.getName(), c.getx(), c.getz(), c.getdfugdt(), c.getdfugdp());
      }
      for (int i = 0; i < p.getNumberOfComponents(); i++) {
        System.out.printf("  %-10s lnphi=%.17g%n", p.getComponent(i).getName(),
            Math.log(p.getComponent(i).getFugacityCoefficient()));
      }
      for (int i = 0; i < p.getNumberOfComponents(); i++) {
        StringBuilder row = new StringBuilder();
        for (int j = 0; j < p.getNumberOfComponents(); j++) {
          row.append(String.format("%24.17g", p.getComponent(i).getdfugdx(j)));
        }
        System.out.printf("  dfugdx[%d] =%s%n", i, row);
      }
      for (int i = 0; i < p.getNumberOfComponents(); i++) {
        StringBuilder row = new StringBuilder();
        for (int j = 0; j < p.getNumberOfComponents(); j++) {
          row.append(String.format("%24.17g", p.getComponent(i).getdfugdn(j)));
        }
        System.out.printf("  dfugdn[%d] =%s%n", i, row);
      }
    }
  }

  public static void main(String[] args) throws Exception {
    SystemPrEos system = new SystemPrEos(330.0, 25.0);
    system.addComponent("methane", 0.6);
    system.addComponent("n-butane", 0.4);
    // `createDatabase` is what attaches `INTER.csv`'s interaction parameters; without
    // it every `kij` is zero and the state is a different fluid from azoth's.
    // The same fluid `FlashTp.java` builds: `setMixingRule(1)` is a different rule from
    // the string form, and the two do not give the same composition.
    system.setMixingRule("classic");
    system.setAttractiveTerm(1);
    report("methane/n-butane", system);
  }
}
