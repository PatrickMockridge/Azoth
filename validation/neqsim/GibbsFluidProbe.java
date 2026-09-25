// One state of the ammonia fluid, printed by NeqSim itself.
//
//     javac -proc:none -cp neqsim-f0c7436.jar GibbsFluidProbe.java
//     java -cp .:neqsim-f0c7436.jar GibbsFluidProbe
//
// **The question this answers.** azoth's cubic agrees with the capture's `phase0_z` to sixteen
// digits at the state the ammonia row converged to, and its fugacity coefficients disagree by a
// factor of two - which cannot both be true of one equation of state unless one of the two
// numbers is not what it looks like. So this builds that exact state in NeqSim and prints what
// NeqSim thinks, rather than what the reactor last cached.
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

public class GibbsFluidProbe {
  public static void main(String[] args) {
    String[] names = { "hydrogen", "nitrogen", "ammonia" };
    double[] moles = { 0.029289705645501173, 0.009763235215167055, 0.48047353217544875 };

    SystemInterface fluid = new SystemPrEos(450.0, 300.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule(2);
    fluid.init(3);

    System.out.println("phases=" + fluid.getNumberOfPhases());
    System.out.println("phase0_type=" + fluid.getPhase(0).getPhaseTypeName());
    System.out.println("phase0_z=" + fluid.getPhase(0).getZ());
    System.out.println("phase0_moles=" + fluid.getPhase(0).getNumberOfMolesInPhase());
    for (int i = 0; i < names.length; i++) {
      System.out.println("x_" + names[i] + "=" + fluid.getPhase(0).getComponent(i).getx());
      System.out.println("phi_" + names[i] + "="
          + fluid.getPhase(0).getComponent(i).getFugacityCoefficient());
    }
    // The same state, stated as mole fractions rather than moles, to see whether the answer is
    // a function of composition alone.
    SystemInterface fractions = new SystemPrEos(450.0, 300.0);
    for (int i = 0; i < names.length; i++) {
      fractions.addComponent(names[i], fluid.getPhase(0).getComponent(i).getx());
    }
    fractions.setMixingRule(2);
    fractions.init(3);
    System.out.println("by_fraction_phase0_z=" + fractions.getPhase(0).getZ());
    for (int i = 0; i < names.length; i++) {
      System.out.println("by_fraction_phi_" + names[i] + "="
          + fractions.getPhase(0).getComponent(i).getFugacityCoefficient());
    }
  }
}
