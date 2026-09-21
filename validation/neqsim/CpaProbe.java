import neqsim.thermo.component.ComponentSrkCPA;
import neqsim.thermo.phase.PhaseCPAInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseSrkCPA;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkCPA;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Prints the association kernel's internals at one state, so azoth's port of it can be
 * compared with the class it came from rather than with itself.
 *
 * <p>
 * The state is a single liquid phase - water and methanol at 300 K and 100 bar - because a
 * two-phase state would make "which phase" the first question and the association kernel
 * the second.
 *
 * <p>
 * Usage: {@code java -cp .:neqsim-f0c7436.jar CpaProbe [T_K] [P_bara] [n_water]}
 */
public final class CpaProbe {

  private CpaProbe() {}

  private static void print(String label, double value) {
    System.out.printf("%-34s %.15g%n", label, value);
  }

  public static void main(String[] args) {
    double temperature = args.length > 0 ? Double.parseDouble(args[0]) : 300.0;
    double pressure = args.length > 1 ? Double.parseDouble(args[1]) : 100.0; // bara
    double nWater = args.length > 2 ? Double.parseDouble(args[2]) : 0.6;

    SystemInterface system = new SystemSrkCPA(temperature, pressure);
    system.addComponent("water", nWater);
    system.addComponent("methanol", 1.0 - nWater);
    system.setMixingRule(10);

    system.init(0);
    system.init(1);

    System.out.printf("components %d%n", system.getNumberOfComponents());
    System.out.printf("phases     %d%n", system.getNumberOfPhases());
    print("T_K", system.getTemperature());
    print("P_bara", system.getPressure());

    PhaseInterface phase = system.getPhase(0);
    System.out.printf("phase0     %s%n", phase.getClass().getSimpleName());
    print("numberOfMolesInPhase", phase.getNumberOfMolesInPhase());
    print("totalVolume_m3", phase.getTotalVolume());
    print("molarVolume_m3_per_mol", phase.getMolarVolume());
    print("B_covolume", phase.getB());
    print("gcpa_g_at_contact", ((PhaseCPAInterface) phase).getGcpa());
    print("gcpav_dlng_dV", ((PhaseCPAInterface) phase).getGcpav());
    print("totalAssociationSites", ((PhaseCPAInterface) phase).getTotalNumberOfAccociationSites());

    // The reduced parameters the cubic actually solved with. The composition is (1, 0)
    // for a pure component, so this is a direct read rather than a mixture average.
    print("phase_getA", phase.getA());
    print("phase_getB", phase.getB());
    print("phase_getZ", phase.getZ());
    print("molarVolume", phase.getMolarVolume());

    for (int i = 0; i < system.getNumberOfComponents(); i++) {
      ComponentSrkCPA c = (ComponentSrkCPA) phase.getComponent(i);
      System.out.printf("-- component %d %s%n", i, c.getComponentName());
      print("moles", c.getNumberOfMolesInPhase());
      print("b_covolume", c.getb());
      print("sites", c.getNumberOfAssociationSites());
      System.out.printf("%-34s %s%n", "associationScheme", c.getAssociationScheme());
      print("associationEnergy", c.getAssociationEnergy());
      print("associationVolume", c.getAssociationVolume());
      print("fugacityCoefficient", c.getFugacityCoefficient());
      print("calc_lngi", c.calc_lngi(phase));
      print("calca", c.calca());
      print("calcb", c.calcb());
      print("getVolumeCorrection", c.getVolumeCorrection());
      print("getBi", c.getBi());
      print("getAi", c.getAi());
      for (int j = 0; j < c.getNumberOfAssociationSites(); j++) {
        print("xsite[" + j + "]", c.getXsite()[j]);
      }
    }

    // The Helmholtz energy and the fugacity term, which are what the port must reproduce.
    print("FCPA_helmholtz_over_RT", ((PhaseSrkCPA) phase).FCPA());
    print("dFCPAdV", ((PhaseSrkCPA) phase).dFCPAdV());
    print("dFCPAdT", ((PhaseSrkCPA) phase).dFCPAdT());
    for (int i = 0; i < system.getNumberOfComponents(); i++) {
      ComponentSrkCPA c = (ComponentSrkCPA) phase.getComponent(i);
      print("dFCPAdN[" + i + "]", c.dFCPAdN(phase, system.getNumberOfComponents(), temperature, pressure));
    }
    for (int i = 0; i < system.getNumberOfComponents(); i++) {
      ComponentSrkCPA c = (ComponentSrkCPA) phase.getComponent(i);
      StringBuilder row = new StringBuilder();
      for (int j = 0; j < system.getNumberOfComponents(); j++) {
        row.append(String.format("%.15g ", c.dFCPAdNdN(j, phase, system.getNumberOfComponents(), temperature, pressure)));
      }
      System.out.printf("%-33s %s%n", "dFCPAdNdN[" + i + "]", row.toString().trim());
    }
    // The flash, so the phase count and the split can be compared as well as the kernel.
    // A state that is single-phase at one temperature and two at another is the case an
    // oracle is needed for: azoth has no way to know which it should be.
    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    operations.TPflash();
    System.out.printf("flashPhases %d%n", system.getNumberOfPhases());
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      PhaseInterface flashed = system.getPhase(i);
      StringBuilder composition = new StringBuilder();
      for (int j = 0; j < system.getNumberOfComponents(); j++) {
        composition.append(String.format("%.15g ", flashed.getComponent(j).getx()));
      }
      System.out.printf(
          "flashPhase %d %s beta=%.15g Z=%.15g molarVolume=%.15g x=%s%n",
          i,
          flashed.getType(),
          flashed.getBeta(),
          flashed.getZ(),
          flashed.getMolarVolume(),
          composition);
    }


  }
}
