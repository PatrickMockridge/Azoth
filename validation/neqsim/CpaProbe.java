import neqsim.thermo.component.ComponentSrkCPA;
import neqsim.thermo.phase.PhaseCPAInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseSrkCPA;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkCPA;

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
 * Usage: {@code java -cp .:neqsim-3.20.0.jar CpaProbe}
 */
public final class CpaProbe {

  private CpaProbe() {}

  private static void print(String label, double value) {
    System.out.printf("%-34s %.15g%n", label, value);
  }

  public static void main(String[] args) {
    double temperature = 300.0;
    double pressure = 100.0; // bara
    double nWater = args.length > 0 ? Double.parseDouble(args[0]) : 0.6;

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

    for (int i = 0; i < system.getNumberOfComponents(); i++) {
      ComponentSrkCPA c = (ComponentSrkCPA) phase.getComponent(i);
      System.out.printf("-- component %d %s%n", i, c.getComponentName());
      print("moles", c.getNumberOfMolesInPhase());
      print("b_covolume", c.getb());
      print("sites", c.getNumberOfAssociationSites());
      System.out.printf("%-34s %s%n", "associationScheme", c.getAssociationScheme());
      print("associationEnergy", c.getAssociationEnergy());
      print("associationVolume", c.getAssociationVolume());
      print("calc_lngi", c.calc_lngi(phase));
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
  }
}
