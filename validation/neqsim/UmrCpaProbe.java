import neqsim.thermo.component.ComponentGEInterface;
import neqsim.thermo.component.ComponentUMRCPA;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.phase.PhaseCPAInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemUMRCPAEoS;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * NeqSim's UMR-CPA layers at one state, in the order the model is built.
 *
 * <p>
 * {@code PhaseUMRCPA} is a Peng-Robinson cubic whose attraction comes from the UMR
 * universal mixing rule - {@code alpha_mix = sum_i x_i (a_i/(b_i R T) + hwfc ln gamma_i)}
 * with {@code hwfc = -1/0.53} - plus the CPA association term. The layers are the
 * component's own parameters, then the mixture's, then the UNIFAC activity coefficients the
 * mixing rule reads, then the root, then the fugacity, then the departure.
 *
 * <p>
 * <b>The attraction parameters are not the cubic's.</b> {@code Component.java:534-539}
 * overrides {@code aCPA}/{@code bCPA}/{@code associationVolume}/{@code associationEnergy}
 * from the {@code UMRCPA_*} columns when {@code UMRCPA_associating} is one, and those
 * columns are in NeqSim's internal scale ({@code a} in {@code Pa*m^6/mol^2 * 1e5},
 * {@code b} in {@code m^3/mol * 1e5}, energy in J/mol). The keys print the values the phase
 * solved with, so the conversion is checked rather than assumed.
 *
 * <p>
 * <b>The mixing rule is the only place this model differs from {@code eos.pr_cpa_phase}.</b>
 * {@code SystemUMRCPAEoS}'s constructor calls {@code setBmixType(1)} - the co-volume
 * combining rule that is the square of the mean of the square roots rather than the mean of
 * the co-volumes, {@code EosMixingRuleHandler.getbij} - but {@code setMixingRule} installs a
 * fresh handler, and the {@code bmixType} key prints what the phase actually solved with.
 *
 * <p>
 * {@code getMolarVolume("m3/mol")} is {@code getMolarMass()/getDensity("kg/m3")}, and the
 * density is only computed by {@code init(3)}, so the key prints the raw field instead:
 * NeqSim carries volume in {@code m^3/mol * 1e5}, the same scale {@code a} and {@code b} use.
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar UmrCpaProbe.java
 * java -cp .:neqsim-3.20.0.jar UmrCpaProbe [T_K P_bara name:z ...]
 * </pre>
 */
public final class UmrCpaProbe {

  /** NeqSim's gas constant, {@code ThermodynamicConstantsInterface.R}. */
  private static final double R = 8.3144621;

  private UmrCpaProbe() {}

  private static void print(String label, double value) {
    System.out.printf("%-34s %.15g%n", label, value);
  }

  private static void probe(double t, double pBar, String[] names, double[] z) {
    SystemInterface system = new SystemUMRCPAEoS(t, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("HV", "UNIFAC_UMRPRU");

    system.init(0);
    system.init(1);
    system.init(2);

    PhaseInterface phase = system.getPhase(0);
    EosMixingRulesInterface rule = (EosMixingRulesInterface) phase.getMixingRule();
    int n = phase.getNumberOfComponents();

    System.out.printf("%n# %s at T=%g K, P=%g bara%n", String.join("/", names), t, pBar);
    System.out.printf("phase                      %s%n", phase.getClass().getSimpleName());
    System.out.printf("%-34s %d%n", "bmixType", rule.getBmixType());
    System.out.printf("%-34s %d%n", "totalAssociationSites",
        ((PhaseCPAInterface) phase).getTotalNumberOfAccociationSites());

    for (int i = 0; i < n; i++) {
      ComponentUMRCPA c = (ComponentUMRCPA) phase.getComponent(i);
      System.out.printf("-- component %d %s x=%.15g%n", i, c.getComponentName(),
          phase.getComponent(i).getx());
      System.out.printf("%-34s %d%n", "attractiveTermNumber", c.getAttractiveTermNumber());
      print("alpha_T", c.getAttractiveTerm().alpha(t));
      print("aT_reduced", c.getaT());
      print("b_covolume", c.getb());
      // `calca`/`calcb` are `ComponentUMRCPA`'s overrides: `aCPA`/`bCPA` for an
      // associating component, the cubic's own values otherwise.
      print("calca", c.calca());
      print("calcb", c.calcb());
      print("sites", c.getNumberOfAssociationSites());
      System.out.printf("%-34s %s%n", "associationScheme", c.getAssociationScheme());
      print("associationEnergy", c.getAssociationEnergy());
      print("associationVolume", c.getAssociationVolume());
      print("qPure_aT_over_bRT", c.getaT() / (c.getb() * R * t));
      print("qPure_calca_over_calcbRT", c.calca() / (c.calcb() * R * t));
      StringBuilder mc = new StringBuilder();
      double[] params = c.getMatiascopemanParamsUMRCPA();
      if (params != null) {
        for (double value : params) {
          mc.append(String.format("%.15g ", value));
        }
      }
      System.out.printf("%-34s %s%n", "matiascopemanParamsUMRCPA", mc.toString().trim());
    }

    // The UNIFAC activity coefficients the mixing rule reads, and their temperature
    // derivative - which is what an enthalpy departure needs and a fugacity does not.
    PhaseInterface ge = rule.getGEPhase();
    if (ge != null) {
      // The handler's own call, `EosMixingRuleHandler.init`, at the phase's initType -
      // which is what makes the derivatives available: at initType 1 the GE component
      // leaves `dlngammadt` at zero.
      ge.init(phase.getNumberOfMolesInPhase() / phase.getBeta(), n, phase.getInitType(),
          phase.getType(), phase.getBeta());
      for (int i = 0; i < n; i++) {
        ComponentGEInterface g = (ComponentGEInterface) ge.getComponent(i);
        print("lnGamma[" + i + "]", g.getLnGamma());
        print("dlnGammadT[" + i + "]", g.getLnGammadt());
      }
    }

    print("A_attraction", phase.getA());
    print("B_covolume", phase.getB());
    // `calcA` is `n * B * R * T * alpha_mix`, so this inverts it rather than guessing.
    print("alpha_mix", phase.getA() / (phase.getNumberOfMolesInPhase() * phase.getB() * R * t));
    print("Z", phase.getZ());
    print("molarVolume_internal_1e5_m3_per_mol", phase.getMolarVolume());
    print("molarVolume_SI_from_ZRT_over_P", phase.getZ() * R * t / (pBar * 1.0e5));
    print("HresTP_J_per_mol", phase.getHresTP() / phase.getNumberOfMolesInPhase());
    print("SresTP_J_per_molK", phase.getSresTP() / phase.getNumberOfMolesInPhase());
    print("gibbsEnergy_J_per_mol", phase.getGibbsEnergy() / phase.getNumberOfMolesInPhase());
    for (int i = 0; i < n; i++) {
      print("lnPhi[" + i + "]", Math.log(phase.getComponent(i).getFugacityCoefficient()));
    }

    // The flash, so the two-phase state and each phase's own layers can be compared too.
    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    operations.TPflash();
    system.init(3);
    System.out.printf("flashPhases %d%n", system.getNumberOfPhases());
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      PhaseInterface flashed = system.getPhase(i);
      StringBuilder composition = new StringBuilder();
      for (int j = 0; j < system.getNumberOfComponents(); j++) {
        composition.append(String.format("%.15g ", flashed.getComponent(j).getx()));
      }
      System.out.printf("flashPhase %d %s beta=%.15g Z=%.15g v=%.15g x=%s%n", i, flashed.getType(),
          flashed.getBeta(), flashed.getZ(), flashed.getZ() * R * t / (pBar * 1.0e5),
          composition.toString().trim());
      for (int j = 0; j < system.getNumberOfComponents(); j++) {
        System.out.printf("flashPhase %d lnPhi[%d]=%.15g%n", i, j,
            Math.log(flashed.getComponent(j).getFugacityCoefficient()));
      }
    }
  }

  public static void main(String[] args) {
    System.out.println("# azoth UmrCpaProbe - NeqSim 3.20.0's SystemUMRCPAEoS (PhaseUMRCPA).");
    System.out.println("# The mixing rule is the lifecycle test's: HV with UNIFAC_UMRPRU.");
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
    // The lifecycle test's nine states, methane/water 0.98/0.02.
    for (double temperature : new double[] {283.15, 298.15, 313.15}) {
      for (double pressure : new double[] {30.0, 70.0, 120.0}) {
        probe(temperature, pressure, new String[] {"methane", "water"},
            new double[] {0.98, 0.02});
      }
    }
  }
}
