// What the three smaller GE electrolyte models answer, before any of them is ported.
//
//     javac -proc:none -cp neqsim-f0c7436.jar GeElectrolyteProbe.java
//     java -cp .:neqsim-f0c7436.jar GeElectrolyteProbe
//
// `PhaseKentEisenberg`, `PhaseDesmukhMather` and `PhaseDuanSun` are a tranche each about
// the size of one Pitzer branch, and they disagree with each other in ways the Pitzer work
// already named: **each states its own reference state and its own insoluble-ion
// constant**, and the three activity-coefficient expressions have nothing in common.
//
// Four things the probe answers that the source does not make obvious:
//
//   1. **Kent-Eisenberg's activity coefficient is identically one** - `PhaseKentEisenberg`
//      overrides `getActivityCoefficient` to `1.0` - so its whole content is the two
//      reference states and the `1e8` an ion gets.
//   2. **Desmukh-Mather's `getSolventWeight` selects by `referenceStateType`**, not by the
//      name `water` the way `PhasePitzer`'s does. Two models, two rules, one quantity.
//   3. **Desmukh-Mather's `gamma` is a mole-fraction ratio**: `m_i M_solvent exp(lngamma) /
//      x_i`, with `exp(lngamma)` as the fallback when that is not finite or not positive.
//   4. **Duan-Sun's salinity is `sum m_i` over the ions**, not `1/2 sum m_i z_i^2`, and its
//      correlation is for three named gases only - every other component gets `gamma = 1`.
//
// # The one-component GE phase, and where it breaks
//
// `SystemDuanSun`'s constructor adds a mole of CO2 and its `addComponent` refuses every
// other name, so **every state it can hold is one component and water-free**. That is the
// shape `PhaseGE.getActivityCoefficientInfDil` cannot take: it clones the phase, passes
// the clone's own component count to `init`, and then a **literal `2`** to
// `getExcessGibbsEnergy` - which every `PhaseGE` subclass loops on. The clone has one
// component, so the subclass reads a component that does not exist.
//
// This is not Duan-Sun's defect and `oneComponentGe` below shows that: a one-component
// `SystemNRTL` or `SystemGEWilson` holding a `solute`-reference component crashes the same
// way, with a different exception per subclass. Duan-Sun's own defect is separate - the
// system holds one mole of CO2 whatever fluid it is handed, because `setModel` catches the
// guard's exception and returns the half-built system, and with the count fixed the
// correlation still contributes nothing: salinity is zero, so `gamma` is exactly one.

import neqsim.thermo.component.ComponentDesmukhMather;
import neqsim.thermo.component.ComponentGeDuanSun;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseDesmukhMather;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemDesmukhMather;
import neqsim.thermo.system.SystemDuanSun;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemKentEisenberg;

public class GeElectrolyteProbe {

  private static final double T = 313.15;
  private static final double P = 5.0;

  private static SystemInterface build(SystemInterface system, String[] names, double[] z) {
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    return system;
  }

  private static PhaseInterface phase(SystemInterface system) {
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      if (system.getPhase(i).getClass().getSimpleName().startsWith("PhaseKentEisenberg")
          || system.getPhase(i).getClass().getSimpleName().startsWith("PhaseDesmukhMather")
          || system.getPhase(i).getClass().getSimpleName().startsWith("PhaseDuanSun")) {
        return system.getPhase(i);
      }
    }
    return system.getPhase(1);
  }

  private static void report(String label, SystemInterface system, String[] names, double[] z,
      boolean fugacity) {
    build(system, names, z);
    PhaseInterface phase = phase(system);
    System.out.printf("%n=== %s ===%n", label);
    System.out.printf("  phase class = %s%n", phase.getClass().getSimpleName());
    // **The system's own state**, not this class's constants: `report` is called with
    // systems built at other temperatures, and printing the constants labelled those
    // blocks with the wrong temperature while their numbers were right.
    System.out.printf("  T = %.2f K   P = %.2f bara   moles = %.10g%n",
        system.getTemperature(), system.getPressure(), phase.getNumberOfMolesInPhase());
    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      ComponentInterface component = phase.getComponent(i);
      // The gamma `fugcoef` builds on, from the phase's own dispatch - which is what
      // `PhaseKentEisenberg` overrides to `1.0` for every component.
      component.fugcoef(phase);
      double gamma = phase.getActivityCoefficient(i, 0);
      System.out.printf(
          "    %-9s z = %+.0f  x = %.12g  gamma = %.15g  ln gamma = %.15g  molality = "
              + "%.12g  refstate = %s",
          component.getComponentName(), component.getIonicCharge(), component.getx(), gamma,
          Math.log(gamma), component.getMolality(phase), component.getReferenceStateType());
      if (fugacity) {
        System.out.printf("  phi = %.15g", component.getFugacityCoefficient());
      }
      System.out.println();
    }
  }

  private static void oneComponentGe(String label, SystemInterface system) {
    system.addComponent("CO2", 1.0);
    system.setMixingRule("classic");
    System.out.printf("  %-15s ", label);
    try {
      system.init(0);
      system.init(1);
      System.out.printf("ok, gamma = %.12g%n",
          system.getPhase(1).getActivityCoefficient(0, 0));
    } catch (Throwable error) {
      StackTraceElement[] frames = error.getStackTrace();
      String site = frames.length > 0 ? frames[0].toString() : "?";
      System.out.printf("%s: %s%n      at %s%n", error.getClass().getSimpleName(),
          error.getMessage(), site);
    }
  }

  /** The infinite-dilution ratio, or the exception that stood in for it. */
  private static String describeInfiniteDilution(neqsim.thermo.phase.PhaseGE phase, int k,
      int waterNumber) {
    try {
      return String.format("%.15g", phase.getActivityCoefficientInfDilWater(k, waterNumber));
    } catch (Throwable error) {
      return error.getClass().getSimpleName() + ": " + error.getMessage();
    }
  }

  public static void main(String[] args) {
    // ---------------------------------------------------------------------------------
    System.out.println("################ a one-component GE phase ################");
    // ---------------------------------------------------------------------------------
    // CO2 rather than water, because `ComponentGE.fugcoef` only reaches the
    // infinite-dilution path for a `solute`-reference component in a phase with no water.
    oneComponentGe("SystemNRTL", new neqsim.thermo.system.SystemNRTL(313.15, 1.0));
    oneComponentGe("SystemGEWilson", new neqsim.thermo.system.SystemGEWilson(313.15, 1.0));
    oneComponentGe("SystemDuanSun", new SystemDuanSun(313.15, 1.0));
    System.out.println("  (a one-component phase holding water, or a two-component one,");
    System.out.println("   takes the solvent branch and is unaffected)");
    {
      SystemInterface waterOnly = new neqsim.thermo.system.SystemNRTL(313.15, 1.0);
      waterOnly.addComponent("water", 1.0);
      waterOnly.setMixingRule("classic");
      waterOnly.init(0);
      waterOnly.init(1);
      System.out.printf("  %-15s ok, gamma = %.12g%n", "SystemNRTL water",
          waterOnly.getPhase(1).getActivityCoefficient(0, 0));
    }
    System.out.println();

    // ---------------------------------------------------------------------------------
    System.out.println("################ Kent-Eisenberg ################");
    // ---------------------------------------------------------------------------------
    report("water + MDEA + Na+ + Cl- + CO2",
        new SystemKentEisenberg(T, P),
        new String[] {"water", "MDEA", "Na+", "Cl-", "CO2"},
        new double[] {0.80, 0.08, 0.04, 0.04, 0.04}, true);
    report("water + Na+ + Cl-",
        new SystemKentEisenberg(T, P),
        new String[] {"water", "Na+", "Cl-"},
        new double[] {0.90, 0.05, 0.05}, true);
    // **The brine a Kent-Eisenberg case can state**: no MDEA, because MDEA's reference
    // state is `solvent` and it carries no vapour pressure, so its `phi` is `NaN` and a
    // port that refuses rather than returning NaN cannot take it. This one exercises the
    // solvent branch (water), the constant ion branch, and the Henry branch (CO2).
    report("water + Na+ + Cl- + CO2", new SystemKentEisenberg(T, P),
        new String[] {"water", "Na+", "Cl-", "CO2"}, new double[] {0.89, 0.04, 0.04, 0.03},
        true);
    // At two more states, because the Henry branch is the only one that moves with `T`.
    for (double t : new double[] {298.15, 373.15}) {
      report(String.format("water + Na+ + Cl- + CO2 at %.2f K", t), new SystemKentEisenberg(t, P),
          new String[] {"water", "Na+", "Cl-", "CO2"}, new double[] {0.89, 0.04, 0.04, 0.03},
          true);
    }

    // ---------------------------------------------------------------------------------
    System.out.println();
    System.out.println("################ Desmukh-Mather ################");
    // ---------------------------------------------------------------------------------
    SystemInterface desmukh = new SystemDesmukhMather(T, P);
    String[] dmNames = {"water", "MDEA", "Na+", "Cl-", "CO2"};
    double[] dmZ = {0.80, 0.08, 0.04, 0.04, 0.04};
    // **Built once.** `addComponent` accumulates, so building a system twice doubles the
    // composition - which is what the first run of this probe did, and it doubles the
    // ionic strength and the solvent weight with it.
    report("water + MDEA + Na+ + Cl- + CO2", desmukh, dmNames, dmZ, true);
    PhaseDesmukhMather phase = (PhaseDesmukhMather) phase(desmukh);
    System.out.printf("%n  I = %.15g   solvent weight = %.15g   solvent molar mass = %.15g%n",
        phase.getIonicStrength(), phase.getSolventWeight(), phase.getSolventMolarMass());
    String[] names = new String[phase.getNumberOfComponents()];
    for (int i = 0; i < names.length; i++) {
      names[i] = phase.getComponent(i).getComponentName();
    }
    System.out.println("  aij / bij, and beta = aij + bij T:");
    for (int i = 0; i < names.length; i++) {
      for (int j = 0; j < names.length; j++) {
        System.out.printf("    %-9s %-9s aij = %20.15g  bij = %20.15g  beta = %20.15g%n", names[i],
            names[j], phase.getAij(i, j), phase.getBij(i, j), phase.getBetaDesMatij(i, j));
      }
    }
    // The diameters are a private field with no getter, so they are read from the same
    // table the constructor reads - and **from `comptemp` first**, falling back to
    // `comp`, which is the order that constructor uses.
    System.out.println("  ionic diameters, as the component's own constructor reads them:");
    try (neqsim.util.database.NeqSimDataBase database = new neqsim.util.database.NeqSimDataBase()) {
      for (int i = 0; i < names.length; i++) {
        String diameter = null;
        String source = "comptemp";
        try (java.sql.ResultSet rows = database
            .getResultSet("SELECT * FROM comptemp WHERE name='" + names[i] + "'")) {
          rows.next();
          rows.getString("FORMULA");
          diameter = rows.getString("DeshMatIonicDiameter");
        } catch (Exception first) {
          source = "comp";
        }
        String fromComp = "no row";
        try (java.sql.ResultSet rows = database
            .getResultSet("SELECT * FROM comp WHERE name='" + names[i] + "'")) {
          rows.next();
          fromComp = rows.getString("DeshMatIonicDiameter");
        } catch (Exception second) {
          fromComp = "unavailable";
        }
        if (diameter == null) {
          diameter = fromComp;
        }
        System.out.printf("    %-9s %-9s = %-12s comp = %s%n", names[i], source, diameter,
            fromComp);
      }
    } catch (Exception ex) {
      System.out.printf("    the database is unavailable: %s%n", ex.getMessage());
    }

    report("water + Na+ + Cl-", new SystemDesmukhMather(T, P),
        new String[] {"water", "Na+", "Cl-"}, new double[] {0.90, 0.05, 0.05}, true);

    // **The brine whose pair the table actually carries.** `MDEA+`/`CO2` is one of the
    // four nonzero `aijDesMath` rows, so this is the whole model: everywhere else the pair
    // sum is identically zero and only the Debye-Huckel term is left.
    report("water + MDEA+ + Cl- + CO2", new SystemDesmukhMather(T, P),
        new String[] {"water", "MDEA+", "Cl-", "CO2"}, new double[] {0.89, 0.04, 0.04, 0.03},
        true);

    // **A mixed solvent, which is the case the reference-state rule is for.** Methanol is
    // tagged `solvent` and carries a real Antoine row, so it joins water in the solvent
    // weight and the mean molar mass is neither component's. `MDEA` would be the amine
    // case and cannot be used: it is `solvent`-tagged with no vapour pressure at all.
    report("water + methanol + Na+ + Cl-", new SystemDesmukhMather(T, P),
        new String[] {"water", "methanol", "Na+", "Cl-"},
        new double[] {0.85, 0.05, 0.05, 0.05}, true);

    // And the infinite-dilution ratio the solute branch divides by, which is a second
    // evaluation of the same arithmetic at a two-component reference state.
    {
      SystemInterface dm = build(new SystemDesmukhMather(T, P),
          new String[] {"water", "MDEA+", "Cl-", "CO2"}, new double[] {0.89, 0.04, 0.04, 0.03});
      PhaseDesmukhMather p2 = (PhaseDesmukhMather) phase(dm);
      System.out.println("  the infinite-dilution activity coefficient `getActivityCoefficientInfDilWater`:");
      int waterNumber = p2.getComponent("water").getComponentNumber();
      for (int i = 0; i < p2.getNumberOfComponents(); i++) {
        ComponentInterface c = p2.getComponent(i);
        System.out.printf("    %-9s reference phase: %s%n", c.getComponentName(),
            neqsim.thermo.component.ComponentGEInterface.class.isAssignableFrom(c.getClass())
                ? describeInfiniteDilution(p2, i, waterNumber)
                : "(not a GE component)");
      }
    }

    // ---------------------------------------------------------------------------------
    System.out.println();
    System.out.println("################ Duan-Sun ################");
    // ---------------------------------------------------------------------------------
    // **`SystemDuanSun` admits CO2 and nothing else** - `addComponent` throws for any
    // other name - although the component's own correlation carries nitrogen and oxygen
    // too. **So no brine topology can be built at all**: the water this phase's own
    // expressions divide by is the first name refused, and the state below is the one the
    // record is about. The other two gases' rows are reachable only through `PhaseDuanSun`
    // directly.
    //
    // **Every attempt is caught rather than left to end the run**, because a probe that
    // dies here records one stack trace and no rows - and `#3841` made the refusal the
    // answer at the *factory* too, so there is no topology left to report on.
    for (String[] feed : new String[][] {
        {"water", "Na+", "Cl-", "CO2"},
        {"nitrogen"},
        {"CO2"}}) {
      try {
        SystemInterface refused = new SystemDuanSun(T, P);
        for (String name : feed) {
          refused.addComponent(name, 1.0 / feed.length);
        }
        System.out.printf("  %-24s accepted%n", String.join(" + ", feed));
      } catch (RuntimeException ex) {
        System.out.printf("  %-24s throws %s%n", String.join(" + ", feed), ex.getMessage());
      }
    }
    try {
      SystemInterface refused = new SystemDuanSun(T, P);
      refused.addComponent("CO2", 1.0);
      refused.setMixingRule("classic");
      refused.init(0);
      refused.init(1);
      System.out.println("  CO2 alone builds and flashes");
    } catch (Exception ex) {
      System.out.printf("  CO2 alone throws %s: %s%n", ex.getClass().getSimpleName(), ex.getMessage());
    }

    // The correlation is a function of `T`, `P` and the salinity alone, so this is the
    // table a port is checked against: lambda and zeta for each of the three gases at
    // each state.
    System.out.println();
    System.out.println("  the three gases against T and P, salinity = 4.0 mol/kg:");
    System.out.printf("  %8s %8s %22s %22s %22s%n", "T / K", "P / bar", "lambda CO2", "lambda N2",
        "lambda O2");
    for (double t : new double[] {298.15, 313.15, 373.15}) {
      for (double p : new double[] {1.0, 5.0, 50.0}) {
        double lambdaCO2 = -0.411370585 + 0.000607632 * t + 97.5347708 / t
            - 0.023762247 * p / t + 0.017065624 * p / (630.0 - t)
            + 1.41335834e-5 * t * Math.log(p);
        double lambdaN2 = -2.4434074 + 0.0036351795 * t + 447.47364 / t - 0.000013711527 * p
            + 0.0000071037217 * p * p / t;
        System.out.printf("  %8.2f %8.2f %22.15g %22.15g %22.15g%n", t, p, lambdaCO2, lambdaN2,
            0.19997);
      }
    }
    System.out.println("  zeta: CO2 = -0.58071053e-2 N2 = -1.2793e-2 O2 = "
        + "0.00033639 - 1.9829898e-5 T + 0.002122208 P/T - 0.005248733 P/(630 - T)");

    // The salinity the phase reports, against `sum m_i` over the ions - which no brine
    // reaches, because the phase admits CO2 alone. **The refusal is the answer here**, and
    // it is printed rather than thrown: a probe that ends on a stack trace ends its capture
    // there too, and the section after it is lost without anything comparing the two.
    System.out.println();
    try {
      SystemInterface duan = build(new SystemDuanSun(T, P),
          new String[] {"water", "Na+", "Cl-", "CO2"},
          new double[] {0.89, 0.04, 0.04, 0.03});
      PhaseInterface duanPhase = phase(duan);
      double sumMolality = 0.0;
      for (int i = 0; i < duanPhase.getNumberOfComponents(); i++) {
        if (duanPhase.getComponent(i).isIsIon()) {
          sumMolality += duanPhase.getComponent(i).getMolality(duanPhase);
        }
      }
      System.out.printf("  sum m_i over the ions = %.15g%n", sumMolality);
      System.out.printf("  gamma(CO2) = %.15g%n",
          duanPhase.getActivityCoefficient(
              duanPhase.getComponent("CO2").getComponentNumber(), 0));
    } catch (RuntimeException ex) {
      System.out.printf("  an aqueous brine is refused: %s%n", ex.getMessage());
      System.out.println("  (so `calcSalinity`'s `sum m_i` has no brine to be computed over,"
          + " and the correlation above is reachable only through `PhaseDuanSun`)");
    }
  }
}
