// The facade around the equilibrium solve, for `reactions.reactive_phase_equilibrium`.
//
// Four things are printed per fluid, and they are four different measurements.
//
// **The phase listing, before and after a flash.** `getReactivePhaseIndex()` is *private*,
// so it is reached by reflection rather than re-implemented: the capture then holds the
// real method's answer and not this file's reading of it. Every phase's
// `getPhaseTypeName()` is printed beside it in system order, so which branch produced the
// answer is visible. The before-flash listing is the state the method's own comment says
// its fallback exists for, and it takes that branch rather than the aqueous one.
//
// **`A`, `b` and `n`.** The element matrix with its charge row, the conserved element
// amounts, and the composition - all read before the solve, because the solve overwrites
// the phase. `b`'s last entry is the point of the whole method: it is the charge the ions
// *outside* the reactive set leave behind, negated, and it is zero only when every ion in
// the phase is reactive. The bicarbonate brine is the fluid where it is not zero.
//
// **What `solveChemEq` answered**, with the three residuals it certifies against. It
// returns `false` on every fluid here - including the one whose first Newton refinement
// converges, which `ChemicalEquilibriumProbe` measures - so `false` alone does not
// distinguish "no reactive phase" from "did not certify". The reflected index does.
//
// **`reacHeat`**, summed over the reactive components, which is how the facade exposes it.
//
// The warm fluid is also the *ordering* measurement: without `chemicalReactionInit()` a
// CO2-water fluid at 400 K and 1.01325 bar is a single gas phase, and with it the flash
// splits into gas and aqueous. The ionic species the machinery adds are what put an
// aqueous phase there, so a caller that flashes before initialising the reactions
// classifies a different phase.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ReactionOperationsProbe.java
//     java -cp .:neqsim-f0c7436.jar ReactionOperationsProbe > captures/reaction_operations_probe.tsv

import java.lang.reflect.Method;

import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.chemicalreactions.chemicalequilibrium.ChemicalEquilibrium;
import neqsim.thermo.ThermodynamicConstantsInterface;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class ReactionOperationsProbe implements ThermodynamicConstantsInterface {
  public static void main(String[] args) {
    // Aqueous after the flash, which is the first branch of the phase search.
    one("co2-water-aqueous", 298.15, new String[] { "CO2", "water" }, new double[] { 0.01, 10.0 });

    // Two aqueous-reaching acids: the fluid whose solve does not converge.
    one("co2-h2s-water", 298.15, new String[] { "CO2", "H2S", "water" },
        new double[] { 0.01, 0.01, 10.0 });

    // Above water's normal boiling point, where the fluid is one gas phase until the
    // reaction machinery adds its ions.
    one("co2-water-warm", 400.0, new String[] { "CO2", "water" }, new double[] { 0.01, 10.0 });

    // A single gas phase after the flash, so no aqueous and nothing for the fallback to
    // take: the `-1` skip, on a fluid whose reaction set is not empty.
    one("gas-condensate-no-aqueous", 400.0,
        new String[] { "methane", "n-butane", "CO2", "water" },
        new double[] { 5.0, 0.5, 0.01, 1.0 });

    // A bicarbonate brine: `HCO3-` is in the reaction set and `Na+` is in no reaction at
    // all, so the spectators leave a charge behind and the charge row's target is not zero.
    one("bicarbonate-brine", 298.15, new String[] { "CO2", "water", "Na+", "HCO3-" },
        new double[] { 0.01, 10.0, 0.02, 0.02 });
  }

  static SystemSrkEos build(double temperature, String[] names, double[] moles) {
    SystemSrkEos system = new SystemSrkEos(temperature, 1.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule(2);
    system.init(0);
    return system;
  }

  static void one(String label, double temperature, String[] names, double[] moles) {
    System.out.println("fluid=" + label);
    System.out.println("temperature_K=" + temperature);
    System.out.println("pressure_bara=1.01325");

    SystemSrkEos system = build(temperature, names, moles);
    ChemicalReactionOperations operations = system.getChemicalReactionOperations();

    System.out.println("phases_before_flash=" + system.getNumberOfPhases());
    listPhases(system);
    System.out.println("reactive_phase_before_flash=" + reactivePhaseIndex(operations));

    system.setMultiPhaseCheck(true);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(system).TPflash();
    system.init(3);

    System.out.println("phases_after_flash=" + system.getNumberOfPhases());
    listPhases(system);
    int reactive = reactivePhaseIndex(operations);
    System.out.println("reactive_phase_after_flash=" + reactive);
    // The phase the facade's own readers fall back to when the search found nothing.
    int readPhase = reactive < 0 ? 0 : reactive;

    operations.setReactiveComponents(readPhase);
    ComponentInterface[] components = operations.getComponents();
    System.out.println("reactive_components=" + components.length);
    for (int i = 0; i < components.length; i++) {
      System.out.println("  component[" + i + "]=" + components[i].getName()
          + " ionic_charge=" + components[i].getIonicCharge());
    }

    // Read before the solve: the solve rewrites the phase, and the matrix and the element
    // amounts are what it was handed.
    System.out.println("Amatrix_built_by_constructor=" + (operations.getAmatrix() != null));
    double[][] a = operations.getAmatrix();
    double[] b = operations.calcBVector();
    double[] n = new double[components.length];
    for (int i = 0; i < components.length; i++) {
      int compNum = components[i].getComponentNumber();
      n[i] = system.getPhase(readPhase).getComponent(compNum).getNumberOfMolesInPhase();
    }
    // `getAllElements` is a `HashSet` iteration, so this line is what makes the capture's
    // `A` rows readable: they are in this order, and the charge row is last.
    String[] elements = operations.getAllElements();
    StringBuilder elementList = new StringBuilder("elements=");
    for (int i = 0; i < elements.length; i++) {
      elementList.append(elements[i]).append(i + 1 < elements.length ? " " : "");
    }
    System.out.println(elementList);
    System.out.println("phase_moles=" + system.getPhase(readPhase).getNumberOfMolesInPhase());
    // `calcBVector`'s `totalPhaseChargeMoles`: summed over **every** component the phase
    // holds, which is what the spectator ions make different from the reactive set's own
    // charge. The brine is the fluid where the two differ.
    double phaseCharge = 0.0;
    for (int i = 0; i < system.getPhase(readPhase).getNumberOfComponents(); i++) {
      ComponentInterface held = system.getPhase(readPhase).getComponent(i);
      phaseCharge += held.getIonicCharge() * held.getNumberOfMolesInPhase();
    }
    System.out.println("phase_charge=" + phaseCharge);
    // **`calcBVector` reads the field `nVector`, not the phase.** It is refreshed by a
    // solve that gets past the phase search, so on a skipped phase it still holds whatever
    // the constructor left there. `b_from_moles` is `A n` over the composition printed
    // below, which is what a caller holding this phase would compute.
    StringBuilder fromMoles = new StringBuilder("b_from_moles=");
    for (int e = 0; e < a.length; e++) {
      double amount = 0.0;
      for (int i = 0; i < components.length; i++) {
        amount += a[e][i] * n[i];
      }
      fromMoles.append(amount);
      if (e + 1 < a.length) {
        fromMoles.append(" ");
      }
    }
    System.out.println(fromMoles);

    // `calcChemRefPot` is public, so the reference potentials are read from the facade
    // rather than reassembled here. The activity coefficients have no such accessor and
    // are recomputed; they are a level-1 quantity, and reading them before `init(1)`
    // gives `-Infinity` for every species.
    double[] chemRef = operations.calcChemRefPot(readPhase);
    system.init(1, readPhase);
    PhaseInterface phase = system.getPhase(readPhase);
    double[] chemRefReduced = new double[components.length];
    double[] logActivity = new double[components.length];
    int waterNumb = 0;
    for (int i = 0; i < components.length; i++) {
      if (components[i].getComponentName().equals("water")) {
        waterNumb = i;
        break;
      }
    }
    for (int i = 0; i < components.length; i++) {
      chemRefReduced[i] = chemRef[i] / (R * phase.getTemperature());
      logActivity[i] = 0.0;
      if (components[i].calcActivity()) {
        logActivity[i] = phase.getLogActivityCoefficient(components[i].getComponentNumber(),
            components[waterNumb].getComponentNumber());
      }

    }
    printVector("chem_ref", chemRef);
    printVector("chem_ref_reduced", chemRefReduced);
    printVector("log_activity", logActivity);

    System.out.println("solveChemEq_returns=" + operations.solveChemEq(readPhase, 0));
    System.out.println("solveChemEq_charge_mols=" + operations.getReactivePhaseChargeMoles());
    System.out.println("solveChemEq_max_element_residual="
        + operations.getMaximumAbsoluteElementBalanceResidual());
    System.out.println("solveChemEq_max_reaction_log_residual="
        + operations.getMaximumAbsoluteReactionLogResidual());

    if (a == null) {
      System.out.println("A=null");
      System.out.println();
      return;
    }
    printMatrix("A", a);
    printVector("b", b);
    printVector("n_before", n);

    double[] after = new double[components.length];
    for (int i = 0; i < components.length; i++) {
      int compNum = components[i].getComponentNumber();
      after[i] = system.getPhase(readPhase).getComponent(compNum).getNumberOfMolesInPhase();
    }
    printVector("n_after", after);
    System.out.println("reacHeat_sum_over_components=" + reacHeatSum(operations, components));
    System.out.println();

    // A second, identically built fluid, because `ChemicalEquilibrium` reads its starting
    // composition from the phase and the phase above has just been overwritten. This is
    // the *first refinement of the ideal path* on the composition printed above, which is
    // what this library ports; `solveChemEq`'s own answer is the one that ran on the
    // other fluid, and the two are not the same measurement.
    if (reactive < 0) {
      System.out.println("direct=none, no phase to solve in");
      System.out.println();
      return;
    }
    SystemSrkEos second = build(temperature, names, moles);
    second.setMultiPhaseCheck(true);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(second).TPflash();
    second.init(3);
    ChemicalReactionOperations ops2 = operationsOn(second);
    ops2.setReactiveComponents(reactive);
    ComponentInterface[] comps2 = ops2.getComponents();
    double[] n2 = new double[comps2.length];
    for (int i = 0; i < comps2.length; i++) {
      n2[i] = second.getPhase(reactive).getComponent(comps2[i].getComponentNumber())
          .getNumberOfMolesInPhase();
    }
    printVector("direct_n0", n2);
    ChemicalEquilibrium direct =
        new ChemicalEquilibrium(ops2.getAmatrix(), ops2.calcBVector(), second, comps2, reactive);
    System.out.println("direct_converged=" + direct.solve());
    System.out.println("direct_iterations=" + direct.getLastIterationCount());
    System.out.println("direct_error=" + direct.getLastError());
    System.out.println("direct_tolerance=" + direct.getConvergenceTolerance());
    printVector("direct_moles", direct.getMoles());
    System.out.println();
  }

  static ChemicalReactionOperations operationsOn(SystemSrkEos system) {
    return system.getChemicalReactionOperations();
  }

  /// `reactionList.reacHeat(phase, name)` for one component at a time, which is the only
  /// way the facade exposes it.
  static double reacHeatSum(ChemicalReactionOperations operations, ComponentInterface[] components) {
    double total = 0.0;
    for (int i = 0; i < components.length; i++) {
      total += operations.reacHeat(0, components[i].getName());
    }
    return total;
  }

  static void listPhases(SystemInterface system) {
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      PhaseInterface phase = system.getPhase(i);
      System.out.println("  phase[" + i + "]=" + phase.getPhaseTypeName() + " moles="
          + phase.getNumberOfMolesInPhase());
    }
  }

  /// `getReactivePhaseIndex()` is private; this calls it rather than copying it.
  static int reactivePhaseIndex(ChemicalReactionOperations operations) {
    try {
      Method method = ChemicalReactionOperations.class.getDeclaredMethod("getReactivePhaseIndex");
      method.setAccessible(true);
      return (Integer) method.invoke(operations);
    } catch (ReflectiveOperationException ex) {
      throw new RuntimeException(ex);
    }
  }

  static void printMatrix(String name, double[][] m) {
    System.out.println(name + "=" + m.length + "x" + (m.length == 0 ? 0 : m[0].length));
    for (int r = 0; r < m.length; r++) {
      StringBuilder row = new StringBuilder("  " + name + "[" + r + "]=");
      for (int c = 0; c < m[r].length; c++) {
        row.append(m[r][c]).append(c + 1 < m[r].length ? " " : "");
      }
      System.out.println(row);
    }
  }

  static void printVector(String name, double[] v) {
    StringBuilder out = new StringBuilder(name + "=");
    for (int i = 0; i < v.length; i++) {
      out.append(v[i]).append(i + 1 < v.length ? " " : "");
    }
    System.out.println(out);
  }
}
