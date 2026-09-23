// The reactive hybrid flash's **coupled loop, pass by pass**, for the port of
// `reactions.reactive_hybrid_eos_ge_flash`.
//
// `HybridEosGeReactiveProbe` captures where the flash ends and this captures how it gets
// there, because the endpoint alone does not pin the loop: the pass count, the two
// convergence tests, the write-back of the adjusted inventory and the per-pass reaction
// delta are all interior, and a port that reports the same final state by another route
// would look identical without them.
//
// `TPHybridEosGeFlash.run` is private work, so this probe **mirrors its body** rather than
// calling it, driving the class's own private `solveFixedTopologyPhaseEquilibrium` and
// `solveAqueousChemicalEquilibrium` by reflection in the order `run` drives them. The mirror
// is verified by its own output: the endpoint it prints is the one
// `HybridEosGeReactiveProbe` already captured, to the digits that capture holds.
//
// The conservative delta is **not** recomputed by hand. `coupledOverallMoles` accumulates
// `old + reactionDeltas`, so the difference between two passes' inventories *is*
// `getConservativeReactionDeltas`'s answer, floored at `1e-45` - which is the quantity the
// projection produces and the one the port has to reproduce, without a second call to the
// method that would run it against an inventory it has already been applied to.
//
//     javac -proc:none -cp neqsim-f0c7436.jar HybridEosGeReactiveLoopProbe.java
//     java -cp .:neqsim-f0c7436.jar HybridEosGeReactiveLoopProbe > captures/hybrid_eos_ge_reactive_loop_probe.tsv

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.HybridEosGeFlashModel;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPitzer;
import neqsim.thermodynamicoperations.flashops.TPHybridEosGeFlash;

public class HybridEosGeReactiveLoopProbe {

  /** `TPHybridEosGeFlash.MAXIMUM_REACTIVE_ITERATIONS`. */
  private static final int MAXIMUM_REACTIVE_ITERATIONS = 100;

  /** `TPHybridEosGeFlash.MINIMUM_REACTIVE_ITERATIONS`. */
  private static final int MINIMUM_REACTIVE_ITERATIONS = 3;

  /** `TPHybridEosGeFlash.REACTIVE_COMPOSITION_TOLERANCE`. */
  private static final double REACTIVE_COMPOSITION_TOLERANCE = 1.0e-10;

  /** `TPHybridEosGeFlash.HYBRID_SOLVER_TOLERANCE`. */
  private static final double HYBRID_SOLVER_TOLERANCE = 1.0e-10;

  public static void main(String[] args) throws Exception {
    one("reactive-gas-aqueous", false);
    one("reactive-gas-oil-aqueous", true);
  }

  static void one(String label, boolean includeOil) throws Exception {
    System.out.println("fluid=" + label);
    SystemPitzer system = new SystemPitzer(313.15, 50.0);
    system.addComponent("methane", 5.0);
    system.addComponent("CO2", 0.05);
    if (includeOil) {
      system.addComponent("n-heptane", 2.0);
    }
    system.addComponent("water", 55.5);
    system.addComponent("Ca++", 6.0e-4);
    system.addComponent("Cl-", 2.0e-4);
    system.addComponent("HCO3-", 1.0e-3);
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule("classic");
    system.setMultiPhaseCheck(true);

    int components = system.getPhase(0).getNumberOfComponents();
    System.out.println("is_chemical_system=" + system.isChemicalSystem());
    System.out.println("components=" + names(system));

    // The exact overall inventory the coupled loop starts from, which is what
    // `updateCoupledOverallComposition` lazily reads on its first call - and **it has to be read
    // before the flash**: `synchronizeHybridEosGeOverallComposition` writes the adjusted
    // inventory into every role's component objects, so the feed totals this reads afterwards
    // are the coupled ones.
    System.out.println("feed_moles=" + join(amounts(system, true)));
    double[] before = conserved(system, true);
    double[] previousCoupled = new double[components];
    for (int i = 0; i < components; i++) {
      previousCoupled[i] = system.getPhase(0).getComponent(i).getNumberOfmoles();
    }

    // `TPflash.runInternal`'s hybrid branch, and then `TPHybridEosGeFlash.run`'s body.
    system.init(0);
    HybridEosGeFlashModel model = system;
    TPHybridEosGeFlash flash = new TPHybridEosGeFlash(system, model);
    Method reactiveIndex = declared("getHybridAqueousPhaseNumber");
    Method chemicalStep = declared("solveAqueousChemicalEquilibrium", boolean.class);
    Method fractionStep = declared("solveFixedTopologyPhaseEquilibrium");
    Field coupled = TPHybridEosGeFlash.class.getDeclaredField("coupledOverallMoles");
    coupled.setAccessible(true);

    model.prepareHybridEosGeFlash();
    system.init(1);
    int aqueousPhase = (Integer) reactiveIndex.invoke(flash);
    System.out.println("aqueous_role_index=" + aqueousPhase);
    System.out.println(
        "role_types=" + system.getPhase(0).getPhaseTypeName() + " "
            + system.getPhase(1).getPhaseTypeName() + " "
            + (system.getNumberOfPhases() > 2 ? system.getPhase(2).getPhaseTypeName() : "-"));

    // The feed-balanced split before chemistry sees the role seeds.
    double residual = (Double) fractionStep.invoke(flash);
    System.out.println("seed_pass_residual=" + residual);
    System.out.println("seed_pass_beta=" + betas(system));
    System.out.println("seed_pass_aqueous_moles=" + join(amounts(system.getPhase(aqueousPhase))));
    System.out.println("seed_pass_coupled_moles=" + join(previousCoupled));

    int passes = 0;
    boolean coupledConverged = false;
    double chemicalDeviation = Double.POSITIVE_INFINITY;
    for (int outer = 0; outer < MAXIMUM_REACTIVE_ITERATIONS; outer++) {
      passes = outer + 1;
      double[] beforeChemistry = amounts(system.getPhase(aqueousPhase));
      chemicalDeviation = (Double) chemicalStep.invoke(flash, outer == 0);
      double[] after = amounts(system.getPhase(aqueousPhase));
      residual = (Double) fractionStep.invoke(flash);
      coupledConverged = outer + 1 >= MINIMUM_REACTIVE_ITERATIONS
          && chemicalDeviation <= REACTIVE_COMPOSITION_TOLERANCE && Double.isFinite(residual)
          && residual <= HYBRID_SOLVER_TOLERANCE;

      double[] inventory = ((double[]) coupled.get(flash)).clone();
      System.out.println("pass=" + passes + " chemical_deviation=" + chemicalDeviation
          + " residual=" + residual + " converged=" + coupledConverged);
      System.out.println("  raw_delta=" + join(delta(beforeChemistry, after)));
      System.out.println("  raw_conservation_residual="
          + maximumConservationResidual(system, delta(beforeChemistry, after))
          + " raw_delta_is_non_negative="
          + rawDeltaIsNonNegative(system, previousCoupled, delta(beforeChemistry, after))
          + " feasible_step=" + feasibleStep(system, previousCoupled, delta(previousCoupled, inventory)));
      System.out.println("  conservative_delta=" + join(delta(previousCoupled, inventory)));
      System.out.println("  coupled_moles=" + join(inventory));
      System.out.println("  aqueous_moles=" + join(after));
      System.out.println("  aqueous_x=" + join(fractions(system.getPhase(aqueousPhase))));
      previousCoupled = inventory;
      if (coupledConverged) {
        break;
      }
    }

    System.out.println("passes=" + passes);
    System.out.println("coupled_converged=" + coupledConverged);
    System.out.println("phases=" + system.getNumberOfPhases());
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      PhaseInterface phase = system.getPhase(p);
      System.out.println("  phase[" + p + "]_type=" + phase.getPhaseTypeName()
          + " beta=" + system.getBeta(p) + " molecules=" + phase.getNumberOfMolesInPhase());
      for (int i = 0; i < phase.getNumberOfComponents(); i++) {
        ComponentInterface component = phase.getComponent(i);
        System.out.println("    x[" + p + "][" + component.getName() + "]=" + component.getx()
            + " moles=" + component.getNumberOfMolesInPhase());
      }
    }
    System.out.println("endpoint_aqueous_hco3_moles="
        + system.getPhase(aqueousPhase).getComponent("HCO3-").getNumberOfMolesInPhase());
    System.out.println("endpoint_aqueous_co3_moles="
        + system.getPhase(aqueousPhase).getComponent("CO3--").getNumberOfMolesInPhase());

    double[] after = conserved(system, false);
    double worstElement = 0.0;
    for (int i = 0; i < before.length - 1; i++) {
      worstElement = Math.max(worstElement, Math.abs(before[i] - after[i]));
    }
    System.out.println("worst_element_residual=" + worstElement);
    System.out.println(
        "net_charge_residual=" + Math.abs(before[before.length - 1] - after[after.length - 1]));
    System.out.println();
  }

  /** A private method of the flash, made callable. */
  static Method declared(String name, Class<?>... parameters) throws Exception {
    Method method = TPHybridEosGeFlash.class.getDeclaredMethod(name, parameters);
    method.setAccessible(true);
    return method;
  }

  /** The reactive components' amounts, from the system's own feed totals. */
  static double[] amounts(SystemInterface system, boolean useFeedTotals) {
    int components = system.getPhase(0).getNumberOfComponents();
    double[] values = new double[components];
    for (int i = 0; i < components; i++) {
      values[i] = useFeedTotals ? system.getPhase(0).getComponent(i).getNumberOfmoles()
          : system.getPhase(0).getComponent(i).getNumberOfMolesInPhase();
    }
    return values;
  }

  /** One phase's component amounts. */
  static double[] amounts(PhaseInterface phase) {
    double[] values = new double[phase.getNumberOfComponents()];
    for (int i = 0; i < values.length; i++) {
      values[i] = phase.getComponent(i).getNumberOfMolesInPhase();
    }
    return values;
  }

  /** One phase's mole fractions. */
  static double[] fractions(PhaseInterface phase) {
    double[] values = new double[phase.getNumberOfComponents()];
    for (int i = 0; i < values.length; i++) {
      values[i] = phase.getComponent(i).getx();
    }
    return values;
  }

  static double[] delta(double[] from, double[] to) {
    double[] values = new double[from.length];
    for (int i = 0; i < from.length; i++) {
      values[i] = to[i] - from[i];
    }
    return values;
  }

  /**
   * `max |A delta|` over the conservation rows, which is the test that decides the branch.
   *
   * Recomputed here from the printed delta and the same `A` matrix the class holds, so the
   * branch `getConservativeReactionDeltas` takes is a **measured** statement rather than one
   * read off the source: at or below `1e-8` the raw delta is used unprojected.
   */
  static double maximumConservationResidual(SystemInterface system, double[] delta) {
    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    ComponentInterface[] reactive = operations.getComponents();
    double[][] conservation = operations.getAmatrix();
    double worst = 0.0;
    for (double[] row : conservation) {
      double value = 0.0;
      for (int r = 0; r < reactive.length; r++) {
        value += row[r] * delta[reactive[r].getComponentNumber()];
      }
      worst = Math.max(worst, Math.abs(value));
    }
    return worst;
  }

  /** `coupledOverallMoles + delta` staying above `-1e-9`, the second clause of the same test. */
  static boolean rawDeltaIsNonNegative(SystemInterface system, double[] inventory, double[] delta) {
    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    for (ComponentInterface component : operations.getComponents()) {
      int index = component.getComponentNumber();
      if (inventory[index] + delta[index] < -1.0e-9) {
        return false;
      }
    }
    return true;
  }

  /** The largest step a delta may take before an inventory runs out, `feasibleStep`'s own. */
  static double feasibleStep(SystemInterface system, double[] inventory, double[] conservative) {
    double step = 1.0;
    for (ComponentInterface component : system.getChemicalReactionOperations().getComponents()) {
      int index = component.getComponentNumber();
      double change = conservative[index];
      if (change >= 0.0) {
        continue;
      }
      double available = Math.max(0.0, inventory[index] - 1.0e-45);
      step = Math.min(step, available / -change);
    }
    return step;
  }

  /** The reactive component set, in the reaction operations' own order. */
  static String names(SystemInterface system) {
    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    StringBuilder line = new StringBuilder();
    for (ComponentInterface component : operations.getComponents()) {
      line.append(component.getName()).append(" ");
    }
    return line.toString().trim();
  }

  static String betas(SystemInterface system) {
    StringBuilder line = new StringBuilder();
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      if (p > 0) {
        line.append(" ");
      }
      line.append(system.getBeta(p));
    }
    return line.toString();
  }

  static String join(double[] values) {
    StringBuilder line = new StringBuilder();
    for (int i = 0; i < values.length; i++) {
      if (i > 0) {
        line.append(" ");
      }
      line.append(values[i]);
    }
    return line.toString();
  }

  /// The conservation matrix applied to the component amounts, then the net charge.
  ///
  /// The same quantity `SystemHybridEosGeFlashTest.getReactiveConservedQuantities` computes,
  /// so what is printed here is what that test asserts.
  static double[] conserved(SystemInterface system, boolean useFeedTotals) {
    ChemicalReactionOperations reactions = system.getChemicalReactionOperations();
    ComponentInterface[] reactiveComponents = reactions.getComponents();
    double[][] conservation = reactions.getAmatrix();
    double[] quantities = new double[conservation.length];
    int components = system.getPhase(0).getNumberOfComponents();

    for (int reactiveIndex = 0; reactiveIndex < reactiveComponents.length; reactiveIndex++) {
      int componentIndex = reactiveComponents[reactiveIndex].getComponentNumber();
      double moles = useFeedTotals ? system.getPhase(0).getComponent(componentIndex).getNumberOfmoles()
          : summed(system, componentIndex);
      for (int row = 0; row < conservation.length - 1; row++) {
        quantities[row] += conservation[row][reactiveIndex] * moles;
      }
    }
    for (int componentIndex = 0; componentIndex < components; componentIndex++) {
      double moles = useFeedTotals ? system.getPhase(0).getComponent(componentIndex).getNumberOfmoles()
          : summed(system, componentIndex);
      quantities[quantities.length - 1] +=
          system.getPhase(0).getComponent(componentIndex).getIonicCharge() * moles;
    }
    return quantities;
  }

  static double summed(SystemInterface system, int componentIndex) {
    double total = 0.0;
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      total += system.getPhase(p).getComponent(componentIndex).getNumberOfMolesInPhase();
    }
    return total;
  }
}
