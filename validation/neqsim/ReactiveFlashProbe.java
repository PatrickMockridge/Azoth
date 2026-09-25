// `ReactiveMultiphaseTPflash`, for `reactions.reactive_tp_flash`.
//
// The package's own stack, not a flash that carries a chemical branch: `reactiveflash/` builds
// an element formula matrix, runs a reactive stability analysis, adds trial phases, drives a
// modified-RAND outer loop and removes the phases that turn out negligible. It is reached from
// `ThermodynamicOperations.reactiveTPflash()` and from nowhere else, and it uses no
// `chemicalreactions` class at all - so its prerequisites are the formation properties and a
// phase model that can return `ln phi` at a trial composition, which is what this capture is
// for.
//
// **The states are the ones the class's own tests use**, because they are the states its
// authors checked: the water-gas shift at 600 K, where the equilibrium constant is about 14
// and the products are strongly favoured, and a four-component methane/water/CO2/hydrogen
// mixture at 1000 K. Both are neutral fluids - no ion is present - which is the branch this
// port can reach, the ionic one being refused by azoth's own component rules.
//
// Everything the flash exposes is printed: the two residuals, the iteration count, the
// reaction count, the Gibbs energy, the Lagrange multipliers and the phase states, because
// the port has to reproduce all of them and a capture that printed only the composition would
// leave the interior unmeasured.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ReactiveFlashProbe.java
//     java -cp .:neqsim-f0c7436.jar ReactiveFlashProbe > captures/reactive_flash_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.flashops.reactiveflash.FormulaMatrix;
import neqsim.thermodynamicoperations.flashops.reactiveflash.ReactiveMultiphaseTPflash;
import neqsim.thermodynamicoperations.flashops.reactiveflash.ReactiveStabilityAnalysis;

public class ReactiveFlashProbe {

  public static void main(String[] args) {
    // The water-gas shift: CO + H2O = CO2 + H2, Kp about 14 at 600 K, so the products are
    // strongly favoured and the direction of the answer is not in doubt.
    one("wgs-600K", 600.0, 1.0,
        new String[] { "CO", "water", "CO2", "hydrogen" },
        new double[] { 0.25, 0.25, 0.25, 0.25 });
    // **The same two states on PR**, which is the cubic the process layer's streams carry:
    // a reactive tray has to flash the fluid it is actually handed, so this is the row its
    // own oracle has to be.
    one("WGS-600K-pr", 600.0, 1.0,
        new String[] { "CO", "water", "CO2", "hydrogen" },
        new double[] { 0.25, 0.25, 0.25, 0.25 }, false);
    // The four-component mixture the package's own multiphase test uses, hot enough that
    // several reactions run at once.
    one("methane-water-co2-hydrogen-1000K", 1000.0, 1.0,
        new String[] { "methane", "water", "CO2", "hydrogen" },
        new double[] { 0.4, 0.2, 0.2, 0.2 });
    one("methane-water-co2-hydrogen-1000K-pr", 1000.0, 1.0,
        new String[] { "methane", "water", "CO2", "hydrogen" },
        new double[] { 0.4, 0.2, 0.2, 0.2 }, false);
    // The ionic fluids: the matrix only, because the flash's ionic branch is one this port
    // refuses and the matrix is what decides how many reactions the fluid has.
    matrixOnly("co2-water-ions", 298.15,
        new String[] { "water", "CO2", "OH-", "H3O+" }, new double[] { 10.0, 0.01, 1e-10, 1e-10 });
    matrixOnly("water-meg", 298.15,
        new String[] { "water", "MEG" }, new double[] { 1.0, 1.0 });
    stabilityOnly("wgs-600K", 600.0, 1.0,
        new String[] { "CO", "water", "CO2", "hydrogen" },
        new double[] { 0.25, 0.25, 0.25, 0.25 });
    // The same chemistry at 50 bar, where the water the shift leaves behind has a liquid to
    // form - the case a stability analysis exists for.
    stabilityOnly("wgs-600K-50bar", 600.0, 50.0,
        new String[] { "CO", "water", "CO2", "hydrogen" },
        new double[] { 0.25, 0.25, 0.25, 0.25 });
    stabilityOnly("methane-water-co2-hydrogen-1000K", 1000.0, 1.0,
        new String[] { "methane", "water", "CO2", "hydrogen" },
        new double[] { 0.4, 0.2, 0.2, 0.2 });
    // A component the feed carries at 1e-40 and the equilibrium leaves there: nitrogen is inert
    // and its mole fraction stays under the class's `MIN_MOLES`, so its reference potential
    // takes the class's sentinel rather than a logarithm.
    stabilityOnly("wgs-600K-trace-nitrogen", 600.0, 1.0,
        new String[] { "CO", "water", "CO2", "hydrogen", "nitrogen" },
        new double[] { 0.25, 0.25, 0.25, 0.25, 1.0e-40 });
    // The phase bookkeeping: what the system holds before the driver touches it, and the
    // second initialisation path - the one where `initializeWithVLEFlash` is reachable.
    forcedOnePhase("wgs-600K", 600.0, 1.0,
        new String[] { "CO", "water", "CO2", "hydrogen" },
        new double[] { 0.25, 0.25, 0.25, 0.25 });
    forcedOnePhase("wgs-300K", 300.0, 1.0,
        new String[] { "CO", "water", "CO2", "hydrogen" },
        new double[] { 0.25, 0.25, 0.25, 0.25 });
    diis();
    // `NR = 0` with more than one phase and no ion: the driver falls back to a conventional
    // VLE flash, which is a path none of the reactive fluids above reaches.
    nonReactive("methane-water-300K-50bar", 300.0, 50.0,
        new String[] { "methane", "water" }, new double[] { 0.5, 0.5 });
  }

  /// `runNonReactiveFlash`: the Wilson-K successive substitution a fluid with no independent
  /// reaction gets instead of the RAND solve.
  static void nonReactive(String label, double temperature, double pressure, String[] names,
      double[] moles) {
    System.out.println("fluid=" + label);
    System.out.println("temperature_K=" + temperature);
    System.out.println("pressure_bara=" + pressure);

    SystemInterface system = new SystemSrkEos(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      System.out.println("  feed[" + names[i] + "]=" + moles[i]);
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    counts("before_flash", system);

    FormulaMatrix matrix = new FormulaMatrix(system);
    System.out.println("independent_reactions=" + matrix.getNumberOfIndependentReactions());
    System.out.println("phase0_type=" + system.getPhase(0).getType());

    ReactiveMultiphaseTPflash flash = new ReactiveMultiphaseTPflash(system);
    flash.run();

    System.out.println("converged=" + flash.isConverged());
    System.out.println("total_iterations=" + flash.getTotalIterations());
    counts("after_flash", system);
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      StringBuilder composition = new StringBuilder("  phase_x[" + phase + "]=");
      for (int i = 0; i < names.length; i++) {
        composition.append(system.getPhase(phase).getComponent(i).getx())
            .append(i + 1 < names.length ? " " : "");
      }
      System.out.println(composition);
    }
    System.out.println();
  }

  /// The DIIS accelerator's own answers, on a sequence it is fed rather than one a solve
  /// produces.
  ///
  /// `DIISAccelerator` is reached only from inside the RAND solve, so its extrapolation is
  /// not visible in any flash's numbers - and it is **load-bearing** there: the 300 K
  /// forced-one-phase run accepts 19 extrapolated steps out of its 35 iterations, which the
  /// capture's `diis_steps_accepted` records. Feeding it directly is what makes the port
  /// checkable: the Pulay system it solves, the rolling buffer's wrap and the `null` it
  /// returns on a singular matrix all become readable.
  static void diis() {
    System.out.println("fluid=diis-accelerator");
    int vectorLength = 3;
    int maxHistory = 4;
    neqsim.thermodynamicoperations.flashops.reactiveflash.DIISAccelerator diis =
        new neqsim.thermodynamicoperations.flashops.reactiveflash.DIISAccelerator(vectorLength,
            maxHistory);
    System.out.println("vector_length=" + vectorLength);
    System.out.println("max_history=" + maxHistory);
    System.out.println("count_at_construction=" + diis.getCount());
    System.out.println("can_extrapolate_at_construction=" + diis.canExtrapolate());
    // Seven entries against a history of four, so the circular buffer wraps and `bufferIndex`
    // has to be right for the extrapolation to be.
    for (int step = 0; step < 7; step++) {
      double[] iterate = new double[vectorLength];
      double[] residual = new double[vectorLength];
      for (int i = 0; i < vectorLength; i++) {
        iterate[i] = 0.5 * (step + 1) * (i + 1);
        residual[i] = 1.0 / (step + 1) * (i + 1) + 0.1 * step;
      }
      printVector("  iterate[" + step + "]", iterate);
      printVector("  residual[" + step + "]", residual);
      diis.addEntry(iterate, residual);
      System.out.println("  count[" + step + "]=" + diis.getCount());
      System.out.println("  can_extrapolate[" + step + "]=" + diis.canExtrapolate());
      double[] extrapolated = diis.extrapolate();
      System.out.println("  extrapolated[" + step + "]="
          + (extrapolated == null ? "null" : java.util.Arrays.toString(extrapolated)));
    }
    diis.reset();
    System.out.println("count_after_reset=" + diis.getCount());
    System.out.println("can_extrapolate_after_reset=" + diis.canExtrapolate());
    System.out.println("extrapolated_after_reset=" + (diis.extrapolate() == null ? "null" : "vector"));

    // Two entries whose residuals are identical: the overlap matrix's rows are then equal, the
    // elimination leaves a zero pivot, and the class returns null rather than a combination.
    neqsim.thermodynamicoperations.flashops.reactiveflash.DIISAccelerator singular =
        new neqsim.thermodynamicoperations.flashops.reactiveflash.DIISAccelerator(3, 4);
    double[] same = { 1.0, 2.0, 3.0 };
    singular.addEntry(new double[] { 1.0, 1.0, 1.0 }, same);
    singular.addEntry(new double[] { 2.0, 2.0, 2.0 }, same);
    System.out.println("singular_count=" + singular.getCount());
    System.out.println("singular_can_extrapolate=" + singular.canExtrapolate());
    System.out.println("singular_extrapolated="
        + (singular.extrapolate() == null ? "null" : "vector"));
    System.out.println();
  }

  /// The driver reached from a forced one-phase state, with the phase bookkeeping printed.
  ///
  /// `SystemSrkEos` **constructs with two phase objects**, each at `beta = 1.0` and each
  /// holding the whole feed, so a flash run on a system that has only been initialised takes
  /// the multiphase outer loop without any trial phase ever being added - which is where the
  /// two identical phases and the doubled Gibbs energy in the blocks above come from.
  /// `setNumberOfPhases(1)` after `init(1)` is the one thing that reaches
  /// `initializeWithVLEFlash`, and what that method then does - add a phase or return - is
  /// what decides which of the driver's two branches runs.
  static void forcedOnePhase(String label, double temperature, double pressure, String[] names,
      double[] moles) {
    System.out.println("fluid=" + label + "-forced-one-phase");
    System.out.println("temperature_K=" + temperature);
    System.out.println("pressure_bara=" + pressure);

    SystemInterface system = new SystemSrkEos(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      System.out.println("  feed[" + names[i] + "]=" + moles[i]);
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule("classic");
    counts("phases_at_construction", system);
    system.init(0);
    counts("after_init0", system);
    system.init(1);
    counts("after_init1", system);
    system.setNumberOfPhases(1);
    counts("after_force_one_phase", system);

    // `initializeWithVLEFlash`'s own Wilson ln K and Rachford-Rice, recomputed because the
    // class keeps neither.
    int nc = names.length;
    double[] lnK = new double[nc];
    boolean hasVolatile = false;
    for (int i = 0; i < nc; i++) {
      ComponentInterface component = system.getPhase(0).getComponent(i);
      double tc = component.getTC();
      double pc = component.getPC();
      double omega = component.getAcentricFactor();
      if (pc > 0 && tc > 0) {
        lnK[i] = Math.log(pc / pressure) + 5.373 * (1.0 + omega) * (1.0 - tc / temperature);
        if (Math.abs(lnK[i]) > 0.1) {
          hasVolatile = true;
        }
      }
    }
    printVector("wilson_lnk", lnK);
    System.out.println("has_volatile=" + hasVolatile);
    double vapour = 0.5;
    if (hasVolatile) {
      for (int iter = 0; iter < 100; iter++) {
        double f = 0.0;
        double derivative = 0.0;
        for (int i = 0; i < nc; i++) {
          double ki = Math.exp(lnK[i]);
          double denominator = 1.0 + vapour * (ki - 1.0);
          if (Math.abs(denominator) < 1e-30) {
            continue;
          }
          f += moles[i] * (ki - 1.0) / denominator;
          derivative -= moles[i] * (ki - 1.0) * (ki - 1.0) / (denominator * denominator);
        }
        if (Math.abs(f) < 1e-10) {
          break;
        }
        if (Math.abs(derivative) > 1e-30) {
          vapour -= f / derivative;
        }
        vapour = Math.max(0.0, Math.min(1.0, vapour));
      }
    }
    System.out.println("rachford_rice_v=" + vapour);

    // **The class's own `initializeWithVLEFlash`, on a fresh one-phase system.** It is private
    // and it mutates the system, so calling it is the only way to see the compositions it sets
    // - the flash below overwrites them, and recomputing its arithmetic here would be a
    // different implementation's answer to the same question.
    SystemInterface fresh = new SystemSrkEos(temperature, pressure);
    for (int i = 0; i < nc; i++) {
      fresh.addComponent(names[i], moles[i]);
    }
    fresh.setMixingRule("classic");
    fresh.init(0);
    fresh.init(1);
    fresh.setNumberOfPhases(1);
    try {
      java.lang.reflect.Method initialize =
          ReactiveMultiphaseTPflash.class.getDeclaredMethod("initializeWithVLEFlash");
      initialize.setAccessible(true);
      initialize.invoke(new ReactiveMultiphaseTPflash(fresh));
      counts("after_vle_init", fresh);
      for (int phase = 0; phase < fresh.getNumberOfPhases(); phase++) {
        StringBuilder composition = new StringBuilder("  vle_x[" + phase + "]=");
        for (int i = 0; i < nc; i++) {
          composition.append(fresh.getPhase(phase).getComponent(i).getx())
              .append(i + 1 < nc ? " " : "");
        }
        System.out.println(composition);
      }
    } catch (ReflectiveOperationException ex) {
      System.out.println("vle_init_failed=" + ex);
    }

    FormulaMatrix matrix = new FormulaMatrix(system);
    System.out.println("independent_reactions=" + matrix.getNumberOfIndependentReactions());

    ReactiveMultiphaseTPflash flash = new ReactiveMultiphaseTPflash(system);
    flash.run();

    System.out.println("converged=" + flash.isConverged());
    System.out.println("total_iterations=" + flash.getTotalIterations());
    System.out.println("final_residual=" + flash.getFinalResidual());
    System.out.println("final_element_residual=" + flash.getFinalElementResidual());
    System.out.println("diis_steps_accepted=" + flash.getDiisStepsAccepted());
    System.out.println("equilibrium_total_moles=" + flash.getEquilibriumTotalMoles());
    System.out.println("final_gibbs_energy=" + flash.getFinalGibbsEnergy());
    printVector("lagrange_multipliers", flash.getLagrangeMultipliers());
    counts("after_flash", system);
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      StringBuilder composition = new StringBuilder("  phase_x[" + phase + "]=");
      for (int i = 0; i < nc; i++) {
        composition.append(system.getPhase(phase).getComponent(i).getx())
            .append(i + 1 < nc ? " " : "");
      }
      System.out.println(composition);
    }
    double[][] equilibriumMoles = flash.getEquilibriumMoles();
    for (int phase = 0; phase < equilibriumMoles.length; phase++) {
      printVector("equilibrium_moles[" + phase + "]", equilibriumMoles[phase]);
    }
    // **The solver's own phase amounts against the phase objects' betas.** `updateSystem`
    // writes each phase's beta and then calls `system.init(1)`, which re-initialises the
    // phase from the *system's* beta array - and the system's copy is only refreshed on the
    // ionic branch. So the two can disagree, and `computeGibbsEnergy` weighs the phases by
    // the phase object's number.
    try {
      java.lang.reflect.Field held = ReactiveMultiphaseTPflash.class.getDeclaredField("solver");
      held.setAccessible(true);
      Object solver = held.get(flash);
      if (solver != null) {
        java.lang.reflect.Method amounts =
            solver.getClass().getDeclaredMethod("getPhaseAmounts");
        amounts.setAccessible(true);
        printVector("solver_phase_amounts", (double[]) amounts.invoke(solver));
        java.lang.reflect.Method total =
            solver.getClass().getDeclaredMethod("getTotalMoles");
        total.setAccessible(true);
        System.out.println("solver_total_moles=" + total.invoke(solver));
      } else {
        System.out.println("solver=null");
      }
      printVector("system_beta", new double[] { system.getBeta(0), system.getBeta(1) });
    } catch (ReflectiveOperationException ex) {
      System.out.println("solver_reflection_failed=" + ex);
    }
    System.out.println();
  }

  static void counts(String label, SystemInterface system) {
    StringBuilder out = new StringBuilder(label + "_phases=" + system.getNumberOfPhases()
        + " max_phases=" + system.getMaxNumberOfPhases());
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      // The phase *type* as well as the fraction: the driver's own code branches on it
      // (`liqIdx` is the index that is not `GAS`), and no azoth result carries a type.
      out.append(" type[").append(phase).append("]=").append(system.getPhase(phase).getPhaseTypeName())
          .append(" beta[").append(phase).append("]=").append(system.getPhase(phase).getBeta());
    }
    System.out.println(out);
  }

  /// The matrix alone, for a fluid whose flash needs the ionic branch this port refuses.
  ///
  /// `FormulaMatrix` is constructed from the phase's components and the element table, so it
  /// answers without a flash being run - and an ionic fluid is where its `Charge` row and the
  /// rank that decides `NR` can be read off the real class rather than reasoned about.
  static void matrixOnly(String label, double temperature, String[] names, double[] moles) {
    System.out.println("fluid=" + label);
    SystemInterface system = new SystemSrkEos(temperature, 1.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    FormulaMatrix matrix = new FormulaMatrix(system);
    System.out.println("elements=" + matrix.getNumberOfElements()
        + " components=" + matrix.getNumberOfComponents()
        + " rank=" + matrix.getRank()
        + " independent_reactions=" + matrix.getNumberOfIndependentReactions());
    StringBuilder elementList = new StringBuilder("  element_names=");
    for (String element : matrix.getElementNames()) {
      elementList.append(element).append(" ");
    }
    System.out.println(elementList.toString().trim());
    for (int row = 0; row < matrix.getNumberOfElements(); row++) {
      StringBuilder line = new StringBuilder("  A[" + row + "]=");
      for (int column = 0; column < matrix.getNumberOfComponents(); column++) {
        line.append(matrix.get(row, column))
            .append(column + 1 < matrix.getNumberOfComponents() ? " " : "");
      }
      System.out.println(line);
    }
    System.out.println();
  }

  /// The stability analysis's own answers, which the flash's final state cannot show.
  ///
  /// `ReactiveStabilityAnalysis` decides whether a second phase forms, and it does that by
  /// bringing the feed to *homogeneous* chemical equilibrium first - so the reference
  /// potentials it tests against are the equilibrated feed's, not the feed's - and then
  /// running a tangent-plane trial from each Wilson and pure-component seed. What it found is
  /// printed: every seed's TPD, which of them were unstable, and the trial compositions they
  /// equilibrated to.
  static void stabilityOnly(String label, double temperature, double pressure, String[] names,
      double[] moles) {
    System.out.println("fluid=" + label);
    SystemInterface system = new SystemSrkEos(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule("classic");
    system.setMaxNumberOfPhases(1);
    system.setNumberOfPhases(1);
    system.init(0);
    system.init(1);

    FormulaMatrix matrix = new FormulaMatrix(system);
    ReactiveStabilityAnalysis stability = new ReactiveStabilityAnalysis(system, matrix);

    // What the class's own `generateTrialPhases` will read, and the Wilson K it builds from
    // them - recomputed here because both are local to the method.
    double[] x = new double[system.getPhase(0).getNumberOfComponents()];
    double[] k = new double[x.length];
    StringBuilder xLine = new StringBuilder("phase_x_after_ce=");
    StringBuilder kLine = new StringBuilder("wilson_k=");
    for (int i = 0; i < x.length; i++) {
      x[i] = system.getPhase(0).getComponent(i).getx();
      double tc = system.getPhase(0).getComponent(i).getTC();
      double pc = system.getPhase(0).getComponent(i).getPC();
      double omega = system.getPhase(0).getComponent(i).getAcentricFactor();
      k[i] = (pc / system.getPressure()) * Math.exp(5.373 * (1.0 + omega) * (1.0 - tc / system.getTemperature()));
      xLine.append(x[i]).append(i + 1 < x.length ? " " : "");
      kLine.append(k[i]).append(i + 1 < k.length ? " " : "");
    }
    System.out.println(xLine);
    System.out.println(kLine);

    boolean unstable = stability.run();

    // **The seeds and their distances, which the class does not expose.** `generateTrialPhases`
    // and `runStabilityTrial` are private, so they are reached by reflection rather than
    // re-derived: a port has to reproduce every seed's TPD, and the class only keeps the ones
    // that came in under its threshold.
    try {
      java.lang.reflect.Method generate =
          ReactiveStabilityAnalysis.class.getDeclaredMethod("generateTrialPhases");
      generate.setAccessible(true);
      java.lang.reflect.Method trial =
          ReactiveStabilityAnalysis.class.getDeclaredMethod("runStabilityTrial", double[].class);
      trial.setAccessible(true);
      // What the same formula and the *live* composition would give, beside what the class
      // returns, so the two are compared in one run rather than across two.
      double xSum = 0.0;
      for (int i = 0; i < x.length; i++) {
        xSum += system.getPhase(0).getComponent(i).getx();
      }
      System.out.println("live_x_sum=" + xSum);
      StringBuilder mineLine = new StringBuilder("recomputed_liquid_seed=");
      StringBuilder vapourLine = new StringBuilder("recomputed_vapour_seed=");
      for (int i = 0; i < x.length; i++) {
        double zi = system.getPhase(0).getComponent(i).getx();
        mineLine.append(zi / k[i]).append(i + 1 < x.length ? " " : "");
        vapourLine.append(k[i] * zi).append(i + 1 < x.length ? " " : "");
      }
      System.out.println(mineLine);
      System.out.println(vapourLine);

      @SuppressWarnings("unchecked")
      java.util.List<double[]> seeds = (java.util.List<double[]>) generate.invoke(stability);
      System.out.println("trial_seeds=" + seeds.size());
      for (int i = 0; i < seeds.size(); i++) {
        printVector("seed[" + i + "]", seeds.get(i));
        System.out.println("  seed_tpd[" + i + "]=" + trial.invoke(stability, (Object) seeds.get(i)));
      }
    } catch (ReflectiveOperationException ex) {
      System.out.println("reflection_failed=" + ex);
    }

    // **The reference potentials the trials are run against**, which the class keeps in a
    // private field and never prints: `d_i = ln x_i + ln phi_i` at the equilibrated feed, with
    // `-100.0` for a component at or below its `MIN_MOLES` floor and `-1000.0` for an ion.
    try {
      java.lang.reflect.Field potentials = ReactiveStabilityAnalysis.class.getDeclaredField("d");
      potentials.setAccessible(true);
      printVector("reference_potentials", (double[]) potentials.get(stability));
    } catch (ReflectiveOperationException ex) {
      System.out.println("reference_potentials_failed=" + ex);
    }

    // What the class leaves in the phase: if its CE step wrote the equilibrated composition
    // with `setx`, this is it, and the seeds are built from that rather than from the feed.
    StringBuilder afterRun = new StringBuilder("phase_x_after_run=");
    for (int i = 0; i < x.length; i++) {
      afterRun.append(system.getPhase(0).getComponent(i).getx()).append(i + 1 < x.length ? " " : "");
    }
    System.out.println(afterRun);
    System.out.println("unstable=" + unstable);
    System.out.println("number_of_unstable_trials=" + stability.getNumberOfUnstableTrials());
    double[] worst = stability.getMostUnstableTrial();
    System.out.println("most_unstable_trial=" + (worst == null ? "null" : java.util.Arrays.toString(worst)));
    java.util.List<double[]> trials = stability.getUnstableTrialCompositions();
    java.util.List<Double> tpds = stability.getTpdValues();
    for (int i = 0; i < trials.size(); i++) {
      printVector("unstable_trial[" + i + "]", trials.get(i));
      System.out.println("  tpd[" + i + "]=" + tpds.get(i));
    }
    System.out.println();
  }

  static void one(String label, double temperature, double pressure, String[] names, double[] moles) {
    one(label, temperature, pressure, names, moles, true);
  }

  /// **The same row on the other cubic.** The class flashes whatever system it is handed, so
  /// `pr` is a state NeqSim answers as readily as `srk` - it is the port's process layer that
  /// is PR-only, which is why a reactive tray needs this row to be oracled against.
  static void one(String label, double temperature, double pressure, String[] names,
      double[] moles, boolean srk) {
    System.out.println("fluid=" + label);
    System.out.println("temperature_K=" + temperature);
    System.out.println("pressure_bara=" + pressure);

    System.out.println("cubic=" + (srk ? "srk" : "pr"));
    SystemInterface system = srk ? new SystemSrkEos(temperature, pressure)
        : new SystemPrEos(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      System.out.println("  feed[" + names[i] + "]=" + moles[i]);
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule("classic");
    system.setMaxNumberOfPhases(1);
    system.setNumberOfPhases(1);
    system.init(0);
    system.init(1);

    FormulaMatrix matrix = new FormulaMatrix(system);
    System.out.println("elements=" + matrix.getNumberOfElements()
        + " components=" + matrix.getNumberOfComponents()
        + " rank=" + matrix.getRank()
        + " independent_reactions=" + matrix.getNumberOfIndependentReactions());
    for (int row = 0; row < matrix.getNumberOfElements(); row++) {
      StringBuilder line = new StringBuilder("  A[" + row + "]=");
      for (int column = 0; column < matrix.getNumberOfComponents(); column++) {
        line.append(matrix.get(row, column)).append(column + 1 < matrix.getNumberOfComponents() ? " " : "");
      }
      System.out.println(line);
    }

    ReactiveMultiphaseTPflash flash = new ReactiveMultiphaseTPflash(system);
    flash.run();

    System.out.println("converged=" + flash.isConverged());
    System.out.println("total_iterations=" + flash.getTotalIterations());
    System.out.println("final_residual=" + flash.getFinalResidual());
    System.out.println("final_element_residual=" + flash.getFinalElementResidual());
    System.out.println("number_of_reactions=" + flash.getNumberOfReactions());
    System.out.println("equilibrium_total_moles=" + flash.getEquilibriumTotalMoles());
    System.out.println("final_gibbs_energy=" + flash.getFinalGibbsEnergy());
    printVector("lagrange_multipliers", flash.getLagrangeMultipliers());

    System.out.println("phases=" + system.getNumberOfPhases());
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      PhaseInterface held = system.getPhase(phase);
      System.out.println("  phase[" + phase + "]=" + held.getPhaseTypeName()
          + " beta=" + held.getBeta()
          + " moles=" + held.getNumberOfMolesInPhase());
      StringBuilder composition = new StringBuilder("    x=");
      for (int i = 0; i < held.getNumberOfComponents(); i++) {
        ComponentInterface component = held.getComponent(i);
        composition.append(component.getName()).append(":").append(component.getx()).append(" ");
      }
      System.out.println(composition.toString().trim());
      System.out.println("    total_moles=" + held.getNumberOfMolesInPhase());
    }
    printVector("equilibrium_moles[0]", flash.getEquilibriumMoles()[0]);
    System.out.println();
  }

  static void printVector(String name, double[] v) {
    StringBuilder out = new StringBuilder(name + "=");
    for (int i = 0; i < v.length; i++) {
      out.append(v[i]).append(i + 1 < v.length ? " " : "");
    }
    System.out.println(out);
  }
}
