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
    // The four-component mixture the package's own multiphase test uses, hot enough that
    // several reactions run at once.
    one("methane-water-co2-hydrogen-1000K", 1000.0, 1.0,
        new String[] { "methane", "water", "CO2", "hydrogen" },
        new double[] { 0.4, 0.2, 0.2, 0.2 });
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
    System.out.println("fluid=" + label);
    System.out.println("temperature_K=" + temperature);
    System.out.println("pressure_bara=" + pressure);

    SystemInterface system = new SystemSrkEos(temperature, pressure);
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
