// The tear's two accelerators, printed by NeqSim itself.
//
//     javac -proc:none -cp neqsim-f0c7436.jar AccelerationProbe.java
//     java -cp .:neqsim-f0c7436.jar AccelerationProbe > captures/process_acceleration.tsv
//
// **Why this probe exists, and why it is not a `ProcessProbe` subcommand.** `Recycle` carries
// three acceleration methods and azoth's port carried none of them: the declaration was accepted,
// validated, and then ignored. What is being ported is arithmetic over a sequence of states - not
// a unit operation and not a flowsheet - so it gets its own file, the way `GibbsFluidProbe` does.
//
// **The two methods are reachable by different doors, which is what fixes the shape.**
// `BroydenAccelerator.accelerate` is public and is driven directly on a scripted sequence, so its
// delay boundary and each step's answer are read off its own return values. `applyWegstein-
// Acceleration` is *private* - the only door is `Recycle.run` - so Wegstein is driven by running a
// `Recycle` pass by pass with a rewritten inlet, and read back through `getWegsteinQFactors()`,
// which is public.
//
// **What this is for.** *Where* a step lands matters more here than whether the loop converges:
// both methods substitute directly for their first two calls, so a port that stepped one call
// early would still converge, to a slightly different place. Printing the first four calls makes
// the delay a row in a table rather than a property someone has to trust.
import java.util.UUID;

import neqsim.process.equipment.stream.Stream;
import neqsim.process.equipment.util.AccelerationMethod;
import neqsim.process.equipment.util.BroydenAccelerator;
import neqsim.process.equipment.util.Recycle;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

public class AccelerationProbe {
  public static void main(String[] args) {
    broyden();
    wegstein();
  }

  /** The linear map Broyden is driven on: `g(x) = A x + b`, with a contraction. */
  private static final double[][] A = {
      { 0.10, 0.05, 0.00 },
      { 0.00, 0.10, 0.05 },
      { 0.05, 0.00, 0.10 },
  };
  private static final double[] B = { 1.0, 2.0, 3.0 };

  private static double[] map(double[] x) {
    double[] out = new double[x.length];
    for (int i = 0; i < x.length; i++) {
      out[i] = B[i];
      for (int j = 0; j < x.length; j++) {
        out[i] += A[i][j] * x[j];
      }
    }
    return out;
  }

  private static String vector(double[] v) {
    StringBuilder text = new StringBuilder();
    for (int i = 0; i < v.length; i++) {
      if (i > 0) {
        text.append(' ');
      }
      text.append(v[i]);
    }
    return text.toString();
  }

  /**
   * `BroydenAccelerator` on a coupled map, driven the way `Recycle` drives it.
   *
   * Each step feeds the accelerator's own answer back as the next input, which is the fixed-point
   * iteration it accelerates - so the sequence is the method's and not a chosen one.
   */
  private static void broyden() {
    // A scalar map first: `g(x) = 0.5x + 1`, whose fixed point is 2 and whose secant slope is
    // exactly 0.5, so the answer is checkable by hand as well as against the port.
    BroydenAccelerator scalar = new BroydenAccelerator();
    double value = 0.0;
    for (int step = 1; step <= 5; step++) {
      double g = 0.5 * value + 1.0;
      double accelerated = scalar.accelerate(new double[] { value }, new double[] { g })[0];
      System.out.println("broyden_scalar_x_" + step + "=" + value);
      System.out.println("broyden_scalar_g_" + step + "=" + g);
      System.out.println("broyden_scalar_acc_" + step + "=" + accelerated);
      value = accelerated;
    }

    System.out.println("broyden_dim=" + A.length);
    // `new BroydenAccelerator()` with no dimension, exactly as `getBroydenAccelerator` makes it.
    BroydenAccelerator accelerator = new BroydenAccelerator();
    double[] x = { 0.0, 0.0, 0.0 };
    for (int step = 1; step <= 8; step++) {
      double[] g = map(x);
      double[] accelerated = accelerator.accelerate(x, g);
      System.out.println("broyden_x_" + step + "=" + vector(x));
      System.out.println("broyden_g_" + step + "=" + vector(g));
      System.out.println("broyden_acc_" + step + "=" + vector(accelerated));
      System.out.println("broyden_count_" + step + "=" + accelerator.getIterationCount());
      x = accelerated;
    }
    System.out.println();
  }

  /**
   * Wegstein through a `Recycle`, whose private accelerator is the only implementation.
   *
   * The inlet's composition is rewritten before each pass, so the maps `g` is asked for are a
   * scripted sequence rather than one that settles. A single inlet is the branch
   * `deactivateOnLowFlow` and the acceleration share, and it is the branch a flowsheet's tear
   * takes.
   */
  private static void wegstein() {
    System.out.println("wegstein_passes=6");
    String[] names = { "methane", "n-butane" };
    SystemInterface fluid = new SystemPrEos(300.0, 5.0);
    for (String name : names) {
      fluid.addComponent(name, 1.0);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("inlet", fluid);
    inlet.setFlowRate(1.0, "mol/sec");
    inlet.run();

    Recycle recycle = new Recycle("probe");
    recycle.addStream(inlet);
    // **`outletStream` starts as null and only `run` would set it** (`Recycle.java:51`, `:503`),
    // and `run` dereferences it before that - so the same trap the `flowsheet` probe's header
    // records has to be worked around here too: make the outlet stream first.
    Stream tear = new Stream("probe_out", fluid.clone());
    tear.setFlowRate(0.0, "mol/sec");
    tear.run();
    recycle.setOutletStream(tear);
    recycle.setAccelerationMethod(AccelerationMethod.WEGSTEIN);

    for (int pass = 1; pass <= 6; pass++) {
      // **A fresh fluid each pass, because a written fraction is not a written amount.** Setting
      // `setx` on one phase and re-running left the overall composition at 0.5/0.5 whatever was
      // written - measured - so the sequence a scripted composition produces is the constant one.
      // Rebuilding the system moves the overall moles, which is what the flash then splits.
      double first = 0.90 - 0.10 * (pass - 1);
      SystemInterface scripted = new SystemPrEos(300.0, 5.0);
      scripted.addComponent("methane", first);
      scripted.addComponent("n-butane", 1.0 - first);
      scripted.setMixingRule(2);
      inlet.setThermoSystem(scripted);
      inlet.run();

      // **`g(x)` is the flash of the scripted inlet**, and the recycle's own clone flashes the
      // same state - so printing the inlet here gives a port a `currentOutput` it can reconstruct,
      // which the outlet alone does not (the outlet is post-acceleration).
      System.out.println("wegstein_in_" + pass + "=" + record(inlet.getThermoSystem(), names));

      recycle.run(UUID.fromString("00000000-0000-0000-0000-0000000000" + String.format("%02d", pass)));

      SystemInterface out = recycle.getOutletStream().getThermoSystem();
      System.out.println("wegstein_out_" + pass + "=" + record(out, names));
      System.out.println("wegstein_iterations_" + pass + "=" + recycle.getIterations());
      double[] q = recycle.getWegsteinQFactors();
      System.out.println("wegstein_q_" + pass + "=" + (q == null ? "null" : vector(q)));
    }
    System.out.println();
  }

  /**
   * `Recycle.extractStreamValues`: `[T, P, n_mol_per_s, x_0 … x_{n-1}]`.
   *
   * The private method itself cannot be called, so it is rebuilt here from the same three
   * readers it uses - which is what makes the port's own `extract` checkable against it.
   */
  private static String record(SystemInterface fluid, String[] names) {
    StringBuilder text = new StringBuilder();
    text.append(fluid.getTemperature()).append(' ').append(fluid.getPressure()).append(' ')
        .append(fluid.getFlowRate("mole/sec"));
    for (String name : names) {
      text.append(' ').append(fluid.getPhase(0).getComponent(name).getx());
    }
    return text.toString();
  }
}
