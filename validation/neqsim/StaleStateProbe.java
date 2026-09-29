// Which call repairs a stale outlet state.
//
// The `process.*` captures show a unit operation publishing an outlet whose own `(T, P, z)`
// implies a different enthalpy. This probe takes that outlet and asks, one call at a time,
// which operation lands it on the value a freshly constructed fluid at the same `(T, P, z)`
// already has.
//
// The candidate list is the whole point. `TPflash` re-solves the split; `init(0)`/`init(1)`/
// `init(2)` rebuild the different property levels; `initProperties` and
// `initThermoProperties` are the named wrappers. A repair that is a *re-flash* and a repair
// that is a *caloric rebuild* are different defects, and only the second one is a fix rather
// than a workaround.
//
//     javac -proc:none -cp neqsim-f0c7436.jar StaleStateProbe.java
//     java -cp .:neqsim-f0c7436.jar StaleStateProbe > captures/stale_state.tsv

import neqsim.process.equipment.splitter.ComponentSplitter;
import neqsim.process.equipment.stream.Stream;
import neqsim.process.equipment.stream.StreamInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

public class StaleStateProbe {

  static Stream feed(String[] names, double[] z, double temperatureK, double pressureBara,
      double molPerSecond) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream stream = new Stream("feed", fluid);
    stream.setFlowRate(molPerSecond, "mol/sec");
    stream.run();
    return stream;
  }

  /** Molar enthalpy and entropy, the way the capture reads them: phase 0, per mole. */
  static double[] hs(SystemInterface fluid) {
    double n = fluid.getTotalNumberOfMoles();
    return new double[] { fluid.getPhase(0).getEnthalpy() / fluid.getPhase(0).getNumberOfMolesInPhase(),
        fluid.getPhase(0).getEntropy() / fluid.getPhase(0).getNumberOfMolesInPhase() };
  }

  static void row(String label, SystemInterface fluid) {
    double[] v = hs(fluid);
    double n = fluid.getTotalNumberOfMoles();
    System.out.printf("%-34s h=%.12f s=%.12f total_h=%.12f phases=%d types=%s betas=%s%n", label,
        v[0], v[1], fluid.getEnthalpy() / n, fluid.getNumberOfPhases(), phaseTypes(fluid),
        betas(fluid));
  }

  static String betas(SystemInterface fluid) {
    StringBuilder b = new StringBuilder();
    for (int i = 0; i < fluid.getNumberOfPhases(); i++) {
      if (i > 0) {
        b.append(",");
      }
      b.append(String.format("%.6f", fluid.getPhase(i).getBeta()));
    }
    return b.toString();
  }

  static String phaseTypes(SystemInterface fluid) {
    StringBuilder b = new StringBuilder();
    for (int i = 0; i < fluid.getNumberOfPhases(); i++) {
      if (i > 0) {
        b.append(",");
      }
      b.append(fluid.getPhase(i).getPhaseTypeName());
    }
    return b.toString();
  }

  /** Every candidate repair, each on its own clone so nothing carries over. */
  static void repairs(SystemInterface reported, String[] names, double[] z, double tK, double pBara,
      double flow) {
    // **The comparator only works built through `Stream.run()`.** `feed()` is the route the
    // published captures use, and it reproduces the feed row's `h` to all printed digits;
    // the obvious `init(0); TPflash()` route does not, because it is *itself* a
    // stale-caloric construction (see the control row below). An acceptance test written on
    // the obvious route compares against a third, equally wrong number.
    row("fresh_fluid (via Stream)", feed(names, z, tK, pBara, flow).getThermoSystem());

    // The control: what the obvious construction returns. Kept because it is the trap.
    SystemInterface naive = new SystemPrEos(tK, pBara);
    for (int i = 0; i < names.length; i++) {
      naive.addComponent(names[i], z[i]);
    }
    naive.setMixingRule(2);
    naive.setTotalFlowRate(flow, "mol/sec");
    naive.init(0);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(naive).TPflash();
    row("control: init(0)+TPflash", naive);

    SystemInterface f1 = reported.clone();
    new neqsim.thermodynamicoperations.ThermodynamicOperations(f1).TPflash();
    row("clone + TPflash", f1);

    SystemInterface f2 = reported.clone();
    f2.init(0);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(f2).TPflash();
    row("clone + init(0) + TPflash", f2);

    SystemInterface f3 = reported.clone();
    f3.init(1);
    row("clone + init(1)", f3);

    SystemInterface f4 = reported.clone();
    f4.init(2);
    row("clone + init(2)", f4);

    SystemInterface f5 = reported.clone();
    f5.initProperties();
    row("clone + initProperties()", f5);

    SystemInterface f6 = reported.clone();
    f6.initThermoProperties();
    row("clone + initThermoProperties()", f6);

    SystemInterface f7 = reported.clone();
    f7.setTemperature(tK);
    f7.init(2);
    row("clone + setT + init(2)", f7);
  }

  public static void main(String[] args) {
    String[] names = new String[] { "methane", "n-butane", "n-pentane" };
    double[] z = new double[] { 0.5, 0.3, 0.2 };

    System.out.println("=== component_splitter near_total_separation ===");
    System.out.println("split_factors=0.98 0.05 0.02");
    Stream inlet = feed(names, z, 300.0, 20.0, 1.0);
    ComponentSplitter splitter = new ComponentSplitter("cs1", inlet);
    splitter.setSplitFactors(new double[] { 0.98, 0.05, 0.02 });
    splitter.run();

    StreamInterface overhead = splitter.getSplitStream(0);
    StreamInterface bottoms = splitter.getSplitStream(1);
    row("reported_overhead", overhead.getThermoSystem());
    System.out.println("overhead_flow_mol_per_sec=" + overhead.getFlowRate("mol/sec"));
    System.out.println("-- overhead repairs --");
    repairs(overhead.getThermoSystem(), new String[] { "methane", "n-butane", "n-pentane" },
        new double[] { 0.962671905697446, 0.029469548133595282, 0.007858546168958742 }, 300.0, 20.0,
        overhead.getFlowRate("mol/sec"));

    System.out.println();
    row("reported_bottoms", bottoms.getThermoSystem());
    System.out.println("bottoms_flow_mol_per_sec=" + bottoms.getFlowRate("mol/sec"));
    System.out.println("-- bottoms repairs --");
    repairs(bottoms.getThermoSystem(), new String[] { "methane", "n-butane", "n-pentane" },
        new double[] { 0.02036659877800409, 0.5804480651731161, 0.39918533604887985 }, 300.0, 20.0,
        bottoms.getFlowRate("mol/sec"));

    heatExchanger();
  }

  /// The `out_temperature_pins_the_cold_side` row of the `heat_exchanger` capture: the cold
  /// side is pinned at 320 K, so the hot side is the one `runSpecifiedStream` moves.
  static void heatExchanger() {
    System.out.println();
    System.out.println("=== heat_exchanger out_temperature_pins_the_cold_side ===");
    String[] hotNames = new String[] { "methane", "n-butane" };
    double[] hotZ = new double[] { 0.8, 0.2 };
    String[] coldNames = new String[] { "n-butane", "n-pentane" };
    double[] coldZ = new double[] { 0.5, 0.5 };

    Stream hotIn = feed(hotNames, hotZ, 400.0, 20.0, 1.0);
    Stream coldIn = feed(coldNames, coldZ, 300.0, 5.0, 1.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger hx =
        new neqsim.process.equipment.heatexchanger.HeatExchanger("hx1", hotIn, coldIn);
    hx.setUAvalue(100.0);
    hx.setOutStreamSpecificationNumber(1);
    hx.setOutTemperature(320.0, "K");
    hx.run();

    row("reported_hot_out", hx.getOutStream(0).getThermoSystem());
    row("reported_cold_out", hx.getOutStream(1).getThermoSystem());

    SystemInterface hot = hx.getOutStream(0).getThermoSystem();
    System.out.println("-- hot_out repairs --");
    repairs(hot, hotNames, hotZ, hot.getTemperature(), hot.getPressure(),
        hx.getOutStream(0).getFlowRate("mol/sec"));

    SystemInterface cold = hx.getOutStream(1).getThermoSystem();
    System.out.println("-- cold_out repairs --");
    repairs(cold, coldNames, coldZ, cold.getTemperature(), cold.getPressure(),
        hx.getOutStream(1).getFlowRate("mol/sec"));
  }
}
