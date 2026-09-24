// The unit operations, for the `process.*` ids.
//
// **The first probe that drives `neqsim.process` rather than `neqsim.thermo`.** Every
// other probe here builds a `ThermoSystem` and calls a flash; a unit operation is a
// `Stream` and an equipment class, and this drives those. One subcommand per unit
// operation, so each gets its own capture and its own recipe.
//
// # What is printed, and why this shape
//
// Every row prints the inlet and the outlet on the palette's own record - `n`, `z`, `P`,
// `T`, `h` - because the oracle is for the *port declaration* as much as for the answer: a
// quantity NeqSim's outlet carries and the declaration has no field for shows up here as a
// line with no counterpart. The equipment's own reported quantities (a duty, a shaft work)
// are printed beside them.
//
// **The state is not on `StreamInterface`.** That interface exposes pressure, temperature
// and flow rate; the molar enthalpy and entropy come through
// `stream.getThermoSystem().getPhase(0)`, divided by the phase's mole count. Reading the
// record off the stream interface alone silently has no `h`.
//
// **`setMixingRule(2)` is the classic rule with kij from NeqSim's database**, which is what
// `azoth_process::Stream::mixture()` resolves: PR with the databank's kij. Type 1 would be
// the same rule with every kij zeroed, and is a different fluid.
//
// **No `ProcessSystem` is needed.** `SimulationInterface.run()` defaults to
// `run(UUID.randomUUID())`, so a stream and an equipment item each run standalone - which is
// what NeqSim's own `Ejector.run` relies on. A `ProcessSystem` is required for a recycle
// alone, and that is P12's subject.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ProcessProbe.java
//     java -cp .:neqsim-f0c7436.jar ProcessProbe pump > captures/process_pump.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe splitter > captures/process_splitter.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe mixer > captures/process_mixer.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe separator > captures/process_separator.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe throttling_valve \
//       > captures/process_throttling_valve.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe heat_exchanger \
//       > captures/process_heat_exchanger.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe heater > captures/process_heater.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe cooler > captures/process_cooler.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe filter > captures/process_filter.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe pipe > captures/process_pipe.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe manifold > captures/process_manifold.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe gas_scrubber \
//       > captures/process_gas_scrubber.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe component_splitter \
//       > captures/process_component_splitter.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe tray > captures/process_column_tray.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe column > captures/process_column.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe condenser \
//       > captures/process_column_condenser.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe reboiler \
//       > captures/process_column_reboiler.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe shortcut_column \
//       > captures/process_shortcut_distillation_column.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe ejector \
//       > captures/process_ejector.tsv
//     java -cp .:neqsim-f0c7436.jar ProcessProbe stream \
//       > captures/process_stream_properties.tsv

import neqsim.process.equipment.stream.Stream;
import neqsim.process.equipment.stream.StreamInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

public class ProcessProbe {
  public static void main(String[] args) {
    String which = args.length > 0 ? args[0] : "pump";
    switch (which) {
      case "pump":
        pump();
        break;
      case "splitter":
        splitter();
        break;
      case "mixer":
        mixer();
        break;
      case "separator":
        separator();
        break;
      case "throttling_valve":
        throttlingValve();
        break;
      case "heat_exchanger":
        heatExchanger();
        break;
      case "heater":
        heatRows("Heater");
        break;
      case "cooler":
        heatRows("Cooler");
        break;
      case "filter":
        filterRows();
        break;
      case "compressor":
        isentropicRows("Compressor");
        break;
      case "expander":
        isentropicRows("Expander");
        break;
      case "pipe":
        pipeRows();
        break;
      case "manifold":
        manifoldRows();
        break;
      case "gas_scrubber":
        gasScrubberRows();
        break;
      case "component_splitter":
        componentSplitterRows();
        break;
      case "three_phase_separator":
        threePhaseSeparatorRows();
        break;
      case "tank":
        tankRows();
        break;
      case "ejector":
        ejectorRows();
        break;
      case "column":
        columnRows();
        break;
      case "column_solvers":
        columnSolverRows(args);
        break;
      case "condenser":
        condenserRows();
        break;
      case "reboiler":
        reboilerRows();
        break;
      case "tray":
        trayRows();
        phFlashRows();
        break;
      case "shortcut_column":
        shortcutColumnRows();
        break;
      case "stream":
        stream();
        break;
      default:
        throw new IllegalArgumentException("no such unit operation: " + which);
    }
  }

  /// A PR fluid on the classic mixing rule with NeqSim's own kij.
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

  static void pump() {
    // **A liquid, because a pump on a gas is not a case.** The fluid is pure n-butane at
    // 250 K and 5 bara, which is below its vapour pressure - the same fluid the kernel's
    // own balance test uses, so the two are comparing like with like.
    for (double efficiency : new double[] { 0.75, 1.0 }) {
      Stream inlet = feed(new String[] { "n-butane" }, new double[] { 1.0 }, 250.0, 5.0, 1.0);
      neqsim.process.equipment.pump.Pump pump =
          new neqsim.process.equipment.pump.Pump("p1", inlet);
      pump.setOutletPressure(20.0);
      pump.setIsentropicEfficiency(efficiency);
      pump.run();

      System.out.println("isentropic_efficiency=" + efficiency);
      print("inlet", inlet);
      print("outlet", pump.getOutletStream());
      System.out.println("power_kW=" + pump.getPower("kW"));
      // **The pump's own intermediates, which the record does not carry.** The total
      // entropy change across the machine is the second-law statement of the same head the
      // shaft work reports: zero iff the step were reversible. `Pump` exposes no getter for
      // the isentropic outlet itself, so this is the layer that localises a wrong head.
      System.out.println(
          "entropy_production_kJ_per_molK=" + pump.getEntropyProduction("kJ/molK"));
      System.out.println();
    }
  }

  /// `Heater` and `Cooler` on the **same rows**, because this is where `Cooler extends
  /// Heater` stops being a claim and becomes a measurement: the class overrides
  /// `runTransient` and a handful of getters, and not `run` - so every row below has to come
  /// back identical through either class. A row that did not would mean the single kernel the
  /// port rests on is wrong.
  ///
  /// **Six rows, one per reachable branch of `run` and one for the branch's own subtlety.**
  /// A stated outlet temperature, the same
  /// with a pressure drop, a stated duty, a negative duty, and neither - where `run` falls
  /// through to `T_in + dT` with `dT` defaulting to zero, so the drop is isothermal.
  ///
  /// **A temperature and a duty cannot both be a row.** The two setters clear each other's
  /// flags, so the class's answer to both is the order they were called in; the port refuses
  /// the pair, and a probe cannot measure a state the record cannot express.
  ///
  /// The fluid is a single-phase gas at the inlet, so the temperature and duty rows move a
  /// state rather than a phase boundary - except the negative-duty one, which is allowed to
  /// reach a two-phase answer because a real cooler does.
  static void heatRows(String which) {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.9, 0.1 };
    heatRow(which, "outlet_temperature_380", names, z, 380.0, null, 0.0);
    heatRow(which, "outlet_temperature_380_drop_2_bara", names, z, 380.0, null, 2.0);
    heatRow(which, "duty_5000_W", names, z, null, 5000.0, 0.0);
    heatRow(which, "duty_minus_5000_W", names, z, null, -5000.0, 0.0);
    heatRow(which, "no_specification", names, z, null, null, 0.0);
    // **And with a drop**, which is the row that settles what this branch is. The outlet
    // holds the inlet temperature and its enthalpy moves anyway - a real fluid's `h` depends
    // on `P` at fixed `T` - so the duty the class reports is that movement and not zero. The
    // reading that a no-specification heater is a throttling is the mistake this measures.
    heatRow(which, "no_specification_drop_2_bara", names, z, null, null, 2.0);
  }

  static void heatRow(String which, String label, String[] names, double[] z, Double outletK,
      Double dutyW, double dropBara) {
    Stream inlet = feed(names, z, 320.0, 30.0, 1.0);
    neqsim.process.equipment.heatexchanger.Heater unit =
        which.equals("Cooler") ? new neqsim.process.equipment.heatexchanger.Cooler("c1", inlet)
            : new neqsim.process.equipment.heatexchanger.Heater("h1", inlet);
    if (outletK != null) {
      unit.setOutletTemperature(outletK);
    }
    if (dutyW != null) {
      unit.setDuty(dutyW);
    }
    if (dropBara != 0.0) {
      unit.setPressureDrop(dropBara);
    }
    unit.run();

    System.out.println(label);
    print("inlet", inlet);
    print("outlet", unit.getOutletStream());
    // `getDuty()` is the class's own recomputation - `newH - oldH`, after the flash and in
    // every branch - so on the temperature rows it is the state's number and not the zero
    // the field was initialised to. On the duty rows it is the enthalpy the flash reached,
    // which is the input only if the flash lands exactly on it.
    System.out.println("duty_W=" + unit.getDuty());
    System.out.println();
  }

  /// `Filter`, whose steady state is a pressure drop and a flash - and two things the class
  /// does that a kernel written from its name would not.
  ///
  /// **It holds the temperature.** `Filter.run` reduces the pressure and runs a `TPflash`, so
  /// the outlet is at the feed's temperature and its enthalpy moves with the pressure. A
  /// throttling valve is the isenthalpic reading of the same words, and the no-specification
  /// heater's capture already prices that difference at `49.6` J/mol over two bar.
  ///
  /// **A drop larger than the inlet pressure is clamped, not refused.** `run` takes
  /// `min(max(0, dP), max(0, P_in - 1e-6 bar))` and logs a warning, so the third row lands a
  /// millionth of a bar above vacuum rather than failing.
  ///
  /// The particle-capture curve is deliberately absent from these rows because it is absent
  /// from the outlet: `concentrationLoadingModelEnabled` is false by default, and with it off
  /// `updateParticleCapturePerformance` zeroes two output fields and returns. The class's own
  /// `Cv = sqrt(dP) / massFlow` is printed, because it is the one number `run` computes from
  /// the two states that is not on either record.
  static void filterRows() {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.9, 0.1 };
    filterRow("pressure_drop_1_bar", names, z, 1.0);
    filterRow("no_pressure_drop", names, z, 0.0);
    filterRow("pressure_drop_35_bar_past_the_inlet", names, z, 35.0);
  }

  static void filterRow(String label, String[] names, double[] z, double dropBara) {
    Stream inlet = feed(names, z, 320.0, 30.0, 1.0);
    neqsim.process.equipment.filter.Filter unit =
        new neqsim.process.equipment.filter.Filter("f1", inlet);
    unit.setDeltaP(dropBara);
    unit.run();

    System.out.println(label);
    print("inlet", inlet);
    print("outlet", unit.getOutletStream());
    // The drop the class *applied*, which differs from the one set on the clamped row, and
    // the coefficient it derives from the mass flow it sees.
    System.out.println("applied_deltaP_bara=" + unit.getDeltaP());
    System.out.println("cv_sqrt_bar_per_kg_per_hr=" + unit.getCvFactor());
    System.out.println();
  }

  /// `Compressor` and `Expander`, which are one machine's route with the efficiency on the
  /// other side of the division.
  ///
  /// **The steady state is the isentropic step, and the inlet entropy is *derived*.** Both
  /// classes read `getEntropy()` off the inlet system, run a `PSflash` at the outlet pressure
  /// on it, and take the enthalpy they land on as the reversible outlet:
  /// `Compressor.run`'s no-chart branch then divides the difference by the isentropic
  /// efficiency, and `Expander.run`'s multiplies by it - which is the same rule, because an
  /// expansion's difference is negative and dividing would make the machine beat the
  /// reversible one.
  ///
  /// **Four rows, two per class**: the isentropic limit at an efficiency of one, and a real
  /// machine at `0.75`. At one, the compressor's shaft work *is* the isentropic difference and
  /// the entropy production is zero - so that row is where the oracle and the port can be
  /// compared without an efficiency convention in between.
  ///
  /// `getEntropyProduction` is the second-law statement of the same step, and the only other
  /// number either class exposes about its interior.
  static void isentropicRows(String which) {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.9, 0.1 };
    for (double efficiency : new double[] { 0.75, 1.0 }) {
      boolean compressing = which.equals("Compressor");
      double inletBara = compressing ? 30.0 : 60.0;
      double outletBara = compressing ? 60.0 : 30.0;
      Stream inlet = feed(names, z, 320.0, inletBara, 1.0);
      neqsim.process.equipment.compressor.Compressor unit = compressing
          ? new neqsim.process.equipment.compressor.Compressor("c1", inlet)
          : new neqsim.process.equipment.expander.Expander("e1", inlet);
      unit.setOutletPressure(outletBara);
      unit.setIsentropicEfficiency(efficiency);
      unit.run();

      System.out.println("isentropic_efficiency=" + efficiency);
      print("inlet", inlet);
      print("outlet", unit.getOutletStream());
      System.out.println("power_kW=" + unit.getPower("kW"));
      System.out.println("entropy_production_kJ_per_molK=" + unit.getEntropyProduction("kJ/molK"));
      System.out.println();
    }
  }

  /// `AdiabaticPipe`, which is a different shape of unit operation: **its outlet pressure is
  /// solved from the geometry**, in a loop, rather than stated.
  ///
  /// `run` iterates `calcPressureOut()` against the state it is evaluating at, to `1e-2` bar or
  /// twenty-five passes, and the arithmetic inside is **a different equation per phase**:
  /// a compressible `P1^2 - P2^2` form for a gas, and Darcy-Weisbach for a liquid - where the
  /// liquid branch deliberately recomputes the velocity from the *physical-properties* density
  /// rather than the cubic's volume, because "cubic EOS liquid volumes are inaccurate for polar
  /// fluids (e.g. water)".
  ///
  /// **Three rows, and the third is the asymmetry.** A gas line, a hydrocarbon liquid line, and
  /// a water line. The probe prints both densities and both viscosities the class can reach, so
  /// that the port's choice is compared against the one the arithmetic actually used:
  ///
  /// * `phase.getPhysicalProperties().getKinematicViscosity()` is what `calcPressureOut` reads.
  ///   `GasPhysicalProperties` and `LiquidPhysicalProperties` both default their viscosity to
  ///   `PFCTViscosityMethodHeavyOil` - for every mixture, water included - which is the
  ///   correlation `eos.viscosity` ports.
  /// * `system.getViscosity("kg/msec")` is a **different route** with its own dispatch, and it
  ///   is the one that gives water `8.55e-4` against the correlation's `5.30e-4`. That
  ///   divergence is real and belongs to a call the pipe does not make.
  static void pipeRows() {
    // **The gas is methane/CO2 and not methane/n-butane, and that is a measurement.**
    // 0.9/0.1 methane/n-butane at 320 K and 30 bara is the one state in this tier where
    // azoth's cubic lands on a different root from NeqSim's - `Z = 0.87296` against
    // `0.91976`, recorded in `crates/azoth-process/tests/stream.rs` - and the pipe reads
    // `Z` for its velocity and its `P1^2 - P2^2` term, so a gas row on that fluid would be
    // measuring `eos.pt_flash` rather than the hydraulics. Methane/CO2 at 300 K and 50 bara
    // is a state where the two roots agree.
    pipeRow("gas_methane_co2_1000m", new String[] { "methane", "CO2" }, new double[] { 0.7, 0.3 },
        300.0, 50.0, 1.0, 1000.0, 0.1, 1.0e-5);
    pipeRow("liquid_n_butane_1000m", new String[] { "n-butane" }, new double[] { 1.0 }, 300.0,
        20.0, 1.0, 1000.0, 0.1, 1.0e-5);
    pipeRow("liquid_water_1000m", new String[] { "water" }, new double[] { 1.0 }, 300.0, 5.0, 1.0,
        1000.0, 0.1, 1.0e-5);
  }

  static void pipeRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, double molPerSecond, double length, double diameter,
      double roughness) {
    Stream inlet = feed(names, z, temperatureK, pressureBara, molPerSecond);
    neqsim.process.equipment.pipeline.AdiabaticPipe pipe =
        new neqsim.process.equipment.pipeline.AdiabaticPipe("pipe1", inlet);
    pipe.setLength(length);
    pipe.setDiameter(diameter);
    pipe.setPipeWallRoughness(roughness);
    pipe.run();

    System.out.println(label);
    print("inlet", inlet);
    print("outlet", pipe.getOutletStream());
    System.out.println("velocity_m_per_s=" + pipe.getVelocity());
    System.out.println("reynolds_number=" + pipe.getReynoldsNumber());
    System.out.println("friction_factor=" + pipe.getFrictionFactor());
    System.out.println("flow_regime=" + pipe.getFlowRegime());
    System.out.println("pressure_drop_bara=" + pipe.getPressureDrop());

    // **The state the arithmetic was evaluated at, and the two routes to a viscosity.** The
    // outlet's system because that is the object `calcPressureOut` reads, and both routes
    // because the port has to use the first and the second is where the known water
    // divergence lives.
    SystemInterface f = pipe.getOutletStream().getThermoSystem();
    System.out.println("phase_type=" + f.getPhase(0).getType());
    System.out.println("phase_total_volume=" + f.getPhase(0).getTotalVolume());
    System.out.println("phase_z=" + f.getPhase(0).getZ());
    System.out.println("phase_density_kg_per_m3=" + f.getPhase(0).getPhysicalProperties().getDensity());
    System.out.println("system_density_kg_per_m3=" + f.getDensity("kg/m3"));
    System.out.println("phase_kinematic_viscosity_m2_per_s="
        + f.getPhase(0).getPhysicalProperties().getKinematicViscosity());
    System.out.println("phase_dynamic_viscosity_kg_per_msec="
        + f.getPhase(0).getPhysicalProperties().getViscosity());
    System.out.println("system_viscosity_kg_per_msec=" + f.getViscosity("kg/msec"));
    System.out.println();
  }

  /// `Manifold`: a mixer and a splitter in one, with a minimum-flow rule between them.
  ///
  /// **Its `run` is four statements and the thirds and fourth are the composition**:
  /// `propagateMinimumFlow()` pushes the manifold's threshold onto both, `localmixer.run()`
  /// joins the feeds, `refreshLocalSplitter()` re-attaches the splitter's inlet to the
  /// mixture, and `localsplitter.run()` divides it. Nothing else is the manifold's own.
  ///
  /// **The threshold's default is `1e-20`** - `ProcessEquipmentBaseClass.DEFAULT_MINIMUM_FLOW`
  /// - so the rule is inert for any feed a case states, and a feed at or below it is dropped
  /// by the mixer rather than joined. The palette declares no minimum-flow parameter, so the
  /// rows below exercise the default only.
  ///
  /// **The palette entry declares one outlet for a class whose outlet count is
  /// `splitFactors.length`**, so the entry is corrected with this kernel: `outlet` gains a
  /// `many` multiplicity and the entry a `split_factors` parameter, as `unit_ops.splitter`
  /// has. The rows carry two outlets and three, which is what makes the correction a
  /// measurement.
  static void manifoldRows() {
    manifoldRow("two_feeds_two_outlets", new double[] { 0.25, 0.75 });
    manifoldRow("two_feeds_three_outlets", new double[] { 0.2, 0.3, 0.5 });
    manifoldRow("a_zero_flow_feed_is_dropped", new double[] { 0.5, 0.5 });
  }

  static void manifoldRow(String label, double[] factors) {
    Stream first = feed(new String[] { "methane", "n-butane" }, new double[] { 0.9, 0.1 }, 320.0, 30.0,
        1.0);
    Stream second = feed(new String[] { "methane", "n-butane" }, new double[] { 0.3, 0.7 }, 300.0, 10.0,
        2.0);

    neqsim.process.equipment.manifold.Manifold manifold =
        new neqsim.process.equipment.manifold.Manifold("mf1");
    manifold.addStream(first);
    manifold.addStream(second);
    manifold.setSplitFactors(factors);
    if (label.equals("a_zero_flow_feed_is_dropped")) {
      // **A stated zero, not an underflow.** `Mixer.mixStream` skips an inlet whose
      // `getFlowRate("kg/hr") <= getMinimumFlow()`, and the threshold's default is `1e-20`
      // kg/hr - so the rule is inert for every flow a case can state, and what this row
      // exercises is the *zero* end of it. A flow of `1e-21 kg/s` was tried first and read
      // back as `0.0`, which is NeqSim's own underflow rather than its comparison.
      first.setFlowRate(0.0, "kg/sec");
      first.run();
    }
    manifold.run();

    StringBuilder given = new StringBuilder("split_factors=");
    for (int i = 0; i < factors.length; i++) {
      given.append(factors[i]).append(i + 1 < factors.length ? " " : "");
    }
    System.out.println(label);
    System.out.println(given);
    print("feed0", first);
    print("feed1", second);
    print("product", manifold.getMixedStream());
    for (int i = 0; i < factors.length; i++) {
      print("products" + i, manifold.getSplitStream(i));
    }
    System.out.println("outlet_count=" + manifold.getNumberOfOutputStreams());
    System.out.println();
  }

  /// `GasScrubber`, whose stream side is `Separator.run` and whose own arithmetic is a
  /// mechanical capacity metric.
  ///
  /// **The class does not override `run`** - measured, `grep -c "public void run"` returns
  /// zero - so the first three rows are the separator's own rows through another class, and
  /// the capture records that they are the same numbers.
  ///
  /// **The fourth row is the metric, and it is what the port leaves out.**
  /// `getCapacityUtilization` needs an internal diameter and a design gas load factor; the
  /// palette entry declares neither, and it is a statement about whether the *vessel* is
  /// big enough rather than about the stream.
  static void gasScrubberRows() {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.7, 0.3 };
    scrubberRow("pressure_drop_bara=0 heat_input_W=0 gas_in_liquid=0", names, z, 0.0, null, 0.0,
        false);
    scrubberRow("pressure_drop_bara=2 heat_input_W=0 gas_in_liquid=0", names, z, 2.0, null, 0.0,
        false);
    scrubberRow("pressure_drop_bara=0 heat_input_W=0 gas_in_liquid=0.05", names, z, 0.0, null, 0.05,
        false);
    scrubberRow("capacity_utilization_dn1000_k_007", names, z, 0.0, null, 0.0, true);
  }

  static void scrubberRow(String label, String[] names, double[] z, double dropBara,
      Double heatInputW, double gasInLiquid, boolean withCapacity) {
    Stream inlet = feed(names, z, 300.0, 20.0, 1.0);
    neqsim.process.equipment.separator.GasScrubber scrubber =
        new neqsim.process.equipment.separator.GasScrubber("gs1", inlet);
    if (dropBara != 0.0) {
      scrubber.setPressureDrop(dropBara);
    }
    if (heatInputW != null) {
      scrubber.setHeatInput(heatInputW);
    }
    if (gasInLiquid != 0.0) {
      scrubber.setEntrainment(gasInLiquid, "mole", "feed", "gas", "liquid");
    }
    if (withCapacity) {
      // The two mechanical parameters the metric needs, and neither is on the palette.
      scrubber.setInternalDiameter(1.0);
      scrubber.setDesignGasLoadFactor(0.07);
    }
    scrubber.run();

    System.out.println(label);
    print("feed", inlet);
    print("vapour", scrubber.getGasOutStream());
    print("liquid", scrubber.getLiquidOutStream());
    if (withCapacity) {
      System.out.println("capacity_utilization=" + scrubber.getCapacityUtilization());
    }
    System.out.println();
  }

  /// **The column, and the one state NeqSim's own tests carry.** `NaphtaliSandholmPublishedStateTest`
  /// (issue #3698) is a deethanizer: 8 stages, the feed on tray 6, a C1-C6 mixture at
  /// 10,000 kg/hr, the condenser at 0 C and the reboiler at 80 C, 24 bara at the top and 25 at
  /// the bottom. The class is written on `SystemSrkEos` there; **this drives it on PR**, because
  /// `azoth_process::Stream::mixture()` resolves a Peng-Robinson fluid and has no other route,
  /// so a row on SRK would not be a case azoth could be held to.
  ///
  /// **The solve is `solveSequential` with the class's default `DIRECT_SUBSTITUTION`**, and its
  /// own stopping rule: the mean tray temperature change under `setTemperatureTolerance`, with a
  /// mass and an energy gate behind it. Every residual below is the class's own, because where a
  /// solve stops is a property of the column and not of the port.
  ///
  /// **The profile is the measurement.** A column's answer is not a scalar, so each tray prints
  /// its temperature, its pressure and both its traffic rates - and one tray prints its whole
  /// record, because a profile that agrees in the temperatures and not in the flows has not
  /// agreed.
  static void columnRows() {
    String[] light = new String[] { "methane", "ethane", "propane", "n-butane" };
    // The deethanizer, on PR rather than the test's SRK.
    columnRow("deethanizer_pr_hard_cap_80", new String[] { "methane", "ethane", "propane", "i-butane",
        "n-butane", "i-pentane", "n-pentane", "n-hexane" },
        new double[] { 0.22, 0.34, 0.20, 0.08, 0.08, 0.03, 0.03, 0.02 }, 283.15, 25.0, 10000.0, 8, 6,
        0.0, 80.0, 24.0, 25.0, 1.0e-5, 80, true, false);
    // The same column with a **soft** iteration limit, which is what it takes to converge: on PR
    // the deethanizer's own `setMaxNumberOfIterations(80, true)` hard cap stops it at a mean
    // tray-temperature change of `7.3e-3` K, so the hard-capped rows above are a *measurement of
    // the cap* rather than of the column. `AUTO`'s ladder is what the soft limit reaches for.
    columnRow("deethanizer_pr_soft_limit", new String[] { "methane", "ethane", "propane",
        "i-butane", "n-butane", "i-pentane", "n-pentane", "n-hexane" },
        new double[] { 0.22, 0.34, 0.20, 0.08, 0.08, 0.03, 0.03, 0.02 }, 283.15, 25.0, 10000.0, 8, 6,
        0.0, 80.0, 24.0, 25.0, 1.0e-5, 200, true, false);
    // **A small binary**, the row a reader can follow by hand: four stages, the feed on tray 2.
    // This is the one that converges rigorously, and it is the port's primary oracle.
    columnRow("binary_methane_butane_4_stages", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0, 19.0, 20.0, 1.0e-6, 200,
        true, false);
    // The same binary with a **looser tolerance**, which is a different measurement rather than
    // a sloppier version of the same one: where a solve stops is what its answer is.
    columnRow("binary_methane_butane_loose", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0, 19.0, 20.0, 1.0e-2, 200,
        true, false);
    // **The five specification types, one row each**, on the same binary column. The first three
    // are *adjusted* - `needsAdjustment` is true for everything but the last two - so each is
    // driven by an outer secant on the end's temperature; the last two are applied directly to
    // the end itself. That split is the class's own, and it is why a reflux ratio is a
    // specification at the end rather than a duty parameter.
    specificationRow("spec_top_purity_0_98_methane", "methane", 0.98, null, null);
    // **The one adjustable row that does not converge**, kept as evidence: the secant on the
    // condenser temperature drives this column to `D = 7.25` - the whole feed - and stops.
    specificationRow("spec_top_recovery_0_97_not_converged", "methane", 0.97, true, null);
    // **The flow-rate target is in mol/hr, which is `ColumnSpecification.defaultTargetUnit`'s
    // own unit for this type** - the class's constructor says so, and a target in mol/s would
    // be a column asked for a millionth of its flow.
    specificationRow("spec_top_flow_rate_14000_mol_per_hour", null, 14000.0, null, null);
    specificationRow("spec_top_reflux_ratio_1_5", null, 1.5, null, null);
    // **A duty specification under a temperature pin, which is inert** - and measured. The
    // condenser's own temperature is stated, so its flash is a TP one and `setHeatInput` is
    // never read: the row reports the pinned `-21323.04` W for a specification of `-20000`. The
    // residual is `0.0` because a DUTY is not *adjusted*, so nothing in the class reports this.
    specificationRow("spec_top_duty_under_a_pin_is_inert", null, -20000.0, null, null);
    // **The same duty without a pin**, where it bites: the condenser's temperature is the
    // flash's own answer and the duty lands on the target to `4e-5` W.
    specificationRow("spec_top_duty_minus_20000", null, -20000.0, null, null);
    // **The same machine at the other location.** `setBottomSpecification` drives the
    // *reboiler's* temperature where the top's drives the condenser's, and its product is the
    // bottom's - so this is the row that says a location is a degree of freedom rather than a
    // mirrored spelling of the top's. The condenser stays pinned at -20 C here.
    specificationRow("spec_bottom_purity_0_98_n_butane", "n-butane", 0.98, null, true);
  }

  static void specificationRow(String label, String component, double target, Boolean recovery,
      Boolean bottom) {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.5, 0.5 };
    SystemInterface fluid = new SystemPrEos(300.0, 20.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("column feed", fluid);
    inlet.setFlowRate(1000.0, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn("col1", 4, true, true);
    column.addFeedStream(inlet, 2);
    column.setReboilerTemperature(100.0, "C");
    column.setTopPressure(19.0);
    column.setBottomPressure(20.0);
    column.setTemperatureTolerance(1.0e-6);
    column.setMaxNumberOfIterations(200, true);
    neqsim.process.equipment.distillation.ColumnSpecification.SpecificationType type;
    if (label.contains("reflux")) {
      type = neqsim.process.equipment.distillation.ColumnSpecification.SpecificationType.REFLUX_RATIO;
    } else if (label.contains("duty")) {
      type = neqsim.process.equipment.distillation.ColumnSpecification.SpecificationType.DUTY;
    } else if (component != null) {
      // A purity or a recovery, which the type tells apart by whether it needs a component name.
      type = recovery != null
          ? neqsim.process.equipment.distillation.ColumnSpecification.SpecificationType.COMPONENT_RECOVERY
          : neqsim.process.equipment.distillation.ColumnSpecification.SpecificationType.PRODUCT_PURITY;
    } else {
      type = neqsim.process.equipment.distillation.ColumnSpecification.SpecificationType.PRODUCT_FLOW_RATE;
    }
    // **A location follows its slot.** `validateColumnSpecification` refuses a `TOP` spec handed
    // to `setBottomSpecification` and names the setter to use, so the field is a spelling of the
    // slot rather than a third degree of freedom - which is why the port's inputs name the end.
    neqsim.process.equipment.distillation.ColumnSpecification.ProductLocation location =
        Boolean.TRUE.equals(bottom)
            ? neqsim.process.equipment.distillation.ColumnSpecification.ProductLocation.BOTTOM
            : neqsim.process.equipment.distillation.ColumnSpecification.ProductLocation.TOP;
    neqsim.process.equipment.distillation.ColumnSpecification specification =
        new neqsim.process.equipment.distillation.ColumnSpecification(type, location, target,
            component);
    if (Boolean.TRUE.equals(bottom)) {
      column.setBottomSpecification(specification);
    } else {
      column.setTopSpecification(specification);
    }
    if (!label.contains("reflux") && !label.contains("minus_20000")) {
      column.setCondenserTemperature(-20.0, "C");
    }
    column.run();

    System.out.println(label);
    print("feed", inlet);
    System.out.println("solved=" + column.solved());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("top_spec_residual=" + column.getLastTopSpecificationResidual());
    if (Boolean.TRUE.equals(bottom)) {
      System.out.println("bottom_spec_residual=" + column.getLastBottomSpecificationResidual());
    }
    System.out.println("condenser_temperature_K=" + column.getCondenser().getTemperature());
    // The profile, which is the port's own interior: a specification's effect is a *state*,
    // so a row that printed only its products could not localise a divergence.
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println("tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + column.getTray(i).getPressure());
      System.out.println("tray" + i + "_gas_n=" + column.getTray(i).getGasOutStream().getFlowRate("mol/sec"));
      System.out.println("tray" + i + "_liquid_n=" + column.getTray(i).getLiquidOutStream().getFlowRate("mol/sec"));
    }
    print("distillate", column.getGasOutStream());
    print("bottoms", column.getLiquidOutStream());
    System.out.println("condenser_duty_W=" + column.getCondenser().getDuty());
    System.out.println("reboiler_duty_W=" + column.getReboiler().getDuty());
    System.out.println();
  }

  /// **The ladder, measured rather than described.** `ColumnSolverFactory` carries ten
  /// strategies and `AUTO` is one of them - a probe over candidates rather than a method - so a
  /// caller who states nothing gets whichever of them the column's own conditions select. These
  /// rows state each one explicitly on the same two columns, which is what decides whether the
  /// rungs are different physics or different paths to one answer.
  static void columnSolverRows(String[] only) {
    String[] binary = new String[] { "methane", "n-butane" };
    double[] binaryZ = new double[] { 0.5, 0.5 };
    String[] deethanizer = new String[] { "methane", "ethane", "propane", "i-butane", "n-butane",
        "i-pentane", "n-pentane", "n-hexane" };
    double[] deethanizerZ = new double[] { 0.22, 0.34, 0.20, 0.08, 0.08, 0.03, 0.03, 0.02 };
    String[] solvers = new String[] { "DIRECT_SUBSTITUTION", "DAMPED_SUBSTITUTION", "INSIDE_OUT",
        "MATRIX_INSIDE_OUT", "WEGSTEIN", "SUM_RATES", "NEWTON", "NAPHTALI_SANDHOLM",
        "MESH_RESIDUAL", "AUTO" };

    for (String solver : solvers) {
      // `args[0]` is this subcommand; a solver name given after it narrows the run.
      if (only.length > 1 && !solver.equalsIgnoreCase(only[1])) {
        continue;
      }
      solverRow("binary_methane_butane_" + solver.toLowerCase(), binary, binaryZ, 300.0, 20.0,
          1000.0, 4, 2, -20.0, 100.0, 19.0, 20.0, 1.0e-6, 200, true, false, solver);
    }
    // **The deethanizer is captured under the one strategy that converges it.** Every other
    // strategy falls back to `DAMPED_SUBSTITUTION` and stops at the iteration cap - which the
    // `deethanizer_pr_soft_limit` row of the column capture already records for
    // `DIRECT_SUBSTITUTION` - and two of them cost more than a sweep can pay: measured here,
    // `MESH_RESIDUAL` takes 594 s for the pair and `AUTO` on this column had not returned after
    // 900 s. The failure rows are a measurement of the cap rather than of the column.
    if (only.length < 2 || solvers[7].equalsIgnoreCase(only[1])) {
      solverRow("deethanizer_naphtali_sandholm", deethanizer, deethanizerZ, 283.15, 25.0, 10000.0,
          8, 6, 0.0, 80.0, 24.0, 25.0, 1.0e-5, 200, true, false, "NAPHTALI_SANDHOLM");
    }
  }

  static void columnRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, double condenserC,
      double reboilerC, double topBara, double bottomBara, double tolerance, int maxIterations,
      boolean condenser, boolean totalCondenser) {
    runColumn(label, names, z, feedTemperatureK, feedPressureBara, kgPerHour, stages, feedTray,
        condenserC, reboilerC, topBara, bottomBara, tolerance, maxIterations, condenser,
        totalCondenser, null);
  }

  /// **The same column under a stated `SolverType`.** `setSolverType` is explicit rather than a
  /// fallback, so the row reports which strategy was actually used and what that strategy's own
  /// residuals are - the two things that decide whether a second solve agrees with the first.
  static void solverRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, double condenserC,
      double reboilerC, double topBara, double bottomBara, double tolerance, int maxIterations,
      boolean condenser, boolean totalCondenser, String solver) {
    runColumn(label, names, z, feedTemperatureK, feedPressureBara, kgPerHour, stages, feedTray,
        condenserC, reboilerC, topBara, bottomBara, tolerance, maxIterations, condenser,
        totalCondenser, solver);
  }

  static void runColumn(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, double condenserC,
      double reboilerC, double topBara, double bottomBara, double tolerance, int maxIterations,
      boolean condenser, boolean totalCondenser, String solver) {
    SystemInterface fluid = new SystemPrEos(feedTemperatureK, feedPressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("column feed", fluid);
    inlet.setFlowRate(kgPerHour, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn("col1", stages, true, condenser);
    column.addFeedStream(inlet, feedTray);
    column.setCondenserTemperature(condenserC, "C");
    column.setReboilerTemperature(reboilerC, "C");
    column.setTopPressure(topBara);
    column.setBottomPressure(bottomBara);
    column.setTemperatureTolerance(tolerance);
    if (totalCondenser) {
      column.setCondenserMode(
          neqsim.process.equipment.distillation.DistillationColumn.CondenserMode.TOTAL);
      column.setCondenserRefluxRatio(1.5);
    }
    column.setMaxNumberOfIterations(maxIterations, true);
    if (solver != null) {
      column.setSolverType(
          neqsim.process.equipment.distillation.DistillationColumn.SolverType.valueOf(solver));
    }
    column.run();

    System.out.println(label);
    System.out.println("stages=" + stages);
    System.out.println("feed_tray=" + feedTray);
    System.out.println("temperature_tolerance=" + tolerance);
    System.out.println("total_condenser=" + totalCondenser);
    print("feed", inlet);
    System.out.println("feed_mol_per_sec=" + inlet.getFlowRate("mol/sec"));
    System.out.println("feed_kg_per_hour=" + inlet.getFlowRate("kg/hr"));
    System.out.println("tray_count=" + column.getNumberOfTrays());
    if (solver != null) {
      System.out.println("solver_requested=" + solver);
    }
    System.out.println("solved=" + column.solved());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("solver=" + column.getLastSolverTypeUsed());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    System.out.println("temperature_step_residual=" + column.getLastAppliedTemperatureStepResidual());
    System.out.println("mass_residual=" + column.getLastMassResidual());
    System.out.println("energy_residual=" + column.getLastEnergyResidual());
    System.out.println("internal_traffic_ratio=" + column.getLastInternalTrafficRatio());
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println("tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + column.getTray(i).getPressure());
      System.out.println("tray" + i + "_gas_n=" + column.getTray(i).getGasOutStream().getFlowRate("mol/sec"));
      System.out.println("tray" + i + "_liquid_n=" + column.getTray(i).getLiquidOutStream().getFlowRate("mol/sec"));
    }
    print("distillate", column.getGasOutStream());
    print("bottoms", column.getLiquidOutStream());
    System.out.println("condenser_duty_W=" + column.getCondenser().getDuty());
    System.out.println("reboiler_duty_W=" + column.getReboiler().getDuty());
    System.out.println();
  }

  /// **The condenser, which is a tray with a reflux split.** `Condenser extends SimpleTray`, so
  /// an equilibrium partial condenser is `super.run` and nothing else; the other two modes add
  /// a specification.
  ///
  /// **Three modes, and each reaches its products differently.** With no reflux set the run is
  /// the tray's own flash and the outlets are its phases. With a ratio set the flash is
  /// `PVrefluxflash(refluxRatio, 0)` - a *temperature* search for the state whose vapour
  /// fraction satisfies `R = 1/beta - 1` - and the outlets are still the tray's phases. With
  /// `setTotalCondenser(true)` the flash is a **bubble-point** one and the reflux is a
  /// `Splitter` over the condensate at `R/(1+R)`; there the *distillate* is reached through
  /// `getGasOutStream()` and the reflux through `getLiquidOutStream()`, even though both are
  /// liquid.
  ///
  /// **The duty is the outlet enthalpy less the inlets',** and the condenser overrides
  /// `getMaterialOutletEnthalpy` to add the separate liquid product when its mode has one.
  static void condenserRows() {
    String[] names = new String[] { "methane", "ethane", "propane", "n-butane" };
    double[] z = new double[] { 0.1, 0.3, 0.4, 0.2 };
    condenserRow("equilibrium_partial_no_reflux", names, z, 300.0, 20.0, null, false, null, null);
    condenserRow("partial_reflux_ratio_1_5", names, z, 300.0, 20.0, 1.5, false, null, null);
    condenserRow("total_reflux_ratio_1_5", names, z, 300.0, 20.0, 1.5, true, null, null);
    condenserRow("total_reflux_ratio_zero", names, z, 300.0, 20.0, 0.0, true, null, null);
    condenserRow("liquid_reflux_split", names, z, 300.0, 20.0, null, false, 0.5, "mol/sec");
  }

  static void condenserRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, Double refluxRatio, boolean total, Double fixedReflux,
      String fixedRefluxUnit) {
    Stream inlet = feed(names, z, temperatureK, pressureBara, 1.0);
    neqsim.process.equipment.distillation.Condenser condenser =
        new neqsim.process.equipment.distillation.Condenser("cond1");
    condenser.addStream(inlet);
    if (total) {
      condenser.setTotalCondenser(true);
    }
    if (refluxRatio != null) {
      condenser.setRefluxRatio(refluxRatio);
    }
    if (fixedReflux != null) {
      condenser.setSeparation_with_liquid_reflux(true, fixedReflux, fixedRefluxUnit);
    }
    condenser.run();

    System.out.println(label);
    print("feed", inlet);
    System.out.println("total_condenser=" + condenser.isTotalCondenser());
    System.out.println("reflux_is_set=" + condenser.isRefluxSet());
    System.out.println("reflux_ratio=" + condenser.getRefluxRatio());
    System.out.println("duty_W=" + condenser.getDuty());
    System.out.println("outlet_temperature_K=" + condenser.getOutletStream().getTemperature());
    print("gas_out", condenser.getGasOutStream());
    print("liquid_out", condenser.getLiquidOutStream());
    if (condenser.getLiquidProductStream() != null) {
      print("liquid_product", condenser.getLiquidProductStream());
    }
    System.out.println();
  }

  /// **The reboiler, which is a tray with a boilup ratio.** `Reboiler extends SimpleTray` and
  /// adds one branch: with no ratio set the run is the tray's flash, and with one it is
  /// `PVrefluxflash(refluxRatio, 1)` - the same temperature search as the condenser's, asking
  /// for the *liquid's* fraction rather than the vapour's, which makes the ratio the boilup
  /// `V/B` rather than the reflux `L/D`.
  static void reboilerRows() {
    String[] names = new String[] { "methane", "ethane", "propane", "n-butane" };
    double[] z = new double[] { 0.1, 0.3, 0.4, 0.2 };
    reboilerRow("equilibrium_no_ratio", names, z, 320.0, 25.0, null);
    reboilerRow("vapor_boilup_ratio_2_0", names, z, 320.0, 25.0, 2.0);
    reboilerRow("vapor_boilup_ratio_0_5", names, z, 320.0, 25.0, 0.5);
  }

  static void reboilerRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, Double boilupRatio) {
    Stream inlet = feed(names, z, temperatureK, pressureBara, 1.0);
    neqsim.process.equipment.distillation.Reboiler reboiler =
        new neqsim.process.equipment.distillation.Reboiler("reb1");
    reboiler.addStream(inlet);
    if (boilupRatio != null) {
      reboiler.setRefluxRatio(boilupRatio);
    }
    reboiler.run();

    System.out.println(label);
    print("feed", inlet);
    System.out.println("reflux_is_set=" + reboiler.isRefluxSet());
    System.out.println("boilup_ratio=" + reboiler.getRefluxRatio());
    System.out.println("duty_W=" + reboiler.getDuty());
    System.out.println("outlet_temperature_K=" + reboiler.getOutletStream().getTemperature());
    print("gas_out", reboiler.getGasOutStream());
    print("liquid_out", reboiler.getLiquidOutStream());
    System.out.println();
  }

  /// **The column's stage, driven on its own.** `SimpleTray` is a `Mixer` plus a flash: it
  /// mixes its inlets, adds `heatInput` to their total enthalpy, and flashes - at a stated
  /// outlet temperature when one is set, otherwise to the enthalpy it holds.
  ///
  /// **`calcMixStreamEnthalpy` is overridden here and that is where the heat input enters.**
  /// `Mixer`'s is the inlets' enthalpy alone; the tray's starts from `heatInput`, subtracts
  /// an energy stream's duty when one is attached, and adds each flowing inlet's - so a tray
  /// with a heat input is a tray whose flash is at a raised enthalpy.
  ///
  /// **The outlets are the flash's phases, not two re-flashed streams.** `getGasOutStream`
  /// and `getLiquidOutStream` extract phase 0's gas and the first liquid-like phase from the
  /// mixed system, scale each by its draw fraction, and report a **zero-flow stream with the
  /// tray's composition** where the phase is absent. That last part is what the one-outlet
  /// rows below measure: a subcooled feed gives a zero-flow vapour, not an error and not a
  /// missing row.
  ///
  /// **`trayPressure` is a bara magnitude and negative means "the inlet's"**, which is the
  /// class's own sentinel rather than an absent value.
  static void trayRows() {
    String[] names = new String[] { "methane", "ethane", "propane", "n-butane" };
    double[] z = new double[] { 0.1, 0.3, 0.4, 0.2 };
    // One two-phase feed, the equilibrium split at the tray's own pressure.
    trayRow("one_two_phase_feed", new double[][] { { 300.0, 20.0, 1.0 } }, names, z, -1.0, null,
        0.0);
    // **Two inlets, which is what a column stage receives**: vapour rising from below and
    // liquid falling from above, at different compositions.
    trayRow("vapour_below_liquid_above", new double[][] { { 320.0, 20.0, 1.0 }, { 290.0, 20.0, 1.0 } },
        names, z, -1.0, null, 0.0);
    // A stated outlet temperature takes the class's `TPflash` branch instead of the
    // enthalpy one.
    trayRow("outlet_temperature_320", new double[][] { { 300.0, 20.0, 1.0 } }, names, z, -1.0,
        320.0, 0.0);
    // A heat input, which is the enthalpy branch with a raised target.
    trayRow("heat_input_5000_W", new double[][] { { 300.0, 20.0, 1.0 } }, names, z, -1.0, null,
        5000.0);
    // A tray pressure of its own, which overrides the inlet's.
    trayRow("tray_pressure_15_bara", new double[][] { { 300.0, 20.0, 1.0 } }, names, z, 15.0, null,
        0.0);
    // **A subcooled feed**, where the vapour outlet is the zero-flow template.
    trayRow("subcooled_feed_one_phase", new double[][] { { 220.0, 20.0, 1.0 } }, names, z, -1.0,
        null, 0.0);
    // **A superheated feed**, the mirror: the liquid outlet is the template.
    trayRow("superheated_feed_one_phase", new double[][] { { 400.0, 20.0, 1.0 } }, names, z, -1.0,
        null, 0.0);
    // **The same duty at two mol/s, which settles what `PHflash(h, 0)`'s unit is.** If the
    // flash divides the total enthalpy by the moles, the temperature rise halves against the
    // one-mol row above; if it reads `h` as a molar quantity, the rise is the same. The whole
    // tier's duty handling rests on which of those it is, and no unit-flow row can tell them
    // apart.
    trayRow("heat_input_5000_W_two_mol_per_second", new double[][] { { 300.0, 20.0, 2.0 } }, names,
        z, -1.0, null, 5000.0);
  }

  /// **Whether NeqSim's own `PHflash` honours the enthalpy it is given.**
  ///
  /// The tray rows show a flash whose returned temperature's phase outlets weigh to an
  /// enthalpy that differs from the one the flash was asked for - by `0.0`, `5.2`, `-0.4`,
  /// `-27.3` and `+179.5` J/mol across five rows. A difference that large on a row whose
  /// step is only `+2500` J/mol cannot be a property of either library's equations of state,
  /// so this asks the question directly: a fresh fluid, a stated molar enthalpy, and the
  /// state's own enthalpy read back. Nothing about a tray is involved.
  static void phFlashRows() {
    String[] names = new String[] { "methane", "ethane", "propane", "n-butane" };
    double[] z = new double[] { 0.1, 0.3, 0.4, 0.2 };
    for (double duty : new double[] { 0.0, 2500.0, 5000.0, -2000.0 }) {
      Stream inlet = feed(names, z, 300.0, 20.0, 1.0);
      SystemInterface fluid = inlet.getThermoSystem().clone();
      double moles = fluid.getTotalNumberOfMoles();
      double requested = fluid.getEnthalpy() / moles + duty;
      neqsim.thermodynamicoperations.ThermodynamicOperations ops =
          new neqsim.thermodynamicoperations.ThermodynamicOperations(fluid);
      ops.PHflash(requested * moles, 0);
      fluid.init(2);
      double got = fluid.getEnthalpy() / fluid.getTotalNumberOfMoles();
      System.out.println("duty_J_per_mol=" + duty);
      System.out.println("requested_h=" + requested);
      System.out.println("returned_h=" + got);
      System.out.println("gap=" + (got - requested));
      System.out.println("temperature_K=" + fluid.getTemperature());
      System.out.println("phases=" + fluid.getNumberOfPhases());
      System.out.println("beta=" + fluid.getBeta());
      System.out.println();
    }
  }

  static void trayRow(String label, double[][] inlets, String[] names, double[] z,
      double trayPressureBara, Double outletTemperatureK, double heatInputW) {
    neqsim.process.equipment.distillation.SimpleTray tray =
        new neqsim.process.equipment.distillation.SimpleTray("tray1");
    java.util.List<Stream> added = new java.util.ArrayList<Stream>();
    for (int i = 0; i < inlets.length; i++) {
      Stream inlet = feed(names, z, inlets[i][0], inlets[i][1], inlets[i][2]);
      added.add(inlet);
      tray.addStream(inlet);
    }
    if (trayPressureBara > 0.0) {
      tray.setPressure(trayPressureBara);
    }
    if (outletTemperatureK != null) {
      tray.setOutTemperature(outletTemperatureK, "K");
    }
    if (heatInputW != 0.0) {
      tray.setHeatInput(heatInputW);
    }
    tray.run();

    System.out.println(label);
    System.out.println("inlets=" + inlets.length);
    if (outletTemperatureK != null) {
      System.out.println("outlet_temperature_K=" + outletTemperatureK);
    }
    System.out.println("heat_input_W=" + heatInputW);
    System.out.println("tray_pressure_bara=" + tray.getPressure());
    System.out.println("tray_temperature_K=" + tray.getTemperature());
    System.out.println("mixed_phases=" + tray.getOutletStream().getThermoSystem().getNumberOfPhases());
    System.out.println("mixed_beta=" + tray.getOutletStream().getThermoSystem().getBeta());
    // The inlets' records, so the enthalpy the flash was given can be assembled from the
    // capture rather than inferred - `calcMixStreamEnthalpy` is the tray's own override.
    for (int i = 0; i < added.size(); i++) {
      print("feed" + i, added.get(i));
    }
    System.out.println("inlet_total_enthalpy_W=" + tray.calcMixStreamEnthalpy0());
    System.out.println("mixed_moles="
        + tray.getOutletStream().getThermoSystem().getTotalNumberOfMoles());
    print("gas_out", tray.getGasOutStream());
    print("liquid_out", tray.getLiquidOutStream());
    System.out.println();
  }

  /// **The family's closed-form member, and the one place the two libraries can be expected
  /// to agree tightly.** `ShortcutDistillationColumn` is Fenske-Underwood-Gilliland with a
  /// Kirkbride feed tray: it flashes the feed for K-values, forms `alpha_i = K_i / K_HK`, and
  /// every answer below is a rearrangement of those. **It is not a `DistillationColumn`** -
  /// it extends `ProcessEquipmentBaseClass` - so there is no tray profile and no MESH
  /// residual to print.
  ///
  /// **Its interior is eight scalars and two streams.** `getMinimumNumberOfStages`,
  /// `getMinimumRefluxRatio`, `getActualNumberOfStages`, `getActualRefluxRatio`,
  /// `getFeedTrayNumber`, `getCondenserDuty`, `getReboilerDuty` and `getRelativeVolatility`
  /// are public; the Underwood root, the feed quality and the per-component split fractions
  /// are private and are **not** in this capture. What is printed beside them is the flash
  /// the class itself performs, re-run here on a clone, so a divergence in `alpha` is
  /// localised before it reaches Fenske. `alpha_lk_hk` is printed twice - once from that
  /// flash, once from `getRelativeVolatility` - and the two agreeing is what makes it the
  /// same flash rather than two.
  ///
  /// **`solved` is printed because the class has a refusal.** A light key less volatile than
  /// the heavy key logs an error, sets `solved = false` and returns with every answer left at
  /// its field initialiser, and the swapped-key row is that path.
  ///
  /// **`setNumberOfTrays` is not settable here.** The class's own override logs
  /// `calculates stages; setNumberOfTrays ignored` and returns, so a row that called it would
  /// measure nothing.
  static void shortcutColumnRows() {
    // A C1/C2/C3/nC4 mixture split propane/n-butane, so every branch of the class's own
    // split-fraction estimate is reached: methane and ethane are lighter than the light key,
    // propane is it, and n-butane is the heavy key.
    String[] names = new String[] { "methane", "ethane", "propane", "n-butane" };
    double[] z = new double[] { 0.1, 0.3, 0.4, 0.2 };
    shortcutColumnRow("propane_nbutane_300K_20bara", names, z, 300.0, 20.0, 1.0, "propane",
        "n-butane", 0.98, 0.98, 1.2, null, null);
    // **The pressures stated**, which move both outlet streams' flash and nothing in the
    // FUG arithmetic - the class reads them only where it creates the products.
    shortcutColumnRow("pressures_stated_18_21_bara", names, z, 300.0, 20.0, 1.0, "propane",
        "n-butane", 0.98, 0.98, 1.2, 18.0, 21.0);
    // **A reflux multiplier of one**, so the actual reflux is the minimum and Gilliland's `X`
    // is zero - the branch where the correlation's own denominator has `sqrt(X)` in it.
    shortcutColumnRow("reflux_multiplier_one", names, z, 300.0, 20.0, 1.0, "propane",
        "n-butane", 0.98, 0.98, 1.0, null, null);
    // **A single-phase gas feed**, which is the Wilson branch: the class falls back to
    // `(Pc/P) exp(5.373 (1 + omega) (1 - Tc/T))` when its flash finds no second phase.
    shortcutColumnRow("wilson_fallback_superheated_gas", names, z, 450.0, 20.0, 1.0, "propane",
        "n-butane", 0.98, 0.98, 1.2, null, null);
    // **A binary**, which is the row a reader can check by hand.
    shortcutColumnRow("binary_methane_nbutane", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, 1.0, "methane", "n-butane", 0.99, 0.99, 1.2,
        null, null);
    // **The refusal**: the light key is the heavy component, so `alpha_LK/HK` is below one.
    shortcutColumnRow("light_key_heavier_than_heavy_key", names, z, 300.0, 20.0, 1.0,
        "n-butane", "propane", 0.98, 0.98, 1.2, null, null);
  }

  static void shortcutColumnRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, double molPerSec, String lightKey, String heavyKey,
      double lightKeyRecovery, double heavyKeyRecovery, double refluxMultiplier,
      Double condenserPressureBara, Double reboilerPressureBara) {
    Stream inlet = feed(names, z, temperatureK, pressureBara, molPerSec);

    // The flash the class performs, on a clone, so the K-values and the relative volatilities
    // are measurable rather than private.
    SystemInterface flashed = inlet.getThermoSystem().clone();
    new neqsim.thermodynamicoperations.ThermodynamicOperations(flashed).TPflash();
    flashed.init(2);
    int gas = shortcutPhaseIndex(flashed, true);
    int liquid = shortcutPhaseIndex(flashed, false);
    double[] kValues = new double[names.length];
    for (int i = 0; i < names.length; i++) {
      if (gas >= 0 && liquid >= 0 && gas != liquid) {
        double yi = flashed.getPhase(gas).getComponent(i).getx();
        double xi = flashed.getPhase(liquid).getComponent(i).getx();
        kValues[i] = xi > 1.0e-20 ? yi / xi : 1.0e10;
      } else {
        double tc = flashed.getPhase(0).getComponent(i).getTC();
        double pc = flashed.getPhase(0).getComponent(i).getPC();
        double omega = flashed.getPhase(0).getComponent(i).getAcentricFactor();
        kValues[i] = (pc / flashed.getPressure())
            * Math.exp(5.373 * (1.0 + omega) * (1.0 - tc / flashed.getTemperature()));
      }
    }
    int lk = flashed.getPhase(0).getComponent(lightKey).getComponentNumber();
    int hk = flashed.getPhase(0).getComponent(heavyKey).getComponentNumber();

    neqsim.process.equipment.distillation.ShortcutDistillationColumn column =
        new neqsim.process.equipment.distillation.ShortcutDistillationColumn("sc1", inlet);
    column.setLightKey(lightKey);
    column.setHeavyKey(heavyKey);
    column.setLightKeyRecoveryDistillate(lightKeyRecovery);
    column.setHeavyKeyRecoveryBottoms(heavyKeyRecovery);
    column.setRefluxRatioMultiplier(refluxMultiplier);
    if (condenserPressureBara != null) {
      column.setCondenserPressure(condenserPressureBara);
    }
    if (reboilerPressureBara != null) {
      column.setReboilerPressure(reboilerPressureBara);
    }
    column.run();

    System.out.println(label);
    print("feed", inlet);
    System.out.println("feed_phases=" + flashed.getNumberOfPhases());
    System.out.println("feed_beta=" + flashed.getBeta());
    System.out.println("alpha_lk_hk_from_flash=" + (kValues[lk] / kValues[hk]));
    System.out.println("alpha_lk_hk_reported=" + column.getRelativeVolatility());
    for (int i = 0; i < names.length; i++) {
      System.out.println("k_" + names[i] + "=" + kValues[i]);
      System.out.println("alpha_" + names[i] + "=" + (kValues[i] / kValues[hk]));
    }
    System.out.println("solved=" + column.isSolved());
    System.out.println("minimum_stages=" + column.getMinimumNumberOfStages());
    System.out.println("minimum_reflux_ratio=" + column.getMinimumRefluxRatio());
    System.out.println("actual_stages=" + column.getActualNumberOfStages());
    System.out.println("actual_reflux_ratio=" + column.getActualRefluxRatio());
    System.out.println("feed_tray_number=" + column.getFeedTrayNumber());
    System.out.println("condenser_duty_W=" + column.getCondenserDuty());
    System.out.println("reboiler_duty_W=" + column.getReboilerDuty());
    if (column.isSolved()) {
      print("distillate", column.getDistillateStream());
      print("bottoms", column.getBottomsStream());
    }
    System.out.println();
  }

  /// The class's own phase lookup, reproduced: `type` asks for the gas phase, otherwise the
  /// first liquid-like one.
  static int shortcutPhaseIndex(SystemInterface system, boolean wantGas) {
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      String type = system.getPhase(i).getPhaseTypeName();
      if (wantGas ? "gas".equals(type)
          : ("liquid".equals(type) || "oil".equals(type) || "aqueous".equals(type))) {
        return i;
      }
    }
    if (!wantGas) {
      for (int i = 0; i < system.getNumberOfPhases(); i++) {
        if (!"gas".equals(system.getPhase(i).getPhaseTypeName())) {
          return i;
        }
      }
    }
    return -1;
  }

  /// `ComponentSplitter`: a per-component routing with each outlet flashed.
  ///
  /// **The factor is per component and the outlet count is two.** `run` loops `for i in 0..2`
  /// and reads `splitFactor[k]` as the fraction of component *k* to the overhead, so the two
  /// outlets carry different *compositions* - which is what separates this from `Splitter`,
  /// whose outlets carry one state at different flows. `run` then sets each outlet's component
  /// moles on an empty fluid and runs a `TPflash`, so each outlet is a state of its own.
  ///
  /// **Three rows, and the third is the edge.** A near-total separation, an even routing
  /// whose outlets are the feed, and one that sends *all* of a component one way.
  static void componentSplitterRows() {
    String[] names = new String[] { "methane", "n-butane", "n-pentane" };
    double[] z = new double[] { 0.5, 0.3, 0.2 };
    componentRow("near_total_separation", names, z, new double[] { 0.98, 0.05, 0.02 });
    componentRow("an_even_routing", names, z, new double[] { 0.5, 0.5, 0.5 });
    componentRow("all_of_one_component", names, z, new double[] { 1.0, 0.5, 0.0 });
  }

  static void componentRow(String label, String[] names, double[] z, double[] factors) {
    Stream inlet = feed(names, z, 300.0, 20.0, 1.0);
    neqsim.process.equipment.splitter.ComponentSplitter splitter =
        new neqsim.process.equipment.splitter.ComponentSplitter("cs1", inlet);
    splitter.setSplitFactors(factors);
    splitter.run();

    StringBuilder given = new StringBuilder("split_factors=");
    for (int i = 0; i < factors.length; i++) {
      given.append(factors[i]).append(i + 1 < factors.length ? " " : "");
    }
    System.out.println(label);
    System.out.println(given);
    print("feed", inlet);
    print("overhead", splitter.getSplitStream(0));
    print("bottoms", splitter.getSplitStream(1));
    // The phase count each outlet settles on, because a material balance can move a stream
    // across a phase boundary and the enthalpy is not a function of the record alone.
    System.out.println("overhead_phases="
        + splitter.getSplitStream(0).getThermoSystem().getNumberOfPhases());
    System.out.println("bottoms_phases="
        + splitter.getSplitStream(1).getThermoSystem().getNumberOfPhases());
    System.out.println("feed_phases=" + inlet.getThermoSystem().getNumberOfPhases());
    System.out.println();
  }

  /// **A feed that splits three ways, so the aqueous outlet is a record rather than the
  /// empty 1e-20 kg/hr system `run` falls back to.** Methane and n-butane give a vapour and
  /// a hydrocarbon liquid and water gives the third, which is the only shape in which
  /// `hasPhaseType("aqueous")` is asked a question it can answer yes to.
  ///
  /// The rows walk `run`'s own sequence: the plain flash at the feed's temperature, the
  /// pressure drop applied before it, a heat input that moves the flash to an enthalpy, and
  /// two of the six entrainment transfers - one out of the vapour and one into it, because
  /// `addPhaseFractionToPhase` reads the phase it moves *from* and the order it runs in
  /// decides what the second sees.
  static void threePhaseSeparatorRows() {
    String[] names = new String[] { "methane", "n-butane", "water" };
    double[] z = new double[] { 0.5, 0.3, 0.2 };

    runThreePhase("pressure_drop_bara=0 heat_input_W=0 entrainment=none", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 0.0, null, null, null, 0.0);

    runThreePhase("pressure_drop_bara=2 heat_input_W=0 entrainment=none", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 2.0, null, null, null, 0.0);

    runThreePhase("pressure_drop_bara=0 heat_input_W=1000 entrainment=none", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 0.0, 1000.0, null, null, 0.0);

    runThreePhase("pressure_drop_bara=0 heat_input_W=0 entrainment=gas_to_oil_0.05", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 0.0, null, "gas", "oil", 0.05);

    runThreePhase("pressure_drop_bara=0 heat_input_W=0 entrainment=oil_to_gas_0.05", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 0.0, null, "oil", "gas", 0.05);

    // **A two-phase feed, for the other side of `run`'s phase test.** No aqueous phase
    // exists, so the water outlet is the empty system - and what that prints is what a port
    // has to reproduce to be reachable on the same inputs.
    String[] dry = new String[] { "methane", "n-butane" };
    double[] dryZ = new double[] { 0.7, 0.3 };
    runThreePhase("dry_feed pressure_drop_bara=0 heat_input_W=0 entrainment=none", dry, dryZ,
        feed(dry, dryZ, 300.0, 20.0, 1.0), 0.0, null, null, null, 0.0);
  }

  static void runThreePhase(String label, String[] names, double[] z, Stream inlet, double dropBara,
      Double heatInputW, String fromPhase, String toPhase, double entrainment) {
    neqsim.process.equipment.separator.ThreePhaseSeparator separator =
        new neqsim.process.equipment.separator.ThreePhaseSeparator("sep3", inlet);
    if (dropBara != 0.0) {
      separator.setPressureDrop(dropBara);
    }
    if (heatInputW != null) {
      separator.setHeatInput(heatInputW);
    }
    if (fromPhase != null) {
      // "mole" and "feed", the pair `Separator`'s own port uses: the fraction is of the
      // from-phase's mole count. The class's `setEntrainment` takes the same two strings
      // per direction, and a direction it does not name is silently ignored upstream.
      separator.setEntrainment(entrainment, "mole", "feed", fromPhase, toPhase);
    }
    separator.run();
    System.out.println(label);
    print("feed", inlet);
    printOrEmpty("vapour", separator.getGasOutStream());
    printOrEmpty("light_liquid", separator.getOilOutStream());
    printOrEmpty("heavy_liquid", separator.getWaterOutStream());
    System.out.println();
  }

  /// `run`'s absent-phase outlet is a 1e-20 kg/hr empty clone, whose enthalpy is 0/0 and
  /// whose molar volume the cubic cannot evaluate. Printing that as a `NaN` is the record;
  /// throwing on it is not, so the columns it cannot produce are named instead.
  static void printOrEmpty(String port, StreamInterface stream) {
    try {
      print(port, stream);
    } catch (Exception error) {
      System.out.println(port + "_unreadable=" + error.getClass().getSimpleName());
    }
  }

  /// **The design volume is the row that matters.** `setVolume` is read by the mechanical
  /// design, the capacity report and the JSON dump, and by nothing in `run` - which flashes
  /// at the *fluid's own* volume and internal energy, `VUflash(thermoSystem2.getVolume(),
  /// thermoSystem2.getInternalEnergy())`. Two rows that differ only in `setVolume` are the
  /// measurement that says whether a stated volume reaches the steady state at all.
  ///
  /// The last two rows ask what `run` does when a phase is absent, because its `else` branch
  /// names the wrong stream: `gasOutStream.setThermoSystemFromPhase(..., "oil")` is written
  /// where the liquid outlet was meant, so the liquid outlet is never refreshed.
  static void tankRows() {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.7, 0.3 };

    runTank("volume_default", names, z, feed(names, z, 300.0, 20.0, 1.0), null);
    runTank("volume_50", names, z, feed(names, z, 300.0, 20.0, 1.0), 50.0);

    String[] one = new String[] { "methane" };
    double[] oneZ = new double[] { 1.0 };
    runTank("single_phase_gas", one, oneZ, feed(one, oneZ, 300.0, 20.0, 1.0), null);

    String[] liquid = new String[] { "n-butane" };
    double[] liquidZ = new double[] { 1.0 };
    runTank("single_phase_liquid", liquid, liquidZ, feed(liquid, liquidZ, 250.0, 5.0, 1.0), null);

    // **A feed that flashes three ways, which is where the second branch bites.** The VU
    // flash finds a gas, an oil and an aqueous; `run` asks for "gas" and "oil" by name and
    // never mentions water, so the aqueous phase is not an outlet and what happens to its
    // material is this row's question.
    String[] wet = new String[] { "methane", "n-butane", "water" };
    double[] wetZ = new double[] { 0.5, 0.3, 0.2 };
    runTank("three_phase_feed", wet, wetZ, feed(wet, wetZ, 300.0, 20.0, 1.0), null);

    // **The same single-phase feed through the sibling class.** The tank's behaviour there
    // is only a defect if the class it shares its outlet construction with disagrees with
    // it, and `Separator.run`'s two phase tests name the stream the test asked about.
    singlePhaseThroughSeparator();
  }

  static void singlePhaseThroughSeparator() {
    String[] names = new String[] { "methane" };
    double[] z = new double[] { 1.0 };
    neqsim.process.equipment.separator.Separator separator =
        new neqsim.process.equipment.separator.Separator("sep1", feed(names, z, 300.0, 20.0, 1.0));
    separator.run();
    System.out.println("separator_on_the_same_single_phase_gas_feed");
    printOrEmpty("gas", separator.getGasOutStream());
    printOrEmpty("liquid", separator.getLiquidOutStream());
    System.out.println();
  }

  static void runTank(String label, String[] names, double[] z, Stream inlet, Double volumeM3) {
    neqsim.process.equipment.tank.Tank tank =
        new neqsim.process.equipment.tank.Tank("tank1", inlet);
    if (volumeM3 != null) {
      tank.setVolume(volumeM3);
    }
    tank.run();
    System.out.println(label);
    print("inlet", inlet);
    printOrEmpty("gas", tank.getGasOutStream());
    printOrEmpty("liquid", tank.getLiquidOutStream());
    // **What `run`'s two phase tests see.** `run` flashes a clone of the inlet at the
    // clone's *own* volume and internal energy and then asks `hasPhaseType("gas")` and
    // `hasPhaseType("oil")`; reproducing that flash here is what says which branch each
    // outlet took, and the type is not deducible from the state.
    SystemInterface vu = inlet.getThermoSystem().clone();
    new neqsim.thermodynamicoperations.ThermodynamicOperations(vu)
        .VUflash(vu.getVolume(), vu.getInternalEnergy());
    System.out.println("vu_phases=" + vu.getNumberOfPhases());
    for (int i = 0; i < vu.getNumberOfPhases(); i++) {
      System.out.println("vu_phase[" + i + "].type=" + vu.getPhase(i).getType() + " beta="
          + vu.getPhase(i).getBeta());
    }
    System.out.println();
  }

  static void splitter() {
    // **A single-phase liquid, so the branch flash is a formality.** `Splitter.run` clones
    // the inlet fluid, subtracts `(1 - f) n` from every component and runs a `TPflash` on
    // each branch. The composition is untouched by that subtraction, so in a single phase
    // the branch's enthalpy is the inlet's as a state function - and the two libraries can
    // be compared without a phase split sitting between them. A two-phase feed would put
    // the flash's own arithmetic into the answer instead of the splitter's.
    // **Two rows: the factors as written, and the same split written unnormalised.**
    // `setSplitFactors` divides by the sum, so the second row is what makes the
    // normalisation an oracled claim rather than one this library asserts about itself.
    double[][] rows = new double[][] { { 0.3, 0.7 }, { 3.0, 7.0 } };
    for (double[] factors : rows) {
      Stream inlet = feed(new String[] { "n-butane", "n-pentane" }, new double[] { 0.6, 0.4 },
          300.0, 10.0, 1.0);

      neqsim.process.equipment.splitter.Splitter splitter =
          new neqsim.process.equipment.splitter.Splitter("sp1", inlet);
      splitter.setSplitFactors(factors);
      splitter.run();

      StringBuilder given = new StringBuilder("split_factors=");
      for (int i = 0; i < factors.length; i++) {
        given.append(factors[i]).append(i + 1 < factors.length ? " " : "");
      }
      System.out.println(given);
      print("feed", inlet);
      double[] applied = splitter.getSplitFactors();
      double weighted = 0.0;
      for (int i = 0; i < applied.length; i++) {
        StreamInterface branch = splitter.getSplitStream(i);
        print("products" + i, branch);
        SystemInterface fluid = branch.getThermoSystem();
        weighted += applied[i] * fluid.getEnthalpy() / fluid.getTotalNumberOfMoles();
      }
      StringBuilder normalised = new StringBuilder("normalised=");
      for (int i = 0; i < applied.length; i++) {
        normalised.append(applied[i]).append(i + 1 < applied.length ? " " : "");
      }
      System.out.println(normalised);
      // **The balance the split owes, printed so the capture carries its own evidence.**
      // A splitter has no duty and no shaft work, so the fraction-weighted outlet enthalpy
      // must equal the inlet's. It is printed rather than left to a reader to derive,
      // because this is the row where NeqSim and this library part.
      SystemInterface feed = inlet.getThermoSystem();
      System.out.println("molar_enthalpy_in=" + feed.getEnthalpy() / feed.getTotalNumberOfMoles());
      System.out.println("molar_enthalpy_out_weighted=" + weighted);
      System.out.println();
    }
  }

  static void mixer() {
    // **Two inlets at different pressures and different compositions.** `Mixer.run` sets
    // the outlet to the lowest active inlet pressure and flashes the joined fluid to the
    // flow-weighted enthalpy, so both rules are exercised: a composition that is not
    // either feed's, and a pressure that is neither inlet's when the outlet is stated.
    // The palette declares the outlet pressure optional, so the second row states one.
    for (Double specified : new Double[] { null, 8.0 }) {
      Stream first = feed(new String[] { "n-butane", "n-pentane" }, new double[] { 0.6, 0.4 },
          300.0, 10.0, 1.0);
      Stream second = feed(new String[] { "n-butane", "n-pentane" }, new double[] { 0.3, 0.7 },
          300.0, 6.0, 2.0);

      neqsim.process.equipment.mixer.Mixer mixer =
          new neqsim.process.equipment.mixer.Mixer("mx1");
      mixer.addStream(first);
      mixer.addStream(second);
      if (specified != null) {
        mixer.setOutletPressure(specified);
      }
      mixer.run();

      System.out.println("specified_outlet_pressure_bara=" + specified);
      print("feed0", first);
      print("feed1", second);
      print("product", mixer.getOutletStream());
      // **The mixer's own number, and the one the record carries only in divided form.**
      // `calcMixStreamEnthalpy` sums the inlets' *totals* - J, for a process fluid whose
      // moles are one second of flow - and the outlet record carries that over the total
      // flow. A wrong weighting would be invisible in the quotient, which is the reason to
      // print the sum.
      //
      // **`getMinInletPressure` is named for the lowest inlet and does not return it.**
      // Its javadoc says "lowest active inlet pressure ... (equal to the outlet pressure)",
      // and the measurement agrees with the parenthesis: on the row that *specifies* 8 bara
      // it reads 8.0, not the 6.0 the lowest inlet is at. So it is a restatement of the
      // outlet pressure and is printed under that name rather than the getter's - it is
      // evidence about the class, not a layer of the mix.
      System.out.println("mix_pressure_bara=" + mixer.getMinInletPressure());
      System.out.println("mixed_enthalpy_W=" + mixer.calcMixStreamEnthalpy());
      System.out.println();
    }
  }

  static void separator() {
    // **Four rows, one per thing `Separator.run` does between the inlet and the outlets.**
    // A pressure drop, an optional heat input before the flash, and an entrainment
    // fraction that carries part of the vapour into the liquid; the first row has none of
    // the three, which is the plain flash at the *feed's* temperature - the state a port
    // that took a temperature of its own would never reach.
    //
    // The fluid is a methane/n-butane mixture at 300 K and 20 bara, chosen so a vapour
    // phase exists and the entrainment row has something to move. Both outlets carry a
    // record rather than one being an empty placeholder.
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.7, 0.3 };

    Stream plain = feed(names, z, 300.0, 20.0, 1.0);
    runSeparator("pressure_drop_bara=0 heat_input_W=0 gas_in_liquid=0", names, z, plain, 0.0, null,
        0.0);

    runSeparator("pressure_drop_bara=2 heat_input_W=0 gas_in_liquid=0", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 2.0, null, 0.0);

    runSeparator("pressure_drop_bara=0 heat_input_W=1000 gas_in_liquid=0", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 0.0, 1000.0, 0.0);

    runSeparator("pressure_drop_bara=0 heat_input_W=0 gas_in_liquid=0.05", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 0.0, null, 0.05);
  }

  static void runSeparator(String label, String[] names, double[] z, Stream inlet, double dropBara,
      Double heatInputW, double gasInLiquid) {
    neqsim.process.equipment.separator.Separator separator =
        new neqsim.process.equipment.separator.Separator("sep1", inlet);
    if (dropBara != 0.0) {
      separator.setPressureDrop(dropBara);
    }
    if (heatInputW != null) {
      separator.setHeatInput(heatInputW);
    }
    if (gasInLiquid != 0.0) {
      // `specType` is "mole": the fraction is of the phase's mole count, which is the
      // basis a five-field port can carry. "volume" and "mass" are the other two and
      // neither is a mole fraction.
      separator.setEntrainment(gasInLiquid, "mole", "feed", "gas", "liquid");
    }
    separator.run();
    System.out.println(label);
    print("feed", inlet);
    print("vapour", separator.getGasOutStream());
    print("liquid", separator.getLiquidOutStream());
    System.out.println();
  }

  static void throttlingValve() {
    // **Three rows, and the third is the one a port gets wrong - in the other direction
    // than it first looks.** A pressure *above* the inlet is not clamped: `acceptNegativeDP`
    // is `true` by default, so the outlet is set to the stated 40 bara and flashed at the
    // inlet's enthalpy, which raises its temperature 5 K. Describing this as a clamp is the
    // natural reading of `run`'s `isAcceptNegativeDP` branch and it is wrong - the branch
    // that clamps is the one the flag is *cleared* for. **A port that refused an outlet
    // above the inlet would be inventing a bound the class does not have**, which is what
    // the case records.
    //
    // The fluid is a single-phase gas, so the isenthalpic drop's temperature fall is a
    // real part of the answer rather than a rounding. `run` chooses the flash by
    // specification - `TPflash` when the pressure did not move or `isIsoThermal`, a
    // `PHflash` otherwise - so row three exercises the first branch and rows one and two
    // the second.
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.9, 0.1 };

    valveRow("outlet_pressure_bara=10", names, z, 10.0, null);
    valveRow("delta_pressure_bara=10", names, z, null, 10.0);
    valveRow("outlet_pressure_bara=40_above_the_inlet", names, z, 40.0, null);
  }

  static void valveRow(String label, String[] names, double[] z, Double outletBara,
      Double deltaBara) {
    Stream inlet = feed(names, z, 320.0, 30.0, 1.0);
    neqsim.process.equipment.valve.ThrottlingValve valve =
        new neqsim.process.equipment.valve.ThrottlingValve("v1", inlet);
    if (outletBara != null) {
      valve.setOutletPressure(outletBara);
    }
    if (deltaBara != null) {
      valve.setDeltaPressure(deltaBara, "bara");
    }
    valve.run();

    System.out.println(label);
    print("inlet", inlet);
    print("outlet", valve.getOutletStream());
    System.out.println();
  }

  static void heatExchanger() {
    // **`HeatExchanger.run`'s default branch is an effectiveness-NTU rating, not a
    // duty.** It seeds each side's outlet by flashing it at the *other* side's inlet
    // temperature to estimate a heat capacity, takes `Cmin`/`Cmax`, forms
    // `NTU = UA / Cmin` and applies `calcThermalEffectivenes(NTU, Cr)` to that swing.
    // `energyInput` is never read: `duty` is this class's *output*.
    //
    // **The first two rows exist because of a dimensional question**, and the `UA` is
    // deliberately small enough that the effectiveness does not saturate at one.
    // `Cmin` is a *molar* heat capacity (J/mol/K) and `UA` is a total conductance (W/K),
    // so `NTU = UA / Cmin` is dimensionless only while the capacity-limiting side flows
    // at one mole per second. Row two doubles *that* side's flow and changes nothing
    // else; a correct `NTU` would fall, and a molar one will not.
    String[] hotNames = new String[] { "methane", "n-butane" };
    double[] hotZ = new double[] { 0.8, 0.2 };
    String[] coldNames = new String[] { "n-butane", "n-pentane" };
    double[] coldZ = new double[] { 0.5, 0.5 };

    exchangerRow("ua_rating_counterflow", hotNames, hotZ, 1.0, coldNames, coldZ, 1.0, 100.0, null,
        0, null);
    exchangerRow("ua_rating_hot_flow_2", hotNames, hotZ, 2.0, coldNames, coldZ, 1.0, 100.0, null, 0,
        null);
    // The arrangement is a parameter, so it gets a row: the same exchanger in parallel flow
    // moves less heat than the counterflow one above, which is the whole reason it is one.
    exchangerRow("ua_rating_parallelflow", hotNames, hotZ, 1.0, coldNames, coldZ, 1.0, 100.0, null,
        0, "concentric tube paralellflow");
    exchangerRow("out_temperature_pins_the_hot_side", hotNames, hotZ, 1.0, coldNames, coldZ, 1.0,
        100.0, 350.0, 0, null);
    // **The pin is a mild one on purpose.** 360 K asked the hot side for more than it had:
    // it came back at 29.7 K, which is far below butane's freezing point, and neither
    // library's cubic is meaningful down there. A pin that leaves both sides in a state the
    // equation of state is for is worth more than one that tests the flash's extrapolation.
    exchangerRow("out_temperature_pins_the_cold_side", hotNames, hotZ, 1.0, coldNames, coldZ, 1.0,
        100.0, 320.0, 1, null);

    // **The pinned rows are not oracle-able, and these two are why.** `runSpecifiedStream`
    // reaches the pinned state by cloning the stream's *already flashed* fluid, setting the
    // temperature and flashing again, and that lands on an enthalpy that belongs to another
    // temperature: a PHflash back from it returns 317.73 K, not the 320 K the stream
    // reports. So these rows give the same two states reached the way that is self
    // consistent - a fresh fluid at the pinned temperature, and a PHflash for the side the
    // balance moves - and the case is pinned to these rather than to the rows above.
    pinReferenceRow("reference_pin_hot_350", hotNames, hotZ, 1.0, coldNames, coldZ, 1.0, 350.0, 0);
    pinReferenceRow("reference_pin_cold_320", hotNames, hotZ, 1.0, coldNames, coldZ, 1.0, 320.0, 1);
  }

  static void pinReferenceRow(String label, String[] hotNames, double[] hotZ, double hotFlow,
      String[] coldNames, double[] coldZ, double coldFlow, double pinnedK, int pinnedSide) {
    Stream hotIn = feed(hotNames, hotZ, 400.0, 20.0, hotFlow);
    Stream coldIn = feed(coldNames, coldZ, 300.0, 5.0, coldFlow);
    boolean hotIsPinned = pinnedSide == 0;

    // The pinned side: a fresh fluid at the pinned temperature, its own inlet pressure.
    Stream pinned = hotIsPinned
        ? feed(hotNames, hotZ, pinnedK, 20.0, hotFlow)
        : feed(coldNames, coldZ, pinnedK, 5.0, coldFlow);
    Stream pinnedIn = hotIsPinned ? hotIn : coldIn;
    double pinnedDuty = duty(pinnedIn, pinned);

    // The other side: a fresh fluid at its own inlet state, moved to the enthalpy the
    // balance leaves it. **Printed from the fluid**, because wrapping it back in a `Stream`
    // re-runs the stream's own machinery over a state that is already the answer.
    // `PHflash` takes the **absolute** molar enthalpy, not a shift, so the target is the
    // other side's own inlet value less the duty it is being handed.
    Stream otherIn = hotIsPinned ? coldIn : hotIn;
    double otherFlow = hotIsPinned ? coldFlow : hotFlow;
    SystemInterface otherFluid = otherIn.getThermoSystem();
    double otherInPerMole = otherFluid.getEnthalpy() / otherFluid.getTotalNumberOfMoles();
    SystemInterface other = hotIsPinned
        ? phFlash(coldNames, coldZ, 300.0, 5.0, coldFlow, otherInPerMole - pinnedDuty / coldFlow)
        : phFlash(hotNames, hotZ, 400.0, 20.0, hotFlow, otherInPerMole - pinnedDuty / hotFlow);

    System.out.println(label);
    print("hot_in", hotIn);
    print("cold_in", coldIn);
    if (hotIsPinned) {
      print("hot_out", pinned);
      printFluid("cold_out", other, otherFlow);
    } else {
      printFluid("hot_out", other, otherFlow);
      print("cold_out", pinned);
    }
    System.out.println();
  }

  /// One fluid's record, for a state that is not carried by a `Stream`.
  static void printFluid(String port, SystemInterface fluid, double molPerSec) {
    System.out.println(port + "_n=" + molPerSec);
    System.out.println(port + "_P=" + fluid.getPressure("bara"));
    System.out.println(port + "_T=" + fluid.getTemperature("K"));
    System.out.println(port + "_h=" + fluid.getEnthalpy() / fluid.getTotalNumberOfMoles());
    double[] overall = fluid.getMolarComposition();
    StringBuilder composition = new StringBuilder(port + "_z=");
    for (int i = 0; i < fluid.getPhase(0).getNumberOfComponents(); i++) {
      composition.append(fluid.getPhase(0).getComponent(i).getName()).append(":")
          .append(overall[i]);
      if (i + 1 < fluid.getPhase(0).getNumberOfComponents()) {
        composition.append(" ");
      }
    }
    System.out.println(composition);
  }

  /// The state a fresh fluid at `(names, z, T, P)` reaches when its enthalpy is moved by
  /// `perMole` joules per mole.
  ///
  /// **A fresh fluid, not a clone of the inlet's.** A cloned stream fluid has already been
  /// through `init(0)` and a flash, and a second flash on it can land on a stale root - the
  /// same reason `runSpecifiedStream`'s own pin does not reproduce a fresh flash's enthalpy.
  static SystemInterface phFlash(String[] names, double[] z, double temperatureK,
      double pressureBara, double molPerSecond, double perMole) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    fluid.setTotalFlowRate(molPerSecond, "mol/sec");
    fluid.init(0);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(fluid).PHflash(perMole, "J/mol");
    return fluid;
  }

  static void exchangerRow(String label, String[] hotNames, double[] hotZ, double hotFlow,
      String[] coldNames, double[] coldZ, double coldFlow, Double ua, Double pinnedK,
      int pinnedSide, String arrangement) {
    Stream hot = feed(hotNames, hotZ, 400.0, 20.0, hotFlow);
    Stream cold = feed(coldNames, coldZ, 300.0, 5.0, coldFlow);

    neqsim.process.equipment.heatexchanger.HeatExchanger hx =
        new neqsim.process.equipment.heatexchanger.HeatExchanger("hx1", hot, cold);
    if (ua != null) {
      hx.setUAvalue(ua);
    }
    if (pinnedK != null) {
      hx.setOutStreamSpecificationNumber(pinnedSide);
      hx.setOutTemperature(pinnedK, "K");
    }
    if (arrangement != null) {
      hx.setFlowArrangement(arrangement);
    }
    hx.run();

    System.out.println(label);
    print("hot_in", hot);
    print("cold_in", cold);
    print("hot_out", hx.getOutStream(0));
    print("cold_out", hx.getOutStream(1));
    // **The balance the exchanger owes, printed so the capture carries its own
    // evidence.** Total enthalpies, because the two sides may carry different flows.
    double hotDuty = duty(hot, hx.getOutStream(0));
    double coldDuty = duty(cold, hx.getOutStream(1));
    System.out.println("hot_duty_W=" + hotDuty);
    System.out.println("cold_duty_W=" + coldDuty);
    // **The rating's own two numbers, and the two the rating does not compute.** `run`
    // keeps `duty` and `thermalEffectiveness` as fields with getters, and those two pin the
    // whole interior: `duty = effectiveness * swing` and the kept side is the one whose
    // swing is larger, so `C_max = duty / (effectiveness * span)` and the effectiveness
    // relation then leaves `C_min` as its only unknown. `NTU` is a package-private field
    // with no getter, so it is reached through them rather than read.
    //
    // The pinned branch sets neither: `runSpecifiedStream` returns without touching `duty`
    // or `thermalEffectiveness`, so both read back as their initialisers. Printing a zero
    // there would say "this branch's duty is zero" when it says only that the branch does
    // not compute one.
    if (pinnedK == null) {
      System.out.println("duty_W=" + hx.getDuty());
      System.out.println("effectiveness=" + hx.getThermalEffectiveness());
    }
    System.out.println();
  }

  /// The total enthalpy a side gained, in W: molar change times molar flow.
  static double duty(StreamInterface inlet, StreamInterface outlet) {
    SystemInterface a = inlet.getThermoSystem();
    SystemInterface b = outlet.getThermoSystem();
    double in = a.getEnthalpy() / a.getTotalNumberOfMoles() * inlet.getFlowRate("mol/sec");
    double out = b.getEnthalpy() / b.getTotalNumberOfMoles() * outlet.getFlowRate("mol/sec");
    return out - in;
  }

  static void stream() {
    // **The stream's own properties, which no port field carries.** `s` is derived from
    // `(T, P, z)` exactly as `h` is, and density and viscosity come from the state rather
    // than from the record - so all three are functions a kernel reaches for rather than
    // quantities a connection passes. This is where they are measured.
    //
    // Four single-phase states and one two-phase one. The two-phase row is the point of
    // the last block: NeqSim's `getDensity()` averages the phases without saying so, and a
    // stream that is two phases has no one density.
    streamRow("butane_liquid", new String[] { "n-butane" }, new double[] { 1.0 }, 300.0, 10.0, 1.0);
    streamRow("methane_butane_gas", new String[] { "methane", "n-butane" },
        new double[] { 0.9, 0.1 }, 320.0, 30.0, 1.0);
    streamRow("water_liquid", new String[] { "water" }, new double[] { 1.0 }, 300.0, 1.0, 1.0);
    streamRow("methane_co2_gas", new String[] { "methane", "CO2" },
        new double[] { 0.7, 0.3 }, 300.0, 50.0, 1.0);
    streamRow("methane_butane_two_phase", new String[] { "methane", "n-butane" },
        new double[] { 0.7, 0.3 }, 300.0, 20.0, 1.0);
  }

  static void streamRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, double molPerSecond) {
    Stream inlet = feed(names, z, temperatureK, pressureBara, molPerSecond);
    SystemInterface fluid = inlet.getThermoSystem();
    double moles = fluid.getTotalNumberOfMoles();
    System.out.println(label);
    System.out.println("phases=" + fluid.getNumberOfPhases());
    System.out.println("molar_mass=" + fluid.getMolarMass());
    System.out.println("mass_flow=" + inlet.getFlowRate("kg/sec"));
    System.out.println("molar_entropy=" + fluid.getEntropy() / moles);
    // **Three densities, and they are three different numbers.** `getPhase(0).getDensity()`
    // with no unit is the cubic's own `M / (Z R T / P)`; the same call *with* a unit adds
    // NeqSim's volume correction, which is 16% on liquid water; and
    // `SystemInterface.getDensity()` averages over the phases. A port reports one of the
    // three, so the capture has to say which - and azoth's `Stream::density()` is the first.
    System.out.println("cubic_density=" + fluid.getPhase(0).getDensity());
    System.out.println("corrected_density=" + fluid.getPhase(0).getDensity("kg/m3"));
    System.out.println("system_density=" + fluid.getDensity("kg/m3"));
    System.out.println("viscosity=" + fluid.getViscosity("kg/msec"));
    System.out.println();
  }

  /// **A gas ejector: a motive stream entrains a suction stream and discharges above the
  /// suction pressure.**
  ///
  /// Four rows. The first is the class at its own defaults - every efficiency the field
  /// initialiser sets, and the mixing pressure the class *estimates* rather than being
  /// given - and the rest move one thing each: a lower motive-nozzle efficiency, a higher
  /// discharge pressure, and a heavier suction fluid. The class's own reported quantities
  /// (the entrainment ratio, the compression ratio, the area ratio, the critical back
  /// pressure and the three Mach numbers) are printed beside the record, because the
  /// machine's whole subject is what those say about the state.
  static void ejectorRows() {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.9, 0.1 };
    runEjector("defaults", names, z, 400.0, 30.0, 300.0, 5.0, 10.0, null, null, null, null);
    runEjector("motive_nozzle_efficiency_0.5", names, z, 400.0, 30.0, 300.0, 5.0, 10.0, 0.5, null, null, null);
    runEjector("discharge_14_bara", names, z, 400.0, 30.0, 300.0, 5.0, 14.0, null, null, null, null);
    runEjector("heavier_suction", new String[] { "methane", "n-butane" }, new double[] { 0.5, 0.5 }, 400.0, 30.0,
        300.0, 5.0, 10.0, null, null, null, null);
  }

  static void runEjector(String label, String[] names, double[] z, double motiveTemperatureK,
      double motivePressureBara, double suctionTemperatureK, double suctionPressureBara, double dischargePressureBara,
      Double motiveEfficiency, Double diffuserEfficiency, Double suctionEfficiency, Double mixingEfficiency) {
    Stream motive = feed(names, z, motiveTemperatureK, motivePressureBara, 1.0);
    Stream suction = feed(names, z, suctionTemperatureK, suctionPressureBara, 0.5);
    neqsim.process.equipment.ejector.Ejector ejector =
        new neqsim.process.equipment.ejector.Ejector("ej", motive, suction);
    ejector.setDischargePressure(dischargePressureBara);
    if (motiveEfficiency != null) {
      ejector.setEfficiencyIsentropic(motiveEfficiency);
    }
    if (diffuserEfficiency != null) {
      ejector.setDiffuserEfficiency(diffuserEfficiency);
    }
    if (suctionEfficiency != null) {
      ejector.setSuctionNozzleEfficiency(suctionEfficiency);
    }
    if (mixingEfficiency != null) {
      ejector.setMixingEfficiency(mixingEfficiency);
    }
    ejector.run();
    System.out.println(label);
    print("motive", motive);
    print("suction", suction);
    printOrEmpty("outlet", ejector.getOutStream());
    System.out.println("entrainment_ratio=" + ejector.getEntrainmentRatio());
    System.out.println("compression_ratio=" + ejector.getCompressionRatio());
    System.out.println("expansion_ratio=" + ejector.getExpansionRatio());
    System.out.println("area_ratio=" + ejector.getAreaRatio());
    System.out.println("critical_back_pressure=" + ejector.getCriticalBackPressure());
    System.out.println("motive_nozzle_mach=" + ejector.getMotiveNozzleMach());
    System.out.println("suction_mach=" + ejector.getSuctionMach());
    System.out.println("mixing_mach=" + ejector.getMixingMach());
    System.out.println("motive_choked=" + ejector.isMotiveChoked());
    System.out.println();
  }

  /// One stream's record, which is what a port carries.
  ///
  /// **Both the enthalpy and the composition are the system's, not phase 0's.** A port's
  /// record describes the *stream*, and a stream that has flashed is two phases with one
  /// overall composition; phase 0's `x` is the vapour's composition and phase 0's enthalpy
  /// over phase 0's moles is the vapour's, neither of which is what the port carries.
  ///
  /// The two agree on a single-phase fluid, which is why the pump's and splitter's
  /// captures cannot tell them apart - but a separator's whole subject is a fluid that is
  /// *not* single-phase, and reading phase 0 there reports the vapour twice.
  static void print(String port, StreamInterface stream) {
    SystemInterface fluid = stream.getThermoSystem();
    System.out.println(port + "_n=" + stream.getFlowRate("mol/sec"));
    System.out.println(port + "_P=" + stream.getPressure("bara"));
    System.out.println(port + "_T=" + stream.getTemperature("K"));
    System.out.println(port + "_h=" + fluid.getEnthalpy() / fluid.getTotalNumberOfMoles());
    // **The entropy at every port, which the record does not carry.** `s` is a function of
    // `(T, P, z)` exactly as `h` is, so a port that reports `h` can report `s` - and a
    // library that derives rather than carries it can be held to the oracle at every port
    // of every unit operation for one line here.
    System.out.println(port + "_s=" + fluid.getEntropy() / fluid.getTotalNumberOfMoles());
    double[] overall = fluid.getMolarComposition();
    StringBuilder composition = new StringBuilder(port + "_z=");
    for (int i = 0; i < fluid.getPhase(0).getNumberOfComponents(); i++) {
      composition.append(fluid.getPhase(0).getComponent(i).getName()).append(":")
          .append(overall[i]);
      if (i + 1 < fluid.getPhase(0).getNumberOfComponents()) {
        composition.append(" ");
      }
    }
    System.out.println(composition);
  }
}
