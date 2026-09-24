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
