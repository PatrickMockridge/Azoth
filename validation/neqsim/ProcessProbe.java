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
      System.out.println();
    }
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
    // **Three rows, and the third is the one a port gets wrong.** `ThrottlingValve.run`
    // clamps a pressure *above* the inlet back to the inlet - `isAcceptNegativeDP` is
    // false by default - so a valve asked to raise the pressure passes the stream
    // through rather than compressing it. A port that took the stated pressure at face
    // value would quietly become a compressor.
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
