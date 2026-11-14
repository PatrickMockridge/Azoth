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
//     javac -proc:none -cp neqsim-118d9af.jar ProcessProbe.java
//     java -cp .:neqsim-118d9af.jar ProcessProbe pump > captures/process_pump.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe splitter > captures/process_splitter.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe splitter_capacity \
//       > captures/process_splitter_capacity.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe splitter_mechanical_design \
//       > captures/process_splitter_mechanical_design.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe mixer > captures/process_mixer.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe separator > captures/process_separator.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe separator_capacity \
//       > captures/process_separator_capacity.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe separator_capacity_assessment \
//       > captures/process_separator_capacity_assessment.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe separator_sizing \
//       > captures/process_separator_sizing.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe separator_mechanical \
//       > captures/process_separator_mechanical_design.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe gas_scrubber_mechanical \
//       > captures/process_gas_scrubber_mechanical_design.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe three_phase_settling \
//       > captures/process_three_phase_settling.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe separator_performance_enhanced \
//       > captures/process_separator_performance_enhanced.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe anti_surge_recycle \
//       > captures/process_anti_surge_recycle.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe throttling_valve \
//       > captures/process_throttling_valve.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe heat_exchanger \
//       > captures/process_heat_exchanger.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe heater > captures/process_heater.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe cooler > captures/process_cooler.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe filter > captures/process_filter.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe pipe > captures/process_pipe.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe manifold > captures/process_manifold.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe gas_scrubber \
//       > captures/process_gas_scrubber.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe component_splitter \
//       > captures/process_component_splitter.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe tray > captures/process_column_tray.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe side_draw \
//       > captures/process_side_draw.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe column_tear \
//       > captures/process_column_tear.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe rate_based \
//       > captures/process_rate_based_packed_column.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe rate_based_billet \
//       > captures/process_rate_based_billet.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe rate_based_solvers \
//       > captures/process_rate_based_solvers.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe column_efficiency \
//       > captures/process_column_efficiency.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe absorber_efficiency \
//       > captures/process_absorber_efficiency.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe capacity \
//       > captures/process_column_capacity.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe designer \
//       > captures/process_internals_designer.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe coupling \
//       > captures/process_hydraulic_coupling.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe mechanical \
//       > captures/process_column_mechanical_design.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe column > captures/process_column.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe condenser \
//       > captures/process_column_condenser.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe reboiler \
//       > captures/process_column_reboiler.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe shortcut_column \
//       > captures/process_shortcut_distillation_column.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe ejector \
//       > captures/process_ejector.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe stream \
//       > captures/process_stream_properties.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe gibbs_reactor \
//       > captures/process_gibbs_reactor.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe gibbs_steps \
//       > captures/process_gibbs_reactor_steps.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe flowsheet \
//       > captures/process_flowsheet.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe flowsheet_accelerated \
//       > captures/process_flowsheet_accelerated.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe gas_dispersion \
//       > captures/process_gas_dispersion.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe gaussian_plume \
//       > captures/process_gaussian_plume.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe probit > captures/process_probit.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe jet_fire \
//       > captures/process_jet_fire.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe pool_fire \
//       > captures/process_pool_fire.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe vce > captures/process_vce.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe bleve > captures/process_bleve.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe event_tree \
//       > captures/process_event_tree.tsv
//     java -cp .:neqsim-118d9af.jar ProcessProbe consequence_analysis \
//       > captures/process_consequence_analysis.tsv

import neqsim.process.equipment.stream.Stream;
import neqsim.process.equipment.util.AccelerationMethod;
import neqsim.process.equipment.stream.StreamInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.process.safety.BoundaryConditions;
import neqsim.process.safety.dispersion.GasDispersionAnalyzer;
import neqsim.process.safety.dispersion.GasDispersionResult;
import neqsim.process.safety.dispersion.GaussianPlume;
import neqsim.process.safety.dispersion.HeavyGasDispersion;
import neqsim.process.safety.dispersion.ProbitModel;
import neqsim.process.safety.fire.BLEVECalculator;
import neqsim.process.safety.fire.JetFireModel;
import neqsim.process.safety.fire.PoolFireModel;
import neqsim.process.safety.fire.VCEModel;
import neqsim.process.safety.qra.ConsequenceAnalysisEngine;
import neqsim.process.safety.release.HomogeneousEquilibriumReleaseModel;
import neqsim.process.safety.release.LeakModel;
import neqsim.process.safety.release.ReleaseFlowRequest;
import neqsim.process.safety.release.ReleaseFlowResult;
import neqsim.process.safety.release.ReleaseState;
import neqsim.process.safety.release.SourceTermResult;
import neqsim.process.safety.risk.eta.EventTreeAnalyzer;
import java.util.Map;

public class ProcessProbe {
  public static void main(String[] args) {
    String which = args.length > 0 ? args[0] : "pump";
    switch (which) {
      case "pump":
        pump();
        break;
      case "pump_mechanical":
        pumpMechanicalRows();
        break;
      case "splitter":
        splitter();
        break;
      case "splitter_capacity":
        splitterCapacity();
        break;
      case "splitter_mechanical_design":
        splitterMechanicalDesign();
        break;
      case "mixer":
        mixer();
        break;
      case "separator":
        separator();
        break;
      case "separator_capacity":
        separatorCapacityRows();
        break;
      case "separator_capacity_assessment":
        separatorCapacityAssessmentRows();
        break;
      case "separator_sizing":
        separatorSizingRows();
        break;
      case "separator_mechanical":
        separatorMechanicalRows();
        break;
      case "separator_autosize":
        separatorAutosizeRows();
        break;
      case "gas_scrubber_mechanical":
        gasScrubberMechanicalRows();
        break;
      case "three_phase_settling":
        threePhaseSettlingRows();
        break;
      case "separator_performance_enhanced":
        separatorPerformanceEnhancedRows();
        break;
      case "separator_performance":
        separatorPerformanceRows();
        break;
      case "throttling_valve":
        throttlingValve();
        break;
      case "heat_exchanger":
        heatExchanger();
        break;
      case "heat_exchanger_delta_t":
        heatExchangerDeltaT();
        break;
      case "heat_exchanger_rating":
        heatExchangerRating();
        break;
      case "heat_exchanger_mechanical":
        heatExchangerMechanicalRows();
        break;
      case "heat_exchanger_transient":
        heatExchangerTransientRows();
        break;
      case "cooler_transient":
        coolerTransientRows();
        break;
      case "heater":
        heatRows("Heater");
        break;
      case "heater_capacity":
        heaterCapacity();
        break;
      case "heater_instruments":
        heaterInstrumentRows();
        break;
      case "heater_electrical":
        heaterElectricalRows();
        break;
      case "heater_vfd":
        heaterVfdRows();
        break;
      case "heater_utility":
        heaterUtilityRows();
        break;
      case "heater_mechanical":
        heaterMechanicalRows();
        break;
      case "cooler":
        heatRows("Cooler");
        break;
      case "filter":
        filterRows();
        break;
      case "filter_models":
        filterModelRows();
        break;
      case "filter_loading":
        filterLoadingRows();
        break;
      case "filter_dynamic":
        filterDynamicRows();
        break;
      case "filter_mechanical_design":
        filterMechanicalDesignRows();
        break;
      case "compressor":
        isentropicRows("Compressor");
        break;
      case "compressor_polytropic":
        compressorPolytropicRows();
        break;
      case "compressor_mechanical":
        compressorMechanicalRows();
        break;
      case "compressor_modes":
        compressorModeRows();
        break;
      case "compressor_deposit":
        compressorDepositRows();
        break;
      case "compressor_curve":
        compressorCurveRows();
        break;
      case "compressor_speed_solve":
        compressorSpeedSolveRows();
        break;
      case "equilibrium_sound_speed":
        equilibriumSoundSpeedRows();
        break;
      case "hem_release":
        hemReleaseRows();
        break;
      case "leak_blowdown":
        leakBlowdownRows();
        break;
      case "gas_dispersion":
        gasDispersionRows();
        break;
      case "gaussian_plume":
        gaussianPlumeRows();
        break;
      case "probit":
        probitRows();
        break;
      case "jet_fire":
        jetFireRows();
        break;
      case "pool_fire":
        poolFireRows();
        break;
      case "vce":
        vceRows();
        break;
      case "bleve":
        bleveRows();
        break;
      case "event_tree":
        eventTreeRows();
        break;
      case "consequence_analysis":
        consequenceAnalysisRows();
        break;
      case "expander_steps":
        expanderStepRows();
        break;
      case "expander_chart":
        expanderChartRows();
        break;
      case "turbo_expander":
        turboExpanderRows();
        break;
      case "turbo_expander_mechanical":
        turboExpanderMechanicalRows();
        break;
      case "anti_surge":
        antiSurgeRows();
        break;
      case "anti_surge_recycle":
        antiSurgeRecycleRows();
        break;
      case "compressor_map_lookup":
        compressorMapLookupRows();
        break;
      case "sound_speed":
        soundSpeedRows();
        break;
      case "compressor_khader":
        compressorKhaderRows();
        break;
      case "map_turbo_expander":
        mapTurboExpanderRows();
        break;
      case "expander":
        isentropicRows("Expander");
        break;
      case "expander_mechanical":
        expanderMechanicalRows();
        break;
      case "expander_capacity":
        expanderCapacityRows();
        break;
      case "expander_power_modes":
        expanderPowerModeRows();
        break;
      case "pipe":
        pipeRows();
        break;
      case "pipe_modes":
        pipeModeRows();
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
      case "flare":
        flareRows();
        break;
      case "stirred_tank_reactor":
        stirredTankReactorRows();
        break;
      case "plug_flow_reactor":
        plugFlowReactorRows();
        break;
      case "gibbs_reactor":
        gibbsReactorRows();
        break;
      case "gibbs_steps":
        gibbsReactorSteps();
        break;
      case "gibbs_algorithms":
        gibbsAlgorithmRows();
        break;
      case "valve_sizing":
        valveSizingRows();
        break;
      case "pump_curve":
        pumpCurveRows();
        break;
      case "flowsheet":
        flowsheetRows();
        break;
      case "flowsheet_accelerated":
        flowsheetAcceleratedRows();
        break;
      case "column":
        columnRows();
        break;
      case "reactive_column":
        reactiveColumnRows();
        break;
      case "column_solvers":
        columnSolverRows(args);
        break;
      case "absorber":
        absorberRows();
        break;
      case "absorber_efficiency":
        absorberEfficiencyRows();
        break;
      case "capacity":
        capacityRows();
        break;
      case "designer":
        designerRows();
        break;
      case "coupling":
        couplingRows();
        break;
      case "mechanical":
        mechanicalRows();
        break;
      case "reactor_mechanical":
        reactorMechanicalRows();
        break;
      case "packed_column":
        packedColumnRows();
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
      case "side_draw":
        sideDrawRows();
        break;
      case "column_tear":
        columnTearRows();
        break;
      case "rate_based":
        rateBasedRows();
        break;
      case "rate_based_billet":
        rateBasedBilletRows();
        break;
      case "rate_based_solvers":
        rateBasedSolverRows();
        break;
      case "column_efficiency":
        murphreeStageRows();
        break;
      case "column_divergence":
        divergenceRows();
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
  /// A feed stated as a **mass** flow, which is how the absorber's own tests state theirs.
  static Stream feedMass(String name, String[] names, double[] z, double temperatureK,
      double pressureBara, double kgPerHour) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream stream = new Stream(name, fluid);
    stream.setFlowRate(kgPerHour, "kg/hr");
    stream.run();
    return stream;
  }

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

  /// `PumpMechanicalDesign.calcDesign`, with the API 610 assessment it runs.
  ///
  /// **The rows are the states the class distinguishes**: the plain isentropic branch at two
  /// efficiencies, the chart branch where the vendor BEP and NPSHr are read, a purchaser-fixed
  /// API 610 type with a stated suction pressure, casing MAWP and shutoff head, and a zero-flow
  /// row. The last is evidence only - the port's model refuses a zero `inlet_n`, so the class's
  /// guard is stated rather than cased.
  ///
  /// **The base class's fields are not one unit.** `getInnerDiameter`, `getOuterDiameter` and
  /// `getTantanLength` are metres here, `getWallThickness` is millimetres, the weights are
  /// kilograms and the module dimensions are metres - the row prints the getters as the class
  /// returns them, and the spec converts each at its own boundary.
  static void pumpMechanicalRows() {
    pumpMechanicalRow("plain_eta_0_75", new String[] {"n-butane"}, new double[] {1.0},
        250.0, 5.0, 1.0, 20.0, 0.75, 1000.0, false, false);
    pumpMechanicalRow("plain_eta_1_0", new String[] {"n-butane"}, new double[] {1.0},
        250.0, 5.0, 1.0, 20.0, 1.0, 1000.0, false, false);
    // The 1 -> 60 bar water state puts the head over the class's 250 m threshold, so
    // `selectPumpType` reaches BETWEEN_BEARINGS and the multi-stage head division.
    pumpMechanicalRow("multistage", new String[] {"water"}, new double[] {1.0},
        293.15, 1.0, 2775.9, 60.0, 0.75, 3000.0, false, false);
    pumpMechanicalRow("purchaser", new String[] {"n-butane"}, new double[] {1.0},
        250.0, 5.0, 1.0, 20.0, 0.75, 1000.0, true, false);
    pumpMechanicalRow("sealless", new String[] {"n-butane"}, new double[] {1.0},
        250.0, 5.0, 1.0, 20.0, 0.75, 1000.0, false, true);
    pumpMechanicalRow("zero_flow", new String[] {"n-butane"}, new double[] {1.0},
        250.0, 5.0, 0.0, 20.0, 0.75, 1000.0, false, false);
    pumpMechanicalChartRow("chart", false);
    pumpMechanicalNpshRow();
  }

  static void pumpMechanicalRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, double molPerSecond, double outletBara, double efficiency,
      double speed, boolean purchaser, boolean sealless) {
    Stream inlet = feed(names, z, temperatureK, pressureBara, molPerSecond);
    neqsim.process.equipment.pump.Pump pump =
        new neqsim.process.equipment.pump.Pump(label, inlet);
    pump.setOutletPressure(outletBara);
    pump.setIsentropicEfficiency(efficiency);
    pump.setSpeed(speed);
    neqsim.process.mechanicaldesign.pump.PumpMechanicalDesign design = pump.getMechanicalDesign();
    if (purchaser) {
      design.setApi610PumpType(
          neqsim.process.mechanicaldesign.pump.PumpApi610DesignCalculator.Api610PumpType.OH2);
      design.setMaximumSuctionPressure(6.0);
      design.setFurnishedCasingMawp(25.0);
      design.setShutoffHead(400.0);
    }
    if (sealless) {
      design.setSealType(neqsim.process.mechanicaldesign.pump.PumpMechanicalDesign.SealType.MAGNETIC_DRIVE);
    }
    pump.run();
    design.calcDesign();
    printPumpDesign(label, inlet, pump, design);
  }

  static void pumpMechanicalChartRow(String label, boolean npshCurve) {
    Stream inlet = feedMass("inlet", new String[] {"water"}, new double[] {1.0}, 293.15, 1.0, 4.0e6);
    neqsim.process.equipment.pump.Pump pump =
        new neqsim.process.equipment.pump.Pump(label, inlet);
    setChartTest((neqsim.process.equipment.pump.PumpChart) pump.getPumpChart(),
        new double[] {0.3, 1.0, 1.0, 1.0});
    pump.setSpeed(12913.0);
    pump.run();
    neqsim.process.mechanicaldesign.pump.PumpMechanicalDesign design = pump.getMechanicalDesign();
    design.calcDesign();
    printPumpDesign(label, inlet, pump, design);
  }

  /// The two-speed chart with a vendor NPSH curve, so `npshrDataSource` is `VENDOR_CURVE`
  /// rather than the screening estimate. The NPSH table is the one the pump's own NPSH row
  /// uses; its outer dimension is one row per speed, which is why it cannot ride the
  /// eight-speed performance chart.
  static void pumpMechanicalNpshRow() {
    Stream inlet = feedMass("inlet", new String[] {"water"}, new double[] {1.0}, 298.15, 2.0,
        39900.0);
    neqsim.process.equipment.pump.Pump pump =
        new neqsim.process.equipment.pump.Pump("chart_npsh", inlet);
    pump.getPumpChart().setCurves(new double[] {}, new double[] {1000.0, 1500.0},
        new double[][] {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0},
            {15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0}},
        new double[][] {{120.0, 118.0, 115.0, 110.0, 103.0, 94.0, 83.0, 70.0},
            {270.0, 265.5, 258.8, 247.5, 231.8, 211.5, 186.8, 157.5}},
        new double[][] {{60.0, 70.0, 78.0, 82.0, 81.0, 76.0, 68.0, 55.0},
            {62.0, 71.0, 79.0, 83.0, 82.0, 77.0, 69.0, 56.0}});
    pump.getPumpChart().setHeadUnit("meter");
    pump.getPumpChart().setNPSHCurve(new double[][] {
        {2.0, 2.2, 2.5, 3.0, 3.8, 4.8, 6.2, 8.0},
        {4.5, 4.95, 5.625, 6.75, 8.55, 10.8, 13.95, 18.0}});
    pump.setSpeed(1000.0);
    pump.setCheckNPSH(true);
    pump.setNPSHMargin(1.3);
    pump.run();
    neqsim.process.mechanicaldesign.pump.PumpMechanicalDesign design = pump.getMechanicalDesign();
    design.calcDesign();
    printPumpDesign("chart_npsh", inlet, pump, design);
  }

  static void printPumpDesign(String label, Stream inlet,
      neqsim.process.equipment.pump.Pump pump,
      neqsim.process.mechanicaldesign.pump.PumpMechanicalDesign design) {
    System.out.println(label);
    System.out.println("pump_type=" + design.getPumpType());
    System.out.println("seal_type=" + design.getSealType());
    System.out.println("api610_pump_type=" + design.getApi610PumpType());
    System.out.println("rated_flow_m3_per_hr=" + design.getRatedFlow());
    System.out.println("rated_head_m=" + design.getRatedHead());
    System.out.println("pump_efficiency=" + design.getPumpEfficiency());
    System.out.println("hydraulic_power_kW=" + design.getHydraulicPower());
    System.out.println("absorbed_power_kW=" + design.getAbsorbedPower());
    System.out.println("rated_speed_rpm=" + design.getRatedSpeed());
    System.out.println("specific_speed=" + design.getSpecificSpeed());
    System.out.println("suction_specific_speed=" + design.getSuctionSpecificSpeed());
    System.out.println("number_of_stages=" + design.getNumberOfStages());
    System.out.println("bep_flow_m3_per_hr=" + design.getBepFlow());
    System.out.println("bep_head_m=" + design.getBepHead());
    System.out.println("npsh_required_m=" + design.getNpshRequired());
    System.out.println("npsh_available_m=" + design.getNpshAvailable());
    System.out.println("design_npsh_margin_m=" + design.getNpshMargin());
    System.out.println("design_npsh_margin_ratio=" + design.getNpshMarginRatio());
    System.out.println("driver_power_kW=" + design.getDriverPower());
    System.out.println("driver_margin=" + design.getDriverMargin());
    System.out.println("design_pressure_bara=" + design.getDesignPressure());
    System.out.println("design_temperature_C=" + design.getDesignTemperature());
    System.out.println("impeller_diameter_mm=" + design.getImpellerDiameter());
    System.out.println("impeller_width_mm=" + design.getImpellerWidth());
    System.out.println("shaft_diameter_mm=" + design.getShaftDiameter());
    System.out.println("casing_wall_thickness_mm=" + design.getCasingWallThickness());
    System.out.println("casing_mawp_bara=" + design.getMaxAllowableWorkingPressure());
    System.out.println("suction_nozzle_size_in=" + design.getSuctionNozzleSize());
    System.out.println("discharge_nozzle_size_in=" + design.getDischargeNozzleSize());
    System.out.println("minimum_flow_m3_per_hr=" + design.getMinimumFlow());
    System.out.println("maximum_flow_m3_per_hr=" + design.getMaximumFlow());
    System.out.println("shutoff_head_m=" + privateDouble(design, "shutoffHead"));
    System.out.println("effective_isentropic_efficiency=" + pump.getIsentropicEfficiency());
    System.out.println("npsh_margin_factor=" + design.getNpshMarginFactor());
    System.out.println("minimum_npsh_margin_m=" + design.getMinimumNpshMarginM());
    System.out.println("por_low_fraction=" + design.getPorLowFraction());
    System.out.println("por_high_fraction=" + design.getPorHighFraction());
    System.out.println("aor_low_fraction=" + design.getAorLowFraction());
    System.out.println("aor_high_fraction=" + design.getAorHighFraction());
    System.out.println(
        "rated_point_low_fraction=" + privateDouble(design, "ratedPointLowFraction"));
    System.out.println(
        "rated_point_high_fraction=" + privateDouble(design, "ratedPointHighFraction"));
    System.out.println(
        "min_continuous_flow_fraction=" + privateDouble(design, "minContinuousFlowFraction"));
    System.out.println(
        "max_continuous_flow_fraction=" + privateDouble(design, "maxContinuousFlowFraction"));
    System.out.println("max_suction_specific_speed=" + design.getMaxSuctionSpecificSpeed());
    System.out.println("hydraulic_power_margin=" + design.getHydraulicPowerMargin());
    System.out.println("head_margin_factor=" + design.getHeadMarginFactor());
    System.out.println("casing_weight_kg=" + design.getCasingWeight());
    System.out.println("motor_weight_kg=" + design.getMotorWeight());
    System.out.println("baseplate_weight_kg=" + design.getBaseplateWeight());
    System.out.println("inner_diameter_m=" + design.getInnerDiameter());
    System.out.println("outer_diameter_m=" + design.getOuterDiameter());
    System.out.println("wall_thickness_mm=" + design.getWallThickness());
    System.out.println("tantan_length_m=" + design.getTantanLength());
    System.out.println("weight_piping_kg=" + design.getWeightPiping());
    System.out.println("weight_electro_instrument_kg=" + design.getWeightElectroInstrument());
    System.out.println("weight_structural_steel_kg=" + design.getWeightStructualSteel());
    System.out.println("weight_total_kg=" + design.getWeightTotal());
    System.out.println("module_length_m=" + design.getModuleLength());
    System.out.println("module_width_m=" + design.getModuleWidth());
    System.out.println("module_height_m=" + design.getModuleHeight());

    neqsim.process.mechanicaldesign.pump.PumpApi610DesignCalculator assessment =
        design.getApi610Assessment();
    System.out.println("assessment_status=" + assessment.getAssessmentStatus());
    System.out.println("selected_driver_power_kW=" + assessment.getSelectedDriverPowerKw());
    System.out.println("required_driver_power_kW=" + assessment.getRequiredDriverPowerKw());
    System.out.println("operating_region=" + assessment.getOperatingRegion());
    System.out.println("operating_flow_ratio=" + assessment.getOperatingFlowRatio());
    System.out.println("assessment_npsh_margin_ratio=" + assessment.getNpshMarginRatio());
    System.out.println("required_casing_pressure_bara=" + assessment.getRequiredCasingPressureBara());
    System.out.println("hydrostatic_test_pressure_bara=" + assessment.getHydrostaticTestPressureBara());
    System.out.println("head_rise_to_shutoff_fraction=" + assessment.getHeadRiseToShutoffFraction());
    System.out.println("bearing_l10_life_hours=" + assessment.getBearingL10LifeHours());
    System.out.println("critical_speed_ratio=" + assessment.getCriticalSpeedRatio());
    System.out.println("assessment_basis=" + assessment.getAssessmentBasis());
    for (neqsim.process.mechanicaldesign.pump.PumpApi610DesignCalculator.Check check :
        assessment.getChecks()) {
      System.out.println("check " + check.getId() + " status=" + check.getStatus()
          + " actual=" + check.getActualValue() + " limit=" + check.getLimitValue()
          + " unit=" + check.getUnit() + " requirement=" + check.getRequirement()
          + " message=" + check.getMessage());
    }

    System.out.println("inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("inlet_density_kg_per_m3=" + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("mass_flow_kg_per_hr=" + inlet.getFlowRate("kg/hr"));
    System.out.println("inlet_number_of_phases=" + inlet.getThermoSystem().getNumberOfPhases());
    print("inlet", inlet);
    print("outlet", pump.getOutletStream());
    System.out.println("power_kW=" + pump.getPower("kW"));
    System.out.println();
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

  /// **`Heater`'s bounded outlet temperatures**, the two HARD constraints its limits create.
  ///
  /// The maximum is a plain `value / design`; the minimum is the inverse, `min / value`, and
  /// the class's `isHardLimitExceeded` reads `maxValue` for the maximum constraint - which
  /// `setMaxOutletTemperature` never sets - so only the minimum can trip it. The mechanical
  /// design's duty and pressure-drop constraints are absent while `maxDesignDuty` and
  /// `maxDesignPressureDrop` are zero.
  static void heaterCapacity() {
    String[] names = new String[] {"methane", "n-butane"};
    double[] z = new double[] {0.9, 0.1};
    heaterCapacityRow("no_limits", names, z, 380.0, null, null);
    heaterCapacityRow("max_above", names, z, 380.0, 400.0, null);
    heaterCapacityRow("max_below", names, z, 380.0, 350.0, null);
    heaterCapacityRow("min_below", names, z, 380.0, null, 300.0);
    heaterCapacityRow("min_above", names, z, 380.0, null, 400.0);
    heaterCapacityRow("both", names, z, 380.0, 400.0, 300.0);
    heaterCapacityRow("min_zero", names, z, 380.0, null, 0.0);
    heaterCapacityRow("min_capped", names, z, 380.0, null, 1.0e6);
    heaterCapacityRow("duty_specified", names, z, null, 400.0, 300.0);
  }

  static void heaterCapacityRow(String label, String[] names, double[] z, Double outletK,
      Double maxLimitK, Double minLimitK) {
    Stream inlet = feed(names, z, 320.0, 30.0, 1.0);
    neqsim.process.equipment.heatexchanger.Heater unit =
        new neqsim.process.equipment.heatexchanger.Heater("h1", inlet);
    if (outletK != null) {
      unit.setOutletTemperature(outletK);
    } else {
      unit.setDuty(5000.0);
    }
    if (maxLimitK != null) {
      unit.setMaxOutletTemperature(maxLimitK);
    }
    if (minLimitK != null) {
      unit.setMinOutletTemperature(minLimitK);
    }
    unit.run();

    System.out.println("case=" + label);
    print("inlet", inlet);
    print("outlet", unit.getOutletStream());
    System.out.println("duty_W=" + unit.getDuty());
    System.out.println("outlet_temperature_k=" + unit.getOutletStream().getTemperature());
    System.out.println("max_outlet_temperature_limit_k=" + unit.getMaxOutletTemperature());
    System.out.println("min_outlet_temperature_limit_k=" + unit.getMinOutletTemperature());
    System.out.println("has_max=" + unit.hasMaxOutletTemperatureLimit());
    System.out.println("has_min=" + unit.hasMinOutletTemperatureLimit());

    Map<String, neqsim.process.equipment.capacity.CapacityConstraint> constraints =
        unit.getCapacityConstraints();
    System.out.println("constraint_count=" + constraints.size());
    for (Map.Entry<String, neqsim.process.equipment.capacity.CapacityConstraint> entry :
        constraints.entrySet()) {
      neqsim.process.equipment.capacity.CapacityConstraint constraint = entry.getValue();
      String name = entry.getKey();
      System.out.println("constraint_" + name + "_type=" + constraint.getType());
      System.out.println("constraint_" + name + "_current=" + constraint.getCurrentValue());
      System.out.println("constraint_" + name + "_design=" + constraint.getDesignValue());
      System.out.println("constraint_" + name + "_min=" + constraint.getMinValue());
      System.out.println("constraint_" + name + "_max=" + constraint.getMaxValue());
      System.out.println("constraint_" + name + "_utilization=" + constraint.getUtilization());
      System.out.println("constraint_" + name + "_violated=" + constraint.isViolated());
      System.out.println(
          "constraint_" + name + "_hard_exceeded=" + constraint.isHardLimitExceeded());
    }
    System.out.println("max_utilization=" + unit.getMaxUtilization());
    System.out.println("capacity_exceeded=" + unit.isCapacityExceeded());
    System.out.println("hard_limit_exceeded=" + unit.isHardLimitExceeded());
    neqsim.process.equipment.capacity.CapacityConstraint bottleneck = unit
        .getBottleneckConstraint();
    System.out.println(
        "bottleneck_constraint=" + (bottleneck == null ? "none" : bottleneck.getName()));
    System.out.println();
  }

  /// `Heater`'s instrument design: the static list, its I/O counts and its three types.
  ///
  /// **The list is not physics.** A heater's design is eight instruments by default - process
  /// inlet/outlet temperature and pressure, a differential pressure, two SIL trips and the
  /// element's own temperature - and the only branch is the type. The estimated costs are USD
  /// and are printed as evidence, not carried.
  static void heaterInstrumentRows() {
    String[] names = new String[] {"methane", "n-butane"};
    double[] z = new double[] {0.9, 0.1};
    for (boolean safety : new boolean[] {true, false}) {
      Stream inlet = feed(names, z, 320.0, 30.0, 1.0);
      neqsim.process.equipment.heatexchanger.Heater heater =
          new neqsim.process.equipment.heatexchanger.Heater("h1", inlet);
      heater.setOutletTemperature(340.0);
      heater.run();
      neqsim.process.instrumentdesign.heatexchanger.HeatExchangerInstrumentDesign design =
          heater.getInstrumentDesign();
      design.setIncludeSafetyInstruments(safety);
      design.calcDesign();
      System.out.println("include_safety=" + safety);
      printInstrumentDesign(design);
    }

    for (neqsim.process.instrumentdesign.heatexchanger.HeatExchangerInstrumentDesign.HeatExchangerType type :
        new neqsim.process.instrumentdesign.heatexchanger.HeatExchangerInstrumentDesign.HeatExchangerType[] {
            neqsim.process.instrumentdesign.heatexchanger.HeatExchangerInstrumentDesign.HeatExchangerType.SHELL_AND_TUBE,
            neqsim.process.instrumentdesign.heatexchanger.HeatExchangerInstrumentDesign.HeatExchangerType.AIR_COOLER,
            neqsim.process.instrumentdesign.heatexchanger.HeatExchangerInstrumentDesign.HeatExchangerType.ELECTRIC_HEATER}) {
      Stream inlet = feed(names, z, 320.0, 30.0, 1.0);
      neqsim.process.equipment.heatexchanger.Heater heater =
          new neqsim.process.equipment.heatexchanger.Heater("h1", inlet);
      heater.setOutletTemperature(340.0);
      heater.run();
      neqsim.process.instrumentdesign.heatexchanger.HeatExchangerInstrumentDesign design =
          heater.getInstrumentDesign();
      design.setHeatExchangerType(type);
      design.calcDesign();
      System.out.println("heat_exchanger_type=" + type.name());
      printInstrumentDesign(design);
    }
  }

  static void printInstrumentDesign(
      neqsim.process.instrumentdesign.InstrumentDesign design) {
    System.out.println("instrument_count=" + design.getInstrumentList().size());
    System.out.println("analog_input_count=" + design.getInstrumentList().getAnalogInputCount());
    System.out.println(
        "analog_output_count=" + design.getInstrumentList().getAnalogOutputCount());
    System.out.println(
        "digital_input_count=" + design.getInstrumentList().getDigitalInputCount());
    System.out.println(
        "digital_output_count=" + design.getInstrumentList().getDigitalOutputCount());
    System.out.println("total_io_count=" + design.getTotalIOCount());
    System.out.println(
        "safety_instrument_count=" + design.getInstrumentList().getSafetyInstrumentCount());
    System.out.println("estimated_cost_usd=" + design.getEstimatedCostUSD());
    System.out.println("hazardous_area_zone=" + design.getHazardousAreaZone());
    System.out.println("protection_concept=" + design.getProtectionConcept());
    System.out.println("instrument_standard=" + design.getInstrumentStandard());
    System.out.println("include_safety_instruments=" + design.isIncludeSafetyInstruments());
    System.out.println("default_sil_level=" + design.getDefaultSilLevel());
    int index = 0;
    for (neqsim.process.instrumentdesign.InstrumentSpecification instrument :
        design.getInstrumentList().getAll()) {
      System.out.println("instrument_" + index + "_tag=" + instrument.getTagNumber());
      System.out.println("instrument_" + index + "_isa_symbol=" + instrument.getIsaSymbol());
      System.out.println("instrument_" + index + "_service=" + instrument.getService());
      System.out.println("instrument_" + index + "_range_min=" + instrument.getRangeMin());
      System.out.println("instrument_" + index + "_range_max=" + instrument.getRangeMax());
      System.out.println("instrument_" + index + "_range_unit=" + instrument.getRangeUnit());
      System.out.println("instrument_" + index + "_io_type=" + instrument.getIoType());
      System.out.println("instrument_" + index + "_sil_rating=" + instrument.getSilRating());
      System.out.println(
          "instrument_" + index + "_safety_related=" + instrument.isSafetyRelated());
      System.out.println(
          "instrument_" + index + "_device_class=" + instrument.getDeviceClassName());
      System.out.println(
          "instrument_" + index + "_estimated_cost_usd=" + instrument.getEstimatedCostUSD());
      index++;
    }
    System.out.println();
  }

  /// `Heater`'s electrical design: the duty as shaft power, the motor table, and the type branches.
  ///
  /// **The electric heater's "motor" is the element bank.** `getProcessShaftPowerKW` returns the
  /// duty in kW, the base class sizes a standard IEC motor on it, and the motor's tabulated
  /// efficiency and power factor are what the electrical input and current come from. A
  /// shell-and-tube heater has no motor: its shaft power is zero and only the auxiliary load is
  /// left. The air-cooler branch sizes fan power at `max(2 x fans, 1% duty)`.
  static void heaterElectricalRows() {
    String[] names = new String[] {"methane", "n-butane"};
    double[] z = new double[] {0.9, 0.1};
    for (double dutyKW : new double[] {5.0, 100.0, 500.0}) {
      Stream inlet = feed(names, z, 320.0, 30.0, dutyKW <= 5.0 ? 1.0 : 100.0);
      neqsim.process.equipment.heatexchanger.Heater heater =
          new neqsim.process.equipment.heatexchanger.Heater("h1", inlet);
      heater.setDuty(dutyKW * 1000.0);
      heater.run();
      neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign design =
          heater.getElectricalDesign();
      design.calcDesign();
      printHeaterElectrical("electric_" + (int) dutyKW + "kw", design);
    }

    Stream shellInlet = feed(names, z, 320.0, 30.0, 100.0);
    neqsim.process.equipment.heatexchanger.Heater shell =
        new neqsim.process.equipment.heatexchanger.Heater("h1", shellInlet);
    shell.setDuty(100_000.0);
    shell.run();
    neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign shellDesign =
        shell.getElectricalDesign();
    shellDesign.setHeatExchangerType(
        neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign.HeatExchangerType.SHELL_AND_TUBE);
    shellDesign.setCoolingWaterPumpKW(3.0);
    shellDesign.calcDesign();
    printHeaterElectrical("shell_and_tube_aux", shellDesign);

    Stream airInlet = feed(names, z, 320.0, 30.0, 100.0);
    neqsim.process.equipment.heatexchanger.Heater air =
        new neqsim.process.equipment.heatexchanger.Heater("h1", airInlet);
    air.setDuty(100_000.0);
    air.run();
    neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign airDesign =
        air.getElectricalDesign();
    airDesign.setHeatExchangerType(
        neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign.HeatExchangerType.AIR_COOLER);
    airDesign.calcDesign();
    printHeaterElectrical("air_cooler_fans", airDesign);

    Stream ie4Inlet = feed(names, z, 320.0, 30.0, 100.0);
    neqsim.process.equipment.heatexchanger.Heater ie4 =
        new neqsim.process.equipment.heatexchanger.Heater("h1", ie4Inlet);
    ie4.setDuty(100_000.0);
    ie4.run();
    neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign ie4Design =
        ie4.getElectricalDesign();
    ie4Design.getMotor().setEfficiencyClass("IE4");
    ie4Design.calcDesign();
    printHeaterElectrical("ie4_100kw", ie4Design);

    Stream polesInlet = feed(names, z, 320.0, 30.0, 100.0);
    neqsim.process.equipment.heatexchanger.Heater poles =
        new neqsim.process.equipment.heatexchanger.Heater("h1", polesInlet);
    poles.setDuty(100_000.0);
    poles.run();
    neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign polesDesign =
        poles.getElectricalDesign();
    polesDesign.getMotor().setPoles(2);
    polesDesign.calcDesign();
    printHeaterElectrical("two_pole_100kw", polesDesign);

    Stream singleInlet = feed(names, z, 320.0, 30.0, 100.0);
    neqsim.process.equipment.heatexchanger.Heater single =
        new neqsim.process.equipment.heatexchanger.Heater("h1", singleInlet);
    single.setDuty(100_000.0);
    single.run();
    neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign singleDesign =
        single.getElectricalDesign();
    singleDesign.setPhases(1);
    singleDesign.calcDesign();
    printHeaterElectrical("single_phase_100kw", singleDesign);

    Stream hzInlet = feed(names, z, 320.0, 30.0, 100.0);
    neqsim.process.equipment.heatexchanger.Heater hz =
        new neqsim.process.equipment.heatexchanger.Heater("h1", hzInlet);
    hz.setDuty(100_000.0);
    hz.run();
    neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign hzDesign =
        hz.getElectricalDesign();
    hzDesign.setFrequencyHz(60.0);
    hzDesign.calcDesign();
    printHeaterElectrical("sixty_hz_100kw", hzDesign);

    for (double lengthM : new double[] {500.0, 2000.0}) {
      Stream inlet = feed(names, z, 320.0, 30.0, 100.0);
      neqsim.process.equipment.heatexchanger.Heater cable =
          new neqsim.process.equipment.heatexchanger.Heater("h1", inlet);
      cable.setDuty(100_000.0);
      cable.run();
      neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign cableDesign =
          cable.getElectricalDesign();
      cableDesign.getPowerCable().setLengthM(lengthM);
      cableDesign.calcDesign();
      printHeaterElectrical("cable_" + (int) lengthM + "m", cableDesign);
    }
  }

  /// `useVFD = true` on the heater's electrical design: the topology arms and their tables.
  ///
  /// The electric-heater voltage taps make the 2-level AFE arm (P > 250 kW at 3300 V or less)
  /// unreachable from a heater, so that block is an air cooler sized to 300 kW of fan power at
  /// 690 V. The rated current is the motor's, and `calcDesign` sets the motor's rated voltage
  /// after it sizes the motor, so it is computed at the motor's own 400 V default.
  static void heaterVfdRows() {
    String[] names = new String[] {"methane", "n-butane"};
    double[] z = new double[] {0.9, 0.1};
    double[] dutiesKW = new double[] {5.0, 100.0, 500.0, 2500.0};
    for (double dutyKW : dutiesKW) {
      double flow = dutyKW <= 5.0 ? 1.0 : (dutyKW >= 2500.0 ? 1000.0 : 100.0);
      Stream inlet = feed(names, z, 320.0, 30.0, flow);
      neqsim.process.equipment.heatexchanger.Heater heater =
          new neqsim.process.equipment.heatexchanger.Heater("h1", inlet);
      heater.setDuty(dutyKW * 1000.0);
      heater.run();
      neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign design =
          heater.getElectricalDesign();
      design.setUseVFD(true);
      design.calcDesign();
      printHeaterElectrical("electric_" + (int) dutyKW + "kw_vfd", design);
    }

    Stream airInlet = feed(names, z, 320.0, 30.0, 10000.0);
    neqsim.process.equipment.heatexchanger.Heater air =
        new neqsim.process.equipment.heatexchanger.Heater("h1", airInlet);
    air.setDuty(30_000_000.0);
    air.run();
    neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign airDesign =
        air.getElectricalDesign();
    airDesign.setHeatExchangerType(
        neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign.HeatExchangerType.AIR_COOLER);
    airDesign.setUseVFD(true);
    airDesign.calcDesign();
    printHeaterElectrical("air_cooler_30mw_vfd", airDesign);
  }

  static void printHeaterElectrical(String label,
      neqsim.process.electricaldesign.heatexchanger.HeatExchangerElectricalDesign design) {
    System.out.println("label=" + label);
    System.out.println("electrical_type=" + design.getHeatExchangerType().name());
    System.out.println("shaft_power_kw=" + design.getShaftPowerKW());
    System.out.println("electrical_input_kw=" + design.getElectricalInputKW());
    System.out.println("apparent_power_kva=" + design.getApparentPowerKVA());
    System.out.println("reactive_power_kvar=" + design.getReactivePowerKVAR());
    System.out.println("power_factor=" + design.getPowerFactor());
    System.out.println("rated_voltage_v=" + design.getRatedVoltageV());
    System.out.println("frequency_hz=" + design.getFrequencyHz());
    System.out.println("phases=" + design.getPhases());
    System.out.println("full_load_current_a=" + design.getFullLoadCurrentA());
    System.out.println("starting_current_a=" + design.getStartingCurrentA());
    System.out.println("total_electrical_losses_kw=" + design.getTotalElectricalLossesKW());
    System.out.println("motor_rated_power_kw=" + design.getMotor().getRatedPowerKW());
    System.out.println("motor_efficiency_percent=" + design.getMotor().getEfficiencyPercent());
    System.out.println("motor_power_factor_fl=" + design.getMotor().getPowerFactorFL());
    System.out.println("motor_rated_speed_rpm=" + design.getMotor().getRatedSpeedRPM());
    System.out.println("motor_poles=" + design.getMotor().getPoles());
    System.out.println("motor_efficiency_class=" + design.getMotor().getEfficiencyClass());
    System.out.println(
        "motor_locked_rotor_multiplier=" + design.getMotor().getLockedRotorCurrentMultiplier());
    System.out.println("motor_frame=" + design.getMotor().getFrameSize());
    System.out.println("total_auxiliary_kw=" + design.getTotalAuxiliaryKW());
    System.out.println("total_connected_load_kw=" + design.getTotalConnectedLoadKW());
    System.out.println("cable_length_m=" + design.getPowerCable().getLengthM());
    System.out.println("cable_cross_section_mm2=" + design.getPowerCable().getCrossSectionMM2());
    System.out.println("cable_ampacity_a=" + design.getPowerCable().getAmpacityA());
    System.out.println("cable_voltage_drop_percent="
        + design.getPowerCable().getVoltageDropPercent());
    System.out.println("cable_short_circuit_withstand_ka="
        + design.getPowerCable().getShortCircuitWithstandKA());
    System.out.println("switchgear_starter=" + design.getSwitchgear().getStarterType());
    System.out.println("switchgear_type=" + design.getSwitchgear().getSwitchgearType());
    System.out.println("switchgear_circuit_breaker_a="
        + design.getSwitchgear().getCircuitBreakerRatingA());
    System.out.println("switchgear_fuse_a=" + design.getSwitchgear().getFuseRatingA());
    System.out.println("switchgear_rated_current_a=" + design.getSwitchgear().getRatedCurrentA());
    System.out.println("switchgear_short_circuit_ka="
        + design.getSwitchgear().getShortCircuitCurrentKA());
    System.out.println("switchgear_weight_kg=" + design.getSwitchgear().getWeightKg());
    if (design.isUseVFD() && design.getVfd() != null) {
      neqsim.process.electricaldesign.components.VariableFrequencyDrive vfd = design.getVfd();
      System.out.println("vfd_rated_power_kw=" + vfd.getRatedPowerKW());
      System.out.println("vfd_rated_current_a=" + vfd.getRatedCurrentA());
      System.out.println("vfd_input_voltage_v=" + vfd.getInputVoltageV());
      System.out.println("vfd_topology=" + vfd.getTopologyType());
      System.out.println("vfd_pulse_configuration=" + vfd.getPulseConfiguration());
      System.out.println("vfd_efficiency_percent=" + vfd.getEfficiencyPercent());
      System.out.println("vfd_input_power_factor=" + vfd.getInputPowerFactor());
      System.out.println("vfd_thd_current_percent=" + vfd.getThdCurrentPercent());
      System.out.println("vfd_has_active_rectifier=" + vfd.isHasActiveRectifier());
      System.out.println("vfd_requires_input_filter=" + vfd.isRequiresInputFilter());
      System.out.println("vfd_heat_dissipation_kw=" + vfd.getHeatDissipationKW());
      System.out.println("vfd_cooling_method=" + vfd.getCoolingMethod());
      System.out.println("vfd_weight_kg=" + vfd.getWeightKg());
      System.out.println("vfd_enclosure_rating=" + vfd.getEnclosureRating());
    }
    System.out.println();
  }

  /// `Heater`'s utility specification: five optional values and their `has*` predicates.
  ///
  /// **No arithmetic.** The class stores a supply and a return temperature, an approach, a
  /// heat-capacity rate and an overall heat-transfer coefficient, each `NaN` until set. The
  /// temperatures and the approach are "specified" when not `NaN`; the two rates are specified
  /// only when they are also positive.
  static void heaterUtilityRows() {
    String[] names = new String[] {"methane", "n-butane"};
    double[] z = new double[] {0.9, 0.1};
    String[] labels = new String[] {
        "default", "supply_450", "return_350", "approach_20", "heat_capacity_500",
        "heat_capacity_0", "u_value_100", "u_value_negative"};
    for (String label : labels) {
      Stream inlet = feed(names, z, 320.0, 30.0, 1.0);
      neqsim.process.equipment.heatexchanger.Heater heater =
          new neqsim.process.equipment.heatexchanger.Heater("h1", inlet);
      heater.run();
      neqsim.process.equipment.heatexchanger.UtilityStreamSpecification utility =
          heater.getUtilitySpecification();
      switch (label) {
        case "supply_450":
          utility.setSupplyTemperature(450.0);
          break;
        case "return_350":
          utility.setReturnTemperature(350.0);
          break;
        case "approach_20":
          utility.setApproachTemperature(20.0);
          break;
        case "heat_capacity_500":
          utility.setHeatCapacityRate(500.0);
          break;
        case "heat_capacity_0":
          utility.setHeatCapacityRate(0.0);
          break;
        case "u_value_100":
          utility.setOverallHeatTransferCoefficient(100.0);
          break;
        case "u_value_negative":
          utility.setOverallHeatTransferCoefficient(-1.0);
          break;
        default:
          break;
      }
      System.out.println("case=" + label);
      System.out.println("supply_temperature_k=" + utility.getSupplyTemperature());
      System.out.println("return_temperature_k=" + utility.getReturnTemperature());
      System.out.println("approach_temperature_k=" + utility.getApproachTemperature());
      System.out.println("heat_capacity_rate_w_per_k=" + utility.getHeatCapacityRate());
      System.out.println("overall_heat_transfer_coefficient="
          + utility.getOverallHeatTransferCoefficient());
      System.out.println("has_supply_temperature=" + utility.hasSupplyTemperature());
      System.out.println("has_return_temperature=" + utility.hasReturnTemperature());
      System.out.println("has_approach_temperature=" + utility.hasApproachTemperature());
      System.out.println("has_heat_capacity_rate=" + utility.hasHeatCapacityRate());
      System.out.println("has_overall_heat_transfer_coefficient="
          + utility.hasOverallHeatTransferCoefficient());
      System.out.println();
    }
  }

  /// `HeatExchangerMechanicalDesign`'s single-stream thermal estimate.
  ///
  /// The heater path: the utility specification's supply and return temperatures and the
  /// process temperatures give an LMTD and an approach; UA is duty/LMTD; the five candidate
  /// types turn that into a required area with their own allowable approach and typical U. The
  /// detailed calculator and the two-stream exchanger branch are separate.
  static void heaterMechanicalRows() {
    heaterMechanicalRow("no_utility", 100_000.0, false, null, null, null, null, null);
    heaterMechanicalRow("supply_450", 100_000.0, true, 450.0, null, null, null, null);
    heaterMechanicalRow("supply_and_capacity", 100_000.0, true, 450.0, 5000.0, null, null, null);
    heaterMechanicalRow("supply_and_u", 100_000.0, true, 450.0, null, 300.0, null, null);
    heaterMechanicalRow("approach_and_capacity", 100_000.0, true, null, 6000.0, null, 20.0, null);
    heaterMechanicalRow("manual_plate", 100_000.0, true, 450.0, null, null, null,
        neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerType.PLATE_AND_FRAME);
    heaterMechanicalRow("cooling", -100_000.0, true, 250.0, null, null, null, null);
  }

  static void heaterMechanicalRow(String label, double dutyW, boolean utility, Double supplyK,
      Double capacityWK, Double uValue, Double approachK,
      neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerType manual) {
    String[] names = new String[] {"methane", "n-butane"};
    double[] z = new double[] {0.9, 0.1};
    Stream inlet = feed(names, z, 320.0, 30.0, 100.0);
    neqsim.process.equipment.heatexchanger.Heater heater =
        new neqsim.process.equipment.heatexchanger.Heater("h1", inlet);
    heater.setDuty(dutyW);
    if (utility) {
      neqsim.process.equipment.heatexchanger.UtilityStreamSpecification spec =
          heater.getUtilitySpecification();
      if (supplyK != null) {
        spec.setSupplyTemperature(supplyK);
      }
      if (capacityWK != null) {
        spec.setHeatCapacityRate(capacityWK);
      }
      if (uValue != null) {
        spec.setOverallHeatTransferCoefficient(uValue);
      }
      if (approachK != null) {
        spec.setApproachTemperature(approachK);
      }
    }
    heater.run();
    neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerMechanicalDesign design =
        heater.getMechanicalDesign();
    if (manual != null) {
      design.setManualSelection(manual);
    }
    design.calcDesign();

    System.out.println("case=" + label);
    System.out.println("duty_w=" + dutyW);
    System.out.println("process_in_k=" + inlet.getTemperature("K"));
    System.out.println("process_out_k=" + heater.getOutletStream().getTemperature("K"));
    System.out.println("utility_supply_k="
        + heater.getUtilitySpecification().getSupplyTemperature());
    System.out.println("utility_return_k="
        + heater.getUtilitySpecification().getReturnTemperature());
    System.out.println("log_mean_temperature_difference_k="
        + design.getLogMeanTemperatureDifference());
    System.out.println("approach_temperature_k=" + design.getApproachTemperature());
    System.out.println("used_u_w_per_m2k=" + design.getUsedOverallHeatTransferCoefficient());
    System.out.println("calculated_ua_w_per_k=" + design.getCalculatedUA());
    for (neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerSizingResult result :
        design.getSizingResults()) {
      String type = result.getType().name().toLowerCase(java.util.Locale.ROOT);
      System.out.println("type_" + type + "_required_area_m2=" + result.getRequiredArea());
      System.out.println("type_" + type + "_required_ua_w_per_k=" + result.getRequiredUA());
      System.out.println("type_" + type + "_u_w_per_m2k="
          + result.getOverallHeatTransferCoefficient());
      System.out.println("type_" + type + "_approach_k=" + result.getApproachTemperature());
      System.out.println("type_" + type + "_inner_diameter_m=" + result.getInnerDiameter());
      System.out.println("type_" + type + "_outer_diameter_m=" + result.getOuterDiameter());
      System.out.println("type_" + type + "_wall_thickness_m=" + result.getWallThickness());
      System.out.println("type_" + type + "_length_m=" + result.getEstimatedLength());
      System.out.println("type_" + type + "_tube_count=" + result.getTubeCount());
      System.out.println("type_" + type + "_tube_passes=" + result.getTubePasses());
      System.out.println("type_" + type + "_fin_surface_area_m2=" + result.getFinSurfaceArea());
      System.out.println("type_" + type + "_estimated_pressure_drop="
          + result.getEstimatedPressureDrop());
      System.out.println("type_" + type + "_estimated_weight_kg=" + result.getEstimatedWeight());
      System.out.println("type_" + type + "_module_length_m=" + result.getModuleLength());
      System.out.println("type_" + type + "_module_width_m=" + result.getModuleWidth());
      System.out.println("type_" + type + "_module_height_m=" + result.getModuleHeight());
    }
    neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerSizingResult selected =
        design.getSelectedSizingResult();
    System.out.println("selected_type="
        + (selected == null ? "none" : selected.getType().name().toLowerCase(java.util.Locale.ROOT)));
    System.out.println("selected_area_m2="
        + (selected == null ? "none" : selected.getRequiredArea()));
    System.out.println();
  }

  /// `HeatExchangerMechanicalDesign.calcDesign`'s two-stream branch.
  ///
  /// **The branch is selected by the equipment, not by an input**: `calcDesign` asks
  /// `instanceof HeatExchanger` first, so a two-stream exchanger takes it and a `Heater`
  /// takes the single-stream one. The temperatures are the exchanger's own four stream
  /// temperatures, ordered by their mean; the UA the estimate starts from is
  /// `getUAvalue()`; and the sizing coefficient is the **type's typical** U rather than the
  /// larger of the type's and a used U - the flag the heater branch calls with `false`.
  ///
  /// **Three things the rows are here to settle.** `calcDesign` resets
  /// `usedOverallHeatTransferCoefficient` to the class's 500 and the two-stream handler
  /// never touches it, so the printed used U is 500 even when the exchanger's UA is 100.
  /// The effectiveness-derived UA branch (`cmin * -ln(1 - effectiveness)`) only replaces
  /// the starting UA when it is **larger**, which the rating that produced it makes
  /// impossible; a pinned or delta-T run leaves the effectiveness field at zero, so the
  /// branch is skipped there. And the pinned rows read the exchanger's `UAvalue` - the
  /// stated one, or the class's 500 default when none was stated - even though `run`
  /// itself never reads it in that mode.
  static void heatExchangerMechanicalRows() {
    hxMechanicalRow("rating_ua_100", 100.0, true, null, null, null, false);
    hxMechanicalRow("rating_ua_300", 300.0, true, null, null, null, false);
    // Near-saturated: the effectiveness is above 0.999 and the approach is 8.5e-5 K, so this
    // row is evidence for the arm's skip rather than a case.
    hxMechanicalRow("rating_ua_1000", 1000.0, true, null, null, null, false);
    // The manual pick is not the area winner here: plate-and-frame is.
    hxMechanicalRow("manual_shell_and_tube", 100.0, true, "SHELL_AND_TUBE", null, null, false);
    hxMechanicalRow("min_weight", 100.0, true, null, "MIN_WEIGHT", null, false);
    hxMechanicalRow("min_pressure_drop", 100.0, true, null, "MIN_PRESSURE_DROP", null, false);
    // A pinned hot outlet: `run` never recomputes `UAvalue`, so the estimate reads the
    // class's 500 default when the caller stated none.
    hxMechanicalRow("pinned_hot_outlet_default_ua", 0.0, false, null, null, 350.0, false);
    hxMechanicalRow("pinned_hot_outlet_stated_ua", 200.0, true, null, null, 350.0, false);
    // The delta-T closure writes `UAvalue = duty / (T_out,set - T_in,set)` and leaves
    // `thermalEffectiveness` at zero; this row is what the estimate reads there. The
    // equal-capacity fluids keep the closure's state non-degenerate.
    hxMechanicalRow("delta_t", 0.0, false, null, null, null, true);
  }

  static void hxMechanicalRow(String label, double ua, boolean setUa, String manual,
      String criterion, Double pinnedK, boolean equalCapacity) {
    Stream hot;
    Stream cold;
    if (equalCapacity) {
      hot = feed(new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5}, 400.0, 20.0, 1.0);
      cold = feed(new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5}, 300.0, 20.0, 1.0);
    } else {
      hot = feed(new String[] {"methane", "n-butane"}, new double[] {0.8, 0.2}, 400.0, 20.0, 1.0);
      cold = feed(new String[] {"n-butane", "n-pentane"}, new double[] {0.5, 0.5}, 300.0, 5.0,
          1.0);
    }
    hot.setName("hot");
    cold.setName("cold");

    neqsim.process.equipment.heatexchanger.HeatExchanger hx =
        new neqsim.process.equipment.heatexchanger.HeatExchanger("hx1", hot, cold);
    if (setUa) {
      hx.setUAvalue(ua);
    }
    if (pinnedK != null) {
      // **The reference state, not `runSpecifiedStream`'s.** That method clones an
      // already-flashed fluid and lands its enthalpy elsewhere, so the port's own pinned solve
      // reaches the fresh-fluid state instead - and this block hands the class exactly the two
      // outlets the port carries: the pin, and the other side at the duty's enthalpy.
      Stream pinned = feed(new String[] {"methane", "n-butane"}, new double[] {0.8, 0.2}, pinnedK,
          20.0, 1.0);
      double hotSwing = duty(hot, pinned);
      SystemInterface coldIn = cold.getThermoSystem();
      SystemInterface coldOut = phFlash(new String[] {"n-butane", "n-pentane"},
          new double[] {0.5, 0.5}, 300.0, 5.0, 1.0,
          coldIn.getEnthalpy() / coldIn.getTotalNumberOfMoles() - hotSwing);
      hx.getOutStream(0).setThermoSystem(pinned.getThermoSystem());
      hx.getOutStream(1).setThermoSystem(coldOut);
    } else {
      if ("delta_t".equals(label)) {
        hx.setDeltaT(1.0);
      }
      hx.run();
    }

    hx.initMechanicalDesign();
    neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerMechanicalDesign design =
        hx.getMechanicalDesign();
    if (manual != null) {
      design.setManualSelection(
          neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerType.valueOf(manual));
    }
    if (criterion != null) {
      design.setSelectionCriterion(
          neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerMechanicalDesign.SelectionCriterion
              .valueOf(criterion));
    }
    design.calcDesign();

    System.out.println("case=" + label);
    System.out.println("duty_w=" + Math.abs(hx.getDuty()));
    System.out.println("hot_in_k=" + hx.getInTemperature(0));
    System.out.println("hot_out_k=" + hx.getOutTemperature(0));
    System.out.println("cold_in_k=" + hx.getInTemperature(1));
    System.out.println("cold_out_k=" + hx.getOutTemperature(1));
    System.out.println("ua_value_w_per_k=" + hx.getUAvalue());
    System.out.println("thermal_effectiveness=" + hx.getThermalEffectiveness());
    System.out.println("log_mean_temperature_difference_k="
        + design.getLogMeanTemperatureDifference());
    System.out.println("approach_temperature_k=" + design.getApproachTemperature());
    System.out.println("used_u_w_per_m2k=" + design.getUsedOverallHeatTransferCoefficient());
    System.out.println("calculated_ua_w_per_k=" + design.getCalculatedUA());
    for (neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerSizingResult result :
        design.getSizingResults()) {
      String type = result.getType().name().toLowerCase(java.util.Locale.ROOT);
      System.out.println("type_" + type + "_required_area_m2=" + result.getRequiredArea());
      System.out.println("type_" + type + "_required_ua_w_per_k=" + result.getRequiredUA());
      System.out.println("type_" + type + "_u_w_per_m2k="
          + result.getOverallHeatTransferCoefficient());
      System.out.println("type_" + type + "_approach_k=" + result.getApproachTemperature());
      System.out.println("type_" + type + "_inner_diameter_m=" + result.getInnerDiameter());
      System.out.println("type_" + type + "_outer_diameter_m=" + result.getOuterDiameter());
      System.out.println("type_" + type + "_wall_thickness_m=" + result.getWallThickness());
      System.out.println("type_" + type + "_length_m=" + result.getEstimatedLength());
      System.out.println("type_" + type + "_tube_count=" + result.getTubeCount());
      System.out.println("type_" + type + "_tube_passes=" + result.getTubePasses());
      System.out.println("type_" + type + "_fin_surface_area_m2=" + result.getFinSurfaceArea());
      System.out.println("type_" + type + "_estimated_pressure_drop="
          + result.getEstimatedPressureDrop());
      System.out.println("type_" + type + "_estimated_weight_kg=" + result.getEstimatedWeight());
      System.out.println("type_" + type + "_module_length_m=" + result.getModuleLength());
      System.out.println("type_" + type + "_module_width_m=" + result.getModuleWidth());
      System.out.println("type_" + type + "_module_height_m=" + result.getModuleHeight());
    }
    neqsim.process.mechanicaldesign.heatexchanger.HeatExchangerSizingResult selected =
        design.getSelectedSizingResult();
    System.out.println("selected_type="
        + (selected == null ? "none" : selected.getType().name().toLowerCase(java.util.Locale.ROOT)));
    System.out.println("selected_area_m2="
        + (selected == null ? "none" : selected.getRequiredArea()));
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

  /// `FilterPressureDropModel`'s four branches and the particle-capture update.
  ///
  /// **The rows are `FilterTest`'s states re-cased on PR methane**, because the port's fluid is
  /// PR. The drop models are fluid-free except the Ergun branch, which reads the density and
  /// viscosity the class's system reports; those are printed and the port takes the viscosity
  /// as an input, because azoth's viscosity dispatch is its own calculation.
  static void filterModelRows() {
    neqsim.process.equipment.filter.Filter fixed = new neqsim.process.equipment.filter.Filter(
        "fixed", filterPrFeed(1000.0));
    fixed.setDeltaP(0.10);
    fixed.run();
    printFilterRow("fixed_methane_1000_kg_per_hr", fixed.getInletStream(), fixed);

    for (double kgPerHour : new double[] {1000.0, 2000.0}) {
      Stream inlet = filterPrFeed(kgPerHour);
      neqsim.process.equipment.filter.Filter filter =
          new neqsim.process.equipment.filter.Filter("flow scaled", inlet);
      filter.setFilterServiceType(neqsim.process.equipment.filter.FilterType.Y_STRAINER);
      filter.setDeltaP(0.10);
      filter.setReferenceFlowRate(flowM3Hr(filterPrFeed(1000.0)));
      filter.run();
      printFilterRow("flow_scaled_methane_" + (int) kgPerHour + "_kg_per_hr",
          filter.getInletStream(), filter);
    }

    Stream tabulatedInlet = filterPrFeed(1000.0);
    double tabulatedFlow = flowM3Hr(tabulatedInlet);
    neqsim.process.equipment.filter.FilterPressureDropCurve curve =
        new neqsim.process.equipment.filter.FilterPressureDropCurve(
            new double[] {0.5 * tabulatedFlow, 1.5 * tabulatedFlow}, new double[] {0.10, 0.50});
    curve.setTestStandard("ISO 3968:2017");
    neqsim.process.equipment.filter.Filter tabulated =
        new neqsim.process.equipment.filter.Filter("tested cartridge", tabulatedInlet);
    tabulated.setPressureDropCurve(curve);
    tabulated.run();
    printFilterRow("tabulated_methane_1000_kg_per_hr", tabulated.getInletStream(), tabulated);

    for (double kgPerHour : new double[] {1000.0, 2000.0}) {
      Stream inlet = filterPrFeed(kgPerHour);
      neqsim.process.equipment.filter.Filter filter =
          new neqsim.process.equipment.filter.Filter("granular guard bed", inlet);
      filter.setFilterServiceType(neqsim.process.equipment.filter.FilterType.GRANULAR_MEDIA);
      filter.setMediaGeometry(1.0, 1.5, 0.002, 0.40);
      filter.run();
      printFilterRow("ergun_methane_" + (int) kgPerHour + "_kg_per_hr",
          filter.getInletStream(), filter);
    }

    Stream captureInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.FilterPerformanceCurve performance =
        new neqsim.process.equipment.filter.FilterPerformanceCurve(
            new double[] {5.0, 10.0, 20.0}, new double[] {2.0, 100.0, 1000.0});
    performance.setTestStandard("ISO 16889:2022");
    neqsim.process.equipment.filter.Filter capture =
        new neqsim.process.equipment.filter.Filter("beta rated filter", captureInlet);
    capture.setPerformanceCurve(performance);
    capture.setParticleSize(10.0);
    capture.setInletParticleConcentration(100.0);
    capture.run();
    System.out.println("capture_beta_ratio_methane_1000_kg_per_hr");
    System.out.println("inlet_m3_per_hr=" + flowM3Hr(captureInlet));
    System.out.println("inlet_density_kg_per_m3="
        + captureInlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("inlet_viscosity_pa_s="
        + captureInlet.getThermoSystem().getViscosity("kg/msec"));
    System.out.println("particle_size_micrometre=" + capture.getParticleSize());
    System.out.println("inlet_particle_concentration_mg_per_kg="
        + capture.getInletParticleConcentration());
    System.out.println("nominal_removal_efficiency=" + capture.getNominalRemovalEfficiency());
    System.out.println("outlet_particle_concentration_mg_per_kg="
        + capture.getOutletParticleConcentration());
    System.out.println("captured_rate_kg_per_hr=" + capture.getCalculatedCapturedRate());
    print("inlet", captureInlet);
    print("outlet", capture.getOutletStream());
    System.out.println();
  }

  /// `Filter`'s steady loaded state: `run()` reading the carried loading, with no time step.
  ///
  /// **The blocks are `FilterTest`'s scenarios re-cased on PR methane**, as `filterModelRows`'
  /// are, because the port's fluid is PR. The loaded block is the test's pre-step state, the
  /// bypass block is `testBypassLimitsPressureDropAndReducesCapture`, and the replacement blocks
  /// reach the two branches of `isReplacementRequired`. The last block is evidence and is not
  /// cased: the clean drop exceeds the inlet, so `run` clamps the bypass-limited drop to the
  /// inlet ceiling and no cubic flash the port carries survives there.
  static void filterLoadingRows() {
    Stream loadedInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter loaded =
        new neqsim.process.equipment.filter.Filter("loaded", loadedInlet);
    loaded.setDeltaP(0.2);
    loaded.setLoadingCapacity(10.0);
    loaded.setPressureDropIncreaseAtCapacity(1.0);
    loaded.setBreakthroughStartFraction(0.5);
    loaded.setSolidsLoading(9.0);
    loaded.run();
    printFilterLoading("loaded_run", loadedInlet, loaded);

    Stream bypassInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter bypass =
        new neqsim.process.equipment.filter.Filter("bypass", bypassInlet);
    bypass.setDeltaP(0.2);
    bypass.setNominalRemovalEfficiency(0.99);
    bypass.setInletParticleConcentration(100.0);
    bypass.setLoadingCapacity(10.0);
    bypass.setPressureDropIncreaseAtCapacity(3.6);
    bypass.setSolidsLoading(5.0);
    bypass.setBypassCrackingDeltaP(0.5);
    bypass.run();
    printFilterLoading("bypass_run", bypassInlet, bypass);

    Stream disabledInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter disabled =
        new neqsim.process.equipment.filter.Filter("disabled", disabledInlet);
    disabled.setDeltaP(0.2);
    disabled.setPressureDropIncreaseAtCapacity(1.0);
    disabled.setBreakthroughStartFraction(0.5);
    disabled.setSolidsLoading(9.0);
    disabled.setLoadingCapacity(0.0);
    disabled.run();
    printFilterLoading("capacity_disabled", disabledInlet, disabled);

    Stream terminalInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter terminal =
        new neqsim.process.equipment.filter.Filter("terminal", terminalInlet);
    terminal.setDeltaP(0.2);
    terminal.setLoadingCapacity(10.0);
    terminal.setPressureDropIncreaseAtCapacity(3.6);
    terminal.setSolidsLoading(5.0);
    terminal.setBypassCrackingDeltaP(0.5);
    terminal.setTerminalDeltaP(0.5);
    terminal.run();
    printFilterLoading("terminal_drop", terminalInlet, terminal);

    Stream overInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter over =
        new neqsim.process.equipment.filter.Filter("over", overInlet);
    over.setDeltaP(0.2);
    over.setLoadingCapacity(10.0);
    over.setPressureDropIncreaseAtCapacity(1.0);
    over.setBreakthroughStartFraction(0.5);
    over.setSolidsLoading(11.0);
    over.run();
    printFilterLoading("loading_beyond_capacity", overInlet, over);

    Stream boundaryInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter boundary =
        new neqsim.process.equipment.filter.Filter("boundary", boundaryInlet);
    boundary.setDeltaP(0.2);
    boundary.setLoadingCapacity(10.0);
    boundary.setPressureDropIncreaseAtCapacity(1.0);
    boundary.setBreakthroughStartFraction(0.5);
    boundary.setSolidsLoading(5.0);
    boundary.run();
    printFilterLoading("breakthrough_boundary", boundaryInlet, boundary);

    Stream clampInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter clamp =
        new neqsim.process.equipment.filter.Filter("clamp", clampInlet);
    clamp.setDeltaP(26.0);
    clamp.setBypassCrackingDeltaP(25.0);
    clamp.run();
    printFilterLoading("clamp_after_bypass_uncased", clampInlet, clamp);
  }

  static void printFilterLoading(String label, StreamInterface inlet,
      neqsim.process.equipment.filter.Filter filter) {
    System.out.println(label);
    System.out.println("inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("inlet_density_kg_per_m3="
        + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("clean_deltaP_pa=" + filter.getCalculatedCleanDeltaP() * 1.0e5);
    System.out.println("unrestricted_deltaP_pa=" + filter.getUnrestrictedDeltaP() * 1.0e5);
    System.out.println("applied_deltaP_pa=" + filter.getDeltaP() * 1.0e5);
    System.out.println("bypass_enabled=" + filter.isBypassEnabled());
    System.out.println("bypass_fraction=" + filter.getBypassFraction());
    System.out.println("loading_fraction=" + filter.getLoadingFraction());
    System.out.println("breakthrough_fraction=" + filter.getBreakthroughFraction());
    System.out.println("current_removal_efficiency=" + filter.getCurrentRemovalEfficiency());
    System.out.println("nominal_removal_efficiency=" + filter.getNominalRemovalEfficiency());
    System.out.println("solids_loading_kg=" + filter.getSolidsLoading());
    System.out.println("differential_pressure_utilization="
        + filter.getDifferentialPressureUtilization());
    System.out.println("replacement_required=" + filter.isReplacementRequired());
    System.out.println("outlet_particle_concentration_mg_per_kg="
        + filter.getOutletParticleConcentration());
    System.out.println("captured_rate_kg_per_hr=" + filter.getCalculatedCapturedRate());
    System.out.println("holdup_residence_time_s=" + filter.getHoldupResidenceTime());
    print("inlet", inlet);
    print("outlet", filter.getOutletStream());
    System.out.println();
  }

  /// `Filter.runTransient`'s `calculateSteadyState=false` branch, one step at a time.
  ///
  /// **The accretion uses the pre-step efficiency and the closing `run` reports the post-step
  /// one**, so each block prints the state before and after its step. The fourth block pins
  /// that sequencing with a concentration model whose breakthrough is live before the step; the
  /// fifth adds a bypass, whose pre-step fraction is in the accretion's efficiency too. The last
  /// block is evidence and is not cased: a zero-flow feed makes the residence time `+Infinity`,
  /// which the port publishes as absent.
  static void filterDynamicRows() {
    Stream firstInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter first =
        new neqsim.process.equipment.filter.Filter("first", firstInlet);
    first.setDeltaP(0.1);
    first.setHoldupVolume(1.0);
    first.setSolidsLoadingRate(5.0);
    first.setLoadingCapacity(10.0);
    first.setPressureDropIncreaseAtCapacity(1.0);
    first.setBreakthroughStartFraction(0.5);
    first.setCalculateSteadyState(false);
    printFilterLoading("dynamic_before_step", firstInlet, first);
    first.runTransient(3600.0, java.util.UUID.randomUUID());
    printFilterLoading("dynamic_first_step_3600", firstInlet, first);
    first.runTransient(1800.0, java.util.UUID.randomUUID());
    printFilterLoading("dynamic_second_step_1800", firstInlet, first);

    Stream regenInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter regen =
        new neqsim.process.equipment.filter.Filter("regen", regenInlet);
    regen.setDeltaP(0.2);
    regen.setSolidsLoading(9.0);
    regen.setLoadingCapacity(10.0);
    regen.setPressureDropIncreaseAtCapacity(1.0);
    regen.setBreakthroughStartFraction(0.5);
    regen.setSolidsLoadingRate(0.0);
    regen.setBackwashRemovalRate(4.0);
    regen.setRegenerationRemovalRate(2.0);
    regen.setCalculateSteadyState(false);
    printFilterLoading("backwash_before_step", regenInlet, regen);
    regen.startBackwash();
    regen.startRegeneration();
    regen.runTransient(1800.0, java.util.UUID.randomUUID());
    printFilterLoading("backwash_after_step", regenInlet, regen);

    Stream concentrationInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter concentration =
        new neqsim.process.equipment.filter.Filter("concentration", concentrationInlet);
    concentration.setDeltaP(0.1);
    concentration.setNominalRemovalEfficiency(0.99);
    concentration.setInletParticleConcentration(100.0);
    concentration.setSolidsLoading(7.9);
    concentration.setLoadingCapacity(10.0);
    concentration.setPressureDropIncreaseAtCapacity(1.0);
    concentration.setBreakthroughStartFraction(0.5);
    concentration.setCalculateSteadyState(false);
    printFilterLoading("concentration_before_step", concentrationInlet, concentration);
    concentration.runTransient(28800.0, java.util.UUID.randomUUID());
    printFilterLoading("concentration_after_8h", concentrationInlet, concentration);

    Stream bypassInlet = filterPrFeed(1000.0);
    neqsim.process.equipment.filter.Filter bypass =
        new neqsim.process.equipment.filter.Filter("bypass dynamic", bypassInlet);
    bypass.setDeltaP(0.2);
    bypass.setNominalRemovalEfficiency(0.99);
    bypass.setInletParticleConcentration(100.0);
    bypass.setSolidsLoading(5.0);
    bypass.setLoadingCapacity(10.0);
    bypass.setPressureDropIncreaseAtCapacity(3.6);
    bypass.setBypassCrackingDeltaP(0.5);
    bypass.setCalculateSteadyState(false);
    printFilterLoading("dynamic_bypass_before_step", bypassInlet, bypass);
    bypass.runTransient(3600.0, java.util.UUID.randomUUID());
    printFilterLoading("dynamic_bypass_after_step", bypassInlet, bypass);

    Stream zeroInlet = filterPrFeed(0.0);
    neqsim.process.equipment.filter.Filter zero =
        new neqsim.process.equipment.filter.Filter("zero", zeroInlet);
    zero.setDeltaP(0.1);
    zero.setHoldupVolume(1.0);
    zero.setSolidsLoadingRate(1.0);
    zero.setLoadingCapacity(10.0);
    zero.setCalculateSteadyState(false);
    zero.runTransient(3600.0, java.util.UUID.randomUUID());
    printFilterLoading("zero_flow_residence_uncased", zeroInlet, zero);
  }

  /// `FilterMechanicalDesign.calcDesign`: the element count, vessel, pressure boundary and weights.
  ///
  /// **The blocks re-case `FilterMechanicalDesignTest` on PR methane**, because the port's fluid
  /// is PR and the design reads the inlet's volumetric flow and density. The overrides block states
  /// every settable margin; the long-bed block exists to make the L/D maximum bind rather than the
  /// minimum; the liquid block crosses the nozzle velocity's density switch at `200` kg/m3.
  static void filterMechanicalDesignRows() {
    Stream coalescerInlet = filterDesignFeed(10000.0, 303.15, 50.0, false);
    neqsim.process.equipment.filter.Filter coalescer =
        new neqsim.process.equipment.filter.Filter("coalescer", coalescerInlet);
    coalescer.setFilterServiceType(neqsim.process.equipment.filter.FilterType.COALESCER);
    coalescer.setDeltaP(0.15);
    coalescer.setTerminalDeltaP(1.0);
    coalescer.setElementCollapsePressure(5.0);
    coalescer.setElementIntegrityVerified(true);
    neqsim.process.equipment.filter.FilterPerformanceCurve performance =
        new neqsim.process.equipment.filter.FilterPerformanceCurve(
            new double[] {1.0, 3.0, 10.0}, new double[] {10.0, 1000.0, 10000.0});
    coalescer.setPerformanceCurve(performance);
    coalescer.run();
    coalescer.getMechanicalDesign().calcDesign();
    printFilterDesign("coalescer_default", coalescerInlet, coalescer,
        (neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign) coalescer
            .getMechanicalDesign());

    Stream basketInlet = filterDesignFeed(10000.0, 303.15, 50.0, false);
    neqsim.process.equipment.filter.Filter basket =
        new neqsim.process.equipment.filter.Filter("basket", basketInlet);
    basket.setFilterServiceType(neqsim.process.equipment.filter.FilterType.BASKET_STRAINER);
    basket.setDeltaP(2.0);
    basket.setTerminalDeltaP(1.0);
    basket.setElementCollapsePressure(1.5);
    basket.run();
    basket.getMechanicalDesign().calcDesign();
    printFilterDesign("basket_strainer_warnings", basketInlet, basket,
        (neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign) basket
            .getMechanicalDesign());

    Stream carbonInlet = filterDesignFeed(10000.0, 303.15, 50.0, false);
    neqsim.process.equipment.filter.Filter carbon =
        new neqsim.process.equipment.filter.Filter("carbon", carbonInlet);
    carbon.setFilterServiceType(neqsim.process.equipment.filter.FilterType.ACTIVATED_CARBON);
    carbon.setElementIntegrityVerified(true);
    carbon.setNominalRemovalEfficiency(0.9);
    carbon.run();
    carbon.getMechanicalDesign().calcDesign();
    printFilterDesign("activated_carbon_granular", carbonInlet, carbon,
        (neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign) carbon
            .getMechanicalDesign());

    Stream highInlet = filterDesignFeed(50000.0, 303.15, 50.0, false);
    neqsim.process.equipment.filter.Filter high =
        new neqsim.process.equipment.filter.Filter("high", highInlet);
    high.setDeltaP(0.1);
    high.setElementIntegrityVerified(true);
    high.run();
    high.getMechanicalDesign().calcDesign();
    printFilterDesign("cartridge_high_flow", highInlet, high,
        (neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign) high
            .getMechanicalDesign());

    Stream lowInlet = filterDesignFeed(1000.0, 303.15, 50.0, false);
    neqsim.process.equipment.filter.Filter low =
        new neqsim.process.equipment.filter.Filter("low", lowInlet);
    low.setFilterServiceType(neqsim.process.equipment.filter.FilterType.Y_STRAINER);
    low.setElementIntegrityVerified(true);
    low.run();
    low.getMechanicalDesign().calcDesign();
    printFilterDesign("y_strainer_low_flow", lowInlet, low,
        (neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign) low
            .getMechanicalDesign());

    Stream twoInlet = filterDesignFeed(1000.0, 303.15, 50.0, false);
    neqsim.process.equipment.filter.Filter two =
        new neqsim.process.equipment.filter.Filter("two", twoInlet);
    two.setFilterServiceType(neqsim.process.equipment.filter.FilterType.Y_STRAINER);
    two.setNumberOfElements(2);
    two.setElementIntegrityVerified(true);
    two.run();
    two.getMechanicalDesign().calcDesign();
    printFilterDesign("two_initial_elements", twoInlet, two,
        (neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign) two
            .getMechanicalDesign());

    Stream overrideInlet = filterDesignFeed(10000.0, 303.15, 50.0, false);
    neqsim.process.equipment.filter.Filter override =
        new neqsim.process.equipment.filter.Filter("override", overrideInlet);
    override.setElementIntegrityVerified(true);
    override.run();
    neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign overrideDesign =
        (neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign) override
            .getMechanicalDesign();
    overrideDesign.setMaxFaceVelocity(0.05);
    overrideDesign.setElementArea(0.5);
    overrideDesign.setMaxNozzleVelocity(10.0);
    overrideDesign.setMinLDRatio(1.5);
    overrideDesign.setMaxLDRatio(3.0);
    overrideDesign.setDesignPressureMargin(1.2);
    overrideDesign.setMinimumDesignPressureMarginBar(5.0);
    overrideDesign.setDesignTemperatureMarginC(15.0);
    overrideDesign.setAllowableStress(100.0);
    overrideDesign.setJointEfficiency(0.7);
    overrideDesign.setCorrosionAllowanceMm(6.0);
    overrideDesign.calcDesign();
    printFilterDesign("overrides", overrideInlet, override, overrideDesign);

    Stream bedInlet = filterDesignFeed(500.0, 303.15, 50.0, false);
    neqsim.process.equipment.filter.Filter bed =
        new neqsim.process.equipment.filter.Filter("bed", bedInlet);
    bed.setFilterServiceType(neqsim.process.equipment.filter.FilterType.GRANULAR_MEDIA);
    bed.setMediaGeometry(1.0, 4.0, 0.002, 0.40);
    bed.setElementIntegrityVerified(true);
    bed.run();
    bed.getMechanicalDesign().calcDesign();
    printFilterDesign("granular_long_bed", bedInlet, bed,
        (neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign) bed
            .getMechanicalDesign());

    Stream liquidInlet = filterDesignFeed(5000.0, 300.0, 5.0, true);
    neqsim.process.equipment.filter.Filter liquid =
        new neqsim.process.equipment.filter.Filter("liquid", liquidInlet);
    liquid.setFilterServiceType(neqsim.process.equipment.filter.FilterType.BAG);
    liquid.setElementIntegrityVerified(true);
    liquid.run();
    liquid.getMechanicalDesign().calcDesign();
    printFilterDesign("liquid_bag_nozzle_density", liquidInlet, liquid,
        (neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign) liquid
            .getMechanicalDesign());
  }

  static Stream filterDesignFeed(double kgPerHour, double kelvin, double bara, boolean liquid) {
    SystemInterface fluid = new SystemPrEos(kelvin, bara);
    fluid.addComponent(liquid ? "n-butane" : "methane", 1.0);
    fluid.setMixingRule(2);
    Stream stream = new Stream("feed", fluid);
    stream.setFlowRate(kgPerHour, "kg/hr");
    stream.run();
    return stream;
  }

  static void printFilterDesign(String label, Stream inlet,
      neqsim.process.equipment.filter.Filter filter,
      neqsim.process.mechanicaldesign.filter.FilterMechanicalDesign design) {
    neqsim.process.equipment.filter.FilterType type = filter.getFilterServiceType();
    System.out.println(label);
    System.out.println("filter_type=" + type.name());
    System.out.println("type_face_velocity_mps=" + type.getDefaultMaximumFaceVelocity());
    System.out.println("type_element_area_m2=" + type.getDefaultElementArea());
    System.out.println("type_element_diameter_m=" + type.getDefaultElementDiameter());
    System.out.println("type_element_length_m=" + type.getDefaultElementLength());
    System.out.println("type_backwash_supported=" + type.isBackwashSupported());
    System.out.println("inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("inlet_density_kg_per_m3="
        + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("design_pressure_bara=" + design.getDesignPressure());
    System.out.println("design_temperature_c=" + design.getDesignTemperatureC());
    System.out.println("required_elements=" + design.getRequiredElements());
    System.out.println("inner_diameter_m=" + design.getInnerDiameter());
    System.out.println("outer_diameter_m=" + design.getOuterDiameter());
    System.out.println("vessel_length_m=" + design.getVesselLength());
    System.out.println("shell_thickness_m=" + design.getShellThickness());
    System.out.println("head_thickness_m=" + design.getHeadThickness());
    System.out.println("element_area_m2=" + design.getElementArea());
    System.out.println("filter_media_area_m2=" + filter.getMediaArea());
    System.out.println("calculated_face_velocity_mps=" + design.getCalculatedFaceVelocity());
    System.out.println("required_nozzle_diameter_m=" + design.getRequiredNozzleDiameter());
    System.out.println("selected_nozzle_diameter_mm=" + design.getSelectedNozzleDiameterMm());
    System.out.println("calculated_nozzle_velocity_mps=" + design.getCalculatedNozzleVelocity());
    System.out.println("empty_vessel_weight_kg=" + design.getEmptyVesselWeight());
    System.out.println("total_equipped_weight_kg=" + design.getTotalEquippedWeight());
    System.out.println("terminal_dp_utilization="
        + design.getTerminalDifferentialPressureUtilization());
    System.out.println("collapse_pressure_utilization="
        + design.getElementCollapsePressureUtilization());
    System.out.println("differential_pressure_design_acceptable="
        + design.isDifferentialPressureDesignAcceptable());
    System.out.println("design_warning_count=" + design.getDesignWarnings().size());
    for (String warning : design.getDesignWarnings()) {
      System.out.println("design_warning=" + warning);
    }
    print("inlet", inlet);
    print("outlet", filter.getOutletStream());
    System.out.println();
  }

  static Stream filterPrFeed(double kgPerHour) {
    SystemInterface fluid = new SystemPrEos(298.15, 20.0);
    fluid.addComponent("methane", 1.0);
    fluid.setMixingRule(2);
    Stream stream = new Stream("feed", fluid);
    stream.setFlowRate(kgPerHour, "kg/hr");
    stream.run();
    return stream;
  }

  static double flowM3Hr(Stream stream) {
    return stream.getFlowRate("m3/hr");
  }

  static void printFilterRow(String label, StreamInterface inlet,
      neqsim.process.equipment.filter.Filter filter) {
    System.out.println(label);
    System.out.println("inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("inlet_density_kg_per_m3="
        + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("inlet_viscosity_pa_s="
        + inlet.getThermoSystem().getViscosity("kg/msec"));
    System.out.println("clean_delta_p_bar=" + filter.getCalculatedCleanDeltaP());
    System.out.println("delta_p_bar=" + filter.getDeltaP());
    print("inlet", inlet);
    print("outlet", filter.getOutletStream());
    System.out.println();
  }

  /// `AdiabaticPipe`'s fittings and `setPressureOut` mode, and `PipeBeggsAndBrills`.
  ///
  /// **The rows are `AdiabaticPipeTest`'s 50 km gas state re-cased on PR** - 95/5
  /// methane/ethane at 100 bara and 25 C - because the port's fluid is PR. The fittings row
  /// adds two standard 90-degree elbows (Crane L/D 30 each), so the effective length grows by
  /// `60 D = 30` m; the pressure-out row states 80 bara and lets `calcFlow` bisect the flow.
  static void pipeModeRows() {
    pipeAdiabaticRow("adiabatic_50km", false, null);
    pipeAdiabaticRow("adiabatic_50km_two_elbows", true, null);
    pipeAdiabaticRow("adiabatic_pressure_out_80_bara", false, 80.0);

    Stream inlet = pipeFeed();
    neqsim.process.equipment.pipeline.PipeBeggsAndBrills pipe =
        new neqsim.process.equipment.pipeline.PipeBeggsAndBrills("bb", inlet);
    pipe.setLength(50000.0);
    pipe.setDiameter(0.5);
    pipe.setPipeWallRoughness(1.0e-5);
    pipe.setElevation(0.0);
    pipe.setNumberOfIncrements(5);
    pipe.run();

    System.out.println("beggs_and_brills_50km");
    System.out.println("inlet_m3_per_day=" + inlet.getFlowRate("m3/day"));
    System.out.println("inlet_mol_per_sec=" + inlet.getFlowRate("mol/sec"));
    print("inlet", inlet);
    print("outlet", pipe.getOutletStream());
    System.out.println("outlet_p_bara=" + pipe.getOutletStream().getPressure("bara"));
    System.out.println("pressure_drop_bar=" + pipe.getPressureDrop());
    System.out.println("pressure_drop_profile=" + pipe.getPressureDropProfile());
    System.out.println();
  }

  static Stream pipeFeed() {
    SystemInterface fluid = new SystemPrEos(278.15, 100.0);
    fluid.addComponent("methane", 95.0);
    fluid.addComponent("ethane", 5.0);
    fluid.setMixingRule(2);
    Stream stream = new Stream("inlet", fluid);
    stream.setPressure(100.0, "bara");
    stream.setTemperature(25.0, "C");
    stream.setFlowRate(10.0, "MSm3/day");
    stream.run();
    return stream;
  }

  static void pipeAdiabaticRow(String label, boolean fittings, Double outletBara) {
    Stream inlet = pipeFeed();
    neqsim.process.equipment.pipeline.AdiabaticPipe pipe =
        new neqsim.process.equipment.pipeline.AdiabaticPipe(label, inlet);
    pipe.setLength(50000.0);
    pipe.setDiameter(0.5);
    pipe.setPipeWallRoughness(1.0e-5);
    if (fittings) {
      pipe.addStandardFittings("elbow_90_standard", 2);
    }
    if (outletBara != null) {
      pipe.setOutletPressure(outletBara);
    }
    pipe.run();

    System.out.println(label);
    System.out.println("inlet_m3_per_day=" + inlet.getFlowRate("m3/day"));
    System.out.println("inlet_kg_per_hr=" + inlet.getFlowRate("kg/hr"));
    System.out.println("inlet_mol_per_sec=" + inlet.getFlowRate("mol/sec"));
    System.out.println("inlet_density_kg_per_m3="
        + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("inlet_viscosity_pa_s="
        + inlet.getThermoSystem().getViscosity("kg/msec"));
    System.out.println("effective_length_m=" + pipe.getEffectiveLength());
    System.out.println("equivalent_length_m=" + pipe.getEquivalentLength());
    System.out.println("total_fittings_ld=" + pipe.getTotalFittingsLdRatio());
    System.out.println("outlet_p_bara=" + pipe.getOutletStream().getPressure("bara"));
    System.out.println("pressure_drop_bar=" + pipe.getPressureDrop());
    System.out.println("velocity_m_per_s=" + pipe.getVelocity());
    System.out.println("reynolds=" + pipe.getReynoldsNumber());
    System.out.println("friction_factor=" + pipe.getFrictionFactor());
    System.out.println("flow_kg_per_hr_after=" + inlet.getFlowRate("kg/hr"));
    print("inlet", inlet);
    print("outlet", pipe.getOutletStream());
    System.out.println();
  }

  /// `ExpanderChartKhader`, the 2-D IGV/U-C map, on the enhancement test's table.
  static void expanderChartRows() {
    double[] igv = {0.6, 1.0};
    double[][] uc = {{0.5, 0.6, 0.7, 0.8, 0.9}, {0.5, 0.6, 0.7, 0.8, 0.9}};
    double[][] eta = {{0.72, 0.80, 0.84, 0.82, 0.76}, {0.78, 0.85, 0.88, 0.86, 0.80}};
    double[][] head = {{30.0, 32.0, 33.0, 32.5, 31.0}, {34.0, 36.0, 37.0, 36.5, 35.0}};

    neqsim.process.equipment.expander.ExpanderChartKhader chart =
        new neqsim.process.equipment.expander.ExpanderChartKhader();
    chart.setCurves(igv, uc, eta, head);

    System.out.println("chart_no_reference");
    System.out.println("map_defined=" + chart.isMapDefined());
    System.out.println("efficiency_07_10=" + chart.getEfficiency(0.7, 1.0));
    System.out.println("efficiency_07_06=" + chart.getEfficiency(0.7, 0.6));
    System.out.println("efficiency_07_08=" + chart.getEfficiency(0.7, 0.8));
    System.out.println("efficiency_04_10=" + chart.getEfficiency(0.4, 1.0));
    System.out.println("efficiency_095_10=" + chart.getEfficiency(0.95, 1.0));
    System.out.println("optimum_uc_10=" + chart.getOptimumVelocityRatio(1.0));
    System.out.println("optimum_uc_06=" + chart.getOptimumVelocityRatio(0.6));
    System.out.println("head_07_10=" + chart.getStageHeadDrop(0.7, 1.0, null));
    System.out.println("head_07_06=" + chart.getStageHeadDrop(0.7, 0.6, null));
    System.out.println("head_07_08=" + chart.getStageHeadDrop(0.7, 0.8, null));
    System.out.println("head_04_10=" + chart.getStageHeadDrop(0.4, 1.0, null));
    System.out.println("head_095_10=" + chart.getStageHeadDrop(0.95, 1.0, null));
    System.out.println();

    SystemInterface reference = new neqsim.thermo.system.SystemSrkEos(298.15, 50.0);
    reference.addComponent("methane", 0.9);
    reference.addComponent("ethane", 0.1);
    reference.setMixingRule(2);
    reference.init(0);
    reference.initThermoProperties();
    neqsim.process.equipment.expander.ExpanderChartKhader withReference =
        new neqsim.process.equipment.expander.ExpanderChartKhader(reference, 0.424);
    withReference.setCurves(igv, uc, eta, head);

    SystemInterface co2 = new neqsim.thermo.system.SystemSrkEos(283.15, 50.0);
    co2.addComponent("CO2", 0.9);
    co2.addComponent("methane", 0.1);
    co2.setMixingRule(2);
    co2.init(0);
    co2.initThermoProperties();

    System.out.println("chart_with_reference");
    System.out.println("reference_sound_speed=" + withReference.getReferenceSoundSpeed());
    System.out.println("reference_system_sound_speed=" + reference.getSoundSpeed());
    System.out.println("reference_phase_sound_speed=" + reference.getPhase(0).getSoundSpeed());
    System.out.println("process_sound_speed=" + co2.getSoundSpeed());
    System.out.println("process_phase_sound_speed=" + co2.getPhase(0).getSoundSpeed());
    System.out.println("head_ref_fluid=" + withReference.getStageHeadDrop(0.7, 1.0, reference));
    System.out.println("head_co2_fluid=" + withReference.getStageHeadDrop(0.7, 1.0, co2));
    System.out.println();
  }

  /// `TurboExpanderCompressor`, the coupled expander-compressor.
  /// `TurboExpanderCompressorMechanicalDesign.calcDesign`, on the coupled machine's states.
  ///
  /// The design reads **both** ends: the expander's inlet and outlet fix its wheel and
  /// casing, and the compressor's suction and discharge fix its impeller and casing. The
  /// third state raises the expander inlet pressure, which moves the expander casing and the
  /// coupled power balance. **The compressor's barrel branch is not reached on these states**:
  /// its discharge is set by the machine's own chart, not by the expander's power, so the
  /// design's `BARREL` arm is stated rather than cased.
  static void turboExpanderMechanicalRows() {
    turboExpanderMechanicalRow("tec_run", 456000.0, 60.95, -23.0, 423448.0, 42.0, 17.0);
    turboExpanderMechanicalRow("tec_low_expander_pressure", 426000.0, 55.95, -23.0, 423448.0,
        42.0, 17.0);
    turboExpanderMechanicalRow("tec_high_expander_pressure", 456000.0, 120.0, -23.0, 423448.0, 42.0, 17.0);
  }

  static void turboExpanderMechanicalRow(String label, double expanderKgHr,
      double expanderBara, double expanderC, double compressorKgHr, double compressorBara,
      double compressorC) {
    Stream expanderFeed = turboFeed(expanderKgHr, expanderBara, expanderC);
    Stream compressorFeed = turboFeed(compressorKgHr, compressorBara, compressorC);
    neqsim.process.equipment.expander.TurboExpanderCompressor machine =
        turboMachine(expanderFeed, compressorFeed);
    machine.run();
    machine.initTECMechanicalDesign();
    neqsim.process.mechanicaldesign.expander.TurboExpanderCompressorMechanicalDesign design =
        machine.getTECMechanicalDesign();
    design.calcDesign();
    System.out.println(label);
    System.out.println("speed_rpm=" + machine.getSpeed());
    System.out.println("gear_ratio=" + machine.getGearRatio());
    System.out.println("expander_power_kW=" + machine.getPowerExpander("kW"));
    System.out.println("compressor_power_kW=" + machine.getPowerCompressor("kW"));
    System.out.println("compressor_polytropic_head_kj_per_kg="
        + machine.getCompressorPolytropicHead());
    System.out.println("expander_wheel_diameter_mm=" + design.getExpanderWheelDiameter());
    System.out.println("expander_tip_speed_m_per_s=" + design.getExpanderTipSpeed());
    System.out.println("expander_type=" + design.getExpanderType());
    System.out.println("expander_stages=" + design.getExpanderStages());
    System.out.println("expander_design_inlet_p_bara=" + design.getExpanderDesignInletPressure());
    System.out.println("expander_design_outlet_p_bara=" + design.getExpanderDesignOutletPressure());
    System.out.println("expander_casing_p_bara=" + design.getExpanderCasingDesignPressure());
    System.out.println("expander_casing_t_C=" + design.getExpanderCasingDesignTemperature());
    System.out.println("expander_power_design_kW=" + design.getExpanderPowerKW());
    System.out.println("expander_enthalpy_drop_kj_per_kg=" + design.getExpanderEnthalpyDropKJkg());
    System.out.println("compressor_impeller_diameter_mm="
        + design.getCompressorImpellerDiameter());
    System.out.println("compressor_tip_speed_m_per_s=" + design.getCompressorTipSpeed());
    System.out.println("compressor_stages=" + design.getCompressorStages());
    System.out.println("compressor_polytropic_head_design=" + design.getCompressorPolytropicHead());
    System.out.println("compressor_head_per_stage_kj_per_kg=" + design.getCompressorHeadPerStage());
    System.out.println("compressor_casing_type=" + design.getCompressorCasingType());
    System.out.println("compressor_design_suction_p_bara="
        + design.getCompressorDesignSuctionPressure());
    System.out.println("compressor_design_discharge_p_bara="
        + design.getCompressorDesignDischargePressure());
    System.out.println("compressor_casing_p_bara=" + design.getCompressorCasingDesignPressure());
    System.out.println("compressor_casing_t_C=" + design.getCompressorCasingDesignTemperature());
    System.out.println("compressor_power_design_kW=" + design.getCompressorPowerKW());
    System.out.println("shaft_diameter_mm=" + design.getShaftDiameter());
    System.out.println("shaft_weight_kg=" + design.getShaftWeight());
    System.out.println("operating_speed_rpm=" + design.getOperatingSpeed());
    System.out.println("max_continuous_speed_rpm=" + design.getMaxContinuousSpeed());
    System.out.println("trip_speed_rpm=" + design.getTripSpeed());
    System.out.println("first_critical_speed_rpm=" + design.getFirstCriticalSpeed());
    System.out.println("second_critical_speed_rpm=" + design.getSecondCriticalSpeed());
    System.out.println("bearing_span_mm=" + design.getBearingSpan());
    System.out.println("bearing_type=" + design.getBearingType());
    System.out.println("expander_seal_type=" + design.getExpanderSealType());
    System.out.println("compressor_seal_type=" + design.getCompressorSealType());
    System.out.println("design_torque_nm=" + design.getDesignTorqueNm());
    System.out.println("design_shaft_stress_mpa=" + design.getDesignShaftStressMPa());
    System.out.println("design_surge_flow_m3_per_hr=" + design.getDesignSurgeFlowM3hr());
    System.out.println("design_stonewall_flow_m3_per_hr=" + design.getDesignStonewallFlowM3hr());
    System.out.println("expander_nozzle_throat_area_m2=" + design.getExpanderNozzleThroatArea());
    System.out.println("expander_wheel_weight_kg=" + design.getExpanderWheelWeight());
    System.out.println("expander_casing_weight_kg=" + design.getExpanderCasingWeight());
    System.out.println("compressor_impeller_weight_kg=" + design.getCompressorImpellerWeight());
    System.out.println("compressor_casing_weight_kg=" + design.getCompressorCasingWeight());
    System.out.println("weight_vessel_kg=" + design.getWeightVessel());
    System.out.println("weight_internals_kg=" + design.getWeigthInternals());
    System.out.println("weight_piping_kg=" + design.getWeightPiping());
    System.out.println("weight_electro_instrument_kg=" + design.getWeightElectroInstrument());
    System.out.println("weight_structural_steel_kg=" + design.getWeightStructualSteel());
    System.out.println("weight_total_kg=" + design.getWeightTotal());
    System.out.println("module_length_m=" + design.getModuleLength());
    System.out.println("module_width_m=" + design.getModuleWidth());
    System.out.println("module_height_m=" + design.getModuleHeight());
    System.out.println("expander_inlet_m3_per_hr=" + expanderFeed.getFlowRate("m3/hr"));
    System.out.println("compressor_inlet_m3_per_hr=" + compressorFeed.getFlowRate("m3/hr"));
    print("expander_inlet", expanderFeed);
    print("compressor_inlet", compressorFeed);
    System.out.println();
  }

  static void turboExpanderRows() {
    turboExpanderRow("turbo_expander_run", 456000.0, 60.95, -23.0, 423448.0, 42.0, 17.0);
    turboExpanderRow("turbo_expander_low_expander_pressure", 426000.0, 55.95, -23.0, 423448.0,
        42.0, 17.0);
    turboExpanderMassRow("turbo_expander_mass_balance");
    turboExpanderInverseRow("turbo_expander_outlet_temperature");
    turboExpanderChartRow("turbo_expander_igv_penalty");
  }

  static Stream turboFeed(double kgPerHour, double bara, double celsius) {
    SystemInterface fluid = new SystemPrEos(273.15 + 42.0, 10.0);
    fluid.addComponent("nitrogen", 0.006);
    fluid.addComponent("CO2", 0.014);
    fluid.addComponent("methane", 0.862);
    fluid.addComponent("ethane", 0.08);
    fluid.addComponent("propane", 0.03);
    fluid.addComponent("i-butane", 0.0024);
    fluid.addComponent("n-butane", 0.004);
    fluid.addComponent("n-hexane", 0.0015);
    fluid.setMixingRule(2);
    fluid.init(0);
    Stream stream = new Stream("feed", fluid);
    stream.setFlowRate(kgPerHour, "kg/hr");
    stream.setTemperature(celsius, "C");
    stream.setPressure(bara, "bara");
    stream.run();
    return stream;
  }

  static neqsim.process.equipment.expander.TurboExpanderCompressor turboMachine(
      Stream expanderFeed, Stream compressorFeed) {
    neqsim.process.equipment.expander.TurboExpanderCompressor machine =
        new neqsim.process.equipment.expander.TurboExpanderCompressor("TurboExpander",
            expanderFeed);
    machine.setCompressorFeedStream(compressorFeed);
    machine.setUCcurve(
        new double[] {0.9964751359624449, 0.7590835113213541, 0.984295619176559,
            0.8827799803397821, 0.9552460269880922, 1.0},
        new double[] {0.984090909090909, 0.796590909090909, 0.9931818181818183,
            0.9363636363636364, 0.9943181818181818, 1.0});
    machine.setQNEfficiencycurve(new double[] {0.5, 0.7, 0.85, 1.0, 1.2, 1.4, 1.6},
        new double[] {0.88, 0.91, 0.95, 1.0, 0.97, 0.85, 0.6});
    machine.setQNHeadcurve(new double[] {0.5, 0.8, 1.0, 1.2, 1.4, 1.6},
        new double[] {1.1, 1.05, 1.0, 0.9, 0.7, 0.4});
    machine.setImpellerDiameter(0.424);
    machine.setDesignSpeed(6850.0);
    machine.setExpanderDesignIsentropicEfficiency(0.88);
    machine.setDesignUC(0.7);
    machine.setDesignQn(0.03328);
    machine.setExpanderOutPressure(42.0);
    machine.setCompressorDesignPolytropicEfficiency(0.81);
    machine.setCompressorDesignPolytropicHead(20.47);
    machine.setMaximumIGVArea(1.637e4);
    return machine;
  }

  static void printTurboRow(String label,
      neqsim.process.equipment.expander.TurboExpanderCompressor machine, Stream expanderFeed,
      Stream compressorFeed) {
    System.out.println(label);
    System.out.println("expander_inlet_m3_per_hr=" + expanderFeed.getFlowRate("m3/hr"));
    System.out.println("expander_inlet_kg_per_hr=" + expanderFeed.getFlowRate("kg/hr"));
    System.out.println("compressor_inlet_m3_per_hr=" + compressorFeed.getFlowRate("m3/hr"));
    System.out.println("compressor_inlet_kg_per_hr=" + compressorFeed.getFlowRate("kg/hr"));
    System.out.println("speed_rpm=" + machine.getSpeed());
    System.out.println("power_expander_w=" + machine.getPowerExpander());
    System.out.println("power_compressor_w=" + machine.getPowerCompressor());
    System.out.println("expander_isentropic_efficiency="
        + machine.getExpanderIsentropicEfficiency());
    System.out.println("compressor_polytropic_efficiency="
        + machine.getCompressorPolytropicEfficiency());
    System.out.println("compressor_polytropic_head_kj_per_kg="
        + machine.getCompressorPolytropicHead());
    System.out.println("uc_ratio_expander=" + machine.getUCratioexpander());
    System.out.println("qn_ratio_compressor=" + machine.getQNratiocompressor());
    System.out.println("igv_opening=" + machine.getIGVopening());
    System.out.println("expander_outlet_p_bara="
        + machine.getExpanderOutletStream().getPressure("bara"));
    System.out.println("expander_outlet_t_k="
        + machine.getExpanderOutletStream().getTemperature("K"));
    System.out.println("compressor_outlet_p_bara="
        + machine.getCompressorOutletStream().getPressure("bara"));
    System.out.println("compressor_outlet_t_k="
        + machine.getCompressorOutletStream().getTemperature("K"));
    print("expander_inlet", expanderFeed);
    print("expander_outlet", machine.getExpanderOutletStream());
    print("compressor_inlet", compressorFeed);
    print("compressor_outlet", machine.getCompressorOutletStream());
    System.out.println();
  }

  static void turboExpanderRow(String label, double expanderKgHr, double expanderBara,
      double expanderC, double compressorKgHr, double compressorBara, double compressorC) {
    Stream expanderFeed = turboFeed(expanderKgHr, expanderBara, expanderC);
    Stream compressorFeed = turboFeed(compressorKgHr, compressorBara, compressorC);
    neqsim.process.equipment.expander.TurboExpanderCompressor machine =
        turboMachine(expanderFeed, compressorFeed);
    machine.run();
    printTurboRow(label, machine, expanderFeed, compressorFeed);
  }

  static void turboExpanderMassRow(String label) {
    SystemInterface fluid = new SystemPrEos(293.15, 50.0);
    fluid.addComponent("methane", 0.95);
    fluid.addComponent("ethane", 0.05);
    fluid.setMixingRule(2);
    fluid.init(0);
    Stream expanderFeed = new Stream("expander feed", fluid);
    expanderFeed.setFlowRate(10000.0, "kg/hr");
    expanderFeed.setTemperature(20.0, "C");
    expanderFeed.setPressure(50.0, "bara");
    expanderFeed.run();
    Stream compressorFeed = new Stream("compressor feed", fluid.clone());
    compressorFeed.setFlowRate(7500.0, "kg/hr");
    compressorFeed.setTemperature(5.0, "C");
    compressorFeed.setPressure(30.0, "bara");
    compressorFeed.run();

    neqsim.process.equipment.expander.TurboExpanderCompressor machine =
        turboMachine(expanderFeed, compressorFeed);
    machine.setExpanderOutPressure(25.0);
    machine.setCompressorDesignPolytropicEfficiency(0.8);
    machine.setCompressorDesignPolytropicHead(10.0);
    machine.run();
    System.out.println("mass_balance_expander_kg_per_hr=" + machine.getMassBalance("kg/hr"));
    printTurboRow(label, machine, expanderFeed, compressorFeed);
  }

  static void turboExpanderInverseRow(String label) {
    Stream expanderFeed = turboFeed(456000.0, 60.95, -23.0);
    Stream compressorFeed = turboFeed(423448.0, 42.0, 17.0);
    neqsim.process.equipment.expander.TurboExpanderCompressor forward =
        turboMachine(expanderFeed, compressorFeed);
    forward.run();
    double targetC = forward.getExpanderOutletStream().getTemperature("C");

    Stream expanderFeed2 = turboFeed(456000.0, 60.95, -23.0);
    Stream compressorFeed2 = turboFeed(423448.0, 42.0, 17.0);
    neqsim.process.equipment.expander.TurboExpanderCompressor inverse =
        turboMachine(expanderFeed2, compressorFeed2);
    inverse.setExpanderOutTemperature(targetC, "C");
    inverse.run();

    System.out.println(label);
    System.out.println("target_outlet_t_c=" + targetC);
    System.out.println("forward_design_efficiency="
        + forward.getExpanderDesignIsentropicEfficiency());
    System.out.println("inverse_design_efficiency="
        + inverse.getExpanderDesignIsentropicEfficiency());
    System.out.println("inverse_isentropic_efficiency="
        + inverse.getExpanderIsentropicEfficiency());
    printTurboRow(label + "_state", inverse, expanderFeed2, compressorFeed2);
  }

  static void turboExpanderChartRow(String label) {
    Stream expanderFeed = turboFeed(456000.0, 60.95, -23.0);
    Stream compressorFeed = turboFeed(423448.0, 42.0, 17.0);
    neqsim.process.equipment.expander.TurboExpanderCompressor machine =
        turboMachine(expanderFeed, compressorFeed);

    double[] igv = {0.6, 1.0};
    double[][] uc = {{0.5, 0.6, 0.7, 0.8, 0.9}, {0.5, 0.6, 0.7, 0.8, 0.9}};
    double[][] eta = {{0.72, 0.80, 0.84, 0.82, 0.76}, {0.78, 0.85, 0.88, 0.86, 0.80}};
    double[][] head = {{30.0, 32.0, 33.0, 32.5, 31.0}, {34.0, 36.0, 37.0, 36.5, 35.0}};
    neqsim.process.equipment.expander.ExpanderChartKhader chart =
        new neqsim.process.equipment.expander.ExpanderChartKhader();
    chart.setCurves(igv, uc, eta, head);
    machine.setExpanderChart(chart);
    machine.setIgvControlMode(true);
    machine.setIGVopening(1.0);
    machine.run();
    double noPenalty = machine.getExpanderIsentropicEfficiency();

    machine.setIgvEfficiencyPenaltyCurve(new double[] {0.4, 1.0}, new double[] {0.85, 0.90});
    machine.setIGVopening(1.0);
    machine.run();

    System.out.println(label);
    System.out.println("efficiency_no_penalty=" + noPenalty);
    System.out.println("efficiency_with_penalty=" + machine.getExpanderIsentropicEfficiency());
    System.out.println("igv_opening=" + machine.getIGVopening());
    System.out.println("igv_penalty_at_10=" + machine.getIgvEfficiencyPenalty(1.0));
    printTurboRow(label + "_state", machine, expanderFeed, compressorFeed);
  }

  /// `MapTurboExpanderCompressor`, the map-based single shaft, re-cased on PR.
  ///
  /// The three states of `MapTurboExpanderCompressorTest` - the balanced-speed mode, the
  /// specified-pressure mode (whose compressor is isentropic, not the map) and the turndown
  /// row whose expander drop is too small to drive the brake compressor. Each prints the
  /// speeds, both powers, the residual and the operating status.
  static void mapTurboExpanderRows() {
    mapTurboExpanderRow("map_turbo_expander_balanced", 25.0, 60.0,
        neqsim.process.equipment.expander.MapTurboExpanderCompressor.ShaftMode.BALANCED_SPEED,
        true, 58.0);
    mapTurboExpanderRow("map_turbo_expander_specified", 25.0, 60.0,
        neqsim.process.equipment.expander.MapTurboExpanderCompressor.ShaftMode.SPECIFIED_PRESSURES,
        false, 58.0);
    mapTurboExpanderRow("map_turbo_expander_turndown", 25.0, 31.0,
        neqsim.process.equipment.expander.MapTurboExpanderCompressor.ShaftMode.BALANCED_SPEED,
        true, 58.0);
  }

  static SystemInterface mapTurboGas(double celsius, double bara) {
    SystemInterface fluid = new SystemPrEos(273.15 + celsius, bara);
    fluid.addComponent("nitrogen", 1.205);
    fluid.addComponent("CO2", 1.340);
    fluid.addComponent("methane", 87.974);
    fluid.addComponent("ethane", 5.258);
    fluid.addComponent("propane", 3.283);
    fluid.addComponent("i-butane", 0.082);
    fluid.addComponent("n-butane", 0.487);
    fluid.setMixingRule(2);
    fluid.setTemperature(celsius, "C");
    fluid.setPressure(bara, "bara");
    return fluid;
  }

  static void mapTurboExpanderRow(String label, double expanderCelsius, double expanderBara,
      neqsim.process.equipment.expander.MapTurboExpanderCompressor.ShaftMode mode,
      boolean attachChart, double compressorOutletBara) {
    Stream expanderFeed = new Stream("expFeed", mapTurboGas(expanderCelsius, expanderBara));
    expanderFeed.setFlowRate(3.0, "MSm3/day");
    expanderFeed.run();
    Stream compressorFeed = new Stream("compFeed", mapTurboGas(20.0, 48.0));
    compressorFeed.setFlowRate(3.0, "MSm3/day");
    compressorFeed.run();

    neqsim.process.equipment.expander.MapTurboExpanderCompressor machine =
        new neqsim.process.equipment.expander.MapTurboExpanderCompressor("TEC", expanderFeed,
            compressorFeed);
    machine.setExpanderOutletPressure(30.0);
    machine.setExpanderIsentropicEfficiency(0.85);
    machine.setMechanicalEfficiency(0.98);
    if (attachChart) {
      machine.getCompressor().setCompressorChart(
          mapTurboChart(compressorFeed));
      machine.getCompressor().getCompressorChart().setHeadUnit("kJ/kg");
      machine.getCompressor().setUsePolytropicCalc(true);
    } else {
      machine.setCompressorOutletPressure(compressorOutletBara);
      machine.getCompressor().setIsentropicEfficiency(0.78);
    }
    machine.setShaftMode(mode);
    machine.run(java.util.UUID.randomUUID());

    System.out.println(label);
    System.out.println("shaft_mode=" + machine.getShaftMode());
    System.out.println("shaft_speed_rpm=" + machine.getShaftSpeed());
    System.out.println("available_shaft_power_w=" + machine.getAvailableShaftPower());
    System.out.println("consumed_compressor_power_w=" + machine.getConsumedCompressorPower());
    System.out.println("power_balance_residual_w=" + machine.getPowerBalanceResidual());
    System.out.println("speed_limited=" + machine.isSpeedLimited());
    System.out.println("operating_status=" + machine.getOperatingStatus());
    System.out.println("feasible=" + machine.isFeasible());
    System.out.println("expander_outlet_p_bara="
        + machine.getExpanderOutStream().getPressure("bara"));
    System.out.println("expander_outlet_t_k="
        + machine.getExpanderOutStream().getTemperature("K"));
    System.out.println("compressor_outlet_p_bara="
        + machine.getCompressorOutStream().getPressure("bara"));
    System.out.println("compressor_outlet_t_k="
        + machine.getCompressorOutStream().getTemperature("K"));
    System.out.println("compressor_flow_m3_per_hr="
        + compressorFeed.getFlowRate("m3/hr"));
    print("expander_inlet", expanderFeed);
    print("expander_outlet", machine.getExpanderOutStream());
    print("compressor_inlet", compressorFeed);
    print("compressor_outlet", machine.getCompressorOutStream());
    System.out.println();
  }

  static neqsim.process.equipment.compressor.CompressorChartKhader2015 mapTurboChart(
      Stream compressorFeed) {
    double[] chartConditions = new double[] {25.0, 50.0, 50.0, 20.0};
    double[] speed = new double[] {12913, 12298, 11683, 11098, 10453, 9224, 8609, 8200};
    double[][] flow = new double[][] {
        {2789.1285, 3174.0375, 3689.2288, 4179.4503, 4570.2768, 4954.7728, 5246.0329, 5661.0331},
        {2571.1753, 2943.7254, 3440.2675, 3837.4448, 4253.0898, 4668.6643, 4997.1926, 5387.4952},
        {2415.3793, 2763.0706, 3141.7095, 3594.7436, 4047.6467, 4494.1889, 4853.7353, 5138.7858},
        {2247.2043, 2799.7342, 3178.3428, 3656.1551, 4102.778, 4394.1591, 4648.3224, 4840.4998},
        {2072.8397, 2463.9483, 2836.4078, 3202.5266, 3599.6333, 3978.0203, 4257.0022, 4517.345},
        {1835.9552, 2208.455, 2618.1322, 2940.8034, 3244.7852, 3530.1279, 3753.3738, 3895.9746},
        {1711.3386, 1965.8848, 2356.9431, 2685.9247, 3008.5154, 3337.2855, 3591.5092},
        {1636.5807, 2002.8708, 2338.0319, 2642.1245, 2896.4894, 3113.6264, 3274.8764, 3411.2977}};
    double[][] head = new double[][] {{80.0375, 78.8934, 76.2142, 71.8678, 67.0062, 60.6061, 53.0499, 39.728},
        {72.2122, 71.8369, 68.9009, 65.8341, 60.7167, 54.702, 47.2749, 35.7471},
        {65.1576, 64.5253, 62.6118, 59.1619, 54.0455, 47.0059, 39.195, 31.6387},
        {58.6154, 56.9627, 54.6647, 50.4462, 44.4322, 38.4144, 32.9084, 28.8109},
        {52.3295, 51.0573, 49.5283, 46.3326, 42.3685, 37.2502, 31.4884, 25.598},
        {40.6578, 39.6416, 37.6008, 34.6603, 30.9503, 27.1116, 23.2713, 20.4546},
        {35.2705, 34.6359, 32.7228, 31.0645, 27.0985, 22.7482, 18.0113},
        {32.192, 31.1756, 29.1329, 26.833, 23.8909, 21.3324, 18.7726, 16.3403}};
    double[][] polyEff = new double[][] {{77.245, 79.415, 80.738, 80.523, 79.221, 75.472, 69.603, 58.732},
        {77.011, 79.307, 80.894, 80.719, 79.531, 75.591, 69.685, 60.004},
        {77.004, 79.169, 80.804, 80.654, 78.853, 73.666, 66.274, 57.672},
        {77.072, 80.463, 81.139, 79.637, 75.381, 69.533, 63.800, 58.812},
        {76.971, 79.834, 80.947, 80.581, 78.046, 73.040, 66.557, 59.862},
        {77.506, 80.206, 81.034, 79.609, 76.381, 70.803, 64.644, 60.530},
        {77.818, 80.065, 81.063, 79.896, 76.198, 69.290, 60.857},
        {78.092, 80.935, 80.790, 78.864, 75.217, 70.311, 65.551, 61.039}};
    neqsim.process.equipment.compressor.CompressorChartKhader2015 chart =
        new neqsim.process.equipment.compressor.CompressorChartKhader2015(compressorFeed, 0.9);
    chart.setCurves(chartConditions, speed, flow, head, flow, polyEff);
    return chart;
  }

  /// `CompressorChartKhader2015`, the Mach-number normalised map, re-cased on PR.
  ///
  /// **`CompressorChartKhader2015Test`'s two states and its numbers**, with the sound speeds the
  /// class uses printed beside them: the reference fluid's *phase* value (which fits the map)
  /// and the process fluid's *system* value (which reads it). The pinned test asserts
  /// `41.56192413` kJ/kg on the first state and `75.11224727` bara on its machine run.
  static void compressorKhaderRows() {
    double[] chartConditions = new double[] {25.0, 50.0, 50.0, 20.0};
    double[] speed = new double[] {12913, 12298, 11683, 11098, 10453, 9224, 8609, 8200};
    double[][] flow = new double[][] {
        {2789.1285, 3174.0375, 3689.2288, 4179.4503, 4570.2768, 4954.7728, 5246.0329, 5661.0331},
        {2571.1753, 2943.7254, 3440.2675, 3837.4448, 4253.0898, 4668.6643, 4997.1926, 5387.4952},
        {2415.3793, 2763.0706, 3141.7095, 3594.7436, 4047.6467, 4494.1889, 4853.7353, 5138.7858},
        {2247.2043, 2799.7342, 3178.3428, 3656.1551, 4102.778, 4394.1591, 4648.3224, 4840.4998},
        {2072.8397, 2463.9483, 2836.4078, 3202.5266, 3599.6333, 3978.0203, 4257.0022, 4517.345},
        {1835.9552, 2208.455, 2618.1322, 2940.8034, 3244.7852, 3530.1279, 3753.3738, 3895.9746},
        {1711.3386, 1965.8848, 2356.9431, 2685.9247, 3008.5154, 3337.2855, 3591.5092},
        {1636.5807, 2002.8708, 2338.0319, 2642.1245, 2896.4894, 3113.6264, 3274.8764, 3411.2977}};
    double[][] head = new double[][] {
        {80.0375, 78.8934, 76.2142, 71.8678, 67.0062, 60.6061, 53.0499, 39.728},
        {72.2122, 71.8369, 68.9009, 65.8341, 60.7167, 54.702, 47.2749, 35.7471},
        {65.1576, 64.5253, 62.6118, 59.1619, 54.0455, 47.0059, 39.195, 31.6387},
        {58.6154, 56.9627, 54.6647, 50.4462, 44.4322, 38.4144, 32.9084, 28.8109},
        {52.3295, 51.0573, 49.5283, 46.3326, 42.3685, 37.2502, 31.4884, 25.598},
        {40.6578, 39.6416, 37.6008, 34.6603, 30.9503, 27.1116, 23.2713, 20.4546},
        {35.2705, 34.6359, 32.7228, 31.0645, 27.0985, 22.7482, 18.0113},
        {32.192, 31.1756, 29.1329, 26.833, 23.8909, 21.3324, 18.7726, 16.3403}};
    double[][] polyEff = new double[][] {
        {77.2452238409573, 79.4154186459363, 80.737960012489, 80.5229826589649, 79.2210931638144,
            75.4719133864634, 69.6034181197298, 58.7322388482707},
        {77.0107837113504, 79.3069974136389, 80.8941189021135, 80.7190194665918, 79.5313242980328,
            75.5912622896367, 69.6846136362097, 60.0043057990909},
        {77.0043065299874, 79.1690958847856, 80.8038169975675, 80.6543975614197, 78.8532389102705,
            73.6664774270613, 66.2735600426727, 57.671664571658},
        {77.0716623789093, 80.4629750233093, 81.1390811169072, 79.6374242667478, 75.380928428817,
            69.5332969549779, 63.7997587622339, 58.8120614497758},
        {76.9705872525642, 79.8335492585324, 80.9468133671171, 80.5806471927835, 78.0462158225426,
            73.0403707523258, 66.5572286338589, 59.8624822515064},
        {77.5063036680357, 80.2056198362559, 81.0339108025933, 79.6085962687939, 76.3814534404405,
            70.8027503005902, 64.6437367160571, 60.5299349982342},
        {77.8175271586685, 80.065165942218, 81.0631362122632, 79.8955051771299, 76.1983240929369,
            69.289982774309, 60.8567149372229},
        {78.0924334304045, 80.9353551568667, 80.7904437766234, 78.8639325223295, 75.2170936751143,
            70.3105081673411, 65.5507568533569, 61.0391468300337}};

    // The nine-component gas, 24 C and 48 bara, at 3 MSm3/day and speed 10000.
    Stream stream = khaderFeed(new SystemPrEos(298.15, 48.0),
        new String[] {"nitrogen", "CO2", "methane", "ethane", "propane", "i-butane", "n-butane",
            "i-pentane", "n-pentane"},
        new double[] {1.205, 1.340, 87.974, 5.258, 3.283, 0.082, 0.487, 0.056, 0.053},
        24.0, 48.0);
    khaderRow("khader_nine_component", stream, chartConditions, speed, flow, head, polyEff, 10000.0);

    // The three-component reference-like gas, 25 C and 50 bara, at 3 MSm3/day.
    Stream reference = khaderFeed(new SystemPrEos(298.15, 50.0),
        new String[] {"methane", "ethane", "propane"},
        new double[] {8.35736E-1, 4.64298E-2, 1.17835E-1}, 25.0, 50.0);
    khaderRow("khader_three_component", reference, chartConditions, speed, flow, head, polyEff,
        10000.0);
  }

  static Stream khaderFeed(SystemInterface fluid, String[] names, double[] z, double celsius,
      double bara) {
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    fluid.setTemperature(celsius, "C");
    fluid.setPressure(bara, "bara");
    fluid.setTotalFlowRate(3.0, "MSm3/day");
    Stream stream = new Stream("stream_1", fluid);
    stream.run();
    return stream;
  }

  static void khaderRow(String label, Stream stream, double[] chartConditions, double[] speed,
      double[][] flow, double[][] head, double[][] polyEff, double compSpeed) {
    neqsim.process.equipment.compressor.Compressor compressor =
        new neqsim.process.equipment.compressor.Compressor("comp1", stream);
    compressor.setSpeed(compSpeed);
    neqsim.process.equipment.compressor.CompressorChartKhader2015 chart =
        new neqsim.process.equipment.compressor.CompressorChartKhader2015(stream, 0.9);
    chart.setCurves(chartConditions, speed, flow, head, flow, polyEff);
    compressor.setCompressorChart(chart);
    compressor.getCompressorChart().setHeadUnit("kJ/kg");

    double inlet = stream.getFlowRate("m3/hr");
    System.out.println(label);
    System.out.println("inlet_m3_per_hr=" + inlet);
    System.out.println("compressor_speed_rpm=" + compSpeed);
    System.out.println("process_sound_speed=" + stream.getFluid().getSoundSpeed());
    System.out.println("process_phase_sound_speed=" + stream.getFluid().getPhase(0).getSoundSpeed());
    System.out.println("reference_sound_speed="
        + chart.getReferenceFluid().getPhase(0).getSoundSpeed());
    System.out.println("reference_temperature_k=" + chart.getReferenceFluid().getTemperature());
    System.out.println("reference_pressure_bara=" + chart.getReferenceFluid().getPressure("bara"));
    System.out.println("reference_molar_mass=" + chart.getReferenceFluid().getMolarMass() * 1000.0);
    System.out.println("head_kj_per_kg=" + chart.getPolytropicHead(inlet, compSpeed));
    System.out.println("efficiency_percent=" + chart.getPolytropicEfficiency(inlet, compSpeed));
    System.out.println("corrected_flow_factor_first="
        + chart.getCorrectedCurves(chartConditions, speed, flow, head, polyEff, polyEff).get(0)
            .correctedFlowFactor[0]);
    System.out.println("real_curve_flow_first=" + chart.getRealCurves().get(0).flow[0]);
    System.out.println("real_curve_head_first=" + chart.getRealCurves().get(0).head[0]);
    compressor.run();
    System.out.println("outlet_p_bara=" + compressor.getOutletStream().getPressure("bara"));
    System.out.println("outlet_t_k=" + compressor.getOutletStream().getTemperature("K"));
    System.out.println("polytropic_head_kj_per_kg=" + compressor.getPolytropicFluidHead());
    System.out.println("polytropic_efficiency_fraction=" + compressor.getPolytropicEfficiency());
    System.out.println("power_w=" + compressor.getPower("W"));
    print("inlet", stream);
    print("outlet", compressor.getOutletStream());
    System.out.println();
  }

  /// `PhaseEos.getSoundSpeed` and `SystemThermo.getSoundSpeed`, on flashed states.
  ///
  /// **The system value is a beta-weighted sum, and an unflashed system doubles it** - two
  /// phase objects with `beta = 1.0` each. Every row here flashes through a `Stream`, so the
  /// system value is the phase value on a single phase and the weighted mean on two. The
  /// phase's own `Cp`, `Cv`, `Z` and density are printed beside it so the port's cubic
  /// derivatives can be held to the same state.
  static void soundSpeedRows() {
    soundSpeedRow("gas_nine_component", new SystemPrEos(298.15, 48.0),
        new String[] {"nitrogen", "CO2", "methane", "ethane", "propane", "i-butane", "n-butane",
            "i-pentane", "n-pentane"},
        new double[] {1.205, 1.340, 87.974, 5.258, 3.283, 0.082, 0.487, 0.056, 0.053},
        24.0, 48.0, 1.0);
    soundSpeedRow("gas_nine_component_hot", new SystemPrEos(298.15, 48.0),
        new String[] {"nitrogen", "CO2", "methane", "ethane", "propane", "i-butane", "n-butane",
            "i-pentane", "n-pentane"},
        new double[] {1.205, 1.340, 87.974, 5.258, 3.283, 0.082, 0.487, 0.056, 0.053},
        -23.0, 60.95, 1.0);
    soundSpeedRow("two_phase_methane_butane", new SystemPrEos(298.15, 50.0),
        new String[] {"methane", "n-butane"}, new double[] {0.5, 0.5}, -23.15, 30.0, 0.0);
    soundSpeedRow("liquid_propane", new SystemPrEos(298.15, 50.0),
        new String[] {"propane"}, new double[] {1.0}, 20.0, 100.0, 0.0);
    soundSpeedRow("gas_methane_low_pressure", new SystemPrEos(298.15, 50.0),
        new String[] {"methane"}, new double[] {1.0}, 27.0, 1.0, 0.0);
  }

  static void soundSpeedRow(String label, SystemInterface fluid, String[] names, double[] z,
      double celsius, double bara, double msm3PerDay) {
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setTemperature(celsius, "C");
    fluid.setPressure(bara, "bara");
    if (msm3PerDay > 0.0) {
      fluid.setTotalFlowRate(msm3PerDay, "MSm3/day");
    }
    Stream stream = new Stream(label, fluid);
    stream.run();
    SystemInterface system = stream.getThermoSystem();
    system.initProperties();

    System.out.println(label);
    System.out.println("phase=" + system.getPhase(0).getType());
    System.out.println("number_of_phases=" + system.getNumberOfPhases());
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      System.out.println("beta_" + i + "=" + system.getBeta(i));
      System.out.println("sound_speed_phase_" + i + "=" + system.getPhase(i).getSoundSpeed());
    }
    System.out.println("sound_speed_system=" + system.getSoundSpeed());
    System.out.println("sound_speed_phase_0=" + system.getPhase(0).getSoundSpeed());
    System.out.println("phase0_z=" + system.getPhase(0).getZ());
    System.out.println("phase0_density_kg_per_m3=" + system.getPhase(0).getDensity("kg/m3"));
    double moles = system.getPhase(0).getNumberOfMolesInPhase();
    System.out.println("phase0_cp_j_per_mol_k=" + system.getPhase(0).getCp() / moles);
    System.out.println("phase0_cv_j_per_mol_k=" + system.getPhase(0).getCv() / moles);
    System.out.println("phase0_molar_mass_kg_per_mol=" + system.getPhase(0).getMolarMass());
    System.out.println("temperature_k=" + system.getTemperature());
    System.out.println("pressure_pa=" + system.getPressure("Pa"));
    System.out.println();
  }

  /// The two alternative compressor map lookups, and a machine run on the extrapolating one.
  ///
  /// **The speed-interpolation block is `CompressorChartSpeedInterpolationTest`'s table**, read
  /// at the points that test asserts (the linear-in-speed head, the zero-speed short circuit,
  /// the two ends scaled by the speed ratio, the single line at `7000/6000`). The base lookup
  /// has no extrapolation: its spline throws outside the flow knots, and the block prints that
  /// exception beside the extrapolating class's answer. The machine row re-cases
  /// `CompressorChartTest.runCurveTest` on PR - twelve components, 29.96 C, 75.73 bara, speed
  /// 8765, the nine-speed table - on both classes.
  static void compressorMapLookupRows() {
    double[] conditions = new double[] {0.3, 1.0, 1.0, 1.0};
    double[] speed = new double[] {6000, 7000, 8000};
    double[][] flow = new double[][] {
        {1000, 2000, 3000}, {1000, 2000, 3000}, {1000, 2000, 3000}};
    double[][] head = new double[][] {
        {65, 60, 50}, {85, 80, 70}, {105, 100, 90}};
    double[][] efficiency = new double[][] {
        {75, 80, 75}, {75, 80, 75}, {75, 80, 75}};

    neqsim.process.equipment.compressor.CompressorChartAlternativeMapLookupExtrapolate extrapolate =
        new neqsim.process.equipment.compressor.CompressorChartAlternativeMapLookupExtrapolate();
    extrapolate.setCurves(conditions, speed, flow, head, efficiency);
    System.out.println("speed_interpolation_extrapolate");
    for (double rpm : new double[] {6001, 6500, 7000, 7500, 8000}) {
      System.out.println("head_2000_at_" + rpm + "=" + extrapolate.getPolytropicHead(2000.0, rpm));
    }
    System.out.println("head_2000_at_0=" + extrapolate.getPolytropicHead(2000.0, 0.0));
    System.out.println("head_2000_at_5000=" + extrapolate.getPolytropicHead(2000.0, 5000.0));
    System.out.println("head_2000_at_9000=" + extrapolate.getPolytropicHead(2000.0, 9000.0));
    System.out.println("efficiency_2000_at_7000="
        + extrapolate.getPolytropicEfficiency(2000.0, 7000.0));
    System.out.println("head_500_at_7000=" + extrapolate.getPolytropicHead(500.0, 7000.0));

    neqsim.process.equipment.compressor.CompressorChartAlternativeMapLookupExtrapolate single =
        new neqsim.process.equipment.compressor.CompressorChartAlternativeMapLookupExtrapolate();
    single.setCurves(conditions, new double[] {6000}, new double[][] {{1000, 2000, 3000}},
        new double[][] {{65, 60, 50}}, new double[][] {{75, 80, 75}});
    System.out.println("single_line_head_2000_at_7000=" + single.getPolytropicHead(2000.0, 7000.0));
    System.out.println();

    neqsim.process.equipment.compressor.CompressorChartAlternativeMapLookup lookup =
        new neqsim.process.equipment.compressor.CompressorChartAlternativeMapLookup();
    lookup.setCurves(conditions, speed, flow, head, efficiency);
    System.out.println("map_lookup");
    for (double rpm : new double[] {6000, 6500, 7000, 8000, 12000}) {
      System.out.println("head_2000_at_" + rpm + "=" + lookup.getPolytropicHead(2000.0, rpm));
    }
    System.out.println("efficiency_2000_at_6500=" + lookup.getPolytropicEfficiency(2000.0, 6500.0));
    try {
      System.out.println("head_500_at_7000=" + lookup.getPolytropicHead(500.0, 7000.0));
    } catch (RuntimeException refused) {
      System.out.println("head_500_at_7000=" + refused.getClass().getSimpleName());
    }
    System.out.println();

    compressorMapRunRow("map_lookup_run_extrapolate", "interpolate and extrapolate");
    compressorMapRunRow("map_lookup_run_lookup", "interpolate");
  }

  static void compressorMapRunRow(String label, String chartType) {
    SystemInterface fluid = new SystemPrEos(273.15 + 29.96, 75.73);
    String[] names = new String[] {"nitrogen", "CO2", "methane", "ethane", "propane", "i-butane",
        "n-butane", "i-pentane", "n-pentane", "n-hexane", "n-heptane", "n-octane", "water"};
    double[] z = new double[] {0.7823, 1.245, 88.4681, 4.7652, 2.3669, 0.3848, 0.873, 0.243,
        0.2947, 0.2455, 0.1735, 0.064, 0.076};
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setTemperature(29.96, "C");
    fluid.setPressure(75.73, "bara");
    fluid.setTotalFlowRate(559401.418270102, "kg/hr");
    Stream inlet = new Stream("Stream1", fluid);
    inlet.run();

    neqsim.process.equipment.compressor.Compressor compressor =
        new neqsim.process.equipment.compressor.Compressor("compressor 1", inlet);
    compressor.setCompressorChartType(chartType);
    compressor.setUsePolytropicCalc(true);
    compressor.setSpeed(8765);
    double[] conditions = new double[] {0.3, 1.0, 1.0, 1.0};
    double[] speed = new double[] {7000, 7500, 8000, 8500, 9000, 9500, 9659, 10000, 10500};
    double[][] flow = new double[][] {{4512.7, 5120.8, 5760.9, 6401, 6868.27},
        {4862.47, 5486.57, 6172.39, 6858.21, 7550.89}, {5237.84, 5852.34, 6583.88, 7315.43, 8046.97, 8266.43},
        {5642.94, 6218.11, 6995.38, 7772.64, 8549.9, 9000.72}, {6221.77, 6583.88, 7406.87, 8229.85, 9052.84, 9768.84},
        {6888.85, 6949.65, 7818.36, 8687.07, 9555.77, 10424.5, 10546.1},
        {7109.83, 7948.87, 8832.08, 9715.29, 10598.5, 10801.6}, {7598.9, 8229.85, 9144.28, 10058.7, 10973.1, 11338.9},
        {8334.1, 8641.35, 9601.5, 10561.6, 11521.8, 11963.5}};
    double[][] head = new double[][] {{61.885, 59.639, 56.433, 52.481, 49.132},
        {71.416, 69.079, 65.589, 61.216, 55.858}, {81.621, 79.311, 75.545, 70.727, 64.867, 62.879,},
        {92.493, 90.312, 86.3, 81.079, 74.658, 70.216}, {103.512, 102.073, 97.83, 92.254, 85.292, 77.638},
        {114.891, 114.632, 110.169, 104.221, 96.727, 87.002, 85.262},
        {118.595, 114.252, 108.203, 100.55, 90.532, 87.54}, {126.747, 123.376, 117.113, 109.056, 98.369, 92.632},
        {139.082, 137.398, 130.867, 122.264, 110.548, 103.247},};
    double[][] polyEff = new double[][] {{78.3, 78.2, 77.2, 75.4, 73.4},
        {78.3, 78.3, 77.5, 75.8, 73}, {78.2, 78.4, 77.7, 76.1, 73.5, 72.5}, {78.2, 78.4, 77.9, 76.4, 74, 71.9},
        {78.3, 78.4, 78, 76.7, 74.5, 71.2}, {78.3, 78.4, 78.1, 77, 74.9, 71.3, 70.5},
        {78.4, 78.1, 77.1, 75, 71.4, 70.2}, {78.3, 78.2, 77.2, 75.2, 71.7, 69.5}, {78.2, 78.2, 77.3, 75.5, 72.2, 69.6}};
    compressor.getCompressorChart().setCurves(conditions, speed, flow, head, polyEff);
    compressor.getCompressorChart().setHeadUnit("kJ/kg");
    compressor.getCompressorChart().getSurgeCurve().setCurve(conditions,
        new double[] {4512.7, 4862.47, 5237.84, 5642.94, 6221.77, 6888.85, 7109.83, 7598.9},
        new double[] {61.885, 71.416, 81.621, 92.493, 103.512, 114.891, 118.595, 126.747});
    compressor.getAntiSurge().setSurgeControlFactor(1.0);
    compressor.run();

    System.out.println(label);
    System.out.println("chart_type=" + chartType);
    System.out.println("inlet_m3_per_hr=" + compressor.getInletStream().getFlowRate("m3/hr"));
    System.out.println("polytropic_head_kj_per_kg=" + compressor.getPolytropicFluidHead());
    System.out.println("polytropic_efficiency_fraction=" + compressor.getPolytropicEfficiency());
    System.out.println("outlet_p_bara=" + compressor.getOutletStream().getPressure("bara"));
    System.out.println("outlet_t_k=" + compressor.getOutletStream().getTemperature("K"));
    System.out.println("power_w=" + compressor.getPower("W"));
    System.out.println("distance_to_surge=" + compressor.getDistanceToSurge());
    print("inlet", inlet);
    print("outlet", compressor.getOutletStream());
    System.out.println();
  }

  /// The anti-surge surface: the boundary splines, the chart's generators, the read-outs, and
  /// `Compressor.run`'s own recycle block.
  ///
  /// **The spline rows are `SafeSplineSurgeCurveTest.testSurgeCurve10250`'s**, read at the
  /// points that test asserts. The chart rows call the class's own `generateSurgeCurve` and
  /// `generateStoneWallCurve` on the eight-speed `CompressorChartTest` table, so the generator
  /// is held as well as the curve. The machine rows re-case the same nine-component gas at
  /// 24 C and 48 bara: one state read out with a 10% control margin, then the 0.30 and
  /// 4.00 MSm3/day states of `AntiSurgeRecycleCalculatorTest` with `AntiSurge` active and the
  /// default 1.05 surge-control factor - one inside surge, one far above it.
  static void antiSurgeRows() {
    double[] flow10250 = new double[] {9758.49, 9578.11, 9397.9, 9248.64, 9006.93, 8749.97,
        8508.5, 8179.81, 7799.81, 7111.75, 6480.26, 6007.91, 5607.45};
    double[] head10250 = new double[] {112.65, 121.13, 127.56, 132.13, 137.29, 140.73, 142.98,
        144.76, 146.14, 148.05, 148.83, 149.54, 150.0};
    neqsim.process.equipment.compressor.SafeSplineSurgeCurve curve =
        new neqsim.process.equipment.compressor.SafeSplineSurgeCurve(flow10250, head10250);
    System.out.println("spline_10250");
    System.out.println("surge_flow_at_130=" + curve.getSurgeFlow(130.0));
    System.out.println("surge_flow_at_190=" + curve.getSurgeFlow(190.0));
    System.out.println("surge_flow_at_10=" + curve.getSurgeFlow(10.0));
    System.out.println("surge_head_at_9321_84315=" + curve.getSurgeHead(9321.84315));
    System.out.println("surge_head_at_6000=" + curve.getSurgeHead(6000.0));
    System.out.println("surge_head_at_9500=" + curve.getSurgeHead(9500.0));
    System.out.println();

    neqsim.process.equipment.compressor.CompressorChart chart = compressorChart();
    chart.generateSurgeCurve();
    chart.generateStoneWallCurve();
    System.out.println("generated_boundaries");
    for (double head : new double[] {31.0, 50.0, 80.0}) {
      System.out.println("surge_flow_at_" + head + "=" + chart.getSurgeCurve().getSurgeFlow(head));
    }
    for (double head : new double[] {31.0, 50.0, 80.0}) {
      System.out.println(
          "stonewall_flow_at_" + head + "=" + chart.getStoneWallCurve().getStoneWallFlow(head));
    }
    System.out.println("surge_flow_at_speed_11683=" + chart.getSurgeFlowAtSpeed(11683.0));
    System.out.println("stonewall_flow_at_speed_11683=" + chart.getStoneWallFlowAtSpeed(11683.0));
    System.out.println();

    antiSurgeReadoutRow("anti_surge_readout");
    antiSurgeRunRow("anti_surge_run_in_surge", 0.30);
    antiSurgeRunRow("anti_surge_run_above_surge", 4.00);
  }

  static neqsim.process.equipment.compressor.Compressor antiSurgeMachine(double msm3PerDay) {
    SystemInterface fluid = new SystemPrEos(298.15, 50.0);
    String[] names = new String[] {"nitrogen", "CO2", "methane", "ethane", "propane", "i-butane",
        "n-butane", "i-pentane", "n-pentane"};
    double[] z = new double[] {1.205, 1.340, 87.974, 5.258, 3.283, 0.082, 0.487, 0.056, 0.053};
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    fluid.setMultiPhaseCheck(true);
    fluid.setTemperature(24.0, "C");
    fluid.setPressure(48.0, "bara");
    fluid.setTotalFlowRate(msm3PerDay, "MSm3/day");
    Stream inlet = new Stream("suction", fluid);
    inlet.run();

    neqsim.process.equipment.compressor.Compressor compressor =
        new neqsim.process.equipment.compressor.Compressor("export", inlet);
    compressor.setUsePolytropicCalc(true);
    compressor.setSpeed(11918);
    compressor.getCompressorChart().setCurves(compressorChart().getChartConditions(),
        compressorChart().getSpeeds(), compressorChart().getFlows(),
        compressorChart().getHeads(), compressorChart().getPolytropicEfficiencies());
    compressor.getCompressorChart().setHeadUnit("kJ/kg");
    compressor.getCompressorChart().setUseCompressorChart(true);
    compressor.getCompressorChart().getSurgeCurve().setCurve(
        compressorChart().getChartConditions(), new double[] {2789.0, 2550.0, 2500.0, 2200.0},
        new double[] {80.0, 72.0, 70.0, 65.0});
    return compressor;
  }

  static void printAntiSurge(String label,
      neqsim.process.equipment.compressor.Compressor compressor) {
    System.out.println(label);
    System.out.println("speed_rpm=" + compressor.getSpeed());
    System.out.println("polytropic_head_kj_per_kg=" + compressor.getPolytropicHead());
    System.out.println("polytropic_fluid_head=" + compressor.getPolytropicFluidHead());
    System.out.println("polytropic_efficiency_fraction=" + compressor.getPolytropicEfficiency());
    System.out.println("power_w=" + compressor.getPower("W"));
    System.out.println("outlet_p_bara=" + compressor.getOutletStream().getPressure("bara"));
    System.out.println("outlet_t_k=" + compressor.getOutletStream().getTemperature("K"));
    System.out.println("inlet_m3_per_hr=" + compressor.getInletStream().getFlowRate("m3/hr"));
    System.out.println("surge_flow_m3_per_hr=" + compressor.getSurgeFlowRate());
    System.out.println("surge_flow_margin_m3_per_hr=" + compressor.getSurgeFlowRateMargin());
    System.out.println("surge_flow_std_m3_per_hr=" + compressor.getSurgeFlowRateStd());
    System.out.println("distance_to_surge=" + compressor.getDistanceToSurge());
    System.out.println("control_line_flow_m3_per_hr=" + compressor.getControlLineFlow());
    System.out.println("distance_to_control_line=" + compressor.getDistanceToControlLine());
    System.out.println(
        "required_recycle_fraction=" + compressor.getRequiredRecycleFractionToControlLine());
    System.out.println("surge_fraction=" + compressor.getAntiSurge().getCurrentSurgeFraction());
    System.out.println("recycle_power_kw="
        + compressor.getAntiSurgeRecyclePower(
            compressor.getAntiSurge().getCurrentSurgeFraction(), "kW"));
    print("inlet", compressor.getInletStream());
    print("outlet", compressor.getOutletStream());
    System.out.println();
  }

  // === process.anti_surge_recycle =============================================

  /// `AntiSurgeRecycleCalculator.solve` on the charted machine `antiSurgeMachine` builds - the
  /// pinned `AntiSurgeRecycleCalculatorTest` states, re-cased on PR.
  ///
  /// The natural suction is the machine's own inlet; `solve` runs it, and when its inlet volume
  /// flow is below the control line it wires a cooled clone of the suction through a mixer and
  /// converges the recycle mass with the calculator's damped fixed point. The class reports
  /// m3/hr for the flows and kg/hr for the recycle mass.
  static void antiSurgeRecycleRows() {
    antiSurgeRecycleRow("low_flow_margin_005", 0.30, 0.05);
    antiSurgeRecycleRow("high_flow_no_recycle", 4.00, 0.05);
    antiSurgeRecycleRow("low_flow_margin_010", 0.30, 0.10);
  }

  static void antiSurgeRecycleRow(String label, double msm3PerDay, double margin) {
    neqsim.process.equipment.compressor.Compressor compressor = antiSurgeMachine(msm3PerDay);
    Stream suction = (Stream) compressor.getInletStream();
    neqsim.process.equipment.compressor.AntiSurgeRecycleCalculator calculator =
        new neqsim.process.equipment.compressor.AntiSurgeRecycleCalculator(compressor, suction);
    calculator.setSurgeControlMargin(margin);
    calculator.setRecycleCoolerTemperature(40.0, "C");
    neqsim.process.equipment.compressor.AntiSurgeRecycleCalculator.Result result =
        calculator.solve();

    System.out.println(label);
    System.out.println("surge_control_margin=" + margin);
    System.out.println("recycle_cooler_temperature_c=40.0");
    System.out.println("recycle_active=" + result.isRecycleActive());
    System.out.println("converged=" + result.isConverged());
    System.out.println("iterations=" + result.getIterations());
    System.out.println("recycle_mass_flow_kg_hr=" + result.getRecycleMassFlow());
    System.out.println("inlet_volume_flow_m3_hr=" + result.getInletVolumeFlow());
    System.out.println("surge_flow_m3_hr=" + result.getSurgeFlow());
    System.out.println("control_flow_m3_hr=" + result.getControlFlow());
    System.out.println("distance_to_surge=" + result.getDistanceToSurge());
    System.out.println("compressor_inlet_m3_hr="
        + compressor.getInletStream().getFlowRate("m3/hr"));
    System.out.println("compressor_inlet_kg_hr="
        + compressor.getInletStream().getFlowRate("kg/hr"));
    System.out.println("compressor_power_kw=" + compressor.getPower("kW"));
    System.out.println("compressor_outlet_p_bara="
        + compressor.getOutletStream().getPressure("bara"));
    print("compressor_inlet", compressor.getInletStream());
    print("compressor_outlet", compressor.getOutletStream());
    System.out.println();
  }

  static void antiSurgeReadoutRow(String label) {
    neqsim.process.equipment.compressor.Compressor compressor = antiSurgeMachine(1.0);
    compressor.setSurgeControlMargin(0.10);
    compressor.run();
    printAntiSurge(label, compressor);
  }

  static void antiSurgeRunRow(String label, double msm3PerDay) {
    neqsim.process.equipment.compressor.Compressor compressor = antiSurgeMachine(msm3PerDay);
    compressor.getAntiSurge().setActive(true);
    compressor.getAntiSurge().setSurgeControlFactor(1.05);
    compressor.run();
    printAntiSurge(label, compressor);
  }

  /// `CompressorChart`, the default fan-law class, and `Compressor.run`'s chart branch.
  ///
  /// **The chart's own rows are fluid-free**, so they compare the reduced fit and the two
  /// inverses alone. The run row re-cases `CompressorCalculationTest.testRun` on PR: its
  /// 9-component gas at 24 C, 48 bara and 1 MSm3/day, the eight-speed `CompressorChartTest`
  /// table in kJ/kg and speed 11918.
  static void compressorCurveRows() {
    neqsim.process.equipment.compressor.CompressorChart chart = compressorChart();

    for (double[] state : new double[][] {{3974.8455, 12913.0}, {3599.6333, 10453.0},
        {1636.5807, 8200.0}, {6000.0, 12913.0}, {100.0, 12913.0}}) {
      System.out.println("chart_row");
      System.out.println("flow_m3_per_hr=" + state[0]);
      System.out.println("speed_rpm=" + state[1]);
      System.out.println("polytropic_head_kj_per_kg=" + chart.getPolytropicHead(state[0], state[1]));
      System.out.println("polytropic_efficiency_percent="
          + chart.getPolytropicEfficiency(state[0], state[1]));
      System.out.println();
    }

    System.out.println("speed_inverse");
    System.out.println("flow_m3_per_hr=3974.8455");
    System.out.println("head_kj_per_kg=70.0");
    System.out.println("speed_rounded=" + chart.getSpeed(3974.8455, 70.0));
    System.out.println("speed_value=" + chart.getSpeedValue(3974.8455, 70.0));
    System.out.println();

    System.out.println("flow_inverse");
    System.out.println("head_kj_per_kg=70.0");
    System.out.println("speed_rpm=12913.0");
    System.out.println("guess_flow_m3_per_hr=4000.0");
    System.out.println("flow_m3_per_hr=" + chart.getFlow(70.0, 12913.0, 4000.0));
    System.out.println();

    double speed = 11683.0;
    double surge = chart.getSurgeFlowAtSpeed(speed);
    double stonewall = chart.getStonewallFlowAtSpeed(speed);
    System.out.println("range_status");
    System.out.println("speed_rpm=" + speed);
    System.out.println("surge_flow_m3_per_hr=" + surge);
    System.out.println("stonewall_flow_m3_per_hr=" + stonewall);
    System.out.println("surge_head_kj_per_kg=" + chart.getSurgeHeadAtSpeed(speed));
    System.out.println("stonewall_head_kj_per_kg=" + chart.getStoneWallHeadAtSpeed(speed));
    System.out.println("status_mid=" + chart.getFlowRangeStatus(0.5 * (surge + stonewall), speed));
    System.out.println("status_low=" + chart.getFlowRangeStatus(surge - 500.0, speed));
    System.out.println("status_high=" + chart.getFlowRangeStatus(stonewall + 500.0, speed));
    System.out.println();

    compressorRunRow("compressor_chart_run", chart);
  }

  /**
   * `EquilibriumSoundSpeedTest`'s states on the pinned class.
   *
   * The test disables `useVolumeCorrection`; this probe keeps the default, because azoth's PR
   * carries the Peneloux shift and the capture is only comparable on the same convention.
   */
  static void equilibriumSoundSpeedRows() {
    eqSpeedRow("dense_co2", eqCo2Feed(null), 1.0e-3);
    eqSpeedRow("dilute_nitrogen", eqIdealGas(), 1.0e-3);

    SystemInterface twoPhase = eqCo2Feed(null);
    double entropy = twoPhase.getEntropy("J/kgK");
    twoPhase.setPressure(40.0);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(twoPhase).PSflash(entropy, "J/kgK");
    twoPhase.init(3);
    eqSpeedRow("two_phase_co2", twoPhase, 1.0e-3);

    for (double pressure : new double[] {50.0, 40.0, 20.0}) {
      SystemInterface mixture = eqReferenceTpRoot(eqCo2Feed("nitrogen"), pressure,
          eqCo2Feed("nitrogen").getEntropy("J/kgK"));
      eqSpeedRow("mixture_nitrogen_" + (int) pressure, mixture, 1.0e-3);
    }
    SystemInterface hydrogen = eqReferenceTpRoot(eqCo2Feed("hydrogen"), 40.0,
        eqCo2Feed("hydrogen").getEntropy("J/kgK"));
    eqSpeedRow("mixture_hydrogen_40", hydrogen, 1.0e-3);

    // **Through a Stream, because a bare setT/setP + TPflash + init(3) leaves `s` stale** -
    // the state the pinned test hands the class reads the feed's entropy at the new T and P.
    // The class's own closure cannot see that; the probe can, by building the state the way a
    // solver does.
    SystemInterface boundaryFluid = eqCo2Feed(null);
    boundaryFluid.setPressure(46.5);
    boundaryFluid.setTemperature(284.3849043789973);
    Stream boundaryStream = new Stream("boundary", boundaryFluid);
    boundaryStream.run();
    eqSpeedRow("phase_entry_46_5", boundaryStream.getThermoSystem(), 0.05);

    SystemInterface steps = eqReferenceTpRoot(eqCo2Feed("nitrogen"), 40.0,
        eqCo2Feed("nitrogen").getEntropy("J/kgK"));
    eqSpeedRow("step_0_0002", steps, 2.0e-4);
    eqSpeedRow("step_0_002", steps, 2.0e-3);
    eqSpeedRow("step_0_01", steps, 1.0e-2);
  }

  static SystemInterface eqCo2Feed(String impurity) {
    SystemInterface fluid = new SystemPrEos(298.15, 150.0);
    fluid.addComponent("CO2", impurity == null ? 1.0 : 0.98);
    if (impurity != null) {
      fluid.addComponent(impurity, 0.02);
    }
    fluid.setMixingRule(2);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(fluid).TPflash();
    fluid.init(3);
    return fluid;
  }

  static SystemInterface eqIdealGas() {
    SystemInterface gas = new SystemPrEos(350.0, 0.001);
    gas.addComponent("nitrogen", 1.0);
    gas.setMixingRule(2);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(gas).TPflash();
    gas.init(3);
    return gas;
  }

  static SystemInterface eqReferenceTpRoot(SystemInterface original, double pressure,
      double entropy) {
    double lower = 200.0;
    double upper = 350.0;
    SystemInterface trial = null;
    for (int i = 0; i < 60; i++) {
      trial = original.clone();
      trial.setPressure(pressure);
      trial.setTemperature((lower + upper) / 2.0);
      new neqsim.thermodynamicoperations.ThermodynamicOperations(trial).TPflash();
      trial.init(3);
      double error = trial.getEntropy("J/kgK") - entropy;
      if (Math.abs(error) < 1.0e-9) {
        break;
      }
      if (error > 0.0) {
        upper = trial.getTemperature();
      } else {
        lower = trial.getTemperature();
      }
    }
    return trial;
  }

  static void eqSpeedRow(String label, SystemInterface fluid, double step) {
    double targetEntropy = fluid.getEntropy("J/kgK");
    double phaseSpeed = fluid.getSoundSpeed();
    double temperature = fluid.getTemperature();
    double pressure = fluid.getPressure();
    int phases = fluid.getNumberOfPhases();
    StringBuilder types = new StringBuilder();
    StringBuilder fractions = new StringBuilder();
    for (int i = 0; i < phases; i++) {
      types.append(i == 0 ? "" : ",").append(fluid.getPhase(i).getType().name());
      fractions.append(i == 0 ? "" : ",").append(fluid.getBeta(i));
    }
    StringBuilder z = new StringBuilder();
    for (int i = 0; i < fluid.getPhase(0).getNumberOfComponents(); i++) {
      z.append(i == 0 ? "" : ",")
          .append(fluid.getPhase(0).getComponent(i).getName()).append(":")
          .append(fluid.getPhase(0).getComponent(i).getz());
    }

    neqsim.thermo.util.EquilibriumSoundSpeed.Result result =
        neqsim.thermo.util.EquilibriumSoundSpeed.calculate(fluid, step);

    System.out.println(label);
    System.out.println("input_p_bara=" + pressure);
    System.out.println("input_t_k=" + temperature);
    System.out.println("input_z=" + z);
    System.out.println("input_phase_count=" + phases);
    System.out.println("input_phase_types=" + types);
    System.out.println("input_phase_fractions=" + fractions);
    System.out.println("specific_entropy_j_per_kgk=" + targetEntropy);
    System.out.println("phase_sound_speed_m_per_s=" + phaseSpeed);
    System.out.println("relative_pressure_step=" + step);
    System.out.println("status=" + result.getStatus());
    System.out.println("message=" + result.getMessage());
    System.out.println("sound_speed_m_per_s=" + result.getSoundSpeed());
    System.out.println("pressure_step_pa=" + result.getPressureStepPa());
    System.out.println("relative_step_error=" + result.getRelativeStepError());
    System.out.println("stencil=" + result.getStencil());
    System.out.println("phase_boundary_encountered=" + result.isPhaseBoundaryEncountered());
    System.out.println("flash_evaluations=" + result.getFlashEvaluations());
    int sample = 0;
    for (neqsim.thermo.util.EquilibriumSoundSpeed.State state : result.getSamples()) {
      StringBuilder phaseTypes = new StringBuilder();
      StringBuilder phaseFractions = new StringBuilder();
      for (int i = 0; i < state.getPhaseTypes().length; i++) {
        phaseTypes.append(i == 0 ? "" : ",").append(state.getPhaseTypes()[i]);
        phaseFractions.append(i == 0 ? "" : ",").append(state.getPhaseFractions()[i]);
      }
      System.out.println("sample_" + sample + "_p_pa=" + state.getPressurePa());
      System.out.println("sample_" + sample + "_t_k=" + state.getTemperatureK());
      System.out.println("sample_" + sample + "_density_kg_per_m3=" + state.getDensity());
      System.out.println("sample_" + sample + "_entropy_residual=" + state.getEntropyResidual());
      System.out.println("sample_" + sample + "_component_residual=" + state.getComponentResidual());
      System.out.println("sample_" + sample + "_fugacity_residual=" + state.getFugacityResidual());
      System.out.println("sample_" + sample + "_phase_types=" + phaseTypes);
      System.out.println("sample_" + sample + "_phase_fractions=" + phaseFractions);
      sample++;
    }
    System.out.println();
  }

  /** `CompressorChartTest.runCurveTestMaxMinSpeed`'s five runs, on the pinned class. */
  static void compressorSpeedSolveRows() {
    double[] chartConditions = new double[] {0.3, 1.0, 1.0, 1.0};
    double[] speed = new double[] {7000, 7500, 8000, 8500, 9000, 9500, 9659, 10000, 10500};
    double[][] flow = new double[][] {{4512.7, 5120.8, 5760.9, 6401, 6868.27},
        {4862.47, 5486.57, 6172.39, 6858.21, 7550.89},
        {5237.84, 5852.34, 6583.88, 7315.43, 8046.97, 8266.43},
        {5642.94, 6218.11, 6995.38, 7772.64, 8549.9, 9000.72},
        {6221.77, 6583.88, 7406.87, 8229.85, 9052.84, 9768.84},
        {6888.85, 6949.65, 7818.36, 8687.07, 9555.77, 10424.5, 10546.1},
        {7109.83, 7948.87, 8832.08, 9715.29, 10598.5, 10801.6},
        {7598.9, 8229.85, 9144.28, 10058.7, 10973.1, 11338.9},
        {8334.1, 8641.35, 9601.5, 10561.6, 11521.8, 11963.5}};
    double[][] head = new double[][] {{61.885, 59.639, 56.433, 52.481, 49.132},
        {71.416, 69.079, 65.589, 61.216, 55.858},
        {81.621, 79.311, 75.545, 70.727, 64.867, 62.879},
        {92.493, 90.312, 86.3, 81.079, 74.658, 70.216},
        {103.512, 102.073, 97.83, 92.254, 85.292, 77.638},
        {114.891, 114.632, 110.169, 104.221, 96.727, 87.002, 85.262},
        {118.595, 114.252, 108.203, 100.55, 90.532, 87.54},
        {126.747, 123.376, 117.113, 109.056, 98.369, 92.632},
        {139.082, 137.398, 130.867, 122.264, 110.548, 103.247}};
    double[][] polyEff = new double[][] {{78.3, 78.2, 77.2, 75.4, 73.4},
        {78.3, 78.3, 77.5, 75.8, 73},
        {78.2, 78.4, 77.7, 76.1, 73.5, 72.5},
        {78.2, 78.4, 77.9, 76.4, 74, 71.9},
        {78.3, 78.4, 78, 76.7, 74.5, 71.2},
        {78.3, 78.4, 78.1, 77, 74.9, 71.3, 70.5},
        {78.4, 78.1, 77.1, 75, 71.4, 70.2},
        {78.3, 78.2, 77.2, 75.2, 71.7, 69.5},
        {78.2, 78.2, 77.3, 75.5, 72.2, 69.6}};
    double[] surgeflow =
        new double[] {4512.7, 4862.47, 5237.84, 5642.94, 6221.77, 6888.85, 7109.83, 7598.9};
    double[] surgehead =
        new double[] {61.885, 71.416, 81.621, 92.493, 103.512, 114.891, 118.595, 126.747};

    SystemInterface fluid = new SystemPrEos(273.15 + 29.96, 75.73);
    fluid.addComponent("nitrogen", 0.7823);
    fluid.addComponent("CO2", 1.245);
    fluid.addComponent("methane", 88.4681);
    fluid.addComponent("ethane", 4.7652);
    fluid.addComponent("propane", 2.3669);
    fluid.addComponent("i-butane", 0.3848);
    fluid.addComponent("n-butane", 0.873);
    fluid.setMixingRule("classic");

    Stream inlet = new Stream("Stream1", fluid);
    inlet.setTemperature(29.96, "C");
    inlet.setPressure(75.73, "bara");
    inlet.setFlowRate(559401.4, "kg/hr");
    inlet.run();

    neqsim.process.equipment.compressor.Compressor compressor =
        new neqsim.process.equipment.compressor.Compressor("compressor 1", inlet);
    compressor.setUsePolytropicCalc(true);
    compressor.setPolytropicEfficiency(0.85);
    compressor.setSpeed(9000);
    compressor.getCompressorChart().setCurves(chartConditions, speed, flow, head, polyEff);
    compressor.getCompressorChart().setHeadUnit("kJ/kg");
    compressor.getCompressorChart().getSurgeCurve().setCurve(chartConditions, surgeflow, surgehead);
    compressor.getAntiSurge().setSurgeControlFactor(1.0);
    compressor.getCompressorChart().setUseCompressorChart(true);
    compressor.setSolveSpeed(true);

    // The pinned test's sequence, run in order on one machine.
    double[][] rows = new double[][] {{240.0, 7000.0, 10500.0, 0.0},
        {240.0, 7000.0, 10500.0, 1.0}, {100.0, 7000.0, 10500.0, 1.0},
        {100.0, 7000.0, 10500.0, 0.0}, {90.0, 0.0, 10500.0, 0.0}};
    double[] initial = new double[] {9000.0, 9000.0, 9000.0, 15000.0, 3000.0};

    for (int i = 0; i < rows.length; i++) {
      if (i == 3) {
        compressor.setSpeed(initial[i]);
      }
      if (i == 4) {
        compressor.setMinimumSpeed(rows[i][1]);
        compressor.setSpeed(initial[i]);
      }
      compressor.setOutletPressure(rows[i][0], "bara");
      compressor.setMaximumSpeed(rows[i][2]);
      if (i != 4) {
        compressor.setMinimumSpeed(rows[i][1]);
      }
      compressor.setLimitSpeed(rows[i][3] > 0.5);
      double ranFrom = compressor.getSpeed();
      compressor.run();

      System.out.println("speed_solve_row_" + (i + 1));
      System.out.println("target_outlet_p_bara=" + rows[i][0]);
      System.out.println("min_speed_rpm=" + rows[i][1]);
      System.out.println("max_speed_rpm=" + rows[i][2]);
      System.out.println("limit_speed=" + (rows[i][3] > 0.5));
      System.out.println("initial_speed_rpm=" + ranFrom);
      System.out.println("inlet_n_mol_per_s=" + inlet.getThermoSystem().getTotalNumberOfMoles());
      System.out.println("speed_rpm=" + compressor.getSpeed());
      System.out.println("outlet_p_bara=" + compressor.getOutletStream().getPressure("bara"));
      System.out.println("outlet_t_k=" + compressor.getOutletStream().getTemperature("K"));
      System.out.println("outlet_h_j_per_mol="
          + compressor.getOutletStream().getThermoSystem().getEnthalpy()
              / compressor.getOutletStream().getThermoSystem().getTotalNumberOfMoles());
      System.out.println("power_w=" + compressor.getPower());
      System.out.println("polytropic_head_kj_per_kg=" + compressor.getPolytropicFluidHead());
      System.out.println("polytropic_head_m=" + compressor.getPolytropicHeadMeter());
      System.out.println("polytropic_efficiency_fraction=" + compressor.getPolytropicEfficiency());
      print("inlet", inlet);
      print("outlet", compressor.getOutletStream());
      System.out.println();
    }
  }

  static neqsim.process.equipment.compressor.CompressorChart compressorChart() {
    double[] chartConditions = new double[] {0.3, 1.0, 1.0, 1.0};
    double[] speed = new double[] {12913, 12298, 11683, 11098, 10453, 9224, 8609, 8200};
    double[][] flow = new double[][] {
        {2789.1285, 3174.0375, 3689.2288, 4179.4503, 4570.2768, 4954.7728, 5246.0329, 5661.0331},
        {2571.1753, 2943.7254, 3440.2675, 3837.4448, 4253.0898, 4668.6643, 4997.1926, 5387.4952},
        {2415.3793, 2763.0706, 3141.7095, 3594.7436, 4047.6467, 4494.1889, 4853.7353, 5138.7858},
        {2247.2043, 2799.7342, 3178.3428, 3656.1551, 4102.778, 4394.1591, 4648.3224, 4840.4998},
        {2072.8397, 2463.9483, 2836.4078, 3202.5266, 3599.6333, 3978.0203, 4257.0022, 4517.345},
        {1835.9552, 2208.455, 2618.1322, 2940.8034, 3244.7852, 3530.1279, 3753.3738, 3895.9746},
        {1711.3386, 1965.8848, 2356.9431, 2685.9247, 3008.5154, 3337.2855, 3591.5092},
        {1636.5807, 2002.8708, 2338.0319, 2642.1245, 2896.4894, 3113.6264, 3274.8764, 3411.2977}};
    double[][] head = new double[][] {
        {80.0375, 78.8934, 76.2142, 71.8678, 67.0062, 60.6061, 53.0499, 39.728},
        {72.2122, 71.8369, 68.9009, 65.8341, 60.7167, 54.702, 47.2749, 35.7471},
        {65.1576, 64.5253, 62.6118, 59.1619, 54.0455, 47.0059, 39.195, 31.6387},
        {58.6154, 56.9627, 54.6647, 50.4462, 44.4322, 38.4144, 32.9084, 28.8109},
        {52.3295, 51.0573, 49.5283, 46.3326, 42.3685, 37.2502, 31.4884, 25.598},
        {40.6578, 39.6416, 37.6008, 34.6603, 30.9503, 27.1116, 23.2713, 20.4546},
        {35.2705, 34.6359, 32.7228, 31.0645, 27.0985, 22.7482, 18.0113},
        {32.192, 31.1756, 29.1329, 26.833, 23.8909, 21.3324, 18.7726, 16.3403}};
    double[][] polyEff = new double[][] {
        {77.2452238409573, 79.4154186459363, 80.737960012489, 80.5229826589649, 79.2210931638144,
            75.4719133864634, 69.6034181197298, 58.7322388482707},
        {77.0107837113504, 79.3069974136389, 80.8941189021135, 80.7190194665918, 79.5313242980328,
            75.5912622896367, 69.6846136362097, 60.0043057990909},
        {77.0043065299874, 79.1690958847856, 80.8038169975675, 80.6543975614197, 78.8532389102705,
            73.6664774270613, 66.2735600426727, 57.671664571658},
        {77.0716623789093, 80.4629750233093, 81.1390811169072, 79.6374242667478, 75.380928428817,
            69.5332969549779, 63.7997587622339, 58.8120614497758},
        {76.9705872525642, 79.8335492585324, 80.9468133671171, 80.5806471927835, 78.0462158225426,
            73.0403707523258, 66.5572286338589, 59.8624822515064},
        {77.5063036680357, 80.2056198362559, 81.0339108025933, 79.6085962687939, 76.3814534404405,
            70.8027503005902, 64.6437367160571, 60.5299349982342},
        {77.8175271586685, 80.065165942218, 81.0631362122632, 79.8955051771299, 76.1983240929369,
            69.289982774309, 60.8567149372229},
        {78.0924334304045, 80.9353551568667, 80.7904437766234, 78.8639325223295, 75.2170936751143,
            70.3105081673411, 65.5507568533569, 61.0391468300337}};
    neqsim.process.equipment.compressor.CompressorChart chart =
        new neqsim.process.equipment.compressor.CompressorChart();
    chart.setCurves(chartConditions, speed, flow, head, polyEff);
    chart.setHeadUnit("kJ/kg");
    return chart;
  }

  static void compressorRunRow(String label,
      neqsim.process.equipment.compressor.CompressorChart chart) {
    SystemInterface fluid = new SystemPrEos(298.15, 50.0);
    String[] names = new String[] {"nitrogen", "CO2", "methane", "ethane", "propane", "i-butane",
        "n-butane", "i-pentane", "n-pentane"};
    double[] z = new double[] {1.205, 1.340, 87.974, 5.258, 3.283, 0.082, 0.487, 0.056, 0.053};
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    fluid.setMultiPhaseCheck(true);
    fluid.setTemperature(24.0, "C");
    fluid.setPressure(48.0, "bara");
    fluid.setTotalFlowRate(1.0, "MSm3/day");
    Stream inlet = new Stream("Stream1", fluid);
    inlet.run();

    neqsim.process.equipment.compressor.Compressor compressor =
        new neqsim.process.equipment.compressor.Compressor("cmp1", inlet);
    compressor.setUsePolytropicCalc(true);
    compressor.setSpeed(11918);
    compressor.getCompressorChart().setCurves(chart.getChartConditions(), chart.getSpeeds(),
        chart.getFlows(), chart.getHeads(), chart.getPolytropicEfficiencies());
    compressor.getCompressorChart().setHeadUnit("kJ/kg");
    compressor.run();

    System.out.println(label);
    System.out.println("inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("inlet_z=" + inlet.getThermoSystem().getZ());
    System.out.println("inlet_gamma2=" + inlet.getThermoSystem().getGamma2());
    System.out.println("inlet_molar_mass_gr_per_mol="
        + inlet.getThermoSystem().getMolarMass("gr/mol"));
    System.out.println("inlet_density_kg_per_m3="
        + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("speed_rpm=" + compressor.getSpeed());
    System.out.println("polytropic_head_chart_kj_per_kg=" + compressor.getPolytropicHead());
    System.out.println("polytropic_efficiency_fraction="
        + compressor.getPolytropicEfficiency());
    System.out.println("power_w=" + compressor.getPower());
    System.out.println("outlet_p_bara=" + compressor.getOutletStream().getPressure("bara"));
    System.out.println("outlet_t_k=" + compressor.getOutletStream().getTemperature("K"));
    print("inlet", inlet);
    print("outlet", compressor.getOutletStream());
    System.out.println();
  }

  /// `Expander`'s stepped polytropic branch, 5 and 40 steps, and the isentropic power.
  static void expanderStepRows() {
    for (int steps : new int[] {5, 40}) {
      Stream inlet = expanderFeed();
      neqsim.process.equipment.expander.Expander expander =
          new neqsim.process.equipment.expander.Expander("expander", inlet);
      expander.setOutletPressure(30.0, "bara");
      expander.setPolytropicEfficiency(0.80);
      expander.setUsePolytropicCalc(true);
      expander.setNumberOfCompressorCalcSteps(steps);
      expander.run();

      System.out.println("polytropic_steps=" + steps);
      System.out.println("inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
      System.out.println("inlet_molar_mass_gr_per_mol="
          + inlet.getThermoSystem().getMolarMass("gr/mol"));
      System.out.println("power_w=" + expander.getPower());
      System.out.println("outlet_p_bara=" + expander.getOutletStream().getPressure("bara"));
      System.out.println("outlet_t_k=" + expander.getOutletStream().getTemperature("K"));
      print("inlet", inlet);
      print("outlet", expander.getOutletStream());
      System.out.println();
    }

    Stream inlet = expanderFeed();
    neqsim.process.equipment.expander.Expander expander =
        new neqsim.process.equipment.expander.Expander("expander", inlet);
    expander.setOutletPressure(30.0, "bara");
    expander.setIsentropicEfficiency(0.80);
    expander.run();
    System.out.println("isentropic_efficiency=0.80");
    System.out.println("power_w=" + expander.getPower());
    System.out.println("outlet_p_bara=" + expander.getOutletStream().getPressure("bara"));
    System.out.println("outlet_t_k=" + expander.getOutletStream().getTemperature("K"));
    print("inlet", inlet);
    print("outlet", expander.getOutletStream());
    System.out.println();
  }

  /// `Expander`'s two remaining steady surfaces: the stated power and the outlet-temperature solve.
  ///
  /// **`setPower` states the total enthalpy change the run applies**, and the class then
  /// back-computes `isentropicEfficiency = dH / (h_isentropic - h_in)` - a value above one is
  /// allowed. **The outlet-temperature surface is the inherited `solveEfficiency`**, a secant
  /// over the efficiency with a one-hundredth first step and a 50-iteration cap; the expander's
  /// own `run` never reads `useOutTemperature`, so the solve is a method rather than a mode.
  static void expanderPowerModeRows() {
    Stream baseInlet = expanderFeed();
    neqsim.process.equipment.expander.Expander base =
        new neqsim.process.equipment.expander.Expander("base", baseInlet);
    base.setOutletPressure(30.0);
    base.setIsentropicEfficiency(0.80);
    base.run();
    double basePower = base.getPower();
    System.out.println("base_power_w=" + basePower);
    System.out.println("base_outlet_t_k=" + base.getOutletStream().getTemperature("K"));
    print("inlet", baseInlet);
    print("outlet", base.getOutletStream());
    System.out.println();

    for (double fraction : new double[] {1.0, 0.5, 0.0, 0.25}) {
      Stream inlet = expanderFeed();
      neqsim.process.equipment.expander.Expander expander =
          new neqsim.process.equipment.expander.Expander("set power", inlet);
      expander.setOutletPressure(30.0);
      expander.setIsentropicEfficiency(0.80);
      expander.setPower(basePower * fraction);
      expander.run();
      System.out.println("power_set_fraction=" + fraction);
      System.out.println("power_w=" + expander.getPower());
      System.out.println("isentropic_efficiency=" + expander.getIsentropicEfficiency());
      System.out.println("outlet_t_k=" + expander.getOutletStream().getTemperature("K"));
      print("inlet", inlet);
      print("outlet", expander.getOutletStream());
      System.out.println();
    }

    Stream warmInlet = expanderFeed();
    neqsim.process.equipment.expander.Expander warm =
        new neqsim.process.equipment.expander.Expander("warm", warmInlet);
    warm.setOutletPressure(30.0);
    warm.setIsentropicEfficiency(0.80);
    warm.setPower(-basePower * 0.1);
    warm.run();
    System.out.println("power_set_fraction=-0.1");
    System.out.println("power_w=" + warm.getPower());
    System.out.println("isentropic_efficiency=" + warm.getIsentropicEfficiency());
    System.out.println("outlet_t_k=" + warm.getOutletStream().getTemperature("K"));
    print("inlet", warmInlet);
    print("outlet", warm.getOutletStream());
    System.out.println();

    for (double target : new double[] {280.0, 250.0, 240.0}) {
      Stream inlet = expanderFeed();
      neqsim.process.equipment.expander.Expander expander =
          new neqsim.process.equipment.expander.Expander("solve", inlet);
      expander.setOutletPressure(30.0);
      expander.setIsentropicEfficiency(0.80);
      double solved = expander.solveEfficiency(target);
      expander.run();
      System.out.println("solve_target_t_k=" + target);
      System.out.println("solved_efficiency=" + solved);
      System.out.println("power_w=" + expander.getPower());
      System.out.println("outlet_t_k=" + expander.getOutletStream().getTemperature("K"));
      print("inlet", inlet);
      print("outlet", expander.getOutletStream());
      System.out.println();
    }
  }

  static Stream expanderFeed() {
    SystemInterface fluid = new SystemPrEos(298.15, 90.0);
    fluid.addComponent("methane", 0.82);
    fluid.addComponent("ethane", 0.10);
    fluid.addComponent("propane", 0.05);
    fluid.addComponent("n-butane", 0.03);
    fluid.setMixingRule(2);
    fluid.setTotalFlowRate(50000.0, "kg/hr");
    Stream stream = new Stream("feed", fluid);
    stream.setTemperature(25.0, "C");
    stream.setPressure(90.0, "bara");
    stream.run();
    return stream;
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
  /// `CompressorDeposit`'s mass-to-multiplier mapping and its run-time effect.
  static void compressorDepositRows() {
    depositRow("clean", 0.0, 0.0);
    depositRow("sulfur5_salt2", 5.0, 2.0);
    depositRow("saturated", 200.0, 0.0);
    // The class's own inverse: the S8 mass that yields a 0.77 efficiency multiplier.
    neqsim.process.equipment.compressor.CompressorDeposit inverse =
        new neqsim.process.equipment.compressor.CompressorDeposit(0.6, 0.02);
    double inverseMass = inverse.depositMassForEfficiencyMultiplier(
        neqsim.process.equipment.compressor.DepositMechanism.SULFUR_S8, 0.77);
    depositRow("inverse_0_77", inverseMass, 0.0);

    // The machine: clean against 3 kg of S8, both on the deposit test's polytropic machine.
    System.out.println("machine_clean");
    printMachineDeposit(null);
    System.out.println("machine_fouled_3kg");
    printMachineDeposit(3.0);
  }

  static void depositRow(String label, double sulfurKg, double saltKg) {
    neqsim.process.equipment.compressor.CompressorDeposit deposit =
        new neqsim.process.equipment.compressor.CompressorDeposit(0.6, 0.02);
    deposit.addDeposit(neqsim.process.equipment.compressor.DepositMechanism.SULFUR_S8, sulfurKg);
    deposit.addDeposit(neqsim.process.equipment.compressor.DepositMechanism.SALT_NACL, saltKg);
    System.out.println(label);
    System.out.println("total_mass_kg=" + deposit.getTotalDepositMass());
    System.out.println("volume_m3=" + deposit.getDepositVolume());
    System.out.println("film_thickness_m=" + deposit.getFilmThickness());
    System.out.println("blockage=" + deposit.getAreaBlockageFraction());
    System.out.println("efficiency_multiplier=" + deposit.getEfficiencyMultiplier());
    System.out.println("head_multiplier=" + deposit.getHeadMultiplier());
    System.out.println();
  }

  static void printMachineDeposit(Double sulfurKg) {
    neqsim.thermo.system.SystemInterface gas =
        new neqsim.thermo.system.SystemPrEos(273.15 + 30.0, 18.8);
    gas.addComponent("methane", 0.86);
    gas.addComponent("ethane", 0.07);
    gas.addComponent("propane", 0.03);
    gas.addComponent("CO2", 0.02);
    gas.addComponent("nitrogen", 0.02);
    gas.setMixingRule(2);
    Stream feed = new Stream("suction", gas);
    feed.setFlowRate(100000.0, "kg/hr");
    feed.setTemperature(30.0, "C");
    feed.setPressure(18.8, "bara");
    feed.run();
    neqsim.process.equipment.compressor.Compressor compressor =
        new neqsim.process.equipment.compressor.Compressor("KA03B", feed);
    compressor.setOutletPressure(42.0);
    compressor.setUsePolytropicCalc(true);
    compressor.setPolytropicEfficiency(0.60);
    if (sulfurKg != null) {
      neqsim.process.equipment.compressor.CompressorDeposit deposit =
          neqsim.process.equipment.compressor.CompressorDeposit.fromCompressor(compressor);
      deposit.addDeposit(neqsim.process.equipment.compressor.DepositMechanism.SULFUR_S8, sulfurKg);
      compressor.setDepositModel(deposit);
      System.out.println("inlet_volume_flow_m3_s=" + feed.getFlowRate("m3/sec"));
      System.out.println("wetted_area_m2=" + deposit.getWettedArea());
      System.out.println("passage_half_height_m=" + deposit.getPassageHalfHeight());
    }
    compressor.run();
    System.out.println("power_kw=" + compressor.getPower() / 1000.0);
    System.out.println("polytropic_efficiency_out=" + compressor.getPolytropicEfficiency());
    System.out.println("degradation_factor=" + compressor.getDegradationFactor());
    System.out.println("fouling_factor=" + compressor.getFoulingFactor());
    System.out.println();
  }

  /// `setCompressionRatio` and `setPower` on the deposit test's machine.
  static void compressorModeRows() {
    compressorModeRow("pr_ratio_2_2", null, 2.2, null, null, null, false, null);
    compressorModeRow("pr_ratio_3_capped_40", null, 3.0, 40.0, null, null, false, null);
    compressorModeRow("pr_power_3mw", 42.0, null, null, 3.0e6, 0.8, false, null);
    compressorModeRow("pr_power_3mw_polytropic", 42.0, null, null, 3.0e6, 0.8, true, 0.75);
  }

  static void compressorModeRow(String label, Double outletBara, Double ratioBara,
      Double maximumBara, Double powerW, Double isentropicEfficiency, boolean polytropic,
      Double polytropicEfficiency) {
    neqsim.thermo.system.SystemInterface gas =
        new neqsim.thermo.system.SystemPrEos(273.15 + 30.0, 18.8);
    gas.addComponent("methane", 0.86);
    gas.addComponent("ethane", 0.07);
    gas.addComponent("propane", 0.03);
    gas.addComponent("CO2", 0.02);
    gas.addComponent("nitrogen", 0.02);
    gas.setMixingRule(2);
    Stream feed = new Stream("suction", gas);
    feed.setFlowRate(100000.0, "kg/hr");
    feed.setTemperature(30.0, "C");
    feed.setPressure(18.8, "bara");
    feed.run();
    neqsim.process.equipment.compressor.Compressor compressor =
        new neqsim.process.equipment.compressor.Compressor(label, feed);
    if (outletBara != null) {
      compressor.setOutletPressure(outletBara);
    }
    if (ratioBara != null) {
      compressor.setCompressionRatio(ratioBara);
    }
    if (maximumBara != null) {
      compressor.setMaxOutletPressure(maximumBara);
    }
    if (powerW != null) {
      compressor.setPower(powerW);
    }
    if (polytropic) {
      compressor.setUsePolytropicCalc(true);
      compressor.setPolytropicEfficiency(polytropicEfficiency);
    } else if (isentropicEfficiency != null) {
      compressor.setIsentropicEfficiency(isentropicEfficiency);
    }
    compressor.run();

    System.out.println(label);
    print("inlet", feed);
    print("outlet", compressor.getOutletStream());
    System.out.println("power_kw=" + compressor.getPower() / 1000.0);
    System.out.println("isentropic_efficiency_out=" + compressor.getIsentropicEfficiency());
    System.out.println("polytropic_efficiency_out=" + compressor.getPolytropicEfficiency());
    System.out.println();
  }

  /// `Compressor.run`'s no-chart polytropic branch on `CompressorDepositTest`'s machine.
  static void compressorPolytropicRows() {
    compressorPolytropicRow("srk_eta_0_60", true, 0.60);
    compressorPolytropicRow("pr_eta_0_60", false, 0.60);
    compressorPolytropicRow("pr_eta_0_75", false, 0.75);
  }

  /// `CompressorMechanicalDesign.calcDesign`, on the state families it branches on.
  ///
  /// **The no-chart polytropic branch is the class's default method**, and its fluid head
  /// is the value the stage count and the impeller screen read. The three outlet pressures
  /// cross the class's casing thresholds - `20` bara stays vertical, `80` is horizontal and
  /// `150` is a barrel - and the NACE row turns the sour-service checks on.
  static void compressorMechanicalRows() {
    compressorMechanicalRow("polytropic_pr_eta_0_60", 18.8, 42.0, 0.60, false);
    compressorMechanicalRow("polytropic_pr_eta_0_75", 18.8, 42.0, 0.75, false);
    compressorMechanicalRow("vertically_split", 10.0, 20.0, 0.75, false);
    compressorMechanicalRow("horizontally_split", 50.0, 80.0, 0.75, false);
    compressorMechanicalRow("barrel", 30.0, 150.0, 0.75, false);
    compressorMechanicalRow("nace", 18.8, 42.0, 0.75, true);
    compressorMechanicalChartRow("chart", 3.0);
    compressorMechanicalChartRow("chart_1msm3", 1.0);
  }

  static void compressorMechanicalRow(String label, double inletBara, double outletBara,
      double efficiency, boolean nace) {
    SystemInterface gas = new SystemPrEos(273.15 + 30.0, inletBara);
    gas.addComponent("methane", 0.86);
    gas.addComponent("ethane", 0.07);
    gas.addComponent("propane", 0.03);
    gas.addComponent("CO2", 0.02);
    gas.addComponent("nitrogen", 0.02);
    gas.setMixingRule(2);
    Stream feed = new Stream("suction", gas);
    feed.setFlowRate(100000.0, "kg/hr");
    feed.setTemperature(30.0, "C");
    feed.setPressure(inletBara, "bara");
    feed.run();
    neqsim.process.equipment.compressor.Compressor compressor =
        new neqsim.process.equipment.compressor.Compressor(label, feed);
    compressor.setOutletPressure(outletBara);
    compressor.setUsePolytropicCalc(true);
    compressor.setPolytropicEfficiency(efficiency);
    compressor.run();
    neqsim.process.mechanicaldesign.compressor.CompressorMechanicalDesign design =
        compressor.getMechanicalDesign();
    if (nace) {
      design.setNaceCompliance(true);
      design.setH2sPartialPressureKPa(10.0);
    }
    design.calcDesign();
    printCompressorDesign(label, feed, compressor, design);
  }

  static void compressorMechanicalChartRow(String label, double msm3PerDay) {
    neqsim.process.equipment.compressor.Compressor compressor = antiSurgeMachine(msm3PerDay);
    compressor.run();
    neqsim.process.mechanicaldesign.compressor.CompressorMechanicalDesign design =
        compressor.getMechanicalDesign();
    design.calcDesign();
    printCompressorDesign(label, compressor.getInletStream(), compressor, design);
  }

  static void printCompressorDesign(String label, StreamInterface feed,
      neqsim.process.equipment.compressor.Compressor compressor,
      neqsim.process.mechanicaldesign.compressor.CompressorMechanicalDesign design) {
    System.out.println(label);
    System.out.println("speed_rpm=" + compressor.getSpeed());
    System.out.println("power_kW=" + compressor.getPower("kW"));
    System.out.println("polytropic_fluid_head_kj_per_kg=" + compressor.getPolytropicFluidHead());
    System.out.println("polytropic_efficiency=" + compressor.getPolytropicEfficiency());
    System.out.println("number_of_stages=" + design.getNumberOfStages());
    System.out.println("impeller_diameter_mm=" + design.getImpellerDiameter());
    System.out.println("shaft_diameter_mm=" + design.getShaftDiameter());
    System.out.println("tip_speed_m_per_s=" + design.getTipSpeed());
    System.out.println("flow_coefficient=" + design.getFlowCoefficient());
    System.out.println("impeller_sizing_speed_rpm=" + design.getImpellerSizingSpeedRPM());
    System.out.println("driver_power_kW=" + design.getDriverPower());
    System.out.println("driver_margin=" + design.getDriverMargin());
    System.out.println("design_pressure_bara=" + design.getDesignPressure());
    System.out.println("design_temperature_C=" + design.getDesignTemperature());
    System.out.println("max_continuous_speed_rpm=" + design.getMaxContinuousSpeed());
    System.out.println("trip_speed_rpm=" + design.getTripSpeed());
    System.out.println("first_critical_speed_rpm=" + design.getFirstCriticalSpeed());
    System.out.println("casing_type=" + design.getCasingType());
    System.out.println("head_per_stage_kj_per_kg=" + design.getHeadPerStage());
    System.out.println("bearing_span_mm=" + design.getBearingSpan());
    System.out.println("casing_weight_kg=" + design.getCasingWeight());
    System.out.println("bundle_weight_kg=" + design.getBundleWeight());
    System.out.println("rotor_weight_kg=" + design.getRotorWeight());
    System.out.println("outer_diameter_m=" + design.getOuterDiameter());
    System.out.println("inner_diameter_m=" + design.getInnerDiameter());
    System.out.println("wall_thickness_m=" + design.getWallThickness());
    System.out.println("tantan_length_m=" + design.getTantanLength());
    System.out.println("weight_vessel_shell_kg=" + design.getWeigthVesselShell());
    System.out.println("weight_internals_kg=" + design.getWeigthInternals());
    System.out.println("weight_nozzle_kg=" + design.getWeightNozzle());
    System.out.println("weight_piping_kg=" + design.getWeightPiping());
    System.out.println("weight_electro_instrument_kg=" + design.getWeightElectroInstrument());
    System.out.println("weight_structural_steel_kg=" + design.getWeightStructualSteel());
    System.out.println("weight_total_kg=" + design.getWeightTotal());
    System.out.println("module_length_m=" + design.getModuleLength());
    System.out.println("module_width_m=" + design.getModuleWidth());
    System.out.println("module_height_m=" + design.getModuleHeight());
    System.out.println("seal_type=" + design.getSealType());
    System.out.println("bearing_type=" + design.getBearingType());
    System.out.println("material_class=" + design.getMaterialClass());
    System.out.println("casing_material_grade=" + design.getCasingMaterialGrade());
    System.out.println("casing_corrosion_allowance_mm=" + design.getCasingCorrosionAllowanceMm());
    System.out.println("h2s_partial_pressure_kpa=" + design.getH2sPartialPressureKPa());
    System.out.println("max_stages_per_casing=" + design.getMaxStagesPerCasing());
    System.out.println("impeller_sizing_issues=" + design.getImpellerSizingIssues());
    System.out.println("inlet_m3_per_hr=" + feed.getFlowRate("m3/hr"));
    System.out.println("inlet_density_kg_per_m3=" + feed.getThermoSystem().getDensity("kg/m3"));
    System.out.println("mass_flow_kg_per_hr=" + feed.getFlowRate("kg/hr"));
    System.out.println("inlet_number_of_phases=" + feed.getThermoSystem().getNumberOfPhases());
    print("inlet", feed);
    print("outlet", compressor.getOutletStream());
    System.out.println();
  }

  static void compressorPolytropicRow(String label, boolean srk, double efficiency) {
    neqsim.thermo.system.SystemInterface gas = srk
        ? new neqsim.thermo.system.SystemSrkEos(273.15 + 30.0, 18.8)
        : new neqsim.thermo.system.SystemPrEos(273.15 + 30.0, 18.8);
    gas.addComponent("methane", 0.86);
    gas.addComponent("ethane", 0.07);
    gas.addComponent("propane", 0.03);
    gas.addComponent("CO2", 0.02);
    gas.addComponent("nitrogen", 0.02);
    gas.setMixingRule(2);
    Stream feed = new Stream("suction", gas);
    feed.setFlowRate(100000.0, "kg/hr");
    feed.setTemperature(30.0, "C");
    feed.setPressure(18.8, "bara");
    feed.run();
    neqsim.process.equipment.compressor.Compressor compressor =
        new neqsim.process.equipment.compressor.Compressor(label, feed);
    compressor.setOutletPressure(42.0);
    compressor.setUsePolytropicCalc(true);
    compressor.setPolytropicEfficiency(efficiency);
    compressor.run();

    System.out.println(label);
    System.out.println("eos=" + (srk ? "srk" : "pr"));
    System.out.println("polytropic_efficiency_in=" + efficiency);
    print("inlet", feed);
    print("outlet", compressor.getOutletStream());
    System.out.println("power_kw=" + compressor.getPower() / 1000.0);
    System.out.println("polytropic_efficiency_out=" + compressor.getPolytropicEfficiency());
    System.out.println("polytropic_exponent=" + compressor.getPolytropicExponent());
    // The branch's own intermediates, reproduced so the capture shows what it divides by.
    neqsim.thermo.system.SystemInterface inletSystem = feed.getThermoSystem();
    double densIn = inletSystem.getDensity();
    double pInBara = inletSystem.getPressure();
    double entropyIn = inletSystem.getEntropy();
    neqsim.thermo.system.SystemInterface isentropic = inletSystem.clone();
    isentropic.setPressure(42.0, "bara");
    new neqsim.thermodynamicoperations.ThermodynamicOperations(isentropic).PSflash(entropyIn);
    isentropic.initProperties();
    double densOut = isentropic.getDensity();
    double hOutIsentropic = isentropic.getEnthalpy();
    double hIn = inletSystem.getEnthalpy();
    double k = Math.log(42.0 / pInBara) / Math.log(densOut / densIn);
    System.out.println("inlet_density_kg_m3=" + densIn);
    System.out.println("isentropic_density_kg_m3=" + densOut);
    System.out.println("inlet_enthalpy_total=" + hIn);
    System.out.println("isentropic_enthalpy_total=" + hOutIsentropic);
    System.out.println("isenthalpic_exponent=" + k);
    System.out.println();
  }

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

  /// `ExpanderMechanicalDesign.calcDesign`, on the states its branches divide.
  ///
  /// **The small gas state reaches the radial branch and the oil-brake load** (recovered power
  /// below the class's 100 kW threshold); the large-flow state reaches `AXIAL`, the multi-stage
  /// head division and the generator. The two temperature states reach the casing's 0.85 and 0.9
  /// allowable-stress derates. `MIXED_FLOW` and `LoadType.COMPRESSOR` are never selected and are
  /// stated rather than cased.
  ///
  /// **The base fields are millimetres here**, unlike the pump's metres: the row prints them in
  /// the class's own units and the spec converts each at its boundary.
  static void expanderMechanicalRows() {
    expanderMechanicalRow("radial_brake", 320.0, 60.0, 1.0, 30.0, 0.75);
    expanderMechanicalRow("radial_eta_1_0", 320.0, 60.0, 1.0, 30.0, 1.0);
    expanderMechanicalRow("axial_multistage", 320.0, 60.0, 250000.0, 30.0, 0.75);
    expanderMechanicalRow("hot_casing", 500.0, 30.0, 10.0, 5.0, 0.8);
    expanderMechanicalRow("cold_casing", 190.0, 30.0, 10.0, 5.0, 0.8);
  }

  static void expanderMechanicalRow(String label, double temperatureK, double inletBara,
      double molPerSecond, double outletBara, double efficiency) {
    Stream inlet = feed(new String[] {"methane", "n-butane"}, new double[] {0.9, 0.1},
        temperatureK, inletBara, molPerSecond);
    neqsim.process.equipment.expander.Expander expander =
        new neqsim.process.equipment.expander.Expander(label, inlet);
    expander.setOutletPressure(outletBara);
    expander.setIsentropicEfficiency(efficiency);
    expander.run();
    neqsim.process.mechanicaldesign.expander.ExpanderMechanicalDesign design =
        expander.getExpanderMechanicalDesign();
    design.calcDesign();

    System.out.println(label);
    System.out.println("expander_type=" + design.getExpanderType());
    System.out.println("load_type=" + design.getLoadType());
    System.out.println("number_of_stages=" + design.getNumberOfStages());
    System.out.println("wheel_diameter_mm=" + design.getWheelDiameter());
    System.out.println("tip_speed_m_per_s=" + design.getTipSpeed());
    System.out.println("shaft_diameter_mm=" + design.getShaftDiameter());
    System.out.println("rated_speed_rpm=" + design.getRatedSpeed());
    System.out.println("max_continuous_speed_rpm=" + privateDouble(design, "maxContinuousSpeed"));
    System.out.println("trip_speed_rpm=" + privateDouble(design, "tripSpeed"));
    System.out.println("first_critical_speed_rpm=" + design.getFirstCriticalSpeed());
    System.out.println("effective_isentropic_efficiency=" + expander.getIsentropicEfficiency());
    System.out.println("casing_pressure_bara=" + design.getCasingDesignPressure());
    System.out.println("casing_temperature_C=" + design.getCasingDesignTemperature());
    System.out.println("casing_wall_thickness_mm=" + privateDouble(design, "casingWallThickness"));
    System.out.println("bearing_span_mm=" + privateDouble(design, "bearingSpan"));
    System.out.println("wheel_weight_kg=" + privateDouble(design, "wheelWeight"));
    System.out.println("shaft_weight_kg=" + privateDouble(design, "shaftWeight"));
    System.out.println("casing_weight_kg=" + privateDouble(design, "casingWeight"));
    System.out.println("generator_weight_kg=" + privateDouble(design, "generatorWeight"));
    System.out.println("seal_type=" + design.getSealType());
    System.out.println("bearing_type=" + design.getBearingType());
    System.out.println("inner_diameter_mm=" + design.getInnerDiameter());
    System.out.println("outer_diameter_mm=" + design.getOuterDiameter());
    System.out.println("wall_thickness_mm=" + design.getWallThickness());
    System.out.println("tantan_length_mm=" + design.getTantanLength());
    System.out.println("weight_vessel_kg=" + design.getWeightVessel());
    System.out.println("weight_internals_kg=" + design.getWeigthInternals());
    System.out.println("weight_piping_kg=" + design.getWeightPiping());
    System.out.println("weight_electro_instrument_kg=" + design.getWeightElectroInstrument());
    System.out.println("weight_structural_steel_kg=" + design.getWeightStructualSteel());
    System.out.println("weight_total_kg=" + design.getWeightTotal());
    System.out.println("module_length_m=" + design.getModuleLength());
    System.out.println("module_width_m=" + design.getModuleWidth());
    System.out.println("module_height_m=" + design.getModuleHeight());

    System.out.println("inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("inlet_density_kg_per_m3=" + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("mass_flow_kg_per_hr=" + inlet.getFlowRate("kg/hr"));
    System.out.println("inlet_number_of_phases=" + inlet.getThermoSystem().getNumberOfPhases());
    print("inlet", inlet);
    print("outlet", expander.getOutletStream());
    System.out.println("power_kW=" + expander.getPower("kW"));
    System.out.println();
  }

  /// **`Expander`'s recovered-power constraint**, the one capacity limit the expander adds.
  ///
  /// The parent's consumed-power constraints are removed and, when a rating is set, one HARD
  /// `recoveredPower` constraint is added: design = the rating, max = 1.1 x the rating, value =
  /// `|getPower("kW")|`. `isCapacityExceeded` is utilization above 1, so the recovered power is
  /// over its rating; `isHardLimitExceeded` is the same value above the 1.1 x rating, a
  /// narrower band than the utilization's own 1.0.
  static void expanderCapacityRows() {
    String[] names = new String[] {"methane", "n-butane"};
    double[] z = new double[] {0.9, 0.1};
    expanderCapacityRow("no_rating", names, z, null);
    expanderCapacityRow("rating_2_kw", names, z, 2.0);
    expanderCapacityRow("rating_1_2_kw", names, z, 1.2);
    expanderCapacityRow("rating_1_kw", names, z, 1.0);
    expanderCapacityRow("rating_zero", names, z, 0.0);
  }

  static void expanderCapacityRow(String label, String[] names, double[] z, Double ratingKW) {
    Stream inlet = feed(names, z, 320.0, 60.0, 1.0);
    neqsim.process.equipment.expander.Expander expander =
        new neqsim.process.equipment.expander.Expander("e1", inlet);
    expander.setOutletPressure(30.0);
    expander.setIsentropicEfficiency(0.75);
    if (ratingKW != null) {
      expander.setRatedRecoveredPower(ratingKW);
    }
    expander.run();

    System.out.println("case=" + label);
    print("inlet", inlet);
    print("outlet", expander.getOutletStream());
    System.out.println("power_kW=" + expander.getPower("kW"));
    System.out.println("recovered_power_kW=" + Math.abs(expander.getPower("kW")));
    System.out.println("rated_recovered_power_kW=" + expander.getRatedRecoveredPower());

    Map<String, neqsim.process.equipment.capacity.CapacityConstraint> constraints =
        expander.getCapacityConstraints();
    System.out.println("constraint_count=" + constraints.size());
    for (Map.Entry<String, neqsim.process.equipment.capacity.CapacityConstraint> entry :
        constraints.entrySet()) {
      String name = entry.getKey();
      neqsim.process.equipment.capacity.CapacityConstraint constraint = entry.getValue();
      System.out.println("constraint_" + name + "_type=" + constraint.getType());
      System.out.println("constraint_" + name + "_enabled=" + constraint.isEnabled());
      System.out.println("constraint_" + name + "_current=" + constraint.getCurrentValue());
      System.out.println("constraint_" + name + "_design=" + constraint.getDesignValue());
      System.out.println("constraint_" + name + "_min=" + constraint.getMinValue());
      System.out.println("constraint_" + name + "_max=" + constraint.getMaxValue());
      System.out.println("constraint_" + name + "_utilization=" + constraint.getUtilization());
      System.out.println("constraint_" + name + "_violated=" + constraint.isViolated());
      System.out.println(
          "constraint_" + name + "_hard_exceeded=" + constraint.isHardLimitExceeded());
    }
    System.out.println("max_utilization=" + expander.getMaxUtilization());
    System.out.println("capacity_exceeded=" + expander.isCapacityExceeded());
    System.out.println("hard_limit_exceeded=" + expander.isHardLimitExceeded());
    neqsim.process.equipment.capacity.CapacityConstraint bottleneck = expander
        .getBottleneckConstraint();
    System.out.println(
        "bottleneck_constraint=" + (bottleneck == null ? "none" : bottleneck.getName()));
    System.out.println();
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
  /// **The last three rows are the metric**, and it is what the port leaves out: `d` and
  /// `K` separately, and then the same vessel with entrainment on, because `getCapacityUtilization`
  /// reads the thermo system `run` already entrained rather than the flash it started from.
  /// `getCapacityUtilization` needs an internal diameter and a design gas load factor; the
  /// palette entry declares neither, and it is a statement about whether the *vessel* is
  /// big enough rather than about the stream.
  static void gasScrubberRows() {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.7, 0.3 };
    scrubberRow("pressure_drop_bara=0 heat_input_W=0 gas_in_liquid=0", names, z, 0.0, null, 0.0,
        null, null);
    scrubberRow("pressure_drop_bara=2 heat_input_W=0 gas_in_liquid=0", names, z, 2.0, null, 0.0,
        null, null);
    scrubberRow("pressure_drop_bara=0 heat_input_W=0 gas_in_liquid=0.05", names, z, 0.0, null, 0.05,
        null, null);
    scrubberRow("capacity_utilization_dn1000_k_007", names, z, 0.0, null, 0.0, 1.0, 0.07);
    scrubberRow("capacity_utilization_dn2000_k_007", names, z, 0.0, null, 0.0, 2.0, 0.07);
    scrubberRow("capacity_utilization_dn1000_k_007_gas_in_liquid_005", names, z, 0.0, null, 0.05,
        1.0, 0.07);
  }

  static void scrubberRow(String label, String[] names, double[] z, double dropBara,
      Double heatInputW, double gasInLiquid, Double diameterM, Double kFactor) {
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
    if (diameterM != null) {
      // The two mechanical parameters the metric needs, and neither is on the palette.
      scrubber.setInternalDiameter(diameterM);
      scrubber.setDesignGasLoadFactor(kFactor);
    }
    scrubber.run();

    System.out.println(label);
    print("feed", inlet);
    print("vapour", scrubber.getGasOutStream());
    print("liquid", scrubber.getLiquidOutStream());
    if (diameterM != null) {
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
    // **The second solve on the same column**, so the column's own capture carries the state
    // the port's mesh solve is held to beside the substitution core's - the two are the pair
    // this tranche pins. The ladder capture holds the same row among the ten strategies.
    solverRow("binary_methane_butane_mesh_solve", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0, 19.0, 20.0, 1.0e-6,
        200, true, false, "NAPHTALI_SANDHOLM");
  
    // The reactive section's rows, which belong with the state they are compared
    // against - see the helper's own note on why they live in this capture.
    reactiveColumnRowsInColumnCapture();

    // The side-draw column rows, which live here for the same reason.
    sideDrawColumnRowsInColumnCapture();

    // The Murphree rows, which live here for the same reason again.
    murphreeRows();

    // **The pumparound's *return*: a liquid draw that comes back to another tray.**
    pumparoundReturnRows();

    // **The draw's *flow* specification, which is a tear loop and not a split.** Every row above
    // states a fraction; `addSideDrawFlowSpecification` states a target *flow* and lets the column
    // move that fraction until the draw delivers it. It is last so that the blocks the cases
    // already name keep their positions.
    sideDrawFlowRows();
  }

  /// **The Murphree efficiency, on the column the port is already held to.** Two efficiencies,
  /// `0.6` and `0.85`, which are the two `DistillationSolverBenchmarkTest` states its own
  /// product-split and multi-solver tests use - and the ideal reference is the
  /// `binary_methane_butane_4_stages` row above, so a reader takes the difference from a state
  /// both libraries already agree on.
  ///
  /// **The correction bites on four of the six trays**, which the profile shows: stage 0 is the
  /// reboiler and stage 5 the condenser, and `applyMurphreeCorrection` returns immediately on
  /// each. The two tray rows that move most are the ones whose vapour is furthest from the
  /// tray below's, which is the definition of the correction rather than a second effect.
  static void murphreeRows() {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.5, 0.5 };
    murphreeRow("binary_murphree_0_6", names, z, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0, 19.0,
        20.0, 1.0e-6, 200, true, false, 0.6);
    murphreeRow("binary_murphree_0_85", names, z, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0, 19.0,
        20.0, 1.0e-6, 200, true, false, 0.85);
    // **The same state at a tighter gate**, which is the row that says whether the `1e-6` one is
    // a *converged* state: the class stops when the mean tray-temperature change falls under the
    // tolerance, and on a corrected column the temperatures settle while the compositions are
    // still moving. If the profile below differs from the one above, the one above was a partial
    // state on the way and not the column's answer.
    murphreeRow("binary_murphree_0_85_tight", names, z, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0,
        19.0, 20.0, 1.0e-9, 2000, true, false, 0.85);
    murphreeRow("binary_murphree_0_6_tight", names, z, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0,
        19.0, 20.0, 1.0e-9, 2000, true, false, 0.6);
  }

  /// **The divergence instrument, and it prints two things the column capture does not.**
  ///
  /// `process_column.tsv` carries each tray's temperature, pressure and two traffic rates, and
  /// **no composition at all** - which is enough to say that the port's `0.85` state differs by
  /// `7` K and not enough to say where. This prints, per tray, the flashed system's *both* phase
  /// compositions (the correction's `y_eq` and its neighbour's `y_in`), the **cached** vapour the
  /// tray hands up, the phase count and the phase-0 type; and per pass, the class's own
  /// `getConvergenceHistory()`.
  ///
  /// **Why those four and no more.** `applyMurphreeCorrection` reads `getThermoSystem().
  /// getPhase(0)` for `y_eq` and the tray below's for `y_in`, and `setCachedGasOutStream` writes
  /// the corrected stream into a *separate* object - so a state's map is fixed by exactly these
  /// three vectors per tray. A warm start built from them, run one pass, says whether NeqSim's
  /// endpoint is a fixed point of this port's arithmetic; the history says whether the two
  /// reached it by the same path.
  static void divergenceRows() {
    divergenceRow("divergence_0_6", 0.6);
    divergenceRow("divergence_0_85", 0.85);
  }

  /// One row of [`divergenceRows`]: [`runColumn`]'s own state with the divergence block turned
  /// on, so the two captures cannot drift apart in how the column was built.
  static void divergenceRow(String label, double efficiency) {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.5, 0.5 };
    runColumn(label, names, z, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0, 19.0, 20.0, 1.0e-6, 200,
        true, false, "DIRECT_SUBSTITUTION", efficiency, -1, 0.0, null, true);
  }

  /// The two things `process_column.tsv` does not carry: the class's per-pass
  /// `getConvergenceHistory()`, and per tray the *flashed* system's two phase compositions
  /// beside the **cached** vapour the tray hands up.
  ///
  /// **Why those and no more.** `applyMurphreeCorrection` reads `getThermoSystem().getPhase(0)`
  /// for `y_eq` and the tray below's for `y_in`, and `setCachedGasOutStream` writes the corrected
  /// stream into a *separate* object - so a state's map is fixed by exactly these vectors, and a
  /// warm start built from them answers whether NeqSim's endpoint is a fixed point of this
  /// port's arithmetic. The history is the cheaper question beside it: whether the two reached
  /// their states by the same path.
  static void divergenceBlock(neqsim.process.equipment.distillation.DistillationColumn column) {
    int pass = 0;
    for (double[] entry : column.getConvergenceHistory()) {
      StringBuilder line = new StringBuilder("pass" + pass + "=");
      for (int k = 0; k < entry.length; k++) {
        if (k > 0) {
          line.append(' ');
        }
        line.append(entry[k]);
      }
      System.out.println(line);
      pass++;
    }
    System.out.println("history_length=" + pass);
    for (int stage = 0; stage < column.getNumberOfTrays(); stage++) {
      String prefix = "tray" + stage + "_";
      SystemInterface system = column.getTray(stage).getThermoSystem();
      System.out.println(prefix + "phases=" + system.getNumberOfPhases());
      System.out.println(prefix + "phase0type=" + system.getPhase(0).getType());
      System.out.println(prefix + "eq_gas_z=" + composition(system.getPhase(0)));
      if (system.getNumberOfPhases() > 1) {
        System.out.println(prefix + "eq_liquid_z=" + composition(system.getPhase(1)));
      }
      StreamInterface handedUp = column.getTray(stage).getGasOutStream();
      if (handedUp != null) {
        System.out.println(prefix + "out_gas_z=" + composition(handedUp.getThermoSystem().getPhase(0)));
        System.out.println(prefix + "out_gas_n=" + handedUp.getFlowRate("mol/sec"));
      }
    }
  }

  /// One phase's composition as `name:value` pairs, in the phase's own component order.
  static String composition(neqsim.thermo.phase.PhaseInterface phase) {
    StringBuilder out = new StringBuilder();
    for (int j = 0; j < phase.getNumberOfComponents(); j++) {
      if (j > 0) {
        out.append(' ');
      }
      out.append(phase.getComponent(j).getName()).append(':').append(phase.getComponent(j).getx());
    }
    return out.toString();
  }

  /// **The per-stage efficiency vector, on that same column.** `getEffectiveMurphreeEfficiency`
  /// resolves two fields - `perStageMurphreeEfficiency[stage]` when it is present and not `NaN`,
  /// `murphreeEfficiency` otherwise - and three rows measure the whole rule: the column-wide
  /// value alone, the same value with one stage overridden through the single-stage setter, and
  /// the same override again spelled as the array. **The second and third have to be one state**
  /// or `NaN` is not the fall-through its javadoc says it is, and both have to differ from the
  /// first or the override did nothing.
  ///
  /// **The array does not carry the column-wide value.** `setMurphreeEfficiencies` writes only
  /// `perStageMurphreeEfficiency`, so both override rows state `0.6` first and then the
  /// override - which is the mixed case: one stage off the default, the rest on it.
  static void murphreeStageRows() {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.5, 0.5 };
    murphreeRow("efficiency_column_wide", names, z, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0,
        19.0, 20.0, 1.0e-6, 200, true, false, 0.6);
    // `setMurphreeEfficiency(3, 0.85)` after the column-wide `0.6`.
    murphreeStageRow("efficiency_one_stage", names, z, 3, 0.85, null);
    // The same override as the array `setMurphreeEfficiencies` takes, `NaN` elsewhere - which is
    // the shape the model's own `tray_murphree_efficiency` declares.
    murphreeStageRow("efficiency_one_stage_array", names, z, 3, 0.85,
        new double[] { Double.NaN, Double.NaN, Double.NaN, 0.85, Double.NaN, Double.NaN });
  }

  /// One row of [`murphreeStageRows`]: `0.6` column-wide, then the override through whichever
  /// setter `array` names - the single-stage one when it is `null`, and `setMurphreeEfficiencies`
  /// with `stage` and `efficiency` left unread when it is not.
  static void murphreeStageRow(String label, String[] names, double[] z, int stage,
      double efficiency, double[] array) {
    runColumn(label, names, z, 300.0, 20.0, 1000.0, 4, 2, -20.0, 100.0, 19.0, 20.0, 1.0e-6, 200,
        true, false, "DIRECT_SUBSTITUTION", 0.6, stage, efficiency, array, false);
  }

  /// **The side-draw column rows, and why they are in this capture rather than one of their own.**
  /// The layer diff pairs one capture with one model - `layers.LAYER_CASES` is grouped by model and
  /// `test_a_capture_has_one_block_per_case` asserts each model reads exactly one capture - so a
  /// side-draw state on `process.distillation_column` has to sit in the capture whose blocks that
  /// model's cases already claim. `ProcessProbe side_draw` writes the same rows to
  /// `captures/process_side_draw.tsv` as the subject's own record, the shape the reactive column's
  /// rows already have.
  static void sideDrawColumnRowsInColumnCapture() {
    sideDrawColumnRows();
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


  /// **The absorber and the stripper, which are one machine with two names.**
  /// `AbsorptionColumn extends DistillationColumn` with `super(name, trays, false, false)`, so
  /// there is no condenser and no reboiler, and its two inlets are fixed: the gas enters stage 0
  /// and the solvent the top stage (`addGasInStream` is `addFeedStream(stream, 0)`,
  /// `addSolventInStream` is `addFeedStream(stream, getNumberOfTrays() - 1)`). `StrippingColumn`
  /// renames the same two inlets. Both rows here are the classes' own tests, re-cased on PR:
  /// each pins **every tray** to one temperature, because an isothermal absorber is what those
  /// tests solve, and each is printed on the class's own product getters.
  /// **The hydrocarbon stripper's two compositions, in one place.** The product capture's
  /// `hydrocarbon_stripper` row and the mechanical capture's stripper row are two readings of one
  /// state, so they share the recipe rather than copying it.
  static final String[] STRIPPER_COMPONENTS =
      new String[] { "methane", "propane", "n-butane", "n-pentane", "n-heptane" };
  static final double[] STRIPPER_GAS_Z =
      new double[] { 0.9900, 0.0080, 0.0015, 0.0005, 0.0 };
  static final double[] STRIPPER_LIQUID_Z =
      new double[] { 0.02, 0.08, 0.12, 0.15, 0.63 };

  static void absorberRows() {
    String[] gas = new String[] { "methane", "ethane", "propane", "n-butane", "n-pentane",
        "n-heptane" };
    double[] gasZ = new double[] { 0.920, 0.040, 0.025, 0.010, 0.005, 0.0 };
    double[] pureHeptane = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
    String[] strip = STRIPPER_COMPONENTS;
    double[] stripZ = STRIPPER_GAS_Z;
    double[] richZ = STRIPPER_LIQUID_Z;

    // **Three rows per machine, and the difference between them is the gate.**
    // `AbsorptionColumnTest` and `StrippingColumnTest` both solve an isothermal column: they
    // pin every stage to one temperature, which makes `solveSequential`'s own gate - the mean
    // tray-temperature change - zero, so the solve stops after its first sweep. Those rows are
    // evidence of that, and the two unpinned rows are the column.
    absorberRow("lean_oil_absorber_pinned_1e2", gas, gasZ, 2000.0, 30.0, gas, pureHeptane, 600.0,
        20.0, 5, 15.0, 298.15, 1.0e-2, "DIRECT_SUBSTITUTION");
    absorberRow("lean_oil_absorber_unpinned_1e2", gas, gasZ, 2000.0, 30.0, gas, pureHeptane,
        600.0, 20.0, 5, 15.0, 0.0, 1.0e-2, "DIRECT_SUBSTITUTION");
    absorberRow("lean_oil_absorber", gas, gasZ, 2000.0, 30.0, gas, pureHeptane, 600.0, 20.0, 5,
        15.0, 0.0, 1.0e-4, "DIRECT_SUBSTITUTION");
    // The stripper's own case states `MESH_RESIDUAL`, which this port refuses, so the row
    // keeps the class's own strategy beside the answer it produces.
    absorberRow("hydrocarbon_stripper_pinned_1e2", strip, stripZ, 150.0, 70.0, strip, richZ,
        900.0, 70.0, 5, 12.0, 343.15, 1.0e-2, "MESH_RESIDUAL");
    absorberRow("hydrocarbon_stripper_unpinned_1e2", strip, stripZ, 150.0, 70.0, strip, richZ,
        900.0, 70.0, 5, 12.0, 0.0, 1.0e-2, "DIRECT_SUBSTITUTION");
    absorberRow("hydrocarbon_stripper", strip, stripZ, 150.0, 70.0, strip, richZ, 900.0, 70.0, 5,
        12.0, 0.0, 1.0e-4, "DIRECT_SUBSTITUTION");

    // **The same two machines with `setReactive(true)`**, which is the sibling ids' reactive
    // section and has no test of its own in the class: no `AbsorptionColumnTest`,
    // `StrippingColumnTest` or `PackedColumnTest` calls `setReactive`. `AbsorptionColumn`
    // overrides no `run`, so this is the base's reactive tray on the base's middle trays - and
    // an absorber has no ends, so `setReactive(true)` flags *every* stage, which the row's own
    // `reactive_tray_flags` prints. On this fluid the carbon-and-hydrogen formula matrix has
    // rank two against six species, so the section is not the NR=0 no-op: the trays really do
    // route to the reactive flash.
    absorberRow("lean_oil_absorber_reactive", gas, gasZ, 2000.0, 30.0, gas, pureHeptane, 600.0,
        20.0, 5, 15.0, 0.0, 1.0e-4, "DIRECT_SUBSTITUTION", true);
    absorberRow("hydrocarbon_stripper_reactive", strip, stripZ, 150.0, 70.0, strip, richZ, 900.0,
        70.0, 5, 12.0, 0.0, 1.0e-4, "DIRECT_SUBSTITUTION", true);
  }

  static void absorberRow(String label, String[] gasNames, double[] gasZ, double gasKgPerHour,
      double gasC, String[] solventNames, double[] solventZ, double solventKgPerHour,
      double solventC, int trays, double bara, double stageC, double tolerance, String solver) {
    absorberRow(label, gasNames, gasZ, gasKgPerHour, gasC, solventNames, solventZ,
        solventKgPerHour, solventC, trays, bara, stageC, tolerance, solver, false);
  }

  /// **The same machine with a Murphree efficiency**, which `AbsorptionColumn` corrects through
  /// its own override rather than the base's: both phases are blended and the flash's vapour
  /// moles are re-allocated across the components. `componentName` is null for the column-wide
  /// row and names one component for the per-component row, which is the level
  /// `setComponentMurphreeEfficiency(String, double)` writes.
  static void absorberEfficiencyRow(String label, String[] gasNames, double[] gasZ,
      double gasKgPerHour, double gasC, String[] solventNames, double[] solventZ,
      double solventKgPerHour, double solventC, int trays, double bara, double stageC,
      double tolerance, double murphree, String componentName, double componentEfficiency) {
    absorberRow(label, gasNames, gasZ, gasKgPerHour, gasC, solventNames, solventZ,
        solventKgPerHour, solventC, trays, bara, stageC, tolerance, "DIRECT_SUBSTITUTION", false,
        murphree, componentName, componentEfficiency);
  }

  /// **The Murphree rows, on the lean-oil absorber the capture already holds.** Three of them:
  /// the column-wide value, one component overridden, and that override reapplied - which is
  /// what separates the two resolution levels from one another.
  static void absorberEfficiencyRows() {
    String[] gas = new String[] { "methane", "ethane", "propane", "n-butane", "n-pentane",
        "n-heptane" };
    double[] gasZ = new double[] { 0.920, 0.040, 0.025, 0.010, 0.005, 0.0 };
    double[] pureHeptane = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
    absorberEfficiencyRow("absorber_murphree_column_wide", gas, gasZ, 2000.0, 30.0, gas,
        pureHeptane, 600.0, 20.0, 5, 15.0, 0.0, 1.0e-2, 0.6, null, 0.0);
    absorberEfficiencyRow("absorber_murphree_one_component", gas, gasZ, 2000.0, 30.0, gas,
        pureHeptane, 600.0, 20.0, 5, 15.0, 0.0, 1.0e-2, 1.0, "methane", 0.6);
  }

  /// **The same machine with `setReactive(true)`.** `AbsorptionColumn` inherits the base's
  /// `setReactive` and overrides no `run`, so its middle trays are the base's middle trays - and
  /// an absorber has no ends at all, so `replaceMiddleTrays` walks *every* stage. The flag's
  /// own javadoc is wrong about the mechanism (the column constructs no `ReactiveTray`), so the
  /// row prints the per-tray flags the class actually set beside its answer.
  static void absorberRow(String label, String[] gasNames, double[] gasZ, double gasKgPerHour,
      double gasC, String[] solventNames, double[] solventZ, double solventKgPerHour,
      double solventC, int trays, double bara, double stageC, double tolerance, String solver,
      boolean reactive) {
    absorberRow(label, gasNames, gasZ, gasKgPerHour, gasC, solventNames, solventZ,
        solventKgPerHour, solventC, trays, bara, stageC, tolerance, solver, reactive, null, null,
        0.0);
  }

  /// One absorber or stripper state, held so two readings of it - the product capture and the
  /// mechanical design - are the same run rather than two copies of one recipe.
  static final class AbsorberCase {
    final Stream gas;
    final Stream solvent;
    final neqsim.process.equipment.absorber.AbsorptionColumn column;

    AbsorberCase(Stream gas, Stream solvent,
        neqsim.process.equipment.absorber.AbsorptionColumn column) {
      this.gas = gas;
      this.solvent = solvent;
      this.column = column;
    }
  }

  /// The construction and the run, separated from the printing.
  static AbsorberCase runAbsorber(String label, String[] gasNames, double[] gasZ,
      double gasKgPerHour, double gasC, String[] solventNames, double[] solventZ,
      double solventKgPerHour, double solventC, int trays, double bara, double stageC,
      double tolerance, String solver, boolean reactive, Double murphree, String componentName,
      double componentEfficiency) {
    Stream gas = feedMass("gas", gasNames, gasZ, gasC + 273.15, bara, gasKgPerHour);
    Stream solvent = feedMass("solvent", solventNames, solventZ, solventC + 273.15, bara,
        solventKgPerHour);

    // **The subclass for the stripper's row**, which is what the class's own test builds: it
    // renames the two inlets and adds no arithmetic, and a row is held to the class it names.
    boolean stripper = label.contains("stripper");
    neqsim.process.equipment.absorber.AbsorptionColumn column = stripper
        ? new neqsim.process.equipment.absorber.StrippingColumn(label, trays)
        : new neqsim.process.equipment.absorber.AbsorptionColumn(label, trays);
    // **The flag is set before the inlets are attached**, which is the order the class's own
    // reactive column tests use.
    column.setReactive(reactive);
    if (stripper) {
      ((neqsim.process.equipment.absorber.StrippingColumn) column).addStrippingGasStream(gas);
      ((neqsim.process.equipment.absorber.StrippingColumn) column).addRichLiquidStream(solvent);
    } else {
      column.addGasInStream(gas);
      column.addSolventInStream(solvent);
    }
    column.setTopPressure(bara);
    column.setBottomPressure(bara);
    // **A stage temperature of zero means no pin at all**, which is the experiment that says
    // whether the class's isothermal rows are a converged counter-current state or the first
    // sweep of one.
    if (stageC > 0.0) {
      // **The setter is kelvin**, which is what the classes' own tests state their stage
      // temperature in; `setTemperature(double)`, the other convention in this library, is
      // Celsius, and the two differ by 273.15.
      for (int tray = 0; tray < column.getNumberOfTrays(); tray++) {
        column.getTray(tray).setOutletTemperature(stageC);
      }
    }
    column.setTemperatureTolerance(tolerance);
    column.setMassBalanceTolerance(5.0e-2);
    column.setEnthalpyBalanceTolerance(5.0e-2);
    column.setMaxNumberOfIterations(80, true);
    column.setSolverType(
        neqsim.process.equipment.distillation.DistillationColumn.SolverType.valueOf(solver));
    // **Before the solve**, which is what `setDoInitializion`'s effect requires - the efficiencies
    // are read by `applyMurphreeCorrection` on every sweep.
    if (murphree != null) {
      column.setMurphreeEfficiency(murphree);
    }
    if (componentName != null) {
      column.setComponentMurphreeEfficiency(componentName, componentEfficiency);
    }
    // **Through a `ProcessSystem`, which is how the class's own tests drive it**: the two feeds
    // are units in the flowsheet and the column runs after them, where a direct `column.run()`
    // leaves the tray inlets in a state the tests never see.
    neqsim.process.processmodel.ProcessSystem process = new neqsim.process.processmodel.ProcessSystem();
    process.add(gas);
    process.add(solvent);
    process.add(column);
    process.run();
    return new AbsorberCase(gas, solvent, column);
  }

  /// The full row, with the Murphree pair the efficiency rows carry and every other row leaves
  /// null.
  static void absorberRow(String label, String[] gasNames, double[] gasZ, double gasKgPerHour,
      double gasC, String[] solventNames, double[] solventZ, double solventKgPerHour,
      double solventC, int trays, double bara, double stageC, double tolerance, String solver,
      boolean reactive, Double murphree, String componentName, double componentEfficiency) {
    AbsorberCase unit = runAbsorber(label, gasNames, gasZ, gasKgPerHour, gasC, solventNames,
        solventZ, solventKgPerHour, solventC, trays, bara, stageC, tolerance, solver, reactive,
        murphree, componentName, componentEfficiency);
    Stream gas = unit.gas;
    Stream solvent = unit.solvent;
    neqsim.process.equipment.absorber.AbsorptionColumn column = unit.column;

    System.out.println(label);
    System.out.println("trays=" + trays);
    System.out.println("pressure_bara=" + bara);
    System.out.println("stage_temperature_K=" + (stageC > 0.0 ? stageC : 0.0));
    System.out.println("temperature_tolerance=" + tolerance);
    System.out.println("solver_requested=" + solver);
    System.out.println("solved=" + column.solved());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("solver=" + column.getLastSolverTypeUsed());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    System.out.println("mass_residual=" + column.getLastMassResidual());
    System.out.println("energy_residual=" + column.getLastEnergyResidual());
    // **Only on a reactive row**, so every row already in the capture keeps its own bytes.
    if (reactive) {
      System.out.println("reactive=true");
      System.out.println("reports_reactive=" + column.isReactive());
      StringBuilder flags = new StringBuilder("reactive_tray_flags=");
      for (int i = 0; i < column.getNumberOfTrays(); i++) {
        flags.append(column.getTray(i).isUseReactiveFlash() ? "1" : "0");
      }
      System.out.println(flags);
    }
    // **The resolved efficiency, per stage and per component.** `getMurphreeEfficiency(stage)`
    // is the base's two-step resolution and `getComponentMurphreeEfficiency(stage, name)` the
    // absorber's third, so these two lines are what the override actually reads - the request
    // the row stated is not, once a component has its own value.
    if (murphree != null || componentName != null) {
      StringBuilder resolved = new StringBuilder("tray_murphree_efficiency=");
      for (int i = 0; i < column.getNumberOfTrays(); i++) {
        if (i > 0) {
          resolved.append(' ');
        }
        resolved.append(column.getMurphreeEfficiency(i));
      }
      System.out.println(resolved);
      if (componentName != null) {
        StringBuilder byComponent = new StringBuilder("component_murphree_efficiency=");
        for (int i = 0; i < column.getNumberOfTrays(); i++) {
          if (i > 0) {
            byComponent.append(' ');
          }
          byComponent.append(column.getComponentMurphreeEfficiency(i, componentName));
        }
        System.out.println(byComponent);
      }
    }
    System.out.println("gas_in_mol_per_sec=" + gas.getFlowRate("mol/sec"));
    System.out.println("solvent_in_mol_per_sec=" + solvent.getFlowRate("mol/sec"));
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println("tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + column.getTray(i).getPressure());
      System.out.println("tray" + i + "_gas_n=" + column.getTray(i).getGasOutStream().getFlowRate("mol/sec"));
      System.out.println("tray" + i + "_liquid_n=" + column.getTray(i).getLiquidOutStream().getFlowRate("mol/sec"));
    }
    print("gas_out", column.getGasOutStream());
    print("liquid_out", column.getLiquidOutStream());
    System.out.println();
  }

  /// **The packed column, whose packing is a report.** `PackedColumn extends DistillationColumn`
  /// and its `run` is `super.run(id)` followed by `calcPackingHydraulics()`, so the separation is
  /// the base's, unchanged, and everything the packing adds is read *after* the column has
  /// converged - into `PackingHydraulicsCalculator` through `ColumnInternalsDesigner`.
  ///
  /// **Two of its claims are measurements rather than readings of the source.**
  ///
  /// **The stage count is a constructor's arithmetic**: `estimateStages(packedHeight, 0.5)` is
  /// `ceil(packedHeight / 0.5)` floored at two, so the stage-count rows build the class at a
  /// stated height and print the trays it made, with no solve at all.
  ///
  /// **And the contactor constructor never reads the height**: `PackedColumn(name, gasInStream)`
  /// fixes ten middle trays with no ends - `AbsorptionColumn(name, 10)`'s own shape -
  /// so `setPackedHeight` afterwards reaches the hydraulics report alone.
  static void packedColumnRows() {
    // **The rule, from the constructor alone.** `middle` is `getNumberOfTrays()` less the two
    // ends, because the base increments its own count once for each end it adds.
    double[] heights = new double[] { -1.0, 0.0, 0.2, 0.4, 0.75, 1.0, 1.2, 2.5, 5.0, 5.3, 10.0 };
    for (double height : heights) {
      neqsim.process.equipment.distillation.PackedColumn column =
          new neqsim.process.equipment.distillation.PackedColumn("stages probe", height,
              "Pall-Ring-50", true, true);
      System.out.println("stage_count_" + height + "m");
      System.out.println("packed_height_m=" + height);
      System.out.println("stage_count_middle=" + (column.getNumberOfTrays() - 2));
      System.out.println("stage_count_trays_with_ends=" + column.getNumberOfTrays());
      System.out.println();
    }

    // **The class's own test state, run exactly as its test runs it**: no condenser or reboiler
    // temperature, so both ends flash at their own enthalpy. Kept as evidence, because what it
    // converges to is not a physical profile - its own test asserts only that the two product
    // streams exist.
    String[] names = new String[] { "methane", "ethane", "propane", "n-butane" };
    double[] z = new double[] { 0.40, 0.30, 0.20, 0.10 };
    packedColumnRow("packed_distillation_class_case_unpinned", names, z, 333.15, 10.0, 10000.0, 5.0,
        "Pall-Ring-25", 5, null, null, 10.0, 10.0, 2.0e-2, 50);

    // **The port's oracle, and it is the base column's own state.** `binary_methane_butane_
    // 4_stages` above is the row this port is held to, so a packed column at the height whose
    // constructor makes four middle trays - 2.0 m - with the same feed, the same feed tray and
    // the same two end pins must land on that row's answer. That is the claim an id for this
    // class makes, and this is the pair of rows that measures it.
    String[] binary = new String[] { "methane", "n-butane" };
    double[] binaryZ = new double[] { 0.5, 0.5 };
    packedColumnRow("packed_distillation_binary_2m", binary, binaryZ, 300.0, 20.0, 1000.0, 2.0,
        "Pall-Ring-50", 2, -20.0, 100.0, 19.0, 20.0, 1.0e-6, 200);
    // **The same height rule, one tray more**: 2.3 m is `ceil(4.6)` = five middle trays, so this
    // row is a different column's answer on the same feed - which is all the packing does.
    packedColumnRow("packed_distillation_binary_2m3", binary, binaryZ, 300.0, 20.0, 1000.0, 2.3,
        "Pall-Ring-50", 2, -20.0, 100.0, 19.0, 20.0, 1.0e-6, 200);

    // The contactor form, `testBasicAbsorber`'s own state re-cased on PR: a twelve-carbon
    // solvent against a light gas, ten middle trays, no ends, the solvent at the top stage.
    packedContactorRow("packed_contactor_6m");

    // **The oracle's own state with a reactive section**, which is the only sibling reactive
    // row that sits on a column this port already reproduces exactly: `packed_distillation_
    // binary_2m` is the base column's own `binary_methane_butane_4_stages`, so the difference
    // here is the section and nothing else. `PackedColumn` inherits `setReactive`, and its
    // middle trays are the four the height derives - the two ends stay equilibrium stages.
    packedColumnRow("packed_distillation_binary_2m_reactive", binary, binaryZ, 300.0, 20.0,
        1000.0, 2.0, "Pall-Ring-50", 2, -20.0, 100.0, 19.0, 20.0, 1.0e-6, 200, true);
  }

  /// The contactor constructor, which fixes ten middle trays and reads its packed height only in
  /// the hydraulics report.
  static void packedContactorRow(String label) {
    String[] gasNames = new String[] { "methane", "ethane", "propane" };
    double[] gasZ = new double[] { 0.90, 0.07, 0.03 };
    String[] solventNames = new String[] { "nC10", "propane" };
    double[] solventZ = new double[] { 0.99, 0.01 };
    Stream gas = feedMass("gas", gasNames, gasZ, 303.15, 50.0, 50000.0);
    Stream solvent = feedMass("solvent", solventNames, solventZ, 298.15, 50.0, 10000.0);

    neqsim.process.equipment.distillation.PackedColumn column =
        new neqsim.process.equipment.distillation.PackedColumn(label, gas);
    column.setPackedHeight(6.0);
    column.setPackingType("Pall-Ring-50");
    column.setStructuredPacking(false);
    column.setDesignFloodFraction(0.70);
    column.addSolventStream(solvent);
    column.setTopPressure(50.0);
    column.setBottomPressure(50.0);
    // **Through a `ProcessSystem`**, which is how the class's own test drives it.
    neqsim.process.processmodel.ProcessSystem process =
        new neqsim.process.processmodel.ProcessSystem();
    process.add(gas);
    process.add(solvent);
    process.add(column);
    process.run();

    System.out.println(label);
    System.out.println("stage_count_middle=" + (column.getNumberOfTrays() - 0));
    System.out.println("packed_height_m=" + column.getPackedHeight());
    System.out.println("tray_count=" + column.getNumberOfTrays());
    System.out.println("solved=" + column.solved());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("solver=" + column.getLastSolverTypeUsed());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    System.out.println("mass_residual=" + column.getLastMassResidual());
    System.out.println("energy_residual=" + column.getLastEnergyResidual());
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println(
          "tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + column.getTray(i).getPressure());
      System.out.println("tray" + i + "_gas_n="
          + column.getTray(i).getGasOutStream().getFlowRate("mol/sec"));
      System.out.println("tray" + i + "_liquid_n="
          + column.getTray(i).getLiquidOutStream().getFlowRate("mol/sec"));
      // **What the designer reads off this tray**, printed per tray so the port's own numbers
      // become a measurement rather than an assumption - and so the middle tray the report
      // actually reads can be identified rather than assumed to be `size() / 2`.
      printTrayProperties(column, i);
    }
    print("gas_out", column.getGasOutStream());
    print("liquid_out", column.getLiquidOutStream());
    printHydraulics(column);
    System.out.println();
  }

  static void packedColumnRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, double packedHeight, String packingType,
      int feedTray, Double condenserC, Double reboilerC, double topBara, double bottomBara,
      double tolerance, int maxIterations) {
    packedColumnRow(label, names, z, feedTemperatureK, feedPressureBara, kgPerHour, packedHeight,
        packingType, feedTray, condenserC, reboilerC, topBara, bottomBara, tolerance,
        maxIterations, false);
  }

  /// **The same column with `setReactive(true)`.** `PackedColumn extends DistillationColumn`
  /// and its `run` is `super.run(id)` before the packing's report, so the reactive section is
  /// the base's own on the stage count the packing's height derives.
  static void packedColumnRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, double packedHeight, String packingType,
      int feedTray, Double condenserC, Double reboilerC, double topBara, double bottomBara,
      double tolerance, int maxIterations, boolean reactive) {
    SystemInterface fluid = new SystemPrEos(feedTemperatureK, feedPressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("feed", fluid);
    inlet.setFlowRate(kgPerHour, "kg/hr");
    inlet.run();

    // **The class's own constructor**, so the stage count is the class's arithmetic and not a
    // number this probe chose. A null end temperature leaves that end unpinned, which is what
    // the class's own test does.
    neqsim.process.equipment.distillation.PackedColumn column =
        new neqsim.process.equipment.distillation.PackedColumn(label, packedHeight, packingType,
            true, true);
    column.setReactive(reactive);
    column.addFeedStream(inlet, feedTray);
    if (condenserC != null) {
      column.setCondenserTemperature(condenserC, "C");
    }
    if (reboilerC != null) {
      column.setReboilerTemperature(reboilerC, "C");
    }
    column.setTopPressure(topBara);
    column.setBottomPressure(bottomBara);
    column.setTemperatureTolerance(tolerance);
    column.setMaxNumberOfIterations(maxIterations, true);
    // **Through a `ProcessSystem`**, which is how the class's own test drives it: the feed is a
    // unit in the flowsheet and the column runs after it, where a direct `column.run()` leaves
    // the tray inlets in a state the tests never see.
    neqsim.process.processmodel.ProcessSystem process =
        new neqsim.process.processmodel.ProcessSystem();
    process.add(inlet);
    process.add(column);
    process.run();

    System.out.println(label);
    if (reactive) {
      System.out.println("reactive=true");
      System.out.println("reports_reactive=" + column.isReactive());
      StringBuilder flags = new StringBuilder("reactive_tray_flags=");
      for (int i = 0; i < column.getNumberOfTrays(); i++) {
        flags.append(column.getTray(i).isUseReactiveFlash() ? "1" : "0");
      }
      System.out.println(flags);
    }
    System.out.println("packed_height_m=" + packedHeight);
    System.out.println("packing_type=" + packingType);
    System.out.println("stage_count_middle=" + (column.getNumberOfTrays() - 2));
    System.out.println("tray_count=" + column.getNumberOfTrays());
    System.out.println("feed_tray=" + feedTray);
    System.out.println("condenser_pin_C=" + (condenserC == null ? "none" : condenserC));
    System.out.println("reboiler_pin_C=" + (reboilerC == null ? "none" : reboilerC));
    System.out.println("temperature_tolerance=" + tolerance);
    System.out.println("feed_mol_per_sec=" + inlet.getFlowRate("mol/sec"));
    System.out.println("feed_kg_per_hour=" + inlet.getFlowRate("kg/hr"));
    System.out.println("solved=" + column.solved());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("solver=" + column.getLastSolverTypeUsed());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    System.out.println("mass_residual=" + column.getLastMassResidual());
    System.out.println("energy_residual=" + column.getLastEnergyResidual());
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println(
          "tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + column.getTray(i).getPressure());
      System.out.println("tray" + i + "_gas_n="
          + column.getTray(i).getGasOutStream().getFlowRate("mol/sec"));
      System.out.println("tray" + i + "_liquid_n="
          + column.getTray(i).getLiquidOutStream().getFlowRate("mol/sec"));
      // **What the designer reads off this tray**, printed per tray so the port's own numbers
      // become a measurement rather than an assumption - and so the middle tray the report
      // actually reads can be identified rather than assumed to be `size() / 2`.
      printTrayProperties(column, i);
    }
    print("distillate", column.getGasOutStream());
    print("bottoms", column.getLiquidOutStream());
    System.out.println("condenser_duty_W=" + column.getCondenser().getDuty());
    System.out.println("reboiler_duty_W=" + column.getReboiler().getDuty());
    printHydraulics(column);
    System.out.println();
  }

  /// **The packing's own report, which this port does not carry.** Printed so the capture holds
  /// what the class reads the packing for, and so a row with `hydraulics_ok=false` or a zero HETP
  /// is visible as the design failure it is rather than as a missing number.
  static void printHydraulics(neqsim.process.equipment.distillation.PackedColumn column) {
    System.out.println("hetp_m=" + column.getHETP());
    System.out.println("theoretical_stages=" + column.getTheoreticalStages());
    System.out.println("percent_flood=" + column.getPercentFlood());
    System.out.println("flooding_velocity_m_per_s=" + column.getFloodingVelocity());
    System.out.println("packing_pressure_drop_Pa=" + column.getPackingPressureDrop());
    System.out.println("hydraulics_ok=" + column.isHydraulicsOk());
    System.out.println("internal_diameter_m=" + column.getInternalDiameter());
    neqsim.process.equipment.distillation.internals.PackingHydraulicsCalculator hydraulics =
        column.getHydraulics();
    if (hydraulics != null) {
      System.out.println("fs_factor=" + hydraulics.getFsFactor());
      System.out.println("kga=" + hydraulics.getKGa());
      System.out.println("kla=" + hydraulics.getKLa());
      System.out.println("wetted_area_m2_per_m3=" + hydraulics.getWettedArea());
      System.out.println("htu_g_m=" + hydraulics.getHtuG());
      System.out.println("htu_l_m=" + hydraulics.getHtuL());
      System.out.println("htu_og_m=" + hydraulics.getHtuOG());
      // **The calculator's own verdicts, printed beside the column's.** `isHydraulicsOk()` above
      // and `isDesignOk()` here are *different predicates* that both read as "the hydraulics are
      // fine" - and every row of this capture has both false, so a port that published one under
      // the other's name would reproduce every row while answering a different question.
      System.out.println("design_ok=" + hydraulics.isDesignOk());
      System.out.println("wetting_ok=" + hydraulics.isWettingOk());
      // **Per metre or in total, which the name does not say.** The column's
      // `getPackingPressureDrop()` is printed above and this is the calculator's own; over a
      // 2.3 m bed the two differ by exactly the height, so the pair settles which is which.
      System.out.println("pressure_drop_per_meter_Pa=" + hydraulics.getPressureDropPerMeter());
      // The geometry the calculator resolved and the diameter it sized, which is what the
      // sizing calc's state has to match.
      System.out.println("resolved_packing=" + hydraulics.getPackingName());
      System.out.println("resolved_packing_factor=" + hydraulics.getPackingFactor());
      System.out.println("calculated_diameter_m=" + hydraulics.getColumnDiameter());
    }
  }

  /// The properties `ColumnInternalsDesigner.getTrayProperties` reads off one tray, with the same
  /// phase selection it makes, plus the interphase call itself.
  ///
  /// **`interphase_surface_tension_N_per_m` is the decisive one.** The designer asks
  /// `fluid.getInterphaseProperties().getSurfaceTension(0, 1)` and answers `0.02` when that is not
  /// finite and positive, so printing what this call returns shows the fallback being taken rather
  /// than leaving it to be inferred from the wetted area. The four properties around it are what
  /// the report's inputs have never been printed as, and they are what turns the port's own
  /// density and viscosity from an assumption into a measured number.
  static void printTrayProperties(neqsim.process.equipment.distillation.PackedColumn column,
      int index) {
    try {
      neqsim.thermo.system.SystemInterface fluid = column.getTray(index).getFluid();
      if (fluid == null) {
        return;
      }
      int phases = fluid.getNumberOfPhases();
      System.out.println("tray" + index + "_phase_count=" + phases);
      // **The order the phases sit in, which `getSurfaceTension(0, 1)` depends on.** The
      // designer resolves its liquid phase by type and then asks by index, so whether those two
      // agree is the whole question.
      for (int p = 0; p < phases; p++) {
        System.out.println(
            "tray" + index + "_phase" + p + "_type=" + fluid.getPhase(p).getType());
      }
      if (phases < 2) {
        return;
      }
      fluid.initProperties();
      double gasDensity = fluid.hasPhaseType("gas")
          ? fluid.getPhase("gas").getDensity("kg/m3") : fluid.getPhase(0).getDensity("kg/m3");
      int liquidPhase = fluid.hasPhaseType("oil") ? fluid.getPhaseNumberOfPhase("oil")
          : (fluid.hasPhaseType("aqueous") ? fluid.getPhaseNumberOfPhase("aqueous") : 1);
      double liquidDensity = fluid.getPhase(liquidPhase).getDensity("kg/m3");
      double liquidViscosity = fluid.getPhase(liquidPhase).getViscosity("kg/msec");
      double sigma = fluid.getInterphaseProperties().getSurfaceTension(0, 1);
      // **The literal pair the designer asks with, against the resolved pair it has already
      // worked out.** `getTrayProperties` resolves its vapour and liquid phases by type and then
      // calls `getSurfaceTension(0, 1)` by index; if the two disagree the resolution is wasted.
      int gasPhase = fluid.hasPhaseType("gas") ? fluid.getPhaseNumberOfPhase("gas") : 0;
      double sigmaResolved = fluid.getInterphaseProperties().getSurfaceTension(gasPhase, liquidPhase);
      System.out.println("tray" + index + "_gas_phase_index=" + gasPhase);
      System.out.println("tray" + index + "_liquid_phase_index=" + liquidPhase);
      System.out.println(
          "tray" + index + "_sigma_by_resolved_indices_N_per_m=" + sigmaResolved);
      System.out.println("tray" + index + "_gas_density_kg_per_m3=" + gasDensity);
      System.out.println("tray" + index + "_liquid_density_kg_per_m3=" + liquidDensity);
      System.out.println("tray" + index + "_liquid_viscosity_Pa_s=" + liquidViscosity);
      System.out.println("tray" + index + "_interphase_surface_tension_N_per_m=" + sigma);
    } catch (Exception ex) {
      System.out.println("tray" + index + "_properties_error=" + ex.getClass().getSimpleName());
    }
  }

  /// **The reactive tray, which the class carries and its own tests never run.** `setReactive`
  /// sets `useReactiveFlash` on the *middle* trays - the ends never - and the tray's flash then
  /// routes to `ThermodynamicOperations.reactiveTPflash`/`reactivePHflash` instead of the
  /// equilibrium one.
  ///
  /// **Both of NeqSim's reactive-column tests assert nothing numeric**: one checks the flags and
  /// never calls `run()`, and the other asserts only that a reactive column on a fluid with no
  /// reaction matches a standard one. So each row here is paired with its standard twin, which
  /// is the measurement the class does make, and the class's own water-gas-shift state - which
  /// its test never runs - is run here to see what it does.
  ///
  /// **On PR**, the cubic the process layer's streams carry: the oracle has to be the fluid the
  /// port flashes.
  static void reactiveColumnRows() {
    String[] wgs = new String[] { "CO", "water", "CO2", "hydrogen" };
    double[] wgsZ = new double[] { 0.25, 0.25, 0.25, 0.25 };
    // The class's own state: three middle trays, both ends, the feed on stage 2, and its own
    // two end temperatures.
    reactiveColumnRow("wgs_standard_column_pr", wgs, wgsZ, 473.15, 10.0, 1000.0, 3, 2, 100.0,
        250.0, 10.0, 10.0, 1.0e-2, 80, false);
    reactiveColumnRow("wgs_reactive_column_pr", wgs, wgsZ, 473.15, 10.0, 1000.0, 3, 2, 100.0,
        250.0, 10.0, 10.0, 1.0e-2, 80, true);
    // The class's own `NR = 0` pair: methane/ethane at -50 C and 15 bara, where its own test
    // asserts the two columns agree. **That state sits at the mixture's pseudocritical
    // temperature** - `0.7 * 190 + 0.3 * 305` is 224 K and the feed is 223.15 - which is the
    // delicacy its own test never exercises numerically.
    String[] light = new String[] { "methane", "ethane" };
    double[] lightZ = new double[] { 0.7, 0.3 };
    reactiveColumnRow("methane_ethane_standard_column_pr", light, lightZ, 223.15, 15.0, 1000.0, 3,
        2, -80.0, -20.0, 15.0, 15.0, 1.0e-2, 100, false);
    reactiveColumnRow("methane_ethane_reactive_column_pr", light, lightZ, 223.15, 15.0, 1000.0, 3,
        2, -80.0, -20.0, 15.0, 15.0, 1.0e-2, 100, true);
  }

  /// **The reactive section on the fluid the column capture's own case is held to.** Methane and
  /// n-butane have no independent reaction either, so the reactive route is the class's `NR = 0`
  /// delegation and the two rows must agree - which its own test asserts to `0.01` kg/hr and
  /// which is measured here to the last digit. They live in that capture because the layer diff
  /// pairs one capture with one model, and this model's case is the binary column there.
  static void reactiveColumnRowsInColumnCapture() {
    String[] binary = new String[] { "methane", "n-butane" };
    double[] binaryZ = new double[] { 0.5, 0.5 };
    reactiveColumnRow("binary_standard_column_pr", binary, binaryZ, 300.0, 20.0, 1000.0, 4, 2,
        -20.0, 100.0, 19.0, 20.0, 1.0e-6, 200, false);
    reactiveColumnRow("binary_reactive_column_pr", binary, binaryZ, 300.0, 20.0, 1000.0, 4, 2,
        -20.0, 100.0, 19.0, 20.0, 1.0e-6, 200, true);
  }

  static void reactiveColumnRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, double condenserC,
      double reboilerC, double topBara, double bottomBara, double tolerance, int maxIterations,
      boolean reactive) {
    SystemInterface fluid = new SystemPrEos(feedTemperatureK, feedPressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("column feed", fluid);
    inlet.setFlowRate(kgPerHour, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn("col1", stages, true, true);
    // **The flag is set before the feed is attached**, which is the order both of the class's
    // own reactive tests use.
    column.setReactive(reactive);
    column.addFeedStream(inlet, feedTray);
    column.setCondenserTemperature(condenserC, "C");
    column.setReboilerTemperature(reboilerC, "C");
    column.setTopPressure(topBara);
    column.setBottomPressure(bottomBara);
    column.setTemperatureTolerance(tolerance);
    column.setMaxNumberOfIterations(maxIterations, true);
    column.run();

    System.out.println(label);
    System.out.println("reactive=" + reactive);
    System.out.println("cubic=pr");
    System.out.println("stages=" + stages);
    System.out.println("feed_tray=" + feedTray);
    System.out.println("temperature_tolerance=" + tolerance);
    System.out.println("tray_count=" + column.getNumberOfTrays());
    System.out.println("reports_reactive=" + column.isReactive());
    StringBuilder flags = new StringBuilder("reactive_tray_flags=");
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      flags.append(column.getTray(i).isUseReactiveFlash() ? "1" : "0");
    }
    System.out.println(flags);
    System.out.println("solved=" + column.solved());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    System.out.println("mass_residual=" + column.getLastMassResidual());
    System.out.println("energy_residual=" + column.getLastEnergyResidual());
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

  static void columnRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, double condenserC,
      double reboilerC, double topBara, double bottomBara, double tolerance, int maxIterations,
      boolean condenser, boolean totalCondenser) {
    runColumn(label, names, z, feedTemperatureK, feedPressureBara, kgPerHour, stages, feedTray,
        condenserC, reboilerC, topBara, bottomBara, tolerance, maxIterations, condenser,
        totalCondenser, null, null);
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
        totalCondenser, solver, null);
  }

  /// **The column under a Murphree efficiency**, which is the same machine with the vapour
  /// leaving each interior tray blended towards its inlet's - `DistillationColumn
  /// .applyMurphreeCorrection`, applied after every tray's run on the sequential path.
  ///
  /// **The column is the port's own primary oracle and not the class's Murphree test column.**
  /// `runBinaryMurphreeColumn` drives its two ends by *reflux ratio* (`getCondenser()
  /// .setRefluxRatio(2.0)`), so its stated end temperatures are inert and its state is a
  /// `PVrefluxflash` one - the route `process.distillation_column` measures as not reproduced.
  /// A row on that column would put the correction and the unported solve in the same number.
  /// So the state here is `binary_methane_butane_4_stages`, which the port already reproduces
  /// exactly, and the *difference* from that row is the correction and nothing else.
  static void murphreeRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, double condenserC,
      double reboilerC, double topBara, double bottomBara, double tolerance, int maxIterations,
      boolean condenser, boolean totalCondenser, double efficiency) {
    runColumn(label, names, z, feedTemperatureK, feedPressureBara, kgPerHour, stages, feedTray,
        condenserC, reboilerC, topBara, bottomBara, tolerance, maxIterations, condenser,
        totalCondenser, "DIRECT_SUBSTITUTION", efficiency);
  }

  static void runColumn(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, double condenserC,
      double reboilerC, double topBara, double bottomBara, double tolerance, int maxIterations,
      boolean condenser, boolean totalCondenser, String solver, Double murphreeEfficiency) {
    runColumn(label, names, z, feedTemperatureK, feedPressureBara, kgPerHour, stages, feedTray,
        condenserC, reboilerC, topBara, bottomBara, tolerance, maxIterations, condenser,
        totalCondenser, solver, murphreeEfficiency, -1, 0.0, null, false);
  }

  /// The same, with the per-stage override: `array` when it is not `null` and the single-stage
  /// setter otherwise, applied after the column-wide value.
  static void runColumn(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, double condenserC,
      double reboilerC, double topBara, double bottomBara, double tolerance, int maxIterations,
      boolean condenser, boolean totalCondenser, String solver, Double murphreeEfficiency,
      int overrideStage, double overrideEfficiency, double[] overrideArray, boolean divergence) {
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
    if (murphreeEfficiency != null) {
      column.setMurphreeEfficiency(murphreeEfficiency);
    }
    if (overrideArray != null) {
      column.setMurphreeEfficiencies(overrideArray);
    } else if (overrideStage >= 0) {
      column.setMurphreeEfficiency(overrideStage, overrideEfficiency);
    }
    column.run();

    System.out.println(label);
    System.out.println("stages=" + stages);
    System.out.println("feed_tray=" + feedTray);
    System.out.println("temperature_tolerance=" + tolerance);
    System.out.println("total_condenser=" + totalCondenser);
    if (murphreeEfficiency != null) {
      System.out.println("murphree_efficiency=" + column.getMurphreeEfficiency());
    }
    // **The resolved per-stage vector is printed, not the request.** `getMurphreeEfficiency(stage)`
    // is `getEffectiveMurphreeEfficiency`, so this line is the class's own answer for every stage
    // - the `murphree_efficiency=` line above is only the column-wide field and says nothing
    // about an override.
    if (overrideArray != null || overrideStage >= 0) {
      StringBuilder resolved = new StringBuilder("tray_murphree_efficiency=");
      for (int stage = 0; stage < column.getNumberOfTrays(); stage++) {
        resolved.append(column.getMurphreeEfficiency(stage));
        if (stage + 1 < column.getNumberOfTrays()) {
          resolved.append(' ');
        }
      }
      System.out.println(resolved);
    }
    print("feed", inlet);
    System.out.println("feed_mol_per_sec=" + inlet.getFlowRate("mol/sec"));
    System.out.println("feed_kg_per_hour=" + inlet.getFlowRate("kg/hr"));
    System.out.println("tray_count=" + column.getNumberOfTrays());
    if (divergence) {
      divergenceBlock(column);
    }
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
    // **The *tray terminals*, which are what the public products were reconciled from.**
    // `DistillationColumn.run` ends in `updateProductsFromExternalComponentBalance`, which
    // rescales the two public products per component so that they carry the feed exactly - so
    // what `distillate` and `bottoms` above print is a *reconciled* record rather than the tray
    // solve's own. The terminals are the solve's, and a reader comparing a port against this
    // capture needs both to see which of the two a difference lives in.
    if (murphreeEfficiency != null) {
      // **What the correction's two phase guards actually see**, per stage, beside the
      // compositions it reads: `getNumberOfPhases() < 2` skips a stage, and a stage whose flash
      // found one phase is a stage this port's own `(Some, Some)` guard also skips - unless the
      // two disagree, which is what a divergence would have to come from.
      for (int i = 0; i < column.getNumberOfTrays(); i++) {
        neqsim.thermo.system.SystemInterface system = column.getTray(i).getThermoSystem();
        System.out.println("tray" + i + "_phases=" + system.getNumberOfPhases() + " phase0type="
            + (system.getNumberOfPhases() > 0 ? system.getPhase(0).getType() : "none") + " beta="
            + system.getBeta());
      }
      print("terminal_top", column.getTray(column.getNumberOfTrays() - 1).getGasOutStream());
      print("terminal_bottom", column.getTray(0).getLiquidOutStream());
      // **What the correction actually produced on each corrected stage, and what it read.**
      // `yEq` is the stage's own equilibrium vapour (`getThermoSystem().getPhase(0)`) and `yIn`
      // the stage below's - so `y_out = y_in + E(y_eq - y_in)` can be checked by hand against a
      // port's own arithmetic rather than inferred from a profile two solvers disagree on.
      for (int i = 0; i < column.getNumberOfTrays(); i++) {
        SystemInterface system = column.getTray(i).getThermoSystem();
        StringBuilder eqZ = new StringBuilder();
        for (int j = 0; j < system.getPhase(0).getNumberOfComponents(); j++) {
          eqZ.append(system.getPhase(0).getComponent(j).getName()).append(':')
              .append(system.getPhase(0).getComponent(j).getx()).append(' ');
        }
        System.out.println("tray" + i + "_eq_gas_z=" + eqZ.toString().trim());
        if (i > 0) {
          SystemInterface below = column.getTray(i - 1).getThermoSystem();
          StringBuilder belowZ = new StringBuilder();
          for (int j = 0; j < below.getPhase(0).getNumberOfComponents(); j++) {
            belowZ.append(below.getPhase(0).getComponent(j).getName()).append(':')
                .append(below.getPhase(0).getComponent(j).getx()).append(' ');
          }
          System.out.println("tray" + i + "_in_gas_z=" + belowZ.toString().trim());
        }
        SystemInterface outSystem = column.getTray(i).getGasOutStream().getThermoSystem();
        // **The stream's *overall* composition, which is what the stage above receives.** The
        // cached corrected system is `phaseToSystem(0)` with the composition replaced and then
        // re-initialised, so its own phase 0 is the *equilibrium* vapour at its T and P again -
        // reading that instead of the overall would report the correction as a no-op.
        StringBuilder corrected = new StringBuilder();
        double[] overall = outSystem.getMolarComposition();
        for (int j = 0; j < overall.length; j++) {
          corrected.append(outSystem.getComponent(j).getName()).append(':').append(overall[j])
              .append(' ');
        }
        System.out.println("tray" + i + "_out_gas_z=" + corrected.toString().trim());
        System.out.println("tray" + i + "_out_gas_phases=" + outSystem.getNumberOfPhases());
        System.out.println("tray" + i + "_out_gas_n=" + outSystem.getTotalNumberOfMoles());
        System.out.println("tray" + i + "_out_gas_T=" + outSystem.getTemperature());
        System.out.println("tray" + i + "_out_gas_P=" + outSystem.getPressure());
        System.out.println("tray" + i + "_out_gas_h=" + outSystem.getPhase(0).getEnthalpy("J/mol"));
      }
    }
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

    // The side-draw tray rows, which belong with the tray's own - see the helper's note.
    sideDrawRowsInTrayCapture();
  }

  /// The side-draw tray rows, in the tray's own capture for the same reason the column rows are in
  /// the column's: the tray has no palette entry, so its rows are the oracle its own Rust tests
  /// read rather than a model case's.
  static void sideDrawRowsInTrayCapture() {
    sideDrawTrayRows();
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

  /// **The side draw, which is a split of the tray's own outlet phase.** From
  /// `SimpleTraySideDrawTest`: `getGasOutStream` is the vapour scaled by `1 - gasSideDrawFraction`,
  /// `getLiquidOutStream` the liquid by `1 - liquidSideDrawFraction - liquidPumparoundDrawFraction`,
  /// and each draw the *same* phase scaled by its own - so a draw carries the tray's composition,
  /// temperature and pressure and a different flow, and a draw taken from the wrong phase is
  /// impossible by construction.
  ///
  /// **Every row is a pure component, and that is what makes the split an oracle rather than an
  /// identity.** The class's own tests run a pure-methane feed (`SystemSrkEos`, 300 K, 10 bara, 100
  /// kg/hr, classic rule) through two standalone trays - a reference and one drawing a quarter of
  /// its vapour - and assert the reference's gas equals the drawing tray's gas plus its draw, and
  /// that the draw is `0.25` of the reference. On a pure fluid the flash cannot move either
  /// composition or the mass flow, so the split is exact arithmetic on `100` kg/hr and every
  /// library agrees on it whatever its equation of state: `25.0` kg/hr is a number, not a band.
  ///
  /// **The fluid is PR, where the class's test uses SRK**, because the process tier is PR-only and
  /// an oracle in a fluid the port cannot carry would be no oracle. On the pure fluids here the two
  /// cubics give the same single phase, so the row is the class's own state in the fluid the port
  /// holds.
  ///
  /// **The pumparound row has no test in the class.** It is here because it is the one state where
  /// *two* liquid fractions compose: `0.20` and `0.30` leave half the liquid in the tray, which is
  /// the composition `validateLiquidSplitFractions` bounds and `SimpleTraySideDrawTest` only
  /// exercises at the setter.
  /// **The coordinated tear: a side-draw flow specification and a pumparound on one column.**
  ///
  /// `solveWithColumnTearVariables` has two fast paths - `hasSingleSideDrawFlowSpecificationOnly`
  /// and `hasPumparoundTearVariablesOnly` - and an outer loop over *every* active tear variable
  /// under one residual test. This port carries the two fast paths (`column::tear` and
  /// `column::pumparound`) and refuses the pair; **this row is what says whether the pair has a
  /// state to be held to at all**, which is what the plan asks before any solver is written.
  ///
  /// The state is the class's own multistage fractionator, re-cased to PR and with the two end
  /// pins stated - the same substitution `side_draw_flow_five_tray_liquid_fractionator` makes and
  /// for the same reason.
  static void columnTearRows() {
    columnTearRow("tear_side_draw_alone", false);
    columnTearRow("tear_side_draw_with_pumparound", true);
  }

  /// One row of [`columnTearRows`], with the pumparound present or absent.
  static void columnTearRow(String label, boolean withPumparound) {
    String[] names = new String[] { "propane", "n-butane", "n-pentane", "n-hexane" };
    double[] z = new double[] { 0.20, 0.35, 0.30, 0.15 };
    SystemInterface fluid = new SystemPrEos(333.15, 8.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("tear feed", fluid);
    inlet.setFlowRate(250.0, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn(label, 5, true, true);
    column.addFeedStream(inlet, 3);
    column.setTopPressure(8.0);
    column.setBottomPressure(8.3);
    column.setCondenserTemperature(30.0, "C");
    column.setReboilerTemperature(110.0, "C");
    neqsim.process.equipment.distillation.DistillationColumn.ColumnSideDrawSpecification spec =
        column.addSideDrawFlowSpecification(
            SIDE_DRAW_TRAY_OF(column, 5, 3),
            neqsim.process.equipment.distillation.DistillationColumn.SideDrawPhase.LIQUID,
            20.160854543137464, "kg/hr");
    spec.setTolerance(1.0e-5);
    if (withPumparound) {
      // **The draw tray and the return tray are the specification's own tray and one below**,
      // which is the shape `DistillationColumnCoordinatedFlowTest` drives.
      column.addLiquidPumparound("tear PA", spec.getTrayNumber(), spec.getTrayNumber() - 1, 0.10,
          5.0);
    }
    column.setMaxColumnTearIterations(30);
    column.run();

    System.out.println(label);
    System.out.println("cubic=pr");
    System.out.println("stages=5");
    System.out.println("feed_tray=3");
    System.out.println("has_condenser=true");
    System.out.println("has_reboiler=true");
    System.out.println("side_draw_tray=" + spec.getTrayNumber());
    System.out.println("side_draw_phase=" + spec.getPhase());
    System.out.println("side_draw_target=" + spec.getTargetFlowRate());
    System.out.println("specification_tolerance=" + spec.getTolerance());
    System.out.println("with_pumparound=" + withPumparound);
    System.out.println("pumparound_count=" + column.getPumparounds().size());
    System.out.println("solved=" + column.solved());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("tear_iteration_limit=30");
    System.out.println("tear_iterations=" + column.getLastColumnTearIterationCount());
    System.out.println("tear_residual=" + column.getLastColumnTearResidual());
    System.out.println("tear_converged=" + column.isLastColumnTearConverged());
    System.out.println("tear_rejected_candidates=" + column.getLastColumnTearRejectedCandidateCount());
    System.out.println("tear_rollbacks=" + column.getLastColumnTearRollbackCount());
    System.out.println("tear_inner_iterations=" + column.getLastColumnTearInnerIterationCount());
    System.out.println("tear_candidate_history=" + column.getLastColumnTearCandidateHistory());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    print("feed", inlet);
    print("gas_out", column.getGasOutStream());
    print("liquid_out", column.getLiquidOutStream());
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println("tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + column.getTray(i).getPressure());
      System.out.println(
          "tray" + i + "_gas_n=" + column.getTray(i).getGasOutStream().getFlowRate("mol/sec"));
      System.out.println(
          "tray" + i + "_liquid_n=" + column.getTray(i).getLiquidOutStream().getFlowRate("mol/sec"));
    }
    System.out.println();
  }

  static void sideDrawRows() {
    sideDrawTrayRows();
    sideDrawColumnRows();
  }

  /// The tray rows, which belong in the tray's own capture - see the helper that calls them there.
  static void sideDrawTrayRows() {
    sideDrawTrayRow("side_draw_reference_tray", "methane", 300.0, 10.0, 100.0, 0.0, 0.0, 0.0);
    sideDrawTrayRow("side_draw_gas_quarter_methane", "methane", 300.0, 10.0, 100.0, 0.25, 0.0,
        0.0);
    sideDrawTrayRow("side_draw_liquid_quarter_pentane", "n-pentane", 300.0, 2.0, 100.0, 0.0, 0.25,
        0.0);
    sideDrawTrayRow("side_draw_liquid_with_pumparound_pentane", "n-pentane", 300.0, 2.0, 100.0, 0.0,
        0.20, 0.30);
  }

  static void sideDrawTrayRow(String label, String name, double temperatureK, double pressureBara,
      double kgPerHour, double gasFraction, double liquidFraction, double pumparoundFraction) {
    Stream inlet = pureFeed(name, temperatureK, pressureBara, kgPerHour);

    neqsim.process.equipment.distillation.SimpleTray tray =
        new neqsim.process.equipment.distillation.SimpleTray("side draw tray");
    tray.addStream(inlet);
    if (gasFraction > 0.0) {
      tray.setGasSideDrawFraction(gasFraction);
    }
    if (liquidFraction > 0.0) {
      tray.setLiquidSideDrawFraction(liquidFraction);
    }
    if (pumparoundFraction > 0.0) {
      tray.setLiquidPumparoundDrawFraction(pumparoundFraction);
    }
    tray.run();

    System.out.println(label);
    System.out.println("component=" + name);
    System.out.println("cubic=pr");
    System.out.println("gas_side_draw_fraction=" + tray.getGasSideDrawFraction());
    System.out.println("liquid_side_draw_fraction=" + tray.getLiquidSideDrawFraction());
    System.out.println("pumparound_fraction=" + tray.getLiquidPumparoundDrawFraction());
    print("feed", inlet);
    System.out.println("feed_mol_per_sec=" + inlet.getFlowRate("mol/sec"));
    System.out.println("feed_kg_per_hour=" + inlet.getFlowRate("kg/hr"));
    System.out.println("tray_temperature_K=" + tray.getTemperature());
    System.out.println("tray_pressure_bara=" + tray.getPressure());
    print("gas_out", tray.getGasOutStream());
    print("liquid_out", tray.getLiquidOutStream());
    System.out.println("gas_out_kg_per_hour=" + tray.getGasOutStream().getFlowRate("kg/hr"));
    System.out.println("liquid_out_kg_per_hour=" + tray.getLiquidOutStream().getFlowRate("kg/hr"));
    if (gasFraction > 0.0) {
      print("gas_side_draw", tray.getGasSideDrawStream());
      System.out.println(
          "gas_side_draw_kg_per_hour=" + tray.getGasSideDrawStream().getFlowRate("kg/hr"));
    }
    if (liquidFraction > 0.0) {
      print("liquid_side_draw", tray.getLiquidSideDrawStream());
      System.out.println(
          "liquid_side_draw_kg_per_hour=" + tray.getLiquidSideDrawStream().getFlowRate("kg/hr"));
    }
    if (pumparoundFraction > 0.0) {
      print("liquid_pumparound_draw", tray.getLiquidPumparoundDrawStream());
      System.out.println("liquid_pumparound_draw_kg_per_hour="
          + tray.getLiquidPumparoundDrawStream().getFlowRate("kg/hr"));
    }
    System.out.println();
  }

  /// The column rows, which belong in the column's own capture - see the helper that calls them
  /// there, and the note on why they live there rather than in a capture of their own.
  static void sideDrawColumnRows() {
    // **The class's own three column states**, from tests 2 and 3: a one-tray column with no ends,
    // fed a pure fluid, drawing a tenth and then a quarter of its vapour - and the class's own
    // second half, the same column on liquid n-pentane drawing a quarter of its liquid.
    sideDrawColumnRow("side_draw_column_gas_tenth_methane", new String[] { "methane" },
        new double[] { 1.0 }, 300.0, 10.0, 100.0, 1, 0, false, false, 0.0, 0.0, 10.0, 10.0, 1.0e-6,
        0, 0.10, 0.0, 0.0);
    sideDrawColumnRow("side_draw_column_gas_quarter_methane", new String[] { "methane" },
        new double[] { 1.0 }, 300.0, 10.0, 100.0, 1, 0, false, false, 0.0, 0.0, 10.0, 10.0, 1.0e-6,
        0, 0.25, 0.0, 0.0);
    sideDrawColumnRow("side_draw_column_liquid_quarter_pentane", new String[] { "n-pentane" },
        new double[] { 1.0 }, 300.0, 2.0, 100.0, 1, 0, false, false, 0.0, 0.0, 2.0, 2.0, 1.0e-6, 0,
        0.0, 0.25, 0.0);
    // **The mechanism on a column that has both products**, which is the state the port's case is
    // held to: the model's own binary column, four stages with both ends, one tray drawing. The
    // draw's placement is tray 3, a middle tray, because a fraction on an end is a draw this port
    // refuses by name - its ends are `column::reboiler` and `column::condenser` rather than stages.
    sideDrawColumnRow("side_draw_column_binary_gas_quarter", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, 1000.0, 4, 2, true, true, -20.0, 100.0, 19.0, 20.0,
        1.0e-6, SIDE_DRAW_TRAY, 0.25, 0.0, 0.0);
    // The same column drawing liquid, and drawing a pumparound off the same liquid beside it: the
    // one state where the pair composes, at a column level rather than a tray's.
    sideDrawColumnRow("side_draw_column_binary_liquid_and_pumparound",
        new String[] { "methane", "n-butane" }, new double[] { 0.5, 0.5 }, 300.0, 20.0, 1000.0, 4,
        2, true, true, -20.0, 100.0, 19.0, 20.0, 1.0e-6, SIDE_DRAW_TRAY, 0.0, 0.10, 0.05);

  }

  /// **The pumparound's return, which is a tray-to-tray recycle rather than a draw.**
  ///
  /// `addLiquidPumparound(name, drawTray, returnTray, fraction, drop)` withdraws liquid from one
  /// tray, cools it by `drop` kelvin, and **returns it to another tray** - where the return joins
  /// the tray's inlets and the whole column is solved again, until the return's own flow stops
  /// changing. The draw alone is what this port already carries; the *return* is what it does
  /// not, so the rows here are on the column the port reproduces exactly, and the difference
  /// from `binary_methane_butane_4_stages` is the recycle and nothing else.
  static void pumparoundReturnRows() {
    pumparoundReturnRow("side_draw_pumparound_return_binary", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, 1000.0, 4, 2, true, true, -20.0, 100.0, 19.0, 20.0,
        1.0e-6, 1, 3, 0.10, 5.0);
  }

  static void pumparoundReturnRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, boolean condenser,
      boolean reboiler, double condenserC, double reboilerC, double topBara, double bottomBara,
      double tolerance, int drawTray, int returnTray, double fraction, double temperatureDrop) {
    SystemInterface fluid = new SystemPrEos(feedTemperatureK, feedPressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("pumparound feed", fluid);
    inlet.setFlowRate(kgPerHour, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn(label, stages, reboiler,
            condenser);
    column.addFeedStream(inlet, feedTray);
    column.setCondenserTemperature(condenserC, "C");
    column.setReboilerTemperature(reboilerC, "C");
    column.setTopPressure(topBara);
    column.setBottomPressure(bottomBara);
    column.setTemperatureTolerance(tolerance);
    column.setMaxNumberOfIterations(200, true);
    column.setMaxPumparoundIterations(12);
    column.setPumparoundTolerance(1.0e-4);
    neqsim.process.equipment.distillation.DistillationColumn.ColumnPumparound pumparound =
        column.addLiquidPumparound(label + " PA", drawTray, returnTray, fraction, temperatureDrop);
    column.run();

    System.out.println(label);
    System.out.println("cubic=pr");
    System.out.println("stages=" + stages);
    System.out.println("feed_tray=" + feedTray);
    System.out.println("has_condenser=" + condenser);
    System.out.println("has_reboiler=" + reboiler);
    System.out.println("pumparound_draw_tray=" + pumparound.getDrawTrayNumber());
    System.out.println("pumparound_return_tray=" + pumparound.getReturnTrayNumber());
    System.out.println("pumparound_fraction=" + pumparound.getDrawFraction());
    System.out.println("pumparound_temperature_drop=" + pumparound.getTemperatureDrop());
    System.out.println("solved=" + column.solved());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    System.out.println("pumparound_relative_change=" + column.getLastPumparoundRelativeChange());
    System.out.println("mass_balance_kg_per_hour=" + column.getMassBalance("kg/hr"));
    System.out.println("energy_balance_error=" + column.getEnergyBalanceError());
    if (pumparound.getDrawStream() != null) {
      print("pumparound_draw", pumparound.getDrawStream());
    }
    if (pumparound.getReturnStream() != null) {
      print("pumparound_return", pumparound.getReturnStream());
    }
    System.out.println("pumparound_duty_W=" + pumparound.getDuty());
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println("tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + column.getTray(i).getPressure());
      System.out.println("tray" + i + "_gas_n=" + column.getTray(i).getGasOutStream().getFlowRate("mol/sec"));
      System.out.println("tray" + i + "_liquid_n=" + column.getTray(i).getLiquidOutStream().getFlowRate("mol/sec"));
    }
    print("distillate", column.getGasOutStream());
    print("bottoms", column.getLiquidOutStream());
    System.out.println();
  }

  /// **The side-draw flow specification, on the class's own two states.**
  ///
  /// NeqSim's `solveSingleSideDrawFlowSpecification` is not an outer loop over the column: it is a
  /// **search over cold candidates**, each a deep copy of the whole column at a proposed fraction
  /// and solved from scratch, with only an accepted inner solve allowed to update the controller.
  /// These rows print enough of that mechanism to hold a port to it: the draw's flow and the
  /// specification's own residual, the tear iteration count and its residual, the accepted and
  /// rejected candidate counts, the inner-iteration total, and the candidate history itself.
  static void sideDrawFlowRows() {
    // The class's own `sideDrawFlowSpecificationClosesProductFlow`: a pure-methane feed into a
    // one-tray column with no ends, drawing 25 kg/hr of vapour from tray 0.
    sideDrawFlowRow("side_draw_flow_one_tray_gas_methane", new String[] { "methane" },
        new double[] { 1.0 }, 300.0, 10.0, 100.0, 1, 0, false, false, 0.0, 0.0, 10.0, 10.0,
        "GAS", 25.0, "kg/hr", 1.0e-5, 12);
    // **The same structure on a fluid that flashes to two phases, which the class's own row
    // does not.** A pure-methane feed at 300 K and 10 bara is all vapour, so its one-tray column
    // has no bottom product at all: NeqSim publishes a *zero-flow* liquid with an enthalpy of
    // `-Infinity`, and this port declines to fabricate that state - a stream's molar enthalpy is
    // finite by construction. So the oracle is this substitution: the same one-tray, no-ends,
    // gas-draw shape on the binary feed this model is already held to, which does have a liquid.
    sideDrawFlowRow("side_draw_flow_one_tray_gas_binary", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, 1000.0, 1, 0, false, false, 0.0, 0.0, 20.0,
        20.0, "GAS", 25.0, "kg/hr", 1.0e-5, 12);
    // **The class's own multistage state, re-cased to PR and with one substitution.** Its
    // `createFractionatorColumn` drives both ends by *reflux ratio*, which this port measures as
    // not reproduced - the end's `PVrefluxflash` searches for a state past the mixture's dew point
    // - so the two end pins are stated instead and the draw's target is the flow the class's own
    // test asks for. **That is a substitution and it is labelled as one**: the tear structure is
    // the class's, the state it runs on is not.
    sideDrawFlowRow("side_draw_flow_five_tray_liquid_fractionator",
        new String[] { "propane", "n-butane", "n-pentane", "n-hexane" },
        new double[] { 0.20, 0.35, 0.30, 0.15 }, 333.15, 8.0, 250.0, 5, 3, true, true, 30.0,
        110.0, 8.0, 8.3, "LIQUID", 20.160854543137464, "kg/hr", 1.0e-5, 30);
  }

  static void sideDrawFlowRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, boolean condenser,
      boolean reboiler, double condenserC, double reboilerC, double topBara, double bottomBara,
      String phase, double targetFlow, String flowUnit, double specificationTolerance,
      int tearIterations) {
    SystemInterface fluid = new SystemPrEos(feedTemperatureK, feedPressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("flow spec feed", fluid);
    inlet.setFlowRate(kgPerHour, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn(label, stages, reboiler,
            condenser);
    column.addFeedStream(inlet, feedTray);
    column.setTopPressure(topBara);
    column.setBottomPressure(bottomBara);
    if (condenser) {
      column.setCondenserTemperature(condenserC, "C");
    }
    if (reboiler) {
      column.setReboilerTemperature(reboilerC, "C");
    }
    neqsim.process.equipment.distillation.DistillationColumn.SideDrawPhase sideDrawPhase =
        neqsim.process.equipment.distillation.DistillationColumn.SideDrawPhase.valueOf(phase);
    neqsim.process.equipment.distillation.DistillationColumn.ColumnSideDrawSpecification spec =
        column.addSideDrawFlowSpecification(SIDE_DRAW_TRAY_OF(column, stages, feedTray), sideDrawPhase,
            targetFlow, flowUnit);
    spec.setTolerance(specificationTolerance);
    column.setMaxColumnTearIterations(tearIterations);
    column.run();

    System.out.println(label);
    System.out.println("cubic=pr");
    System.out.println("stages=" + stages);
    System.out.println("feed_tray=" + feedTray);
    System.out.println("has_condenser=" + condenser);
    System.out.println("has_reboiler=" + reboiler);
    System.out.println("side_draw_tray=" + spec.getTrayNumber());
    System.out.println("side_draw_phase=" + spec.getPhase());
    System.out.println("side_draw_target=" + spec.getTargetFlowRate());
    System.out.println("side_draw_unit=" + spec.getFlowUnit());
    System.out.println("specification_tolerance=" + spec.getTolerance());
    System.out.println("tear_iteration_limit=" + tearIterations);
    System.out.println("solved=" + column.solved());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    System.out.println("tear_iterations=" + column.getLastColumnTearIterationCount());
    System.out.println("tear_residual=" + column.getLastColumnTearResidual());
    System.out.println("tear_converged=" + column.isLastColumnTearConverged());
    System.out.println("tear_rejected_candidates=" + column.getLastColumnTearRejectedCandidateCount());
    System.out.println("tear_rollbacks=" + column.getLastColumnTearRollbackCount());
    System.out.println("tear_inner_iterations=" + column.getLastColumnTearInnerIterationCount());
    System.out.println("tear_candidate_history=" + column.getLastColumnTearCandidateHistory());
    System.out.println("draw_fraction="
        + column.getTray(spec.getTrayNumber()).getGasSideDrawFraction());
    System.out.println("spec_last_actual=" + spec.getLastActualFlowRate());
    System.out.println("spec_last_residual=" + spec.getLastRelativeResidual());
    System.out.println("mass_balance_kg_per_hour=" + column.getMassBalance("kg/hr"));
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println("tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_gas_n=" + column.getTray(i).getGasOutStream().getFlowRate("mol/sec"));
      System.out.println("tray" + i + "_liquid_n=" + column.getTray(i).getLiquidOutStream().getFlowRate("mol/sec"));
    }
    print("distillate", column.getGasOutStream());
    print("bottoms", column.getLiquidOutStream());
    System.out.println();
  }

  /// The tray the flow specification is placed on, as this probe's own rows state it.
  static int SIDE_DRAW_TRAY_OF(
      neqsim.process.equipment.distillation.DistillationColumn column, int stages, int feedTray) {
    return stages > 1 ? 3 : feedTray;
  }

  static void sideDrawColumnRow(String label, String[] names, double[] z, double feedTemperatureK,
      double feedPressureBara, double kgPerHour, int stages, int feedTray, boolean condenser,
      boolean reboiler, double condenserC, double reboilerC, double topBara, double bottomBara,
      double tolerance, int drawTray, double gasFraction, double liquidFraction,
      double pumparoundFraction) {
    SystemInterface fluid = new SystemPrEos(feedTemperatureK, feedPressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("column feed", fluid);
    inlet.setFlowRate(kgPerHour, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn("col1", stages, reboiler,
            condenser);
    column.addFeedStream(inlet, feedTray);
    if (condenser) {
      column.setCondenserTemperature(condenserC, "C");
    }
    if (reboiler) {
      column.setReboilerTemperature(reboilerC, "C");
    }
    column.setTopPressure(topBara);
    column.setBottomPressure(bottomBara);
    column.setTemperatureTolerance(tolerance);
    if (gasFraction > 0.0) {
      column.setGasSideDrawFraction(drawTray, gasFraction);
    }
    if (liquidFraction > 0.0) {
      column.setLiquidSideDrawFraction(drawTray, liquidFraction);
    }
    if (pumparoundFraction > 0.0) {
      column.getTray(drawTray).setLiquidPumparoundDrawFraction(pumparoundFraction);
    }
    column.run();

    System.out.println(label);
    System.out.println("cubic=pr");
    System.out.println("stages=" + stages);
    System.out.println("feed_tray=" + feedTray);
    System.out.println("has_condenser=" + condenser);
    System.out.println("has_reboiler=" + reboiler);
    System.out.println("temperature_tolerance=" + tolerance);
    System.out.println("gas_side_draw_fraction=" + gasFraction);
    System.out.println("liquid_side_draw_fraction=" + liquidFraction);
    System.out.println("pumparound_fraction=" + pumparoundFraction);
    print("feed", inlet);
    System.out.println("feed_mol_per_sec=" + inlet.getFlowRate("mol/sec"));
    System.out.println("feed_kg_per_hour=" + inlet.getFlowRate("kg/hr"));
    System.out.println("tray_count=" + column.getNumberOfTrays());
    System.out.println("solved=" + column.solved());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("mass_balance_kg_per_hour=" + column.getMassBalance("kg/hr"));
    System.out.println("energy_balance_error=" + column.getEnergyBalanceError());
    if (condenser) {
      System.out.println("condenser_duty_W=" + column.getCondenser().getDuty());
    }
    if (reboiler) {
      System.out.println("reboiler_duty_W=" + column.getReboiler().getDuty());
    }
    System.out.println("outlet_count=" + column.getOutletStreams().size());
    System.out.println("side_draw_count=" + column.getSideDrawStreams().size());
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println("tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + column.getTray(i).getPressure());
      System.out.println("tray" + i + "_gas_n="
          + column.getTray(i).getGasOutStream().getFlowRate("mol/sec"));
      System.out.println("tray" + i + "_liquid_n="
          + column.getTray(i).getLiquidOutStream().getFlowRate("mol/sec"));
    }
    print("distillate", column.getGasOutStream());
    print("bottoms", column.getLiquidOutStream());
    if (gasFraction > 0.0) {
      print("gas_side_draw",
          column.getSideDrawStream(drawTray,
              neqsim.process.equipment.distillation.DistillationColumn.SideDrawPhase.GAS));
    }
    if (liquidFraction > 0.0) {
      print("liquid_side_draw",
          column.getSideDrawStream(drawTray,
              neqsim.process.equipment.distillation.DistillationColumn.SideDrawPhase.LIQUID));
    }
    if (pumparoundFraction > 0.0) {
      print("liquid_pumparound_draw",
          column.getTray(drawTray).getLiquidPumparoundDrawStream());
    }
    System.out.println();
  }

  /// **The tray the binary rows draw from**, which is a middle one: a fraction on an end is a draw
  /// this port refuses by name, because its ends are `column::reboiler` and `column::condenser`
  /// rather than stages.
  static final int SIDE_DRAW_TRAY = 3;

  /// The class's own side-draw feed: one pure component, its own temperature and pressure, and a
  /// stated mass flow - `SimpleTraySideDrawTest.createMethaneFeed`, in PR.
  static Stream pureFeed(String name, double temperatureK, double pressureBara, double kgPerHour) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBara);
    fluid.addComponent(name, 1.0);
    fluid.setMixingRule(2);
    Stream stream = new Stream("side draw feed", fluid);
    stream.setFlowRate(kgPerHour, "kg/hr");
    stream.run();
    return stream;
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
    // **The fix's own regression state**: an asymmetric ethane/propane split whose material
    // balance the class's test writes out by hand (D = 32.9, B = 67.1 kmol/h).
    shortcutColumnRow("asymmetric_kirkbride_regression", new String[] { "ethane", "propane" },
        new double[] { 0.3, 0.7 }, 303.15, 20.0, 100.0 / 3.6, "ethane", "propane", 0.98, 0.95,
        1.5, null, null);
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

    // The three outlet valve fractions are stored and clamped; the steady split ignores
    // them, and `getMassBalance` is the report that reads the split's outlets.
    threePhaseValves("valve_fractions_010_020_030", names, z,
        feed(names, z, 300.0, 20.0, 1.0));

    // The vessel rows. **Retention is a mechanical-design read**: NLL and NIL are fractions
    // of the internal diameter (0.5 and 0.2 by default) and the effective liquid length is
    // `0.82 * separatorLength`, so the geometry is what makes the numbers exist.
    threePhaseVessel("vessel_d1000_l3000_retention", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 1.0, 3.0, null, null);
    // The level setters clamp at zero; a stated level does not move the retention times,
    // which read the mechanical design's NLL/NIL rather than these fields.
    threePhaseVessel("levels_water_030_oil_060", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 1.0, 3.0, 0.3, 0.6);
    threePhaseVessel("levels_negative_clamped", names, z,
        feed(names, z, 300.0, 20.0, 1.0), 1.0, 3.0, -1.0, -0.5);
    threePhaseVessel("vessel_dry_feed", dry, dryZ, feed(dry, dryZ, 300.0, 20.0, 1.0), 1.0, 3.0,
        null, null);
  }

  static void threePhaseValves(String label, String[] names, double[] z, Stream inlet) {
    neqsim.process.equipment.separator.ThreePhaseSeparator separator =
        new neqsim.process.equipment.separator.ThreePhaseSeparator("sep3", inlet);
    separator.setGasOutletFlowFraction(0.1);
    separator.setOilOutletFlowFraction(0.2);
    separator.setWaterOutletFlowFraction(0.3);
    separator.run();
    System.out.println(label);
    print("feed", inlet);
    printOrEmpty("vapour", separator.getGasOutStream());
    printOrEmpty("light_liquid", separator.getOilOutStream());
    printOrEmpty("heavy_liquid", separator.getWaterOutStream());
    System.out.println("gas_outlet_fraction=" + separator.getGasOutletFlowFraction());
    System.out.println("oil_outlet_fraction=" + separator.getOilOutletFlowFraction());
    System.out.println("water_outlet_fraction=" + separator.getWaterOutletFlowFraction());
    System.out.println("mass_balance_kg_hr=" + separator.getMassBalance("kg/hr"));
    System.out.println();
  }

  static void runThreePhase(String label, String[] names, double[] z, Stream inlet, double dropBara,
      Double heatInputW, String fromPhase, String toPhase, double entrainment) {
    threePhase(label, names, z, inlet, dropBara, heatInputW, fromPhase, toPhase, entrainment, null,
        null, null, null);
  }

  /// A row that also states the vessel: geometry for the retention times, and the levels the
  /// class's setters clamp.
  static void threePhaseVessel(String label, String[] names, double[] z, Stream inlet,
      double diameterM, double lengthM, Double waterLevelM, Double oilLevelM) {
    threePhase(label, names, z, inlet, 0.0, null, null, null, 0.0, diameterM, lengthM,
        waterLevelM, oilLevelM);
  }

  static void threePhase(String label, String[] names, double[] z, Stream inlet, double dropBara,
      Double heatInputW, String fromPhase, String toPhase, double entrainment, Double diameterM,
      Double lengthM, Double waterLevelM, Double oilLevelM) {
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
    if (diameterM != null) {
      separator.setInternalDiameter(diameterM);
    }
    if (lengthM != null) {
      separator.setSeparatorLength(lengthM);
    }
    if (waterLevelM != null) {
      separator.setWaterLevel(waterLevelM);
    }
    if (oilLevelM != null) {
      separator.setOilLevel(oilLevelM);
    }
    separator.run();
    System.out.println(label);
    print("feed", inlet);
    printOrEmpty("vapour", separator.getGasOutStream());
    printOrEmpty("light_liquid", separator.getOilOutStream());
    printOrEmpty("heavy_liquid", separator.getWaterOutStream());
    System.out.println("mass_balance_kg_hr=" + separator.getMassBalance("kg/hr"));
    if (diameterM != null || waterLevelM != null || oilLevelM != null) {
      System.out.println("water_level_m=" + separator.getWaterLevel());
      System.out.println("oil_level_m=" + separator.getOilLevel());
      System.out.println("oil_thickness_m=" + separator.getOilThickness());
    }
    if (diameterM != null) {
      SystemInterface system = separator.getThermoSystem();
      System.out.println("nll_m=" + separator.getMechanicalDesign().getNLL());
      System.out.println("nil_m=" + separator.getMechanicalDesign().getNIL());
      System.out.println("effective_liquid_length_m="
          + separator.getMechanicalDesign().getEffectiveLengthLiquid());
      for (String phase : new String[] {"oil", "aqueous"}) {
        if (system.hasPhaseType(phase)) {
          System.out.println(phase + "_density_kg_m3="
              + system.getPhase(phase).getDensity("kg/m3"));
          System.out.println(phase + "_flow_m3_hr="
              + system.getPhase(phase).getFlowRate("m3/hr"));
        }
      }
      System.out.println("oil_retention_time_min=" + separator.calcOilRetentionTime());
      System.out.println("water_retention_time_min=" + separator.calcWaterRetentionTime());
      System.out.println("interface_settling_time_min=" + separator.calcInterfaceSettlingTime());
    }
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

  /// **`Tank.autoSize`'s residence-time sizing and the three constraints it initializes.**
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

  // **Nine rows, one per branch of the class's lazy capacity map.** The splitter's map is
  // created by `setDesignPressureDrop`, `setMaxDesignVelocity`, or `setCapacityAnalysisEnabled
  // (true)` on an empty map - not by the constructor and not by `setHeaderDiameter` - so the
  // rows walk each trigger, both non-positive suppressions, the 9.99 cap, and a two-phase feed.
  // The last uses a one-outlet split (`[1.0]`): the class's branch path subtracts zero moles,
  // which keeps the measurement away from its known-broken `clone/init(0)/addComponent` route.
  static void splitterCapacity() {
    splitterCapacityRow("no_setters", new double[] { 0.3, 0.7 }, null, null, null, false,
        10.0, 300.0);
    splitterCapacityRow("header_alone", new double[] { 0.3, 0.7 }, 0.02, null, null, false,
        10.0, 300.0);
    splitterCapacityRow("enabled_defaults", new double[] { 0.3, 0.7 }, null, null, null, true,
        10.0, 300.0);
    splitterCapacityRow("sized_header_design_drop", new double[] { 0.3, 0.7 }, 0.02, 0.02, null,
        true, 10.0, 300.0);
    splitterCapacityRow("tight_velocity_disabled", new double[] { 0.3, 0.7 }, 0.02, null, 0.1,
        false, 10.0, 300.0);
    splitterCapacityRow("zero_design_drop", new double[] { 0.3, 0.7 }, 0.02, 0.0, null, false,
        10.0, 300.0);
    splitterCapacityRow("zero_max_velocity", new double[] { 0.3, 0.7 }, 0.02, null, 0.0, false,
        10.0, 300.0);
    splitterCapacityRow("cap_9_99", new double[] { 0.3, 0.7 }, 0.02, null, 0.01, false, 10.0,
        300.0);
    splitterCapacityRow("two_phase_one_outlet", new double[] { 1.0 }, 0.05, null, null, true,
        1.5, 300.0);
  }

  static void splitterCapacityRow(String name, double[] factors, Double headerDiameter,
      Double designPressureDropBar, Double maxDesignVelocity, boolean analysisEnabled,
      double pressureBara, double temperatureK) {
    Stream inlet = feed(new String[] { "n-butane", "n-pentane" }, new double[] { 0.6, 0.4 },
        temperatureK, pressureBara, 1.0);
    neqsim.process.equipment.splitter.Splitter splitter =
        new neqsim.process.equipment.splitter.Splitter("sp1", inlet);
    splitter.setSplitFactors(factors);

    neqsim.process.mechanicaldesign.splitter.SplitterMechanicalDesign design =
        (neqsim.process.mechanicaldesign.splitter.SplitterMechanicalDesign) splitter
            .getMechanicalDesign();
    if (headerDiameter != null) {
      design.setHeaderDiameter(headerDiameter);
    }
    if (designPressureDropBar != null) {
      splitter.setDesignPressureDrop(designPressureDropBar);
    }
    if (maxDesignVelocity != null) {
      splitter.setMaxDesignVelocity(maxDesignVelocity);
    }
    if (analysisEnabled) {
      splitter.setCapacityAnalysisEnabled(true);
    }
    splitter.run();

    System.out.println("case=" + name);
    print("feed", inlet);
    System.out.println("feed_phases=" + inlet.getThermoSystem().getNumberOfPhases());
    for (int i = 0; i < inlet.getThermoSystem().getNumberOfPhases(); i++) {
      System.out.println("feed_phase" + i + "_type=" + inlet.getThermoSystem().getPhase(i).getType()
          + " feed_phase" + i + "_beta=" + inlet.getThermoSystem().getPhase(i).getBeta());
    }
    System.out.println("feed_volume_flow_m3_s=" + inlet.getFlowRate("m3/sec"));
    System.out.println("header_diameter_m=" + design.getHeaderDiameter());
    System.out.println("design_pressure_drop_pa=" + splitter.getDesignPressureDrop() * 1.0e5);
    System.out.println("max_design_velocity_m_s=" + splitter.getMaxDesignVelocity());
    System.out.println("capacity_analysis_enabled=" + splitter.isCapacityAnalysisEnabled());

    // The map is unmodifiable: iterate it, never write to it. The pressure-drop constraint's
    // values are bar, so they are printed in pascals under the spec's own unit.
    Map<String, neqsim.process.equipment.capacity.CapacityConstraint> constraints =
        splitter.getCapacityConstraints();
    System.out.println("constraint_count=" + constraints.size());
    for (Map.Entry<String, neqsim.process.equipment.capacity.CapacityConstraint> entry :
        constraints.entrySet()) {
      String constraintName = entry.getKey();
      neqsim.process.equipment.capacity.CapacityConstraint constraint = entry.getValue();
      if ("pressureDrop".equals(constraintName)) {
        System.out.println("pressure_drop_value_pa=" + constraint.getCurrentValue() * 1.0e5);
        System.out.println("pressure_drop_design_pa=" + constraint.getDesignValue() * 1.0e5);
        System.out.println("pressure_drop_utilization=" + constraint.getUtilization());
      } else if ("velocity".equals(constraintName)) {
        System.out.println("velocity_value_m_s=" + constraint.getCurrentValue());
        System.out.println("velocity_design_m_s=" + constraint.getDesignValue());
        System.out.println("velocity_utilization=" + constraint.getUtilization());
      }
      System.out.println("constraint_" + constraintName + "_type=" + constraint.getType());
      System.out.println("constraint_" + constraintName + "_enabled=" + constraint.isEnabled());
    }

    System.out.println("max_utilization=" + splitter.getMaxUtilization());
    System.out.println("capacity_exceeded=" + splitter.isCapacityExceeded());
    neqsim.process.equipment.capacity.CapacityConstraint bottleneck = splitter
        .getBottleneckConstraint();
    System.out.println(
        "bottleneck_constraint=" + (bottleneck == null ? "none" : bottleneck.getName()));
    System.out.println();
  }

  // `SplitterMechanicalDesign.calcDesign` is the caller's explicit call: `run` never reaches
  // it, and it creates no capacity constraint. Each block moves one of its inputs so the
  // capture holds the sizing, the wall formula's units, the discarded joint efficiency and
  // the capacity map's read of the sized diameter.
  static void splitterMechanicalDesign() {
    splitterDesignRow("default_design", 1.0, 10.0, 300.0, new double[] {0.3, 0.7}, null, null, null,
        null, null, false);
    splitterDesignRow("twenty_mol", 20.0, 10.0, 300.0, new double[] {0.3, 0.7}, null, null, null,
        null, null, false);
    splitterDesignRow("design_velocity_ten", 69.0, 10.0, 300.0, new double[] {0.3, 0.7}, 10.0, null,
        null, null, null, false);
    splitterDesignRow("high_pressure", 20.0, 10.0, 300.0, new double[] {0.3, 0.7}, null, 300.0, null,
        null, null, false);
    splitterDesignRow("low_tensile", 20.0, 10.0, 300.0, new double[] {0.3, 0.7}, null, null, 300.0,
        null, null, false);
    splitterDesignRow("corrosion", 20.0, 10.0, 300.0, new double[] {0.3, 0.7}, null, null, null, null,
        1.5, false);
    splitterDesignRow("joint_efficiency", 20.0, 10.0, 300.0, new double[] {0.3, 0.7}, null, null, null,
        0.7, null, false);
    splitterDesignRow("two_phase", 1.0, 1.5, 300.0, new double[] {1.0}, null, null, null, null, null,
        false);
    splitterDesignRow("design_with_capacity", 1.0, 10.0, 300.0, new double[] {0.3, 0.7}, null, null,
        null, null, null, true);
  }

  static void splitterDesignRow(String name, double flow, double pressureBara, double temperatureK,
      double[] factors, Double designVelocity, Double maxOperationPressureBara,
      Double tensileStrengthMpa, Double jointEfficiency, Double corrosionAllowanceMm,
      boolean capacityEnabled) {
    Stream inlet = feed(new String[] {"n-butane", "n-pentane"}, new double[] {0.6, 0.4}, temperatureK,
        pressureBara, flow);
    neqsim.process.equipment.splitter.Splitter splitter =
        new neqsim.process.equipment.splitter.Splitter("sp1", inlet);
    splitter.setSplitFactors(factors);
    neqsim.process.mechanicaldesign.splitter.SplitterMechanicalDesign design =
        (neqsim.process.mechanicaldesign.splitter.SplitterMechanicalDesign) splitter
            .getMechanicalDesign();
    if (designVelocity != null) {
      design.setDesignVelocity(designVelocity);
    }
    if (maxOperationPressureBara != null) {
      design.setMaxOperationPressure(maxOperationPressureBara, "bara");
    }
    if (tensileStrengthMpa != null) {
      design.setTensileStrength(tensileStrengthMpa);
    }
    if (jointEfficiency != null) {
      design.setJointEfficiency(jointEfficiency);
    }
    if (corrosionAllowanceMm != null) {
      design.setCorrosionAllowance(corrosionAllowanceMm);
    }
    if (capacityEnabled) {
      splitter.setCapacityAnalysisEnabled(true);
    }
    splitter.run();
    design.calcDesign();

    System.out.println("case=" + name);
    print("feed", inlet);
    System.out.println("feed_phases=" + inlet.getThermoSystem().getNumberOfPhases());
    System.out.println("feed_volume_flow_m3_s=" + inlet.getFlowRate("m3/sec"));
    System.out.println("mechanical_design=true");
    System.out.println("design_velocity_m_s=" + design.getDesignVelocity());
    System.out.println("max_operation_pressure_pa=" + design.getMaxOperationPressure() * 1.0e5);
    System.out.println("tensile_strength_pa=" + design.getTensileStrength() * 1.0e6);
    System.out.println("joint_efficiency=" + design.getJointEfficiency());
    System.out.println("corrosion_allowance_mm=" + design.getCorrosionAllowance());
    System.out.println("design_header_diameter_m=" + design.getHeaderDiameter());
    System.out.println("design_branch_diameter_m=" + design.getBranchDiameter());
    System.out.println("design_header_length_m=" + design.getHeaderLength());
    System.out.println("design_wall_thickness_m=" + design.getHeaderWallThickness() / 1000.0);
    System.out.println("manifold_volume_m3=" + design.getManifoldVolume());
    System.out.println("design_total_pressure_drop_pa=" + design.getTotalPressureDrop() * 1.0e5);
    System.out.println("header_weight_kg=" + design.weigthVesselShell);
    System.out.println("branch_weight_kg=" + design.weightPiping);
    System.out.println("total_weight_kg=" + design.getWeightTotal());

    if (capacityEnabled) {
      // The map reads the sized diameter live, so these are the post-design answers.
      Map<String, neqsim.process.equipment.capacity.CapacityConstraint> constraints =
          splitter.getCapacityConstraints();
      System.out.println("constraint_count=" + constraints.size());
      for (Map.Entry<String, neqsim.process.equipment.capacity.CapacityConstraint> entry :
          constraints.entrySet()) {
        String constraintName = entry.getKey();
        neqsim.process.equipment.capacity.CapacityConstraint constraint = entry.getValue();
        if ("velocity".equals(constraintName)) {
          System.out.println("velocity_value_m_s=" + constraint.getCurrentValue());
          System.out.println("velocity_design_m_s=" + constraint.getDesignValue());
          System.out.println("velocity_utilization=" + constraint.getUtilization());
        }
      }
      System.out.println("max_utilization=" + splitter.getMaxUtilization());
      System.out.println("capacity_exceeded=" + splitter.isCapacityExceeded());
      neqsim.process.equipment.capacity.CapacityConstraint bottleneck = splitter
          .getBottleneckConstraint();
      System.out.println(
          "bottleneck_constraint=" + (bottleneck == null ? "none" : bottleneck.getName()));
    }
    System.out.println();
  }

  static void mixer() {
    // **Two inlets at different pressures and different compositions.** `Mixer.run` sets
    // the outlet to the lowest active inlet pressure and flashes the joined fluid to the
    // flow-weighted enthalpy, so both rules are exercised: a composition that is not
    // either feed's, and a pressure that is neither inlet's when the outlet is stated.
    // The palette declares the outlet pressure optional, so the second row states one.
    double[][] pairs = new double[][] {{10.0, 6.0}, {10.0, 6.0}, {10.0, 9.6}};
    Double[] specifieds = new Double[] { null, 8.0, null };
    for (int row = 0; row < pairs.length; row++) {
      Double specified = specifieds[row];
      Stream first = feed(new String[] { "n-butane", "n-pentane" }, new double[] { 0.6, 0.4 },
          300.0, pairs[row][0], 1.0);
      Stream second = feed(new String[] { "n-butane", "n-pentane" }, new double[] { 0.3, 0.7 },
          300.0, pairs[row][1], 2.0);

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
      System.out.println("inlet_pressure_spread_bara=" + mixer.getInletPressureSpread());
      System.out.println("pressure_mismatch=" + mixer.isPressureMismatch());
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

  // === process.separator capacity ==============================================

  /// The five steady capacity getters on the vessel `Separator.run` leaves behind.
  ///
  /// **`getCapacityUtilization` and `getGasLoadFactor()` are two density rules and not
  /// one.** The sizing family (`getMaxAllowableGasVelocity`, and the utilisation over it)
  /// reads the liquid's own density and falls back to `1000.0` when there is no liquid - no
  /// floor. The load factor reads the same phases and then substitutes `1000.0` whenever
  /// `rho_l - rho_g < 10.0`. A near-dry state therefore makes the two disagree by design,
  /// which the near-dry rows are here to pin.
  ///
  /// **The PR rows are the oracle**, because azoth's databank route resolves
  /// Peng-Robinson. The SRK rows re-case `SeparatorTest`'s own states as evidence, and the
  /// scrubber row prints the defaults the port deliberately does not adopt.
  static void separatorCapacityRows() {
    // The two-phase fluid of `process_separator.tsv`, so the capacity rows sit on a state
    // the port already carries.
    separatorCapacityRow("pr_two_phase_vertical_d1000_k011", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.8, 0.0, 0.0);
    separatorCapacityRow("pr_two_phase_horizontal_d1000_k011", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "horizontal", 0.8, 0.0, 0.0);
    // Issue #3154: vertical capacity is the full cross-section, whatever the design liquid
    // fraction says; horizontal derates it by `1 - Fg`.
    separatorCapacityRow("pr_two_phase_vertical_fg08", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.8, 0.0, 0.0);
    separatorCapacityRow("pr_two_phase_vertical_fg02", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.2, 0.0, 0.0);
    separatorCapacityRow("pr_two_phase_horizontal_fg02", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "horizontal", 0.2, 0.0, 0.0);
    // Dry gas: the sizing rule's `1000.0` fallback, and `getMaxAllowableGasVelocity`'s
    // `50.0` only applies where there is no gas at all - a state utilisation answers `0.0`
    // for before it reaches the velocity.
    separatorCapacityRow("pr_dry_gas_vertical_d1000_k011", false, false, null,
        new String[] {"methane"}, new double[] {1.0}, 300.0, 50.0, 1.0, null, 1.0, 0.11,
        "vertical", 0.8, 0.0, 0.0);
    // The pin's near-dry scrubber state: a trace heavy end at 90 bara, where the two-phase
    // densities sit close enough for the load factor's `10.0` floor to decide.
    separatorCapacityRow("pr_near_dry_96_3_005_vertical_d3000", false, true, null,
        new String[] {"methane", "ethane", "n-heptane"}, new double[] {96.0, 3.0, 0.05}, 298.15,
        90.0, null, 400000.0, 3.0, 0.11, "vertical", 0.8, 0.0, 0.0);
    // A vessel far too small for its flow, so utilisation crosses one and the overload flag
    // is reachable.
    separatorCapacityRow("pr_overloaded_vertical_d020_k011", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 50.0, null,
        0.2, 0.11, "vertical", 0.8, 0.0, 0.0);
    // Entrainment moves moles into the liquid before the getters read, so the liquid density
    // is the entrained composition's.
    separatorCapacityRow("pr_entrained_gas_in_liquid_005_vertical_d1000", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.8, 0.05, 0.0);
    // The metric is at the outlet state, which a pressure drop moves.
    separatorCapacityRow("pr_pressure_drop_outlet_state_vertical_d1000", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.8, 0.0, 2.0);

    // SRK evidence, re-casing the pin's own tests.
    separatorCapacityRow("srk_two_phase_vertical_d1000_k011", true, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.8, 0.0, 0.0);
    separatorCapacityRow("srk_near_dry_bounded_vertical_d3000", true, true, null,
        new String[] {"methane", "ethane", "n-heptane"}, new double[] {96.0, 3.0, 0.05}, 298.15,
        90.0, null, 400000.0, 3.0, 0.11, "vertical", 0.8, 0.0, 0.0);
    separatorCapacityRow("srk_orientation_vertical_d3000", true, true, null,
        new String[] {"methane", "ethane", "n-heptane"}, new double[] {90.0, 6.0, 4.0}, 298.15,
        80.0, null, 300000.0, 3.0, 0.11, "vertical", 0.8, 0.0, 0.0);
    separatorCapacityRow("srk_orientation_horizontal_d3000", true, true, null,
        new String[] {"methane", "ethane", "n-heptane"}, new double[] {90.0, 6.0, 4.0}, 298.15,
        80.0, null, 300000.0, 3.0, 0.11, "horizontal", 0.8, 0.0, 0.0);

    // `GasScrubber` at minimum setup - a diameter and nothing else - which is where the
    // class's own vertical orientation and 0.1 design liquid level come from. The port
    // adopts neither default and refuses a half-stated vessel.
    scrubberCapacityDefaultsRow("srk_gas_scrubber_minimum_setup_d1200", true, true, null,
        new String[] {"methane", "ethane", "n-hexane"}, new double[] {90.0, 5.0, 5.0}, 303.15,
        50.0, null, 50000.0, 1.2);
  }

  @SuppressWarnings({"checkstyle:MethodLength", "checkstyle:ParameterNumber"})
  static void separatorCapacityRow(String label, boolean srk, boolean multiPhase,
      Double heatInputW, String[] names, double[] amounts, double temperature, double pressureBar,
      Double molPerSecond, Double kgPerHour, double diameter, double kFactor, String orientation,
      double designLiquidLevelFraction, double gasInLiquid, double pressureDropBar) {
    Stream inlet = capacityFeed(srk, multiPhase, names, amounts, temperature, pressureBar,
        molPerSecond, kgPerHour);
    neqsim.process.equipment.separator.Separator separator =
        new neqsim.process.equipment.separator.Separator("sep1", inlet);
    if (heatInputW != null) {
      separator.setHeatInput(heatInputW);
    }
    if (pressureDropBar != 0.0) {
      separator.setPressureDrop(pressureDropBar);
    }
    if (gasInLiquid != 0.0) {
      separator.setEntrainment(gasInLiquid, "mole", "feed", "gas", "liquid");
    }
    separator.setInternalDiameter(diameter);
    separator.setDesignGasLoadFactor(kFactor);
    separator.setOrientation(orientation);
    separator.setDesignLiquidLevelFraction(designLiquidLevelFraction);
    separator.run();
    printCapacity(label, srk, inlet, separator, diameter, kFactor, orientation,
        designLiquidLevelFraction, gasInLiquid, pressureDropBar);
  }

  @SuppressWarnings({"checkstyle:MethodLength", "checkstyle:ParameterNumber"})
  static void scrubberCapacityDefaultsRow(String label, boolean srk, boolean multiPhase,
      Double heatInputW, String[] names, double[] amounts, double temperature, double pressureBar,
      Double molPerSecond, Double kgPerHour, double diameter) {
    Stream inlet = capacityFeed(srk, multiPhase, names, amounts, temperature, pressureBar,
        molPerSecond, kgPerHour);
    neqsim.process.equipment.separator.GasScrubber scrubber =
        new neqsim.process.equipment.separator.GasScrubber("gs1", inlet);
    scrubber.setInternalDiameter(diameter);
    scrubber.run();
    printCapacity(label, srk, inlet, scrubber, diameter, scrubber.getDesignGasLoadFactor(),
        scrubber.getOrientation(), scrubber.getDesignLiquidLevelFraction(), 0.0, 0.0);
  }

  static Stream capacityFeed(boolean srk, boolean multiPhase, String[] names, double[] amounts,
      double temperature, double pressureBar, Double molPerSecond, Double kgPerHour) {
    SystemInterface fluid =
        srk ? new SystemSrkEos(temperature, pressureBar) : new SystemPrEos(temperature, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], amounts[i]);
    }
    fluid.setMixingRule(2);
    if (multiPhase) {
      fluid.setMultiPhaseCheck(true);
    }
    Stream stream = new Stream("feed", fluid);
    if (molPerSecond != null) {
      stream.setFlowRate(molPerSecond, "mol/sec");
    } else {
      stream.setFlowRate(kgPerHour, "kg/hr");
    }
    stream.run();
    return stream;
  }

  /// The five getters, the phases they read, and the geometry they divide by.
  ///
  /// `gas_area_m2` is printed as the instrument computes it - full cross-section vertical,
  /// `(1 - Fg)` of it horizontal - rather than read back, because the class has no getter
  /// for it. `max_allowable_gas_flow_m3_s` is the getter itself.
  @SuppressWarnings({"checkstyle:MethodLength", "checkstyle:ParameterNumber"})
  static void printCapacity(String label, boolean srk, Stream inlet,
      neqsim.process.equipment.separator.Separator separator, double diameter,
      double kFactor, String orientation, double designLiquidLevelFraction, double gasInLiquid,
      double pressureDropBar) {
    SystemInterface system = separator.getThermoSystem();
    // **`run` leaves the densities stale where it entrained**: the phase labels are not
    // revisited but `initPhysicalProperties` recomputes the density of the oil that
    // absorbed the carried moles, and every capacity getter calls it. Doing it here means
    // the printed phase densities are the ones the getters below actually read.
    system.initPhysicalProperties();
    System.out.println(label);
    System.out.println("eos=" + (srk ? "srk" : "pr"));
    System.out.println("temperature_k=" + inlet.getTemperature("K"));
    System.out.println("pressure_bara=" + inlet.getPressure("bara"));
    System.out.println("flow_mol_s=" + inlet.getFlowRate("mol/sec"));
    System.out.println("flow_kg_hr=" + inlet.getFlowRate("kg/hr"));
    StringBuilder z = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      z.append(i == 0 ? "" : ",").append(system.getPhase(0).getComponent(i).getName())
          .append(":").append(system.getMolarComposition()[i]);
    }
    System.out.println("z=" + z);
    System.out.println("internal_diameter_m=" + diameter);
    System.out.println("design_gas_load_factor_m_per_s=" + kFactor);
    System.out.println("orientation=" + orientation);
    System.out.println("design_liquid_level_fraction=" + designLiquidLevelFraction);
    System.out.println("gas_in_liquid=" + gasInLiquid);
    System.out.println("pressure_drop_bara=" + pressureDropBar);
    System.out.println("phase_count=" + system.getNumberOfPhases());
    StringBuilder types = new StringBuilder();
    StringBuilder densities = new StringBuilder();
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      types.append(i == 0 ? "" : ",").append(system.getPhase(i).getType().name());
      densities.append(i == 0 ? "" : ",")
          .append(system.getPhase(i).getPhysicalProperties().getDensity());
    }
    System.out.println("phase_types=" + types);
    System.out.println("phase_densities_kg_m3=" + densities);
    System.out.println("has_gas=" + system.hasPhaseType("gas"));
    System.out.println("has_oil=" + system.hasPhaseType("oil"));
    System.out.println("has_aqueous=" + system.hasPhaseType("aqueous"));
    System.out.println("gas_density_kg_m3=" + system.getPhase(0).getPhysicalProperties().getDensity());
    System.out.println("liquid_density_kg_m3="
        + (system.getNumberOfPhases() < 2 ? "none"
            : system.getPhase(1).getPhysicalProperties().getDensity()));
    System.out.println("gas_volumetric_flow_m3_s=" + system.getPhase(0).getFlowRate("m3/sec"));
    double area = Math.PI * diameter * diameter / 4.0;
    double gasArea = orientation.equals("horizontal") ? area * (1.0 - designLiquidLevelFraction)
        : area;
    System.out.println("gas_area_m2=" + gasArea);
    System.out.println("max_allowable_gas_velocity_m_per_s="
        + separator.getMaxAllowableGasVelocity());
    System.out.println("max_allowable_gas_flow_m3_s=" + separator.getMaxAllowableGasFlowRate());
    System.out.println("gas_superficial_velocity_m_per_s=" + separator.getGasSuperficialVelocity());
    System.out.println("gas_load_factor_m_per_s=" + separator.getGasLoadFactor());
    System.out.println("capacity_utilization=" + separator.getCapacityUtilization());
    System.out.println("is_overloaded=" + separator.isOverloaded());
    System.out.println();
  }

  // === process.separator sizing and derating ==================================

  /// `Separator`'s sizing and derating getters, on the state `run` leaves behind.
  ///
  /// **The phase-number overload is not the no-argument one.** It divides the Souders-Brown
  /// term by the horizontal gas-area fraction `1 - liquidVolume / separatorVolume` and reads
  /// the stated phase's density; the no-argument form reads phase 1 and the design liquid
  /// level fraction. Both are printed on every row, so a state where they agree (vertical, or
  /// a level equal to the design fraction) does not hide the difference.
  ///
  /// **The derating is the interphase surface tension's**: `1 - 0.5 (10e-3 - sigma)/10e-3`
  /// below `10e-3` N/m, applied to the superficial velocity. The rows include a near-critical
  /// light gas where the tension falls under the threshold, and print the tension each
  /// overload reads.
  ///
  /// **`sizeFromFlow` is a closed form and it mutates the vessel.** It runs last, after every
  /// getter has read the stated diameter, and the answer is printed beside that diameter.
  static void separatorSizingRows() {
    separatorSizingRow("pr_two_phase_vertical", false, false,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, 5.0,
        0.5, "vertical", 0.8, 0.11, 1.0);
    separatorSizingRow("pr_two_phase_horizontal", false, false,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, 5.0,
        0.5, "horizontal", 0.8, 0.11, 1.0);
    separatorSizingRow("pr_two_phase_horizontal_level_02", false, false,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, 5.0,
        0.2, "horizontal", 0.2, 0.11, 1.0);
    // A near-critical light gas: the interphase tension is the row's measurement, and the
    // derating branch is reachable only where it falls under 10 mN/m.
    separatorSizingRow("pr_near_critical_derating", false, false,
        new String[] {"methane", "ethane"}, new double[] {0.9, 0.1}, 200.0, 40.0, 1.0, 5.0,
        0.5, "horizontal", 0.8, 0.11, 1.0);
    separatorSizingRow("srk_near_critical_derating", true, false,
        new String[] {"methane", "ethane"}, new double[] {0.9, 0.1}, 200.0, 40.0, 1.0, 5.0,
        0.5, "horizontal", 0.8, 0.11, 1.0);
    // Three phases: the indexed overloads read the oil and the aqueous density separately.
    separatorSizingRow("pr_three_phase_horizontal", false, true,
        new String[] {"methane", "n-butane", "water"}, new double[] {0.5, 0.4, 0.1}, 300.0, 20.0,
        1.0, 5.0, 0.5, "horizontal", 0.8, 0.11, 1.0);
    // A dry gas: the sizing rule's 1000 kg/m3 liquid fallback and the derating's own.
    separatorSizingRow("pr_dry_gas_vertical", false, false, new String[] {"methane"},
        new double[] {1.0}, 300.0, 50.0, 1.0, 5.0, 0.5, "vertical", 0.8, 0.11, 1.0);
    // No gas at all: `sizeFromFlow` returns without sizing, which is the class's own branch.
    separatorSizingRow("pr_liquid_only_no_sizing", false, false, new String[] {"n-butane"},
        new double[] {1.0}, 300.0, 20.0, 1.0, 5.0, 0.5, "vertical", 0.8, 0.11, 1.0);
    // The safety factor is the one input `sizeFromFlow` takes; the same state at 1.5.
    separatorSizingRow("pr_sizing_safety_15", false, false,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, 5.0,
        0.5, "horizontal", 0.8, 0.11, 1.5);
  }

  @SuppressWarnings({"checkstyle:MethodLength", "checkstyle:ParameterNumber"})
  static void separatorSizingRow(String label, boolean srk, boolean multiPhase, String[] names,
      double[] amounts, double temperature, double pressureBar, double diameter, double length,
      double liquidLevelFraction, String orientation, double designLiquidLevelFraction,
      double kFactor, double safetyFactor) {
    Stream inlet = capacityFeed(srk, multiPhase, names, amounts, temperature, pressureBar, 1.0,
        null);
    neqsim.process.equipment.separator.Separator separator =
        new neqsim.process.equipment.separator.Separator("sep1", inlet);
    separator.setOrientation(orientation);
    separator.setInternalDiameter(diameter);
    separator.setSeparatorLength(length);
    separator.setLiquidLevel(liquidLevelFraction);
    separator.setDesignLiquidLevelFraction(designLiquidLevelFraction);
    separator.setDesignGasLoadFactor(kFactor);
    separator.run();
    printSizing(label, srk, inlet, separator, diameter, length, liquidLevelFraction, orientation,
        designLiquidLevelFraction, kFactor, safetyFactor);
  }

  @SuppressWarnings({"checkstyle:MethodLength", "checkstyle:ParameterNumber"})
  static void printSizing(String label, boolean srk, Stream inlet,
      neqsim.process.equipment.separator.Separator separator, double diameter, double length,
      double liquidLevelFraction, String orientation, double designLiquidLevelFraction,
      double kFactor, double safetyFactor) {
    SystemInterface system = separator.getThermoSystem();
    system.initPhysicalProperties();
    System.out.println(label);
    System.out.println("eos=" + (srk ? "srk" : "pr"));
    System.out.println("temperature_k=" + inlet.getTemperature("K"));
    System.out.println("pressure_bara=" + inlet.getPressure("bara"));
    System.out.println("flow_mol_s=" + inlet.getFlowRate("mol/sec"));
    StringBuilder z = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      z.append(i == 0 ? "" : ",").append(system.getPhase(0).getComponent(i).getName())
          .append(":").append(system.getMolarComposition()[i]);
    }
    System.out.println("z=" + z);
    System.out.println("orientation=" + orientation);
    System.out.println("internal_diameter_m=" + diameter);
    System.out.println("separator_length_m=" + length);
    System.out.println("liquid_level_fraction=" + liquidLevelFraction);
    System.out.println("design_liquid_level_fraction=" + designLiquidLevelFraction);
    System.out.println("design_gas_load_factor_m_per_s=" + kFactor);
    System.out.println("safety_factor=" + safetyFactor);
    System.out.println("has_geometry=" + separator.hasGeometry());
    System.out.println("phase_count=" + system.getNumberOfPhases());
    StringBuilder types = new StringBuilder();
    StringBuilder densities = new StringBuilder();
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      types.append(i == 0 ? "" : ",").append(system.getPhase(i).getType().name());
      densities.append(i == 0 ? "" : ",")
          .append(system.getPhase(i).getPhysicalProperties().getDensity());
    }
    System.out.println("phase_types=" + types);
    System.out.println("phase_densities_kg_m3=" + densities);
    System.out.println("has_gas=" + system.hasPhaseType("gas"));
    System.out.println("has_oil=" + system.hasPhaseType("oil"));
    System.out.println("has_aqueous=" + system.hasPhaseType("aqueous"));
    // The tension each derating overload reads: the no-argument one reads the (0, 1) pair.
    double classSigma = 0.0;
    for (int i = 0; i + 1 < system.getNumberOfPhases(); i++) {
      double sigma = system.getInterphaseProperties().getSurfaceTension(i, i + 1);
      if (i == 0) {
        classSigma = sigma;
      }
      System.out.println("surface_tension_" + i + "_" + (i + 1) + "_n_per_m=" + sigma);
    }
    // **The class's interphase model is a clone of the pre-flash system.** It carries zero
    // densities, so the Parachor it computes is 0.0 whatever the live state - and the
    // derating is therefore `1 - 0.5` on every row. A fresh model over the live system
    // measures what the formula would answer if it read this state, which is the difference
    // the spec records.
    double liveSigma = Double.NaN;
    if (system.getNumberOfPhases() >= 2) {
      try {
        liveSigma =
            new neqsim.physicalproperties.interfaceproperties.surfacetension.ParachorSurfaceTension(
                system).calcSurfaceTension(0, 1);
      } catch (Exception ex) {
        liveSigma = Double.NaN;
      }
    }
    System.out.println("live_parachor_surface_tension_n_per_m=" + liveSigma);
    System.out.println("class_derating_factor="
        + (classSigma < 10.0e-3 ? 1.0 - 0.5 * (10.0e-3 - classSigma) / 10.0e-3 : 1.0));
    System.out.println("intended_derating_factor="
        + (!Double.isNaN(liveSigma) && liveSigma < 10.0e-3
            ? 1.0 - 0.5 * (10.0e-3 - liveSigma) / 10.0e-3
            : 1.0));
    double area = Math.PI * diameter * diameter / 4.0;
    double vesselVolume = area * length;
    double liquidVolume = separator.calcLiquidVolume();
    System.out.println("separator_volume_m3=" + vesselVolume);
    System.out.println("liquid_volume_m3=" + liquidVolume);
    System.out.println("gas_area_fraction="
        + (orientation.equals("horizontal") ? 1.0 - liquidVolume / vesselVolume : 1.0));
    System.out.println("gas_superficial_velocity_m_per_s="
        + separator.getGasSuperficialVelocity());
    System.out.println("gas_load_factor_m_per_s=" + separator.getGasLoadFactor());
    for (int i = 1; i < system.getNumberOfPhases(); i++) {
      System.out.println("gas_load_factor_phase_" + i + "_m_per_s="
          + separator.getGasLoadFactor(i));
    }
    System.out.println("derated_gas_load_factor_m_per_s="
        + separator.getDeRatedGasLoadFactor());
    for (int i = 1; i < system.getNumberOfPhases(); i++) {
      System.out.println("derated_gas_load_factor_phase_" + i + "_m_per_s="
          + separator.getDeRatedGasLoadFactor(i));
    }
    System.out.println("reported_design_gas_load_factor_m_per_s="
        + separator.getDesignGasLoadFactor());
    // `sizeFromFlow` mutates the vessel, so it is the row's last read.
    separator.sizeFromFlow(safetyFactor);
    System.out.println("sized_diameter_m=" + separator.getInternalDiameter());
    System.out.println();
  }

  // === process.separator mechanical design ====================================

  /// `SeparatorMechanicalDesign.calcDesign` on the class's own test fixture.
  ///
  /// **The fixture is `SeparatorMechanicalDesignTest.setUp`**: an SRK
  /// methane/ethane/propane/n-pentane fluid at 298.15 K and 20 bara, 10,000 kg/hr at
  /// 25 C, run through a `Separator` and then `calcDesign` with no company standard set.
  /// The row prints the design-standard keys the class finds, the response record's own
  /// fields and the getters the response does not carry, and any exception with its stack
  /// - because whether a default standard exists is the first question this slice has to
  /// answer, and a probe that swallowed it would answer nothing.
  static void separatorMechanicalRows() {
    // **Two rows, because the class's test is SRK and the port's separator is PR.** The SRK
    // row is the class's own fixture; the PR row is the fluid `Stream::mixture()` resolves,
    // so it is the oracle the port can be held to. The arithmetic is cubic-agnostic; only
    // the densities and the viscosity move.
    separatorMechanicalRow("class_test_fixture_srk", true);
    separatorMechanicalRow("port_fixture_pr", false);
  }

  /// `Separator.autoSize(safetyFactor)` on the PR fixture, at the class's no-argument
  /// `1.2` and at an overridden `1.5`.
  static void separatorAutosizeRows() {
    separatorAutosizeRow("pr_fixture_safety_1_2", 1.2);
    separatorAutosizeRow("pr_fixture_safety_1_5", 1.5);
  }

  static void separatorAutosizeRow(String label, double safetyFactor) {
    SystemInterface fluid = new SystemPrEos(298.15, 20.0);
    fluid.addComponent("methane", 0.7);
    fluid.addComponent("ethane", 0.1);
    fluid.addComponent("propane", 0.1);
    fluid.addComponent("n-pentane", 0.1);
    fluid.setMixingRule(2);
    Stream feed = new Stream("Feed", fluid);
    feed.setFlowRate(10000.0, "kg/hr");
    feed.setTemperature(25.0, "C");
    feed.setPressure(20.0, "bara");
    feed.run();
    neqsim.process.equipment.separator.Separator separator =
        new neqsim.process.equipment.separator.Separator("TestSeparator", feed);
    separator.run();
    separator.autoSize(safetyFactor);
    System.out.println(label);
    System.out.println("safety_factor=" + safetyFactor);
    System.out.println("auto_sized=" + separator.isAutoSized());
    System.out.println("design_inner_diameter_m=" + separator.getInternalDiameter());
    System.out.println("design_separator_length_m=" + separator.getSeparatorLength());
    System.out.println("design_wall_thickness_m="
        + separator.getMechanicalDesign().getWallThickness());
    System.out.println("design_outer_diameter_m="
        + separator.getMechanicalDesign().getOuterDiameter());
    System.out.println("design_vessel_volume_m3="
        + separator.getMechanicalDesign().getVolumeTotal());
    System.out.println("liquid_retention_time_s="
        + separator.getMechanicalDesign().calcLiquidRetentionTime());
    System.out.println("design_gas_load_factor_m_per_s="
        + separator.getMechanicalDesign().getGasLoadFactor());
    System.out.println("fg=" + separator.getMechanicalDesign().getFg());
    System.out.println("volume_safety_factor="
        + separator.getMechanicalDesign().getVolumeSafetyFactor());
    System.out.println("max_design_volume_flow_m3_hr="
        + separator.getMechanicalDesign().getMaxDesignVolumeFlow());
    System.out.println("max_design_gas_volume_flow_m3_hr="
        + separator.getMechanicalDesign().getMaxDesignGassVolumeFlow());
    System.out.println();
  }

  @SuppressWarnings({"checkstyle:MethodLength"})
  static void separatorMechanicalRow(String label, boolean srk) {
    SystemInterface fluid = srk ? new SystemSrkEos(298.15, 20.0) : new SystemPrEos(298.15, 20.0);
    fluid.addComponent("methane", 0.7);
    fluid.addComponent("ethane", 0.1);
    fluid.addComponent("propane", 0.1);
    fluid.addComponent("n-pentane", 0.1);
    if (srk) {
      fluid.setMixingRule("classic");
    } else {
      fluid.setMixingRule(2);
    }
    Stream feed = new Stream("Feed", fluid);
    feed.setFlowRate(10000.0, "kg/hr");
    feed.setTemperature(25.0, "C");
    feed.setPressure(20.0, "bara");
    feed.run();
    neqsim.process.equipment.separator.Separator separator =
        new neqsim.process.equipment.separator.Separator("TestSeparator", feed);
    feed.setTemperature(25.0, "C");
    separator.run();
    neqsim.process.mechanicaldesign.separator.SeparatorMechanicalDesign design =
        (neqsim.process.mechanicaldesign.separator.SeparatorMechanicalDesign) separator
            .getMechanicalDesign();

    System.out.println(label);
    System.out.println("eos=" + (srk ? "srk" : "pr"));
    System.out.println("standard_keys=" + design.getDesignStandard().keySet());
    for (String key : design.getDesignStandard().keySet()) {
      Object standard = design.getDesignStandard().get(key);
      System.out.println("standard_" + key.replace(' ', '_') + "="
          + (standard == null ? "null" : standard.getClass().getSimpleName()));
    }
    try {
      design.calcDesign();
      System.out.println("calc_design=ok");
    } catch (Exception ex) {
      System.out.println("calc_design_error=" + ex);
      ex.printStackTrace(System.out);
    }
    // **`super.calcDesign()` installs the default company standards when none were set**,
    // so the keys are read again here rather than before the call.
    // **The key names only, never the set's `toString`.** A `DesignStandard` value's string
    // form reaches `MechanicalDesign.hashCode`, which recurses through the cost estimator's
    // own `hashCode` and overflows the stack - the class's defect, and not this probe's.
    System.out.println("standard_keys_after=" + design.getDesignStandard().size());
    for (String key : design.getDesignStandard().keySet()) {
      Object standard = design.getDesignStandard().get(key);
      System.out.println("standard_after_" + key.replace(' ', '_') + "="
          + (standard == null ? "null" : standard.getClass().getSimpleName()));
    }
    // The default standards' own numbers, which is what a port has to declare as inputs.
    neqsim.process.mechanicaldesign.designstandards.SeparatorDesignStandard separatorStandard =
        (neqsim.process.mechanicaldesign.designstandards.SeparatorDesignStandard) design
            .getDesignStandard().get("separator process design");
    neqsim.process.mechanicaldesign.designstandards.MaterialPlateDesignStandard materialStandard =
        (neqsim.process.mechanicaldesign.designstandards.MaterialPlateDesignStandard) design
            .getDesignStandard().get("material plate design codes");
    neqsim.process.mechanicaldesign.designstandards.JointEfficiencyPlateStandard jointStandard =
        (neqsim.process.mechanicaldesign.designstandards.JointEfficiencyPlateStandard) design
            .getDesignStandard().get("plate Joint Efficiency design codes");
    neqsim.process.mechanicaldesign.designstandards.PressureVesselDesignStandard vesselStandard =
        (neqsim.process.mechanicaldesign.designstandards.PressureVesselDesignStandard) design
            .getDesignStandard().get("pressure vessel design code");
    System.out.println("separator_standard_name=" + separatorStandard.getClass().getSimpleName());
    System.out.println("separator_standard_gas_load_factor=" + separatorStandard.getGasLoadFactor());
    System.out.println("separator_standard_fg=" + separatorStandard.getFg());
    System.out.println("separator_standard_volume_safety_factor="
        + separatorStandard.getVolumetricDesignFactor());
    System.out.println("separator_standard_retention_s="
        + separatorStandard.getLiquidRetentionTime("API12J", design));
    System.out.println("material_standard_division_class=" + materialStandard.getDivisionClass());
    System.out.println("joint_efficiency=" + jointStandard.getJEFactor());
    System.out.println("vessel_standard_name=" + vesselStandard.getClass().getSimpleName());
    System.out.println("max_operation_pressure_bara=" + design.getMaxOperationPressure());
    System.out.println("max_operation_temperature_k=" + design.getMaxOperationTemperature());
    System.out.println("corrosion_allowance=" + design.getCorrosionAllowance());
    System.out.println("orientation=" + separator.getOrientation());
    // The sizing intermediates: the two densities, the two raw volumes the class divides by
    // `1e5`, and the two volume flows it forms from them.
    System.out.println("gas_density_kg_m3="
        + separator.getThermoSystem().getPhase(0).getPhysicalProperties().getDensity());
    System.out.println("liquid_density_kg_m3="
        + separator.getThermoSystem().getPhase(1).getPhysicalProperties().getDensity());
    System.out.println("liquid_viscosity_pa_s="
        + separator.getThermoSystem().getPhase(1).getPhysicalProperties().getViscosity());
    System.out.println("phase0_volume_raw="
        + separator.getThermoSystem().getPhase(0).getVolume());
    System.out.println("liquid_volume_raw=" + separator.getThermoSystem().getLiquidVolume());
    System.out.println("max_design_volume_flow_m3_hr=" + design.getMaxDesignVolumeFlow());
    System.out.println("max_design_gas_volume_flow_m3_hr=" + design.getMaxDesignGassVolumeFlow());
    System.out.println("inner_diameter_m=" + design.getInnerDiameter());
    System.out.println("outer_diameter_m=" + design.getOuterDiameter());
    System.out.println("vessel_shell_weight_kg=" + design.getWeigthVesselShell());
    System.out.println("internals_weight_kg=" + design.getWeigthInternals());
    System.out.println("piping_weight_kg=" + design.getWeightPiping());
    System.out.println("structural_weight_kg=" + design.getWeightStructualSteel());
    System.out.println("electrical_weight_kg=" + design.getWeightElectroInstrument());
    System.out.println("total_weight_kg=" + design.getWeightTotal());
    System.out.println("module_width_m=" + design.getModuleWidth());
    System.out.println("module_length_m=" + design.getModuleLength());
    System.out.println("module_height_m=" + design.getModuleHeight());
    System.out.println("vessel_volume_m3=" + design.getVolumeTotal());
    System.out.println("tantan_length_m=" + design.getTantanLength());
    System.out.println("wall_thickness_m=" + design.getWallThickness());
    System.out.println("gas_load_factor_m_per_s=" + design.getGasLoadFactor());
    System.out.println("fg=" + design.getFg());
    System.out.println("volume_safety_factor=" + design.getVolumeSafetyFactor());
    System.out.println("retention_time_s=" + design.getRetentionTime());
    System.out.println("hhll_m=" + design.getHHLL());
    System.out.println("hll_m=" + design.getHLL());
    System.out.println("nll_m=" + design.getNLL());
    System.out.println("lll_m=" + design.getLLL());
    System.out.println("llll_m=" + design.getLLLL());
    System.out.println("weir_height_m=" + design.getWeirHeight());
    System.out.println("hil_m=" + design.getHIL());
    System.out.println("nil_m=" + design.getNIL());
    System.out.println("lil_m=" + design.getLIL());
    System.out.println("effective_length_liquid_m=" + design.getEffectiveLengthLiquid());
    System.out.println("effective_length_gas_m=" + design.getEffectiveLengthGas());
    System.out.println("calc_liquid_retention_time_s=" + design.calcLiquidRetentionTime());
    System.out.println("calc_surge_volume_m3=" + design.calcSurgeVolume());
    System.out.println("calc_surge_time_s=" + design.calcSurgeTime());
    System.out.println("calc_inlet_nozzle_id_m=" + design.calcInletNozzleID());
    System.out.println("calc_gas_outlet_nozzle_id_m=" + design.calcGasOutletNozzleID());
    System.out.println("calc_oil_outlet_nozzle_id_m=" + design.calcOilOutletNozzleID());
    neqsim.process.mechanicaldesign.separator.SeparatorMechanicalDesignResponse response =
        design.getResponse();
    if (response != null) {
      System.out.println("response_orientation=" + response.getOrientation());
      System.out.println("response_separator_type=" + response.getSeparatorType());
      System.out.println("response_gas_load_factor=" + response.getGasLoadFactor());
      System.out.println("response_volume_safety_factor=" + response.getVolumeSafetyFactor());
      System.out.println("response_liquid_level_fraction=" + response.getLiquidLevelFraction());
      System.out.println("response_retention_time=" + response.getRetentionTime());
      System.out.println("response_demister_type=" + response.getDemisterType());
      System.out.println("response_demister_efficiency=" + response.getDemisterEfficiency());
      System.out.println("response_inlet_nozzles=" + response.getNumberOfInletNozzles());
      System.out.println("response_inlet_nozzle_diameter=" + response.getInletNozzleDiameter());
      System.out.println("response_gas_outlet_nozzle_diameter="
          + response.getGasOutletNozzleDiameter());
      System.out.println("response_liquid_outlet_nozzle_diameter="
          + response.getLiquidOutletNozzleDiameter());
      System.out.println("response_water_outlet_nozzle_diameter="
          + response.getWaterOutletNozzleDiameter());
      System.out.println("response_head_type=" + response.getHeadType());
      System.out.println("response_head_thickness=" + response.getHeadThickness());
      System.out.println("response_shell_thickness=" + response.getShellThickness());
      System.out.println("response_design_code=" + response.getDesignCode());
      System.out.println("response_gas_design_velocity=" + response.getGasDesignVelocity());
      System.out.println("response_actual_gas_velocity=" + response.getActualGasVelocity());
      System.out.println("response_allowable_gas_velocity=" + response.getAllowableGasVelocity());
      System.out.println("response_design_gas_flow=" + response.getDesignGasFlow());
      System.out.println("response_design_liquid_flow=" + response.getDesignLiquidFlow());
      System.out.println("response_design_water_flow=" + response.getDesignWaterFlow());
      System.out.println("response_normal_liquid_level=" + response.getNormalLiquidLevel());
      System.out.println("response_high_liquid_level=" + response.getHighLiquidLevel());
      System.out.println("response_low_liquid_level=" + response.getLowLiquidLevel());
      System.out.println("response_interface_level=" + response.getInterfaceLevel());
      System.out.println("response_surge_volume=" + response.getSurgeVolume());
    }
    System.out.println();
  }

  // === process.gas_scrubber mechanical design =================================

  /// `GasScrubberMechanicalDesign.calcDesign` on the same two-cubic fixture.
  ///
  /// **The subclass is a second sizing machine.** It calls the separator's `calcDesign` first
  /// - which installs the default standards and sizes a horizontal-vessel diameter - and then
  /// re-sizes with the scrubber's own vertical `Fg = 1.0` (or `1 - design level fraction`)
  /// and `length = 5 d`, replacing the diameter while the levels the base pass set stay on
  /// the base pass's diameter. The row prints both bases' K and safety factor and the levels
  /// beside the final diameter, so that staleness is visible rather than inferred.
  static void gasScrubberMechanicalRows() {
    gasScrubberMechanicalRow("class_fixture_srk", true, false);
    gasScrubberMechanicalRow("port_fixture_pr", false, false);
    // The mesh pad's configuration is design state, set after `calcDesign` sizes the vessel.
    gasScrubberMechanicalRow("port_fixture_pr_mesh_pad", false, true);
  }

  @SuppressWarnings({"checkstyle:MethodLength"})
  static void gasScrubberMechanicalRow(String label, boolean srk, boolean meshPad) {
    SystemInterface fluid = srk ? new SystemSrkEos(298.15, 20.0) : new SystemPrEos(298.15, 20.0);
    fluid.addComponent("methane", 0.7);
    fluid.addComponent("ethane", 0.1);
    fluid.addComponent("propane", 0.1);
    fluid.addComponent("n-pentane", 0.1);
    if (srk) {
      fluid.setMixingRule("classic");
    } else {
      fluid.setMixingRule(2);
    }
    Stream feed = new Stream("Feed", fluid);
    feed.setFlowRate(10000.0, "kg/hr");
    feed.setTemperature(25.0, "C");
    feed.setPressure(20.0, "bara");
    feed.run();
    neqsim.process.equipment.separator.GasScrubber scrubber =
        new neqsim.process.equipment.separator.GasScrubber("TestScrubber", feed);
    feed.setTemperature(25.0, "C");
    scrubber.run();
    neqsim.process.mechanicaldesign.separator.GasScrubberMechanicalDesign design =
        (neqsim.process.mechanicaldesign.separator.GasScrubberMechanicalDesign) scrubber
            .getMechanicalDesign();

    System.out.println(label);
    System.out.println("eos=" + (srk ? "srk" : "pr"));
    try {
      design.calcDesign();
      System.out.println("calc_design=ok");
    } catch (Exception ex) {
      System.out.println("calc_design_error=" + ex);
      ex.printStackTrace(System.out);
    }
    if (meshPad) {
      double area = Math.PI * Math.pow(design.getInnerDiameter() / 2.0, 2);
      design.setMeshPad(area, 150.0);
      design.setMeshPadElevationM(0.75);
    }
    System.out.println("standard_keys_after=" + design.getDesignStandard().size());
    for (String key : design.getDesignStandard().keySet()) {
      Object standard = design.getDesignStandard().get(key);
      System.out.println("standard_after_" + key.replace(' ', '_') + "="
          + (standard == null ? "null" : standard.getClass().getSimpleName()));
    }
    neqsim.process.mechanicaldesign.designstandards.GasScrubberDesignStandard scrubberStandard =
        (neqsim.process.mechanicaldesign.designstandards.GasScrubberDesignStandard) design
            .getDesignStandard().get("gas scrubber process design");
    neqsim.process.mechanicaldesign.designstandards.SeparatorDesignStandard separatorStandard =
        (neqsim.process.mechanicaldesign.designstandards.SeparatorDesignStandard) design
            .getDesignStandard().get("separator process design");
    System.out.println("scrubber_standard_gas_load_factor="
        + scrubberStandard.getGasLoadFactor());
    System.out.println("scrubber_standard_volume_safety_factor="
        + scrubberStandard.getVolumetricDesignFactor());
    System.out.println("separator_standard_gas_load_factor="
        + separatorStandard.getGasLoadFactor());
    System.out.println("separator_standard_volume_safety_factor="
        + separatorStandard.getVolumetricDesignFactor());
    System.out.println("orientation=" + scrubber.getOrientation());
    System.out.println("design_liquid_level_fraction="
        + scrubber.getDesignLiquidLevelFraction());
    System.out.println("gas_density_kg_m3="
        + scrubber.getThermoSystem().getPhase(0).getPhysicalProperties().getDensity());
    System.out.println("liquid_density_kg_m3="
        + scrubber.getThermoSystem().getPhase(1).getPhysicalProperties().getDensity());
    System.out.println("phase0_volume_raw=" + scrubber.getThermoSystem().getPhase(0).getVolume());
    System.out.println("liquid_volume_raw=" + scrubber.getThermoSystem().getLiquidVolume());
    System.out.println("max_design_volume_flow_m3_hr=" + design.getMaxDesignVolumeFlow());
    System.out.println("max_design_gas_volume_flow_m3_hr=" + design.getMaxDesignGassVolumeFlow());
    System.out.println("inner_diameter_m=" + design.getInnerDiameter());
    System.out.println("tantan_length_m=" + design.getTantanLength());
    System.out.println("outer_diameter_m=" + design.getOuterDiameter());
    System.out.println("wall_thickness_m=" + design.getWallThickness());
    System.out.println("vessel_shell_weight_kg=" + design.getWeigthVesselShell());
    System.out.println("internals_weight_kg=" + design.getWeigthInternals());
    System.out.println("piping_weight_kg=" + design.getWeightPiping());
    System.out.println("structural_weight_kg=" + design.getWeightStructualSteel());
    System.out.println("electrical_weight_kg=" + design.getWeightElectroInstrument());
    System.out.println("total_weight_kg=" + design.getWeightTotal());
    System.out.println("module_width_m=" + design.getModuleWidth());
    System.out.println("module_length_m=" + design.getModuleLength());
    System.out.println("module_height_m=" + design.getModuleHeight());
    System.out.println("vessel_volume_m3=" + design.getVolumeTotal());
    System.out.println("gas_load_factor_m_per_s=" + design.getGasLoadFactor());
    System.out.println("fg=" + design.getFg());
    System.out.println("volume_safety_factor=" + design.getVolumeSafetyFactor());
    System.out.println("retention_time_s=" + design.getRetentionTime());
    System.out.println("hhll_m=" + design.getHHLL());
    System.out.println("hll_m=" + design.getHLL());
    System.out.println("nll_m=" + design.getNLL());
    System.out.println("lll_m=" + design.getLLL());
    System.out.println("llll_m=" + design.getLLLL());
    System.out.println("weir_height_m=" + design.getWeirHeight());
    System.out.println("hil_m=" + design.getHIL());
    System.out.println("nil_m=" + design.getNIL());
    System.out.println("lil_m=" + design.getLIL());
    System.out.println("effective_length_liquid_m=" + design.getEffectiveLengthLiquid());
    System.out.println("effective_length_gas_m=" + design.getEffectiveLengthGas());
    System.out.println("calc_liquid_retention_time_s=" + design.calcLiquidRetentionTime());
    System.out.println("calc_surge_volume_m3=" + design.calcSurgeVolume());
    System.out.println("calc_surge_time_s=" + design.calcSurgeTime());
    System.out.println("calc_inlet_nozzle_id_m=" + design.calcInletNozzleID());
    System.out.println("calc_gas_outlet_nozzle_id_m=" + design.calcGasOutletNozzleID());
    System.out.println("calc_oil_outlet_nozzle_id_m=" + design.calcOilOutletNozzleID());
    System.out.println("has_mesh_pad=" + design.hasMeshPad());
    System.out.println("mesh_pad_area_m2=" + design.getMeshPadAreaM2());
    System.out.println("mesh_pad_thickness_mm=" + design.getMeshPadThicknessMm());
    System.out.println("mesh_pad_elevation_m=" + design.getMeshPadElevationM());
    System.out.println("has_inlet_cyclones=" + design.hasInletCyclones());
    System.out.println("has_demisting_cyclones=" + design.hasDemistingCyclones());
    if (meshPad) {
      // `sizeFromFlow` is inherited from `Separator`; it mutates the vessel, so it is last.
      scrubber.sizeFromFlow(1.5);
      System.out.println("sized_diameter_m=" + scrubber.getInternalDiameter());
    }
    System.out.println();
  }

  // === process.three_phase_separator interface settling ======================

  /// `ThreePhaseSeparator.calcInterfaceSettlingTime` on the two vessel states.
  ///
  /// The class settles a **150 µm water droplet** through the oil layer, Stokes on the oil's
  /// own viscosity, and answers **minutes**; `isWithinAllLimits` requires at least `1.0`.
  /// The levels are the mechanical design's, with `0.20`/`0.50` of the diameter and
  /// `0.82 * length` as fallbacks, which is what a vessel that never ran `calcDesign` reads.
  /// The valve fractions are set on the second state and the mass balance printed: they are
  /// read by that report alone, which is why they are not model inputs.
  static void threePhaseSettlingRows() {
    String[] names = new String[] {"methane", "n-butane", "water"};
    double[] z = new double[] {0.5, 0.3, 0.2};
    threePhaseSettlingRow("vessel_d1000_l3000", names, z, 1.0, 3.0, false);
    threePhaseSettlingRow("vessel_d1200_l4000_valves", names, z, 1.2, 4.0, true);
  }

  static void threePhaseSettlingRow(String label, String[] names, double[] z, double diameterM,
      double lengthM, boolean valves) {
    Stream inlet = feed(names, z, 300.0, 20.0, 1.0);
    neqsim.process.equipment.separator.ThreePhaseSeparator separator =
        new neqsim.process.equipment.separator.ThreePhaseSeparator("sep3", inlet);
    separator.setInternalDiameter(diameterM);
    separator.setSeparatorLength(lengthM);
    if (valves) {
      separator.setGasOutletFlowFraction(0.10);
      separator.setOilOutletFlowFraction(0.20);
      separator.setWaterOutletFlowFraction(0.30);
    }
    separator.run();
    neqsim.process.mechanicaldesign.separator.SeparatorMechanicalDesign design =
        (neqsim.process.mechanicaldesign.separator.SeparatorMechanicalDesign) separator
            .getMechanicalDesign();
    System.out.println(label);
    System.out.println("phase_count=" + separator.getThermoSystem().getNumberOfPhases());
    System.out.println("has_oil=" + separator.getThermoSystem().hasPhaseType("oil"));
    System.out.println("has_aqueous=" + separator.getThermoSystem().hasPhaseType("aqueous"));
    System.out.println("internal_diameter_m=" + separator.getInternalDiameter());
    System.out.println("separator_length_m=" + separator.getSeparatorLength());
    System.out.println("nil_m=" + design.getNIL());
    System.out.println("nll_m=" + design.getNLL());
    System.out.println("effective_length_liquid_m=" + design.getEffectiveLengthLiquid());
    System.out.println("oil_density_kg_m3="
        + separator.getThermoSystem().getPhase("oil").getDensity("kg/m3"));
    System.out.println("water_density_kg_m3="
        + separator.getThermoSystem().getPhase("aqueous").getDensity("kg/m3"));
    System.out.println("oil_viscosity_pa_s="
        + separator.getThermoSystem().getPhase("oil").getPhysicalProperties().getViscosity());
    // **The class's own call reads the aqueous density before the phase has initialised it.**
    // On a freshly run separator that is `0.0`, the settling velocity is negative and the
    // answer is `Infinity`. The same call after `initPhysicalProperties()` is the state the
    // method means; both are printed because the difference is the finding.
    System.out.println("interface_settling_time_uninitialised_min="
        + separator.calcInterfaceSettlingTime());
    System.out.println("within_all_limits_uninitialised=" + separator.isWithinAllLimits());
    separator.getThermoSystem().initPhysicalProperties();
    System.out.println("water_density_initialised_kg_m3="
        + separator.getThermoSystem().getPhase("aqueous").getDensity("kg/m3"));
    System.out.println("oil_viscosity_initialised_pa_s="
        + separator.getThermoSystem().getPhase("oil").getPhysicalProperties().getViscosity());
    System.out.println("interface_settling_time_initialised_min="
        + separator.calcInterfaceSettlingTime());
    System.out.println("within_all_limits_initialised=" + separator.isWithinAllLimits());
    if (valves) {
      System.out.println("gas_outlet_fraction=" + separator.getGasOutletFlowFraction());
      System.out.println("oil_outlet_fraction=" + separator.getOilOutletFlowFraction());
      System.out.println("water_outlet_fraction=" + separator.getWaterOutletFlowFraction());
      System.out.println("mass_balance_kg_hr=" + separator.getMassBalance("kg/hr"));
    }
    System.out.println();
  }

  // === process.separator_performance, the enhanced chain =====================

  /// `SeparatorPerformanceCalculator.calculateEnhanced` on the class's own test states.
  ///
  /// The chain is flow regime -> inlet device -> geometry -> gravity -> mist eliminator ->
  /// liquid-liquid. **Every row prints each stage's intermediates**, because the enhanced
  /// path is a second solver and a capture of the final efficiency cannot say which stage
  /// moved. The constants are `EnhancedPerformanceCalculatorTest`'s own.
  static void separatorPerformanceEnhancedRows() {
    double gasDensity = 50.0;
    double oilDensity = 750.0;
    double waterDensity = 1025.0;
    double gasViscosity = 1.2e-5;
    double oilViscosity = 2.0e-3;
    double waterViscosity = 8.0e-4;
    double velocity = 2.0;
    double diameter = 2.4;
    double length = 7.0;
    double level = 0.5;

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator twoPhase =
        enhancedCalculator();
    twoPhase.setInletPipeDiameter(0.2);
    twoPhase.setSurfaceTension(0.025);
    twoPhase.setGasBubbleDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            500e-6, 2.0));
    runEnhanced("two_phase_horizontal", twoPhase, gasDensity, oilDensity, 0.0, gasViscosity,
        oilViscosity, 0.0, velocity, diameter, length, "horizontal", level);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator threePhase =
        enhancedCalculator();
    threePhase.setInletPipeDiameter(0.25);
    threePhase.setSurfaceTension(0.025);
    threePhase.setOilWaterInterfacialTension(0.030);
    threePhase.setGasBubbleDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            500e-6, 2.0));
    threePhase.setWaterInOilDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.logNormal(
            200e-6, 0.6));
    threePhase.setOilInWaterDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.logNormal(
            150e-6, 0.5));
    runEnhanced("three_phase_horizontal", threePhase, gasDensity, oilDensity, waterDensity,
        gasViscosity, oilViscosity, waterViscosity, velocity, diameter, length, "horizontal",
        level);

    // The two geometry positions the class exposes as setters. The inlet device's is a
    // private constant with no setter, so the gravity-section length is the read-back.
    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator mistShift =
        enhancedCalculator();
    neqsim.process.equipment.separator.entrainment.SeparatorGeometryCalculator mistGeometry =
        new neqsim.process.equipment.separator.entrainment.SeparatorGeometryCalculator();
    mistGeometry.setMistEliminatorPosition(0.80);
    mistShift.setGeometryCalculator(mistGeometry);
    runEnhanced("mist_position_0_80", mistShift, gasDensity, oilDensity, 0.0, gasViscosity,
        oilViscosity, 0.0, velocity, diameter, length, "horizontal", level);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator weirShift =
        enhancedCalculator();
    neqsim.process.equipment.separator.entrainment.SeparatorGeometryCalculator weirGeometry =
        new neqsim.process.equipment.separator.entrainment.SeparatorGeometryCalculator();
    weirGeometry.setWeirPosition(0.50);
    weirShift.setGeometryCalculator(weirGeometry);
    runEnhanced("weir_position_0_50", weirShift, gasDensity, oilDensity, waterDensity,
        gasViscosity, oilViscosity, waterViscosity, velocity, diameter, length, "horizontal",
        level);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator vertical =
        enhancedCalculator();
    vertical.setInletPipeDiameter(0.15);
    vertical.setSurfaceTension(0.025);
    runEnhanced("two_phase_vertical", vertical, gasDensity, oilDensity, 0.0, gasViscosity,
        oilViscosity, 0.0, velocity, 1.5, 5.0, "vertical", 0.4);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator none =
        enhancedCalculator();
    none.setInletDeviceModel(
        new neqsim.process.equipment.separator.entrainment.InletDeviceModel(
            neqsim.process.equipment.separator.entrainment.InletDeviceModel.InletDeviceType.NONE));
    runEnhanced("inlet_device_none", none, gasDensity, oilDensity, 0.0, gasViscosity, oilViscosity,
        0.0, velocity, diameter, length, "horizontal", level);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator cyclone =
        enhancedCalculator();
    cyclone.setInletDeviceModel(
        new neqsim.process.equipment.separator.entrainment.InletDeviceModel(
            neqsim.process.equipment.separator.entrainment.InletDeviceModel.InletDeviceType.INLET_CYCLONE));
    runEnhanced("inlet_device_cyclone", cyclone, gasDensity, oilDensity, 0.0, gasViscosity,
        oilViscosity, 0.0, velocity, diameter, length, "horizontal", level);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator flooding =
        enhancedCalculator();
    runEnhanced("high_velocity_flooding", flooding, gasDensity, oilDensity, 0.0, gasViscosity,
        oilViscosity, 0.0, 15.0, diameter, length, "horizontal", level);

    // The internals database's other two default records, whose design K the calculator reads
    // when none is stated.
    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator vane =
        enhancedCalculator();
    vane.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.vanePackDefault());
    vane.setInletPipeDiameter(0.2);
    vane.setSurfaceTension(0.025);
    vane.setGasBubbleDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            500e-6, 2.0));
    runEnhanced("vane_pack_design_k", vane, gasDensity, oilDensity, 0.0, gasViscosity, oilViscosity,
        0.0, velocity, diameter, length, "horizontal", level);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator overrides =
        enhancedCalculator();
    neqsim.process.equipment.separator.entrainment.InletDeviceModel overrideDevice =
        new neqsim.process.equipment.separator.entrainment.InletDeviceModel(
            neqsim.process.equipment.separator.entrainment.InletDeviceModel.InletDeviceType
                .HALF_PIPE);
    overrideDevice.setOverrideBulkEfficiency(0.5);
    overrideDevice.setOverrideDsdMultiplier(1.5);
    overrides.setInletDeviceModel(overrideDevice);
    overrides.setInletPipeDiameter(0.2);
    overrides.setSurfaceTension(0.025);
    runEnhanced("inlet_device_overrides", overrides, gasDensity, oilDensity, 0.0, gasViscosity,
        oilViscosity, 0.0, velocity, diameter, length, "horizontal", level);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator axial =
        enhancedCalculator();
    axial.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve
            .axialCycloneDefault());
    axial.setInletPipeDiameter(0.2);
    axial.setSurfaceTension(0.025);
    axial.setGasBubbleDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            500e-6, 2.0));
    runEnhanced("axial_cyclone_design_k", axial, gasDensity, oilDensity, 0.0, gasViscosity,
        oilViscosity, 0.0, velocity, diameter, length, "horizontal", level);
  }

  static neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator
      enhancedCalculator() {
    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator calc =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    calc.setUseEnhancedCalculation(true);
    calc.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.wireMeshDefault());
    return calc;
  }

  @SuppressWarnings({"checkstyle:MethodLength", "checkstyle:ParameterNumber"})
  static void runEnhanced(String label,
      neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator calc,
      double gasDensity, double oilDensity, double waterDensity, double gasViscosity,
      double oilViscosity, double waterViscosity, double gasVelocity, double diameter,
      double length, String orientation, double level) {
    calc.calculate(gasDensity, oilDensity, waterDensity, gasViscosity, oilViscosity,
        waterViscosity, gasVelocity, diameter, length, orientation, level);
    System.out.println(label);
    System.out.println("gas_density_kg_m3=" + gasDensity);
    System.out.println("oil_density_kg_m3=" + oilDensity);
    System.out.println("water_density_kg_m3=" + waterDensity);
    System.out.println("gas_velocity_m_per_s=" + gasVelocity);
    System.out.println("diameter_m=" + diameter);
    System.out.println("length_m=" + length);
    System.out.println("orientation=" + orientation);
    System.out.println("liquid_level_fraction=" + level);
    System.out.println("inlet_pipe_diameter_m=" + calc.getInletPipeDiameter());
    System.out.println("surface_tension_n_per_m=" + calc.getSurfaceTension());
    System.out.println("oil_water_interfacial_tension_n_per_m="
        + calc.getOilWaterInterfacialTension());
    System.out.println("overall_gas_liquid_efficiency=" + calc.getOverallGasLiquidEfficiency());
    System.out.println("gravity_section_efficiency=" + calc.getGravitySectionEfficiency());
    System.out.println("mist_eliminator_efficiency=" + calc.getMistEliminatorEfficiency());
    System.out.println("gravity_cut_diameter_m=" + calc.getGravityCutDiameter());
    System.out.println("liquid_liquid_gravity_efficiency="
        + calc.getLiquidLiquidGravityEfficiency());
    System.out.println("oil_in_gas_fraction=" + calc.getOilInGasFraction());
    System.out.println("water_in_gas_fraction=" + calc.getWaterInGasFraction());
    System.out.println("gas_in_oil_fraction=" + calc.getGasInOilFraction());
    System.out.println("gas_in_water_fraction=" + calc.getGasInWaterFraction());
    System.out.println("oil_in_water_fraction=" + calc.getOilInWaterFraction());
    System.out.println("water_in_oil_fraction=" + calc.getWaterInOilFraction());
    System.out.println("k_factor_m_per_s=" + calc.getKFactor());
    System.out.println("k_factor_utilization=" + calc.getKFactorUtilization());
    System.out.println("inlet_flow_regime="
        + (calc.getInletFlowRegime() == null ? "null" : calc.getInletFlowRegime().name()));
    System.out.println("inlet_device_bulk_efficiency=" + calc.getInletDeviceBulkEfficiency());
    neqsim.process.equipment.separator.entrainment.InletDeviceModel device =
        calc.getInletDeviceModel();
    if (device != null) {
      System.out.println("inlet_device_type=" + device.getDeviceType().name());
      System.out.println("inlet_device_bulk_separation_efficiency="
          + device.getBulkSeparationEfficiency());
      System.out.println("inlet_device_pressure_drop_pa=" + device.getPressureDrop());
      System.out.println("inlet_nozzle_velocity_m_per_s=" + device.getNozzleVelocity());
      System.out.println("inlet_momentum_flux_pa=" + device.getMomentumFlux());
      neqsim.process.equipment.separator.entrainment.DropletSizeDistribution downstreamDsd =
          device.getDownstreamDSD();
      if (downstreamDsd != null) {
        System.out.println("post_inlet_d50_m=" + downstreamDsd.getD50());
        System.out.println("post_inlet_d32_m=" + downstreamDsd.getSauterMeanDiameter());
      }
    }
    neqsim.process.equipment.separator.entrainment.SeparatorGeometryCalculator geometry =
        calc.getGeometryCalculator();
    if (geometry != null) {
      System.out.println("geometry_gas_area_m2=" + geometry.getGasArea());
      System.out.println("geometry_liquid_area_m2=" + geometry.getLiquidArea());
      System.out.println("geometry_gas_volume_m3=" + geometry.getGasVolume());
      System.out.println("geometry_liquid_volume_m3=" + geometry.getLiquidVolume());
      System.out.println("geometry_total_volume_m3="
          + (geometry.getGasVolume() + geometry.getLiquidVolume()));
      System.out.println("geometry_gas_settling_height_m="
          + geometry.getEffectiveGasSettlingHeight());
      System.out.println("geometry_liquid_settling_height_m="
          + geometry.getEffectiveLiquidSettlingHeight());
      System.out.println("geometry_gas_residence_time_s=" + geometry.getGasResidenceTime());
      System.out.println("geometry_liquid_residence_time_s=" + geometry.getLiquidResidenceTime());
      System.out.println("geometry_gravity_section_length_m="
          + geometry.getGravitySectionLength());
      System.out.println("geometry_oil_pad_thickness_m=" + geometry.getOilPadThickness());
      System.out.println("geometry_water_layer_height_m=" + geometry.getWaterLayerHeight());
    }
    neqsim.process.equipment.separator.entrainment.DropletSettlingCalculator.ApiComplianceResult
        api = calc.getApiComplianceResult();
    if (api != null) {
      System.out.println("api_gas_liquid_compliant=" + api.gasLiquidSectionCompliant);
      System.out.println("api_liquid_compliant=" + api.liquidSectionCompliant);
      System.out.println("api_fully_compliant=" + api.isFullyCompliant());
      System.out.println("api_gravity_cut_diameter_um=" + api.gravityCutDiameter_um);
      System.out.println("api_k_factor_utilization=" + api.kFactorUtilization);
      System.out.println("api_gas_liquid_comment=" + api.gasLiquidComment);
      System.out.println("api_liquid_comment=" + api.liquidComment);
    }
    System.out.println();
  }

  // === process.separator_performance ==========================================

  /// `SeparatorPerformanceCalculator`'s standard path.
  ///
  /// **The configurations are the class's own test states**, so the capture pins the
  /// arithmetic the test only bounds: the two-phase wire mesh, the vertical vane-pack scrubber,
  /// the three-phase state, gravity alone, the gas-bubble carry-under and the plate-pack
  /// coalescer with no gas-liquid distribution. Two rows follow the calibration test's base and
  /// tuned machines, one turns the gravity section off, and one uses a log-normal distribution
  /// with custom lookup tables on both internals and coalescer.
  static void separatorPerformanceRows() {
    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator twoPhase =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    twoPhase.setGasLiquidDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            100.0e-6, 2.6));
    twoPhase.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.wireMeshDefault());
    runSeparatorPerformance("two_phase_horizontal_wire_mesh", twoPhase, 50.0, 800.0, 0.0, 1.5e-5,
        2.0e-3, 0.0, 1.0, 2.0, 8.0, "horizontal", 0.5, null, null, null);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator vertical =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    vertical.setGasLiquidDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            30.0e-6, 2.6));
    vertical.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.vanePackDefault());
    runSeparatorPerformance("vertical_scrubber_vane_pack", vertical, 60.0, 850.0, 0.0, 1.5e-5,
        3.0e-3, 0.0, 2.0, 1.0, 4.0, "vertical", 0.1, null, null, null);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator threePhase =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    threePhase.setGasLiquidDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            80.0e-6, 2.6));
    threePhase.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.wireMeshDefault());
    threePhase.setWaterInOilDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            300.0e-6, 2.0));
    threePhase.setOilInWaterDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            200.0e-6, 2.0));
    runSeparatorPerformance("three_phase_horizontal", threePhase, 40.0, 800.0, 1020.0, 1.5e-5,
        5.0e-3, 8.0e-4, 1.5, 2.5, 10.0, "horizontal", 0.6, null, null, null);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator gravityOnly =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    gravityOnly.setGasLiquidDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            20.0e-6, 2.6));
    runSeparatorPerformance("gravity_only_no_mist_eliminator", gravityOnly, 30.0, 750.0, 0.0,
        1.2e-5, 2.0e-3, 0.0, 0.5, 3.0, 12.0, "horizontal", 0.5, null, null, null);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator bubble =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    bubble.setGasLiquidDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            100.0e-6, 2.6));
    bubble.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.wireMeshDefault());
    bubble.setGasBubbleDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            500.0e-6, 2.0));
    runSeparatorPerformance("gas_bubble_carry_under", bubble, 50.0, 800.0, 0.0, 1.5e-5, 3.0e-3,
        0.0, 1.0, 2.0, 8.0, "horizontal", 0.5, null, null, null);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator plate =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    plate.setWaterInOilDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            100.0e-6, 2.0));
    plate.setOilInWaterDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            80.0e-6, 2.0));
    plate.setOilWaterCoalescerCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.platePack(30.0e-6,
            0.98));
    runSeparatorPerformance("liquid_liquid_plate_pack_no_gas_liquid", plate, 40.0, 800.0, 1020.0,
        1.5e-5, 5.0e-3, 8.0e-4, 1.5, 2.5, 10.0, "horizontal", 0.6, null, null, null);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator calibrationBase =
        calibrationMachine();
    runSeparatorPerformance("calibration_base", calibrationBase, 50.0, 800.0, 1020.0, 1.5e-5,
        3.0e-3, 8.0e-4, 1.2, 2.0, 8.0, "horizontal", 0.5, null, null, null);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator calibrationTuned =
        calibrationMachine();
    runSeparatorPerformance("calibration_tuned", calibrationTuned, 50.0, 800.0, 1020.0, 1.5e-5,
        3.0e-3, 8.0e-4, 1.2, 2.0, 8.0, "horizontal", 0.5, 1.5, 1.25, 0.8);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator custom =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    custom.setGasLiquidDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.logNormal(50.0e-6,
            0.8));
    custom.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.custom(
            new double[] {1.0e-6, 10.0e-6, 100.0e-6}, new double[] {0.1, 0.5, 0.99}));
    custom.setGasBubbleDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.logNormal(400.0e-6,
            0.7));
    custom.setWaterInOilDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.logNormal(150.0e-6,
            0.7));
    custom.setOilInWaterDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.logNormal(100.0e-6,
            0.9));
    custom.setOilWaterCoalescerCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.custom(
            new double[] {5.0e-6, 50.0e-6, 500.0e-6}, new double[] {0.05, 0.4, 0.95}));
    runSeparatorPerformance("log_normal_custom_curves", custom, 40.0, 820.0, 1010.0, 1.4e-5,
        2.5e-3, 9.0e-4, 0.8, 2.2, 9.0, "horizontal", 0.55, null, null, null);

    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator noGravity =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    noGravity.setGasLiquidDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            100.0e-6, 2.6));
    noGravity.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.wireMeshDefault());
    noGravity.setIncludeGravitySection(false);
    runSeparatorPerformance("gravity_section_off", noGravity, 50.0, 800.0, 0.0, 1.5e-5, 2.0e-3,
        0.0, 1.0, 2.0, 8.0, "horizontal", 0.5, null, null, null);
  }

  /// The calibration test's machine, before either factor is set.
  static neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator
      calibrationMachine() {
    neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator calc =
        new neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator();
    calc.setGasLiquidDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            120.0e-6, 2.2));
    calc.setMistEliminatorCurve(
        neqsim.process.equipment.separator.entrainment.GradeEfficiencyCurve.vanePackDefault());
    calc.setGasBubbleDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            400.0e-6, 2.0));
    calc.setWaterInOilDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            120.0e-6, 2.0));
    calc.setOilInWaterDSD(
        neqsim.process.equipment.separator.entrainment.DropletSizeDistribution.rosinRammler(
            120.0e-6, 2.0));
    return calc;
  }

  static void runSeparatorPerformance(String label,
      neqsim.process.equipment.separator.entrainment.SeparatorPerformanceCalculator calc,
      double gasDensity, double oilDensity, double waterDensity, double gasViscosity,
      double oilViscosity, double waterViscosity, double gasVelocity, double diameter,
      double length, String orientation, double level, Double liquidInGasFactor,
      Double gasCarryUnderFactor, Double liquidLiquidFactor) {
    if (liquidInGasFactor != null) {
      calc.setLiquidInGasCalibrationFactor(liquidInGasFactor);
    }
    if (gasCarryUnderFactor != null) {
      calc.setGasCarryUnderCalibrationFactor(gasCarryUnderFactor);
    }
    if (liquidLiquidFactor != null) {
      calc.setLiquidLiquidCalibrationFactor(liquidLiquidFactor);
    }
    calc.calculate(gasDensity, oilDensity, waterDensity, gasViscosity, oilViscosity,
        waterViscosity, gasVelocity, diameter, length, orientation, level);
    System.out.println(label);
    System.out.println("gas_density_kg_m3=" + gasDensity);
    System.out.println("oil_density_kg_m3=" + oilDensity);
    System.out.println("water_density_kg_m3=" + waterDensity);
    System.out.println("gas_viscosity_Pa_s=" + gasViscosity);
    System.out.println("oil_viscosity_Pa_s=" + oilViscosity);
    System.out.println("water_viscosity_Pa_s=" + waterViscosity);
    System.out.println("gas_velocity_m_per_s=" + gasVelocity);
    System.out.println("vessel_diameter_m=" + diameter);
    System.out.println("vessel_length_m=" + length);
    System.out.println("orientation=" + orientation);
    System.out.println("liquid_level_fraction=" + level);
    System.out.println("oil_in_gas_fraction=" + calc.getOilInGasFraction());
    System.out.println("water_in_gas_fraction=" + calc.getWaterInGasFraction());
    System.out.println("gas_in_oil_fraction=" + calc.getGasInOilFraction());
    System.out.println("gas_in_water_fraction=" + calc.getGasInWaterFraction());
    System.out.println("oil_in_water_fraction=" + calc.getOilInWaterFraction());
    System.out.println("water_in_oil_fraction=" + calc.getWaterInOilFraction());
    System.out.println("gravity_section_efficiency=" + calc.getGravitySectionEfficiency());
    System.out.println("mist_eliminator_efficiency=" + calc.getMistEliminatorEfficiency());
    System.out.println("overall_gas_liquid_efficiency=" + calc.getOverallGasLiquidEfficiency());
    System.out.println(
        "liquid_liquid_gravity_efficiency=" + calc.getLiquidLiquidGravityEfficiency());
    System.out.println("gravity_cut_diameter_m=" + calc.getGravityCutDiameter());
    System.out.println();
  }

  // === process.separator_capacity_assessment ==================================

  /// `Separator.getCapacityAssessment`, the detached fail-closed snapshot.
  ///
  /// **The mechanical limit is stated as a fraction of the Souders-Brown flow**, because it
  /// is the one number the pin's test derives rather than states: its `0.8` limit is
  /// `0.8 * getSoudersBrownMaximumGasFlowM3PerS()`. The probe computes it after the run and
  /// prints the absolute m3/s it set, which is what a case states.
  ///
  /// The failing states - a zero diameter, an orientation that is neither of the two, a
  /// liquid level of one, a liquid-only feed, a blank provenance - are the class's own
  /// `UNASSESSED` answers, not exceptions; each prints its diagnostic.
  static void separatorCapacityAssessmentRows() {
    assessmentRow("pr_feasible_sb_governs", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.8, null, "design basis A");
    assessmentRow("pr_infeasible", false, false, null, new String[] {"methane", "n-butane"},
        new double[] {0.7, 0.3}, 300.0, 20.0, 50.0, null, 0.2, 0.11, "vertical", 0.8, null,
        "an undersized vessel");
    assessmentRow("pr_mechanical_governs", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.8, 0.8, "mechanical design basis");
    assessmentRow("pr_mechanical_above_sb", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.8, 1.2, "a mechanical limit that does not bind");
    assessmentRow("pr_blank_provenance", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 0.8, null, " ");
    assessmentRow("pr_liquid_only", false, false, null, new String[] {"n-butane"},
        new double[] {1.0}, 250.0, 5.0, 1.0, null, 1.0, 0.11, "vertical", 0.8, null,
        "liquid only");
    // The class checks the gas phase before the liquid level, so this answers the geometry
    // sentence rather than the level's.
    assessmentRow("pr_liquid_only_level_one", false, false, null, new String[] {"n-butane"},
        new double[] {1.0}, 250.0, 5.0, 1.0, null, 1.0, 0.11, "vertical", 1.0, null,
        "liquid only");
    assessmentRow("pr_zero_k", false, false, null, new String[] {"methane", "n-butane"},
        new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null, 1.0, 0.0, "vertical", 0.8, null,
        "design basis A");
    // `1e308 * sqrt((rho_l - rho_g)/rho_g)` overflows, which is the class's non-finite
    // evidence branch rather than a geometry refusal.
    assessmentRow("pr_overflowing_k", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 1.0e308, "vertical", 0.8, null, "design basis A");
    assessmentRow("pr_horizontal_feasible", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "horizontal", 0.8, null, "design basis A");
    assessmentRow("pr_zero_diameter", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        0.0, 0.11, "vertical", 0.8, null, "design basis A");
    assessmentRow("pr_invalid_orientation", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "inclined", 0.8, null, "design basis A");
    assessmentRow("pr_level_one", false, false, null,
        new String[] {"methane", "n-butane"}, new double[] {0.7, 0.3}, 300.0, 20.0, 1.0, null,
        1.0, 0.11, "vertical", 1.0, null, "design basis A");

    // The pin's #3154 state, mirrored on PR: 3 m vertical, methane/ethane at 50000 kg/hr.
    assessmentRow("pr_3154_feasible", false, false, null,
        new String[] {"methane", "ethane"}, new double[] {95.0, 5.0}, 303.15, 50.0, null,
        50000.0, 3.0, 0.11, "vertical", 0.8, null, "Issue #3154 design basis");
    assessmentRow("pr_3154_mechanical", false, false, null,
        new String[] {"methane", "ethane"}, new double[] {95.0, 5.0}, 303.15, 50.0, null,
        50000.0, 3.0, 0.11, "vertical", 0.8, 0.8, "Issue #3154 mechanical design basis");
    assessmentRow("pr_3154_blank_provenance", false, false, null,
        new String[] {"methane", "ethane"}, new double[] {95.0, 5.0}, 303.15, 50.0, null,
        50000.0, 3.0, 0.11, "vertical", 0.8, null, " ");

    // SRK evidence, the cubic `SeparatorTest` itself uses.
    assessmentRow("srk_3154_feasible", true, false, null,
        new String[] {"methane", "ethane"}, new double[] {95.0, 5.0}, 303.15, 50.0, null,
        50000.0, 3.0, 0.11, "vertical", 0.8, null, "Issue #3154 design basis");
    assessmentRow("srk_3154_mechanical", true, false, null,
        new String[] {"methane", "ethane"}, new double[] {95.0, 5.0}, 303.15, 50.0, null,
        50000.0, 3.0, 0.11, "vertical", 0.8, 0.8, "Issue #3154 mechanical design basis");
    assessmentRow("srk_3154_blank_provenance", true, false, null,
        new String[] {"methane", "ethane"}, new double[] {95.0, 5.0}, 303.15, 50.0, null,
        50000.0, 3.0, 0.11, "vertical", 0.8, null, " ");
  }

  @SuppressWarnings({"checkstyle:MethodLength", "checkstyle:ParameterNumber"})
  static void assessmentRow(String label, boolean srk, boolean multiPhase, Double heatInputW,
      String[] names, double[] amounts, double temperature, double pressureBar,
      Double molPerSecond, Double kgPerHour, double diameter, double kFactor, String orientation,
      double designLiquidLevelFraction, Double mechanicalFraction, String provenance) {
    Stream inlet = capacityFeed(srk, multiPhase, names, amounts, temperature, pressureBar,
        molPerSecond, kgPerHour);
    neqsim.process.equipment.separator.Separator separator =
        new neqsim.process.equipment.separator.Separator("sep1", inlet);
    separator.setInternalDiameter(diameter);
    separator.setDesignGasLoadFactor(kFactor);
    separator.setOrientation(orientation);
    separator.setDesignLiquidLevelFraction(designLiquidLevelFraction);
    separator.run();

    Double mechanicalLimit = null;
    if (mechanicalFraction != null) {
      mechanicalLimit = mechanicalFraction * separator.getMaxAllowableGasFlowRate();
      // The class stores m3/hr; the snapshot reports m3/s, which is the port's input unit.
      separator.getMechanicalDesign().setMaxDesignGassVolumeFlow(mechanicalLimit * 3600.0);
    }
    neqsim.process.equipment.separator.SeparatorCapacityAssessment assessment =
        separator.getCapacityAssessment(provenance);

    System.out.println(label);
    System.out.println("eos=" + (srk ? "srk" : "pr"));
    System.out.println("temperature_k=" + inlet.getTemperature("K"));
    System.out.println("pressure_bara=" + inlet.getPressure("bara"));
    System.out.println("flow_mol_s=" + inlet.getFlowRate("mol/sec"));
    System.out.println("flow_kg_hr=" + inlet.getFlowRate("kg/hr"));
    StringBuilder z = new StringBuilder();
    for (int i = 0; i < inlet.getThermoSystem().getPhase(0).getNumberOfComponents(); i++) {
      z.append(i == 0 ? "" : ",")
          .append(inlet.getThermoSystem().getPhase(0).getComponent(i).getName()).append(":")
          .append(inlet.getThermoSystem().getMolarComposition()[i]);
    }
    System.out.println("z=" + z);
    System.out.println("internal_diameter_m=" + diameter);
    System.out.println("design_gas_load_factor_m_per_s=" + kFactor);
    System.out.println("orientation=" + orientation);
    System.out.println("design_liquid_level_fraction=" + designLiquidLevelFraction);
    System.out.println("mechanical_design_max_gas_flow_m3_s="
        + (mechanicalLimit == null ? "none" : mechanicalLimit));
    System.out.println("design_basis_provenance=" + provenance);
    System.out.println("status=" + assessment.getStatus());
    System.out.println("is_feasible=" + assessment.isFeasible());
    System.out.println("diagnostic=" + assessment.getDiagnostic());
    System.out.println("design_basis_provenance_echo=" + assessment.getDesignBasisProvenance());
    System.out.println("design_basis_source=" + assessment.getDesignBasisSource());
    System.out.println("design_basis_confidence=" + assessment.getDesignBasisConfidence());
    System.out.println("governing_capacity_basis=" + assessment.getGoverningCapacityBasis());
    System.out.println("assessment_orientation=" + assessment.getOrientation());
    System.out.println("assessment_internal_diameter_m="
        + assessment.getInternalDiameterM());
    System.out.println("assessment_design_liquid_level_fraction="
        + assessment.getDesignLiquidLevelFraction());
    System.out.println("gas_area_m2=" + assessment.getGasAreaM2());
    System.out.println("gas_density_kg_m3=" + assessment.getGasDensityKgM3());
    System.out.println("liquid_density_kg_m3=" + assessment.getLiquidDensityKgM3());
    System.out.println("design_gas_load_factor_assessed_m_per_s="
        + assessment.getDesignGasLoadFactorMPerS());
    System.out.println("maximum_gas_velocity_m_per_s="
        + assessment.getMaximumGasVelocityMPerS());
    System.out.println("current_gas_flow_m3_s=" + assessment.getCurrentGasFlowM3PerS());
    System.out.println("mechanical_design_maximum_gas_flow_m3_s="
        + assessment.getMechanicalDesignMaximumGasFlowM3PerS());
    System.out.println("souders_brown_maximum_gas_flow_m3_s="
        + assessment.getSoudersBrownMaximumGasFlowM3PerS());
    System.out.println("maximum_gas_flow_m3_s=" + assessment.getMaximumGasFlowM3PerS());
    System.out.println("utilization=" + assessment.getUtilization());
    String trimmedProvenance = provenance.trim();
    System.out.println("json_has_provenance="
        + (trimmedProvenance.isEmpty() ? "n/a" : assessment.toJson().contains(trimmedProvenance)));
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
    System.out.println("kv=" + valve.getKv());
    System.out.println("cv=" + valve.getCv());
    System.out.println("cg=" + valve.getCg());
    System.out.println("gas_valve=" + valve.isGasValve());
    neqsim.thermo.system.SystemInterface inFluid = inlet.getThermoSystem();
    System.out.println("inlet_gamma2=" + inFluid.getGamma2());
    System.out.println("inlet_density_kg_m3=" + inFluid.getDensity("kg/m3"));
    System.out.println("inlet_molar_mass_kg_mol=" + inFluid.getMolarMass());
    System.out.println("inlet_z_factor=" + inFluid.getPhase(0).getZ());
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

  /// `HeatExchanger.runDeltaT`, which `setDeltaT` enables.
  static void heatExchangerDeltaT() {
    String[] hotNames = new String[] { "methane", "n-butane" };
    double[] hotZ = new double[] { 0.8, 0.2 };
    String[] coldNames = new String[] { "n-butane", "n-pentane" };
    double[] coldZ = new double[] { 0.5, 0.5 };

    // **The class's own test fluid on PR.** The rating capture's 100 K approach is
    // degenerate here: the closure asks the hot side for more enthalpy than it has and the
    // class answers a 0 K outlet with a NaN enthalpy, which no port should reproduce.
    Stream prHot = prTestFluid(333.15, 20.0, new String[] {"methane", "ethane", "n-heptane"},
        new double[] {120.0, 120.0, 3.0}, "Stream1");
    Stream prCold = prTestFluid(313.15, 20.0, new String[] {"methane", "ethane"},
        new double[] {220.0, 120.0}, "Stream2");
    deltaTRow("pr_test_run2_delta_t_1", prHot, prCold, 1.0);
    // The streams swapped: the class's `sign` and its swap are index-based, not hot/cold.
    deltaTRow("pr_test_run2_swapped_delta_t_1",
        prTestFluid(313.15, 20.0, new String[] {"methane", "ethane"}, new double[] {220.0, 120.0},
            "Stream2"),
        prTestFluid(333.15, 20.0, new String[] {"methane", "ethane", "n-heptane"},
            new double[] {120.0, 120.0, 3.0}, "Stream1"),
        1.0);
    // A near-equal-capacity pair: both sides are the same fluid at the same flow.
    deltaTRow("pr_near_equal_cr_delta_t_1",
        feed(hotNames, new double[] {0.5, 0.5}, 400.0, 20.0, 1.0),
        feed(hotNames, new double[] {0.5, 0.5}, 300.0, 20.0, 1.0), 1.0);

    // The class's own `HeatExchangerTest.testRun2` fluid, SRK. Its expected UA is
    // `15780.77130` to `1e-3`; this row is that number's control and is not a PR case.
    neqsim.thermo.system.SystemInterface hotFluid = new neqsim.thermo.system.SystemSrkEos(
        273.15 + 60.0, 20.00);
    hotFluid.addComponent("methane", 120.00);
    hotFluid.addComponent("ethane", 120.0);
    hotFluid.addComponent("n-heptane", 3.0);
    hotFluid.createDatabase(true);
    hotFluid.setMixingRule(2);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(hotFluid).TPflash();
    Stream srkHot = new Stream("Stream1", hotFluid);

    neqsim.thermo.system.SystemInterface coldFluid = new neqsim.thermo.system.SystemSrkEos(
        273.15 + 40.0, 20.00);
    coldFluid.addComponent("methane", 220.0);
    coldFluid.addComponent("ethane", 120.0);
    coldFluid.setMixingRule(2);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(coldFluid).TPflash();
    Stream srkCold = new Stream("Stream2", coldFluid);

    deltaTRow("srk_test_run2_delta_t_1", srkHot, srkCold, 1.0);
  }

  /// `HeatExchanger.runThermalHydraulicRating`, which `setUseRatingPressureDrop` enables.
  ///
  /// **The class's own fouling-rating test, re-cased on PR.** The class test is SRK
  /// (`SystemSrkEos`); every row here is `SystemPrEos` with the databank's kij, because that
  /// is the cubic `azoth_process::Stream` builds. The SRK row is printed separately as the
  /// class's control and is never presented as a PR case.
  static void heatExchangerRating() {
    String[] names = new String[] {"n-heptane", "n-octane"};
    double[] z = new double[] {0.7, 0.3};

    // The two orientations, the equal-inlet state, the deposit, and Bell-Delaware.
    hxRatingRow("tube_30_shell_110", names, z, 303.15, 383.15, 0.0, 0.0, false);
    hxRatingRow("tube_110_shell_30", names, z, 383.15, 303.15, 0.0, 0.0, false);
    hxRatingRow("equal_inlets_60", names, z, 333.15, 333.15, 0.0, 0.0, false);
    hxRatingRow("deposit_layer", names, z, 303.15, 383.15, 0.001, 0.2, false);
    hxRatingRow("bell_delaware", names, z, 303.15, 383.15, 0.0, 0.0, true);

    // The class's own SRK control row.
    hxRatingSrkControl();

    // The refusals the class's own tests assert.
    hxRatingMultiphaseRefusal();
    hxRatingPhaseCrossingRefusal();
    hxRatingConfigurationRefusal("missing_calculator", "none");
    hxRatingConfigurationRefusal("zero_area", "zero_area");
    hxRatingConfigurationRefusal("conflicting_pin", "pin");
  }

  static void hxRatingRow(String label, String[] names, double[] z, double tubeT, double shellT,
      double depositThickness, double depositConductivity, boolean bellDelaware) {
    Stream tube = feedMass("tube feed", names, z, tubeT, 20.0, 8.0 * 3600.0);
    Stream shell = feedMass("shell feed", names, z, shellT, 20.0, 10.0 * 3600.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger exchanger =
        new neqsim.process.equipment.heatexchanger.HeatExchanger("feed preheater", tube, shell);
    neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator calculator =
        new neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator();
    calculator.setTubeIDm(0.020);
    calculator.setTubeODm(0.025);
    calculator.setTubePitchm(0.032);
    calculator.setTubeCount(100);
    calculator.setTubePasses(2);
    calculator.setTubeLengthm(4.0);
    calculator.setFoulingTube(0.0);
    calculator.setFoulingShell(0.0);
    if (depositThickness > 0.0) {
      calculator.setTubeFoulingLayer(depositThickness, depositConductivity);
    }
    if (bellDelaware) {
      calculator.setShellSideMethod(
          neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator.ShellSideMethod
              .BELL_DELAWARE);
    }
    exchanger.setRatingCalculator(calculator);
    exchanger.setRatingArea(Math.PI * 0.025 * 4.0 * 100);
    exchanger.setUseRatingPressureDrop(true);

    System.out.println(label);
    try {
      exchanger.run();
    } catch (RuntimeException ex) {
      System.out.println("run_exception=" + ex.getClass().getSimpleName());
      System.out.println("run_message=" + ex.getMessage());
      System.out.println();
      return;
    }
    printHxRating(exchanger, calculator);
    System.out.println();
  }

  /// The class's own `HeatExchangerTest` fluid, SRK: the control row, not a PR case.
  static void hxRatingSrkControl() {
    String[] names = new String[] {"n-heptane", "n-octane"};
    double[] z = new double[] {0.7, 0.3};
    Stream tube = srkMassFeed("tube feed", names, z, 303.15, 20.0, 8.0 * 3600.0);
    Stream shell = srkMassFeed("shell feed", names, z, 383.15, 20.0, 10.0 * 3600.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger exchanger =
        new neqsim.process.equipment.heatexchanger.HeatExchanger("srk control", tube, shell);
    neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator calculator =
        new neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator();
    calculator.setTubeIDm(0.020);
    calculator.setTubeODm(0.025);
    calculator.setTubePitchm(0.032);
    calculator.setTubeCount(100);
    calculator.setTubePasses(2);
    calculator.setTubeLengthm(4.0);
    calculator.setFoulingTube(0.0);
    calculator.setFoulingShell(0.0);
    exchanger.setRatingCalculator(calculator);
    exchanger.setRatingArea(Math.PI * 0.025 * 4.0 * 100);
    exchanger.setUseRatingPressureDrop(true);
    System.out.println("srk_control");
    try {
      exchanger.run();
    } catch (RuntimeException ex) {
      System.out.println("run_exception=" + ex.getClass().getSimpleName());
      System.out.println("run_message=" + ex.getMessage());
      System.out.println();
      return;
    }
    printHxRating(exchanger, calculator);
    System.out.println();
  }

  static Stream srkMassFeed(String name, String[] names, double[] z, double temperatureK,
      double pressureBara, double kgPerHour) {
    SystemInterface fluid = new SystemSrkEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream stream = new Stream(name, fluid);
    stream.setFlowRate(kgPerHour, "kg/hr");
    stream.run();
    return stream;
  }

  static void printHxRating(neqsim.process.equipment.heatexchanger.HeatExchanger exchanger,
      neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator calculator) {
    print("tube_in", exchanger.getInStream(0));
    print("shell_in", exchanger.getInStream(1));
    print("tube_out", exchanger.getOutStream(0));
    print("shell_out", exchanger.getOutStream(1));
    System.out.println("duty_w_java=" + exchanger.getDuty());
    System.out.println("rating_u_w_per_m2k=" + exchanger.getRatingU());
    System.out.println("rating_area_m2=" + exchanger.getRatingArea());
    System.out.println("ua_value_w_per_k=" + exchanger.getUAvalue());
    System.out.println("thermal_effectiveness=" + exchanger.getThermalEffectiveness());
    double t0 = exchanger.getInTemperature(0);
    double t1 = exchanger.getInTemperature(1);
    double span = t0 - t1;
    double c0 = span == 0.0 ? 0.0
        : Math.abs((ratingLimitEnthalpy(exchanger.getInStream(0), t1)
            - exchanger.getInStream(0).getThermoSystem().getEnthalpy()) / span);
    double c1 = span == 0.0 ? 0.0
        : Math.abs((ratingLimitEnthalpy(exchanger.getInStream(1), t0)
            - exchanger.getInStream(1).getThermoSystem().getEnthalpy()) / span);
    double cMin = Math.min(c0, c1);
    System.out.println("derived_c_tube_w_per_k=" + c0);
    System.out.println("derived_c_shell_w_per_k=" + c1);
    System.out.println("derived_c_min_w_per_k=" + cMin);
    System.out.println("derived_ntu=" + (cMin > 0.0 ? exchanger.getUAvalue() / cMin : 0.0));
    // The five properties the calculator received, so a divergence is localised to the
    // property layer rather than read as a correlation error.
    for (int side = 0; side < 2; side++) {
      SystemInterface fluid = exchanger.getInStream(side).getThermoSystem();
      String prefix = side == 0 ? "tube_props" : "shell_props";
      System.out.println(prefix + "_density_kg_per_m3=" + fluid.getDensity("kg/m3"));
      System.out.println(prefix + "_viscosity_Pa_s=" + fluid.getViscosity("kg/msec"));
      System.out.println(prefix + "_cp_J_per_kgK=" + fluid.getCp("J/kgK"));
      System.out.println(prefix + "_conductivity_W_per_mK=" + fluid.getThermalConductivity("W/mK"));
      System.out.println(prefix + "_mass_flow_kg_per_s=" + fluid.getFlowRate("kg/sec"));
    }
    // The calculator's own quantities, which the equipment result publishes.
    System.out.println("calc_tube_side_htc=" + calculator.getTubeSideHTC());
    System.out.println("calc_shell_side_htc=" + calculator.getShellSideHTC());
    System.out.println("calc_overall_u=" + calculator.getOverallU());
    System.out.println("calc_tube_side_pressure_drop_Pa=" + calculator.getTubeSidePressureDrop());
    System.out.println("calc_shell_side_pressure_drop_Pa=" + calculator.getShellSidePressureDrop());
    System.out.println("calc_tube_side_velocity=" + calculator.getTubeSideVelocity());
    System.out.println("calc_shell_side_velocity=" + calculator.getShellSideVelocity());
    System.out.println("calc_tube_side_re=" + calculator.getTubeSideRe());
    System.out.println("calc_shell_side_re=" + calculator.getShellSideRe());
  }

  static double ratingLimitEnthalpy(StreamInterface inlet, double temperatureK) {
    SystemInterface limit = inlet.getThermoSystem().clone();
    limit.setTemperature(temperatureK);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(limit).TPflash();
    return limit.getEnthalpy();
  }

  static void hxRatingMultiphaseRefusal() {
    Stream tube = feedMass("tube feed", new String[] {"n-heptane", "n-octane"},
        new double[] {0.7, 0.3}, 303.15, 20.0, 8.0 * 3600.0);
    Stream shell = feedMass("shell feed", new String[] {"n-heptane", "n-octane"},
        new double[] {0.7, 0.3}, 383.15, 20.0, 10.0 * 3600.0);
    SystemInterface wet = new SystemPrEos(303.15, 20.0);
    wet.addComponent("methane", 0.5);
    wet.addComponent("n-heptane", 0.5);
    wet.setMixingRule(2);
    tube.setThermoSystem(wet);
    tube.setFlowRate(8.0 * 3600.0, "kg/hr");
    tube.run();
    hxRatingRefusal("multiphase_feed", tube, shell, null);
  }

  /// The class's own `singlePhaseEndpointsOnOppositeSidesOfBoilingAreRejected`: pure
  /// n-heptane at 1 bara, 30 C against 200 C.
  static void hxRatingPhaseCrossingRefusal() {
    Stream tube = feedMass("tube feed", new String[] {"n-heptane"}, new double[] {1.0}, 303.15,
        1.0, 8.0 * 3600.0);
    Stream shell = feedMass("shell feed", new String[] {"n-heptane"}, new double[] {1.0}, 473.15,
        1.0, 1.0 * 3600.0);
    hxRatingRefusal("phase_crossing", tube, shell, null);
  }

  static void hxRatingConfigurationRefusal(String label, String mode) {
    Stream tube = feedMass("tube feed", new String[] {"n-heptane", "n-octane"},
        new double[] {0.7, 0.3}, 303.15, 20.0, 8.0 * 3600.0);
    Stream shell = feedMass("shell feed", new String[] {"n-heptane", "n-octane"},
        new double[] {0.7, 0.3}, 383.15, 20.0, 10.0 * 3600.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger exchanger =
        new neqsim.process.equipment.heatexchanger.HeatExchanger(label, tube, shell);
    exchanger.setUseRatingPressureDrop(true);
    if (!"none".equals(mode)) {
      neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator calculator =
          new neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator();
      calculator.setTubeIDm(0.020);
      calculator.setTubeODm(0.025);
      calculator.setTubePitchm(0.032);
      calculator.setTubeCount(100);
      calculator.setTubePasses(2);
      calculator.setTubeLengthm(4.0);
      exchanger.setRatingCalculator(calculator);
      exchanger.setRatingArea("zero_area".equals(mode) ? 0.0 : Math.PI * 0.025 * 4.0 * 100);
    }
    if ("pin".equals(mode)) {
      exchanger.setOutTemperature(60.0, "C");
    }
    hxRatingRefusal(label, tube, shell, exchanger);
  }

  static void hxRatingRefusal(String label, Stream tube, Stream shell,
      neqsim.process.equipment.heatexchanger.HeatExchanger prepared) {
    System.out.println(label);
    try {
      neqsim.process.equipment.heatexchanger.HeatExchanger exchanger = prepared;
      if (exchanger == null) {
        exchanger = new neqsim.process.equipment.heatexchanger.HeatExchanger(label, tube, shell);
        neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator calculator =
            new neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator();
        calculator.setTubeIDm(0.020);
        calculator.setTubeODm(0.025);
        calculator.setTubePitchm(0.032);
        calculator.setTubeCount(100);
        calculator.setTubePasses(2);
        calculator.setTubeLengthm(4.0);
        exchanger.setRatingCalculator(calculator);
        exchanger.setRatingArea(Math.PI * 0.025 * 4.0 * 100);
        exchanger.setUseRatingPressureDrop(true);
      }
      exchanger.run();
      System.out.println("run_result=accepted");
    } catch (RuntimeException ex) {
      System.out.println("run_exception=" + ex.getClass().getSimpleName());
      System.out.println("run_message=" + ex.getMessage());
    }
    System.out.println();
  }

  static Stream prTestFluid(double temperature, double pressureBar, String[] names,
      double[] amounts, String streamName) {
    neqsim.thermo.system.SystemInterface fluid = new neqsim.thermo.system.SystemPrEos(temperature,
        pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], amounts[i]);
    }
    fluid.setMixingRule(2);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(fluid).TPflash();
    return new Stream(streamName, fluid);
  }

  static void deltaTRow(String label, Stream firstIn, Stream secondIn, double deltaT) {
    neqsim.process.equipment.heatexchanger.HeatExchanger exchanger =
        new neqsim.process.equipment.heatexchanger.HeatExchanger("hx", firstIn, secondIn);
    exchanger.setDeltaT(deltaT);
    exchanger.run();

    System.out.println(label);
    print("first_in", firstIn);
    print("second_in", secondIn);
    print("first_out", exchanger.getOutStream(0));
    print("second_out", exchanger.getOutStream(1));
    System.out.println("delta_t_k=" + deltaT);
    // `getDuty` carries the class's sign: the enthalpy change of the side it set.
    System.out.println("duty_w_java=" + exchanger.getDuty());
    System.out.println("ua_value_w_per_k=" + exchanger.getUAvalue());
    System.out.println();
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

  // === process.heat_exchanger, runTransient =====================================

  /// `HeatExchanger.runTransient`'s enabled branch, one step at a time.
  ///
  /// **The wall ODE is flow-independent** while no holdup is stated - it reads only the two
  /// inlet temperatures - so the wall trajectory is the same at every flow. The flow moves only
  /// the outlet enthalpy shift, and the class applies that shift as a **power where a joule
  /// belongs**: `outH = inH - q`, not `q * dt`, despite the comment in the source. The upstream
  /// `1000 kg/hr` rows are therefore printed as evidence, and the cased wall rows use
  /// `100 000 kg/hr` so the PH flash stays in the range both libraries can evaluate.
  ///
  /// **The two CSTR temperatures have no getters**, so they are read by reflection; the capture
  /// prints them before and after every step, which is what lets a stateless one-step port
  /// replay the trajectory. The class's stream 0 is the **shell** and stream 1 the **tube** -
  /// the reverse of the thermal-hydraulic rating's stream-0-is-tube convention.
  static void heatExchangerTransientRows() {
    String[] methane = new String[] {"methane"};
    double[] pure = new double[] {1.0};

    // 1. The gate: a disabled dynamic model delegates to Heater.runTransient, the steady run.
    Stream steadyHot = feedMass("steady_hot", methane, pure, 353.15, 10.0, 100000.0);
    Stream steadyCold = feedMass("steady_cold", methane, pure, 293.15, 10.0, 100000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger steady = dynamicExchanger(
        "gate_disabled", steadyHot, steadyCold, 200.0, 20.0, 500.0, 800.0);
    steady.setDynamicModelEnabled(false);
    transientStep("gate_disabled_default", steady, steadyHot, steadyCold, 1.0);

    // The class's other two gates.
    Stream massZeroHot = feedMass("mass_zero_hot", methane, pure, 353.15, 10.0, 100000.0);
    Stream massZeroCold = feedMass("mass_zero_cold", methane, pure, 293.15, 10.0, 100000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger massZero = dynamicExchanger(
        "gate_wall_mass_zero", massZeroHot, massZeroCold, 0.0, 20.0, 500.0, 800.0);
    transientStep("gate_wall_mass_zero", massZero, massZeroHot, massZeroCold, 1.0);

    Stream areaZeroHot = feedMass("area_zero_hot", methane, pure, 353.15, 10.0, 100000.0);
    Stream areaZeroCold = feedMass("area_zero_cold", methane, pure, 293.15, 10.0, 100000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger areaZero = dynamicExchanger(
        "gate_area_zero", areaZeroHot, areaZeroCold, 200.0, 0.0, 500.0, 800.0);
    transientStep("gate_area_zero", areaZero, areaZeroHot, areaZeroCold, 1.0);

    // 2. The first-step seed and ten steps of the wall-only ODE, at the cased flow.
    Stream wallHot = feedMass("wall_hot", methane, pure, 353.15, 10.0, 100000.0);
    Stream wallCold = feedMass("wall_cold", methane, pure, 293.15, 10.0, 100000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger wall = dynamicExchanger(
        "wall_ten_steps", wallHot, wallCold, 200.0, 20.0, 500.0, 800.0);
    java.util.UUID wallId = java.util.UUID.randomUUID();
    for (int step = 1; step <= 10; step++) {
      transientStep("wall_only_step_" + step, wall, wallHot, wallCold, 1.0, wallId);
    }

    // 3. The same trajectory at the upstream test's flow. Evidence only: the rate-as-joules
    //    shift is about -17 kJ/mol here, which neither cubic's PH flash is meaningful over.
    Stream upstreamHot = feedMass("upstream_hot", methane, pure, 353.15, 10.0, 1000.0);
    Stream upstreamCold = feedMass("upstream_cold", methane, pure, 293.15, 10.0, 1000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger upstream = dynamicExchanger(
        "wall_upstream", upstreamHot, upstreamCold, 200.0, 20.0, 500.0, 800.0);
    for (int step = 1; step <= 10; step++) {
      try {
        transientStep("wall_upstream_1000kg_step_" + step, upstream, upstreamHot, upstreamCold, 1.0,
            wallId);
      } catch (RuntimeException failure) {
        System.out.println("wall_upstream_1000kg_step_" + step + "_outcome=threw:"
            + failure.getClass().getSimpleName());
        System.out.println();
        break;
      }
    }

    // 4. A stated wall is carried, not seeded.
    Stream statedHot = feedMass("stated_hot", methane, pure, 353.15, 10.0, 100000.0);
    Stream statedCold = feedMass("stated_cold", methane, pure, 293.15, 10.0, 100000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger stated = dynamicExchanger(
        "stated_wall", statedHot, statedCold, 200.0, 20.0, 500.0, 800.0);
    stated.setWallTemperature(323.15);
    transientStep("first_step_stated_wall", stated, statedHot, statedCold, 1.0);

    // 5. The outlet shift does not scale with dt; the wall does.
    for (double dt : new double[] {0.0, 0.5, 1.0}) {
      Stream quirkHot = feedMass("quirk_hot_" + dt, methane, pure, 353.15, 10.0, 100000.0);
      Stream quirkCold = feedMass("quirk_cold_" + dt, methane, pure, 293.15, 10.0, 100000.0);
      neqsim.process.equipment.heatexchanger.HeatExchanger quirk = dynamicExchanger(
          "dt_quirk_" + dt, quirkHot, quirkCold, 200.0, 20.0, 500.0, 800.0);
      quirk.setWallTemperature(323.15);
      transientStep("dt_quirk_" + dt, quirk, quirkHot, quirkCold, dt);
    }

    // 6. The two CSTR states, and the zero-holdup control.
    for (boolean holdup : new boolean[] {true, false}) {
      String label = holdup ? "holdup" : "holdup_zero";
      Stream cstrHot = feedMass(label + "_hot", methane, pure, 423.15, 30.0, 100000.0);
      Stream cstrCold = feedMass(label + "_cold", methane, pure, 303.15, 30.0, 100000.0);
      neqsim.process.equipment.heatexchanger.HeatExchanger cstr = dynamicExchanger(
          label, cstrHot, cstrCold, 400.0, 2.0, 500.0, 800.0);
      if (holdup) {
        cstr.setShellHoldupVolume(40.0);
        cstr.setTubeHoldupVolume(30.0);
      }
      java.util.UUID cstrId = java.util.UUID.randomUUID();
      for (int step = 1; step <= 5; step++) {
        transientStep(label + "_step_" + step, cstr, cstrHot, cstrCold, 1.0, cstrId);
      }
    }

    // 7. The setters' clamps.
    Stream clampHot = feedMass("clamp_hot", methane, pure, 353.15, 10.0, 1000.0);
    Stream clampCold = feedMass("clamp_cold", methane, pure, 293.15, 10.0, 1000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger clamp = dynamicExchanger(
        "clamp", clampHot, clampCold, 200.0, 20.0, 500.0, 800.0);
    clamp.setWallMass(-10.0);
    clamp.setShellHoldupVolume(-1.0);
    clamp.setTubeHoldupVolume(-1.0);
    clamp.setHeatTransferArea(-5.0);
    System.out.println("clamp_getters");
    System.out.println("wall_mass_kg=" + clamp.getWallMass());
    System.out.println("shell_holdup_volume_m3=" + clamp.getShellHoldupVolume());
    System.out.println("tube_holdup_volume_m3=" + clamp.getTubeHoldupVolume());
    System.out.println("heat_transfer_area_m2=" + clamp.getHeatTransferArea());
    System.out.println();

    // 8. The side mapping: the class's stream 0 is the shell, so constructing it with the cold
    //    stream first reverses the two sides and the class's duty sign follows.
    Stream mappingCold = feedMass("mapping_cold", methane, pure, 293.15, 10.0, 100000.0);
    Stream mappingHot = feedMass("mapping_hot", methane, pure, 353.15, 10.0, 100000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger mapping = dynamicExchanger(
        "tube_hot", mappingCold, mappingHot, 200.0, 20.0, 500.0, 800.0);
    transientStep("tube_hot_mapping", mapping, mappingCold, mappingHot, 1.0);

    // 9. Evidence that is not cased.
    Stream zeroHot = feedMass("zero_hot", methane, pure, 353.15, 10.0, 0.0);
    Stream zeroCold = feedMass("zero_cold", methane, pure, 293.15, 10.0, 100000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger zero = dynamicExchanger(
        "zero_flow", zeroHot, zeroCold, 200.0, 20.0, 500.0, 800.0);
    try {
      zero.runTransient(1.0, java.util.UUID.randomUUID());
      System.out.println("zero_flow_uncased=answered");
      System.out.println("wall_temperature_K=" + zero.getWallTemperature());
    } catch (RuntimeException failure) {
      System.out.println("zero_flow_uncased=threw:" + failure.getClass().getSimpleName());
    }
    System.out.println();

    String[] mix = new String[] {"methane", "n-butane"};
    double[] half = new double[] {0.5, 0.5};
    Stream twoPhaseHot = feedMass("two_phase_hot", mix, half, 250.0, 5.0, 1000.0);
    Stream twoPhaseCold = feedMass("two_phase_cold", methane, pure, 293.15, 5.0, 1000.0);
    System.out.println("two_phase_uncased");
    System.out.println("shell_phases=" + twoPhaseHot.getThermoSystem().getNumberOfPhases());
    neqsim.process.equipment.heatexchanger.HeatExchanger twoPhase = dynamicExchanger(
        "two_phase", twoPhaseHot, twoPhaseCold, 200.0, 20.0, 500.0, 800.0);
    try {
      twoPhase.runTransient(1.0, java.util.UUID.randomUUID());
      System.out.println("two_phase_outcome=answered");
      System.out.println("wall_temperature_K=" + twoPhase.getWallTemperature());
      print("shell_out", twoPhase.getOutStream(0));
    } catch (RuntimeException failure) {
      System.out.println("two_phase_outcome=threw:" + failure.getClass().getSimpleName());
    }
    System.out.println();

    // 10. calculateSteadyState=false: Heater.runTransient copies the outlet pressure to the
    //     inlet before re-running. On the ordinary (pressure-drop-free) rating this is a no-op;
    //     it is printed so the non-port has its evidence.
    Stream cssHot = feedMass("css_hot", methane, pure, 353.15, 10.0, 1000.0);
    Stream cssCold = feedMass("css_cold", methane, pure, 293.15, 10.0, 1000.0);
    neqsim.process.equipment.heatexchanger.HeatExchanger css = dynamicExchanger(
        "steady_state_false", cssHot, cssCold, 200.0, 20.0, 500.0, 800.0);
    css.setDynamicModelEnabled(false);
    css.run();
    System.out.println("calculate_steady_state_false_uncased");
    System.out.println("shell_in_pressure_bara_before=" + cssHot.getPressure("bara"));
    System.out.println("shell_out_pressure_bara_before=" + css.getOutStream(0).getPressure("bara"));
    css.setCalculateSteadyState(false);
    try {
      css.runTransient(1.0, java.util.UUID.randomUUID());
      System.out.println("calculate_steady_state_false_outcome=answered");
    } catch (RuntimeException failure) {
      System.out.println("calculate_steady_state_false_outcome=threw:"
          + failure.getClass().getSimpleName());
    }
    System.out.println("shell_in_pressure_bara_after=" + cssHot.getPressure("bara"));
    System.out.println();
  }

  /// One step of the dynamic exchanger, with the carried state printed before and after it.
  static void transientStep(String label,
      neqsim.process.equipment.heatexchanger.HeatExchanger hx, StreamInterface shellIn,
      StreamInterface tubeIn, double dt) {
    transientStep(label, hx, shellIn, tubeIn, dt, java.util.UUID.randomUUID());
  }

  static void transientStep(String label,
      neqsim.process.equipment.heatexchanger.HeatExchanger hx, StreamInterface shellIn,
      StreamInterface tubeIn, double dt, java.util.UUID id) {
    System.out.println(label);
    System.out.println("time_step_s=" + dt);
    System.out.println("in_wall_temperature_K=" + hx.getWallTemperature());
    System.out.println(
        "in_shell_fluid_temperature_K=" + privateDouble(hx, "shellFluidTemperature"));
    System.out.println(
        "in_tube_fluid_temperature_K=" + privateDouble(hx, "tubeFluidTemperature"));
    hx.runTransient(dt, id);
    System.out.println("out_wall_temperature_K=" + hx.getWallTemperature());
    System.out.println(
        "out_shell_fluid_temperature_K=" + privateDouble(hx, "shellFluidTemperature"));
    System.out.println(
        "out_tube_fluid_temperature_K=" + privateDouble(hx, "tubeFluidTemperature"));
    print("shell_in", shellIn);
    print("tube_in", tubeIn);
    print("shell_out", hx.getOutStream(0));
    print("tube_out", hx.getOutStream(1));
    System.out.println("shell_duty_W=" + duty(shellIn, hx.getOutStream(0)));
    System.out.println("tube_duty_W=" + duty(tubeIn, hx.getOutStream(1)));
    System.out.println("duty_W_java=" + hx.getDuty());
    System.out.println();
  }

  /// A dynamic exchanger with the class's own defaults for the fields the capture varies.
  static neqsim.process.equipment.heatexchanger.HeatExchanger dynamicExchanger(String name,
      StreamInterface shellIn, StreamInterface tubeIn, double wallMass, double area,
      double shellHtc, double tubeHtc) {
    neqsim.process.equipment.heatexchanger.HeatExchanger hx =
        new neqsim.process.equipment.heatexchanger.HeatExchanger(name, shellIn, tubeIn);
    hx.setUAvalue(100.0);
    hx.setWallMass(wallMass);
    hx.setWallCp(500.0);
    hx.setHeatTransferArea(area);
    hx.setShellSideHtc(shellHtc);
    hx.setTubeSideHtc(tubeHtc);
    hx.setDynamicModelEnabled(true);
    return hx;
  }

  /// A private field with no getter, which is how the two CSTR temperatures are reached.
  static double privateDouble(Object target, String name) {
    try {
      java.lang.reflect.Field field = target.getClass().getDeclaredField(name);
      field.setAccessible(true);
      return field.getDouble(target);
    } catch (ReflectiveOperationException failure) {
      throw new IllegalStateException("no private double " + name, failure);
    }
  }

  // === process.cooler, runTransient ============================================

  /// `Cooler.runTransient`'s dynamic temperature-control branch, one step at a time.
  ///
  /// The enabled branch is a first-order valve model: the requested opening moves the valve
  /// by `dt / (tau_actuator + dt)`, an NTU scaled by opening and inverse mass flow gives the
  /// equilibrium outlet, and a second first-order lag moves the outlet toward it. Without a
  /// controller the requested opening is the step-initial one, so the valve does not move and
  /// only the thermal lag acts; the manual-controller rows supply the requested opening.
  static void coolerTransientRows() {
    // 1. The gate: no dynamic configuration -> Heater.runTransient, the steady run.
    neqsim.process.equipment.heatexchanger.Cooler steady =
        new neqsim.process.equipment.heatexchanger.Cooler("cooler_gate_disabled",
            coolerFeed("cooler_gate_disabled", 353.15));
    steady.setOutTemperature(35.0, "C");
    steady.run();
    steady.setCalculateSteadyState(false);
    coolerStep("cooler_gate_disabled", steady, 2.0);

    // 2. Configured, no controller: the requested opening is the step-initial one.
    neqsim.process.equipment.heatexchanger.Cooler configured =
        new neqsim.process.equipment.heatexchanger.Cooler("cooler_configured",
            coolerFeed("cooler_configured", 353.15));
    configured.setOutTemperature(35.0, "C");
    configured.configureDynamicTemperatureControl(10000.0, 80.0, 35.0, 20.0, "C");
    configured.setDynamicTimeConstants(5.0, 10.0);
    configured.run();
    configured.setCalculateSteadyState(false);
    coolerStep("cooler_configured_no_controller", configured, 2.0);

    // 2b. A hotter inlet: the design point no longer matches, and the thermal lag moves the
    //     outlet toward the equilibrium the current opening reaches.
    neqsim.process.equipment.heatexchanger.Cooler hotter =
        new neqsim.process.equipment.heatexchanger.Cooler("cooler_hotter_inlet",
            coolerFeed("cooler_hotter_inlet", 353.15));
    hotter.setOutTemperature(35.0, "C");
    hotter.configureDynamicTemperatureControl(10000.0, 80.0, 35.0, 20.0, "C");
    hotter.setDynamicTimeConstants(5.0, 10.0);
    hotter.run();
    hotter.setCalculateSteadyState(false);
    hotter.getInletStream().setTemperature(100.0, "C");
    hotter.getInletStream().run();
    coolerStep("cooler_configured_hotter_inlet", hotter, 2.0);

    // 2c. Two steps: the second call's inputs are the first call's outputs.
    neqsim.process.equipment.heatexchanger.Cooler twoSteps =
        coolerController("cooler_two_steps", 353.15, 100.0);
    coolerStep("cooler_manual_100_step_1", twoSteps, 2.0);
    coolerStep("cooler_manual_100_step_2", twoSteps, 2.0);

    // 3. A manual controller at 100%: the actuator moves the valve before the thermal lag.
    coolerStep("cooler_manual_100", coolerController("cooler_manual_100", 353.15, 100.0), 2.0);

    // 4. Manual at 0%: the valve closes and the equilibrium is the inlet.
    coolerStep("cooler_manual_0", coolerController("cooler_manual_0", 353.15, 0.0), 2.0);

    // 5. Zero time constants and a feed colder than the medium: the outlet is the inlet.
    neqsim.process.equipment.heatexchanger.Cooler cold =
        new neqsim.process.equipment.heatexchanger.Cooler("cooler_cold",
            coolerFeed("cooler_cold", 283.15));
    cold.setOutTemperature(35.0, "C");
    cold.configureDynamicTemperatureControl(10000.0, 80.0, 35.0, 20.0, "C");
    cold.setDynamicTimeConstants(0.0, 0.0);
    cold.run();
    cold.setCalculateSteadyState(false);
    cold.getInletStream().setTemperature(10.0, "C");
    cold.getInletStream().run();
    coolerStep("cooler_cold_feed_zero_taus", cold, 1.0);

    // 6. The opening clamp: a requested 150 or -10 is clamped to [0, 100].
    coolerStep(
        "cooler_manual_150_clamped", coolerController("cooler_manual_150", 353.15, 150.0), 2.0);
    coolerStep(
        "cooler_manual_negative_clamped", coolerController("cooler_manual_neg", 353.15, -10.0), 2.0);

    // 7. The NTU clamp: a low flow with the valve open reaches the 20 cap.
    neqsim.process.equipment.heatexchanger.Cooler lowFlow =
        new neqsim.process.equipment.heatexchanger.Cooler("cooler_low_flow",
            coolerFeed("cooler_low_flow", 353.15, 100.0));
    lowFlow.setOutTemperature(35.0, "C");
    lowFlow.configureDynamicTemperatureControl(10000.0, 80.0, 35.0, 20.0, "C");
    lowFlow.setDynamicTimeConstants(5.0, 10.0);
    lowFlow.run();
    lowFlow.setCalculateSteadyState(false);
    lowFlow.setCoolingValveOpening(100.0);
    coolerStep("cooler_ntu_clamp_low_flow", lowFlow, 2.0);

    // 8. Evidence: the same UUID does not advance the state twice. That cache is session
    //    state, and this one-step id has no step identity; the capture records it uncased.
    neqsim.process.equipment.heatexchanger.Cooler cached =
        coolerController("cooler_cached", 353.15, 100.0);
    java.util.UUID shared = java.util.UUID.randomUUID();
    System.out.println("cooler_same_uuid_first");
    cached.runTransient(2.0, shared);
    System.out.println("out_cooling_valve_opening_percent=" + cached.getCoolingValveOpening());
    System.out.println("out_outlet_temperature_K=" + cached.getOutletStream().getTemperature("K"));
    System.out.println("cooler_same_uuid_second");
    cached.runTransient(2.0, shared);
    System.out.println("out_cooling_valve_opening_percent=" + cached.getCoolingValveOpening());
    System.out.println("out_outlet_temperature_K=" + cached.getOutletStream().getTemperature("K"));
    System.out.println();

    // 9. An invalid design point is the class's own IllegalArgumentException.
    try {
      neqsim.process.equipment.heatexchanger.Cooler invalid =
          new neqsim.process.equipment.heatexchanger.Cooler("cooler_invalid",
              coolerFeed("cooler_invalid", 353.15));
      invalid.configureDynamicTemperatureControl(10000.0, 80.0, 15.0, 20.0, "C");
      System.out.println("cooler_invalid_design=answered");
    } catch (RuntimeException failure) {
      System.out.println("cooler_invalid_design=threw:" + failure.getClass().getSimpleName());
    }
    System.out.println();
  }

  static Stream coolerFeed(String name, double kelvin) {
    return coolerFeed(name, kelvin, 10000.0);
  }

  static Stream coolerFeed(String name, double kelvin, double kgPerHour) {
    return feedMass(name, new String[] {"methane"}, new double[] {1.0}, kelvin, 50.0, kgPerHour);
  }

  /// A configured cooler with a manual controller, which is how a requested opening differs
  /// from the step-initial one.
  static neqsim.process.equipment.heatexchanger.Cooler coolerController(
      String name, double kelvin, double manualOutput) {
    neqsim.process.equipment.heatexchanger.Cooler cooler =
        new neqsim.process.equipment.heatexchanger.Cooler(name, coolerFeed(name, kelvin));
    cooler.setOutTemperature(35.0, "C");
    cooler.configureDynamicTemperatureControl(10000.0, 80.0, 35.0, 20.0, "C");
    cooler.setDynamicTimeConstants(5.0, 10.0);
    neqsim.process.controllerdevice.ControllerDeviceBaseClass controller =
        new neqsim.process.controllerdevice.ControllerDeviceBaseClass(name + " controller");
    controller.setTransmitter(
        new neqsim.process.measurementdevice.TemperatureTransmitter(name + " outlet",
            cooler.getOutletStream()));
    controller.setMode(
        neqsim.process.controllerdevice.ControllerDeviceInterface.ControllerMode.MANUAL);
    // The limits are wider than the class's own, so the clamp rows exercise the class's.
    controller.setOutputLimits(-20.0, 200.0);
    controller.setManualOutput(manualOutput);
    cooler.setController(controller);
    cooler.run();
    cooler.setCalculateSteadyState(false);
    return cooler;
  }

  /// One dynamic-cooler step, with the carried opening and outlet printed before and after it.
  static void coolerStep(String label, neqsim.process.equipment.heatexchanger.Cooler cooler,
      double dt) {
    System.out.println(label);
    System.out.println("time_step_s=" + dt);
    System.out.println("in_cooling_valve_opening_percent=" + cooler.getCoolingValveOpening());
    System.out.println("in_outlet_temperature_K=" + cooler.getOutletStream().getTemperature("K"));
    cooler.runTransient(dt, java.util.UUID.randomUUID());
    System.out.println("out_cooling_valve_opening_percent=" + cooler.getCoolingValveOpening());
    System.out.println("out_outlet_temperature_K=" + cooler.getOutletStream().getTemperature("K"));
    System.out.println("out_duty_W_java=" + cooler.getDuty());
    print("outlet", cooler.getOutletStream());
    System.out.println();
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

  /// **The flare, whose whole steady state is its report.**
  ///
  /// `Flare.run` clones the inlet into the outlet and computes two numbers beside it: the
  /// heat duty, `inStream.LCV() * flowSm3sec`, and the CO2 emission, the carbon the gas
  /// carries times `44.01e-3` kg/mol. Both are printed, and so are the inlet and the outlet,
  /// because the pass-through is half of what the class does - a reader has to see that the
  /// record does not move.
  ///
  /// **The two volumetric references are not the same one, and the duty carries both.**
  /// `Stream.LCV()` builds `new Standard_ISO6976(fluid, 0, 15.55, "volume")` - joules per
  /// normal cubic metre at 0 C - while `getFlowRate("Sm3/sec")` is
  /// `n * R * standardStateTemperature / atm` with `standardStateTemperature` 288.15 K, a
  /// volume at 15 C. The product is therefore an energy density at one reference times a
  /// volumetric flow at another, and the duty it reports is 5.5 per cent above the same gas
  /// measured consistently. Reproduced here rather than corrected: it is what the class
  /// **The stirred-tank reactor, whose reaction is a stoichiometry map and a conversion.**
  ///
  /// `StoichiometricReaction.react` takes the limiting reactant's moles times the conversion,
  /// scales every coefficient by `coeff / |coeff_limiting|`, and adds the result to the fluid -
  /// so the reaction is a *movement of moles* and the enthalpy that follows is the flash's.
  /// The class then sets the outlet pressure (a stated one, or the inlet less the drop), holds
  /// the temperature when it is isothermal, and flashes: `TPflash` when isothermal, and
  /// `PHflash(inlet enthalpy)` when adiabatic.
  ///
  /// The reaction is `methanecombustion` - CH4 + 2 O2 -> CO2 + 2 H2O, the one entry the
  /// reaction data carries whose substances are all components here - and the three rows are
  /// the two flash branches plus a stated reactor pressure, on the same feed.
  static void stirredTankReactorRows() {
    reactorRow("adiabatic", false, null, null);
    reactorRow("isothermal_800K", true, 800.0, null);
    reactorRow("isothermal_800K_3bara", true, 800.0, 3.0);
    // The three stored fields are unread by `run`; this row moves all three so the capture
    // shows `toMap`'s echo and `getAgitatorPower`'s product rather than only the defaults.
    reactorRow("stored_fields_set", false, null, null, 20.0, 0.5, 2.0);
  }

  static void reactorRow(String label, boolean isothermal, Double reactorT, Double reactorP) {
    reactorRow(label, isothermal, reactorT, reactorP, null, null, null);
  }

  static void reactorRow(String label, boolean isothermal, Double reactorT, Double reactorP,
      Double vesselVolumeM3, Double residenceHours, Double agitatorPowerPerVolume) {
    String[] names = new String[] { "methane", "oxygen", "CO2", "water", "nitrogen" };
    double[] z = new double[] { 0.05, 0.10, 0.02, 0.03, 0.80 };
    Stream inlet = feed(names, z, 500.0, 5.0, 1.0);
    neqsim.process.equipment.reactor.StirredTankReactor reactor =
        new neqsim.process.equipment.reactor.StirredTankReactor("cstr", inlet);
    reactor.setIsothermal(isothermal);
    if (reactorT != null) {
      reactor.setReactorTemperature(reactorT);
    }
    if (reactorP != null) {
      reactor.setReactorPressure(reactorP);
    }
    if (vesselVolumeM3 != null) {
      reactor.setVesselVolume(vesselVolumeM3);
    }
    if (residenceHours != null) {
      reactor.setResidenceTime(residenceHours, "hr");
    }
    if (agitatorPowerPerVolume != null) {
      reactor.setAgitatorPowerPerVolume(agitatorPowerPerVolume);
    }
    neqsim.process.equipment.reactor.StoichiometricReaction reaction =
        new neqsim.process.equipment.reactor.StoichiometricReaction("methanecombustion");
    reaction.addReactant("methane", 1.0);
    reaction.addReactant("oxygen", 2.0);
    reaction.addProduct("CO2", 1.0);
    reaction.addProduct("water", 2.0);
    reaction.setLimitingReactant("oxygen");
    reaction.setConversion(0.5);
    reactor.addReaction(reaction);
    reactor.run();
    System.out.println(label);
    print("feed", inlet);
    printOrEmpty("product", reactor.getOutletStream());
    System.out.println("heat_duty_W=" + reactor.getHeatDuty());
    System.out.println("vessel_volume_m3=" + reactor.getVesselVolume());
    System.out.println("residence_time_hr=" + reactor.getResidenceTime());
    System.out.println("agitator_power_per_volume_kW_per_m3=" + reactor.getAgitatorPowerPerVolume());
    System.out.println("agitator_power_kW=" + reactor.getAgitatorPower());
    System.out.println();
  }

  /// **A plug-flow reactor is a profile, so the profile is the capture.**
  ///
  /// `PlugFlowReactor.run` marches `[F1..Fn, T, P]` along the tube with a scheme it writes out
  /// inline - `IntegrationMethod.EULER` or `RK4`, RK4 by default, a fixed step
  /// `length / numberOfSteps`, no step-size control - and stores every station in a
  /// `ReactorAxialProfile` whose `toCSV()` is printed here verbatim. A single outlet row would
  /// pin one number where the answer is a curve.
  ///
  /// The reaction is a power-law methane combustion, so the products are species the feed does
  /// not carry: `ensureProductComponentsExist` adds each at `1.0e-20` and re-inits, which is the
  /// behaviour that made this class a tier rather than a kernel.
  ///
  /// **Properties are frozen by default.** `thermodynamicCoupling` defaults to
  /// `FROZEN_PROPERTIES`, which re-flashes only every `propertyUpdateFrequency` (10) steps and
  /// never for the RK4 sub-states - so the rate law reads a *stale* temperature and composition
  /// between updates. The last two rows measure that: the same reactor at a frequency of 1 and
  /// the default 10 are different answers, and the difference is the coupling and not the
  /// chemistry.
  ///
  /// `concentration_probe` is the other measurement a port needs and cannot guess:
  /// `getPhase(0).getDensity("mol/m3")` against `getDensity("kg/m3") / molarMass`, which says
  /// whether the unit-carrying molar density carries the Peneloux volume translation.
  static void plugFlowReactorRows() {
    plugFlowReactorRow("rk4_default", "RK4", 100, 10, false, 1.0);
    plugFlowReactorRow("euler_default", "EULER", 100, 10, false, 1.0);
    plugFlowReactorRow("rk4_coarse", "RK4", 50, 10, false, 1.0);
    plugFlowReactorRow("rk4_every_step", "RK4", 100, 1, false, 1.0);
    plugFlowReactorRow("isothermal_rk4", "RK4", 100, 10, false, 1.0);
    plugFlowReactorRow("bed_rk4", "RK4", 100, 10, true, 1.0);
    plugFlowReactorRow("bed_half_activity", "RK4", 100, 10, true, 0.5);
    // **The ten-step rows are the case's, and they are the setter's own floor.**
    // `setNumberOfSteps` clamps to `max(10, n)`, and at ten steps with the default frequency
    // of ten the refresh never fires inside the loop - `s % 10 == 0 && s < 10` is false for
    // every `s` - so these rows are a *purely frozen* march, and the eleven stations are short
    // enough for a case file to carry.
    plugFlowReactorRow("rk4_ten_steps", "RK4", 10, 10, false, 1.0);
    plugFlowReactorRow("euler_ten_steps", "EULER", 10, 10, false, 1.0);
    concentrationProbe();
  }

  static void plugFlowReactorRow(String label, String method, int steps, int updateEvery,
      boolean withBed, double activity) {
    String[] names = new String[] { "methane", "oxygen", "nitrogen" };
    double[] z = new double[] { 0.05, 0.10, 0.85 };
    Stream inlet = feed(names, z, 600.0, 5.0, 1.0);

    neqsim.process.equipment.reactor.PlugFlowReactor reactor =
        new neqsim.process.equipment.reactor.PlugFlowReactor("pfr", inlet);
    reactor.setLength(5.0, "m");
    reactor.setDiameter(0.10, "m");
    reactor.setNumberOfSteps(steps);
    reactor.setIntegrationMethod(method);
    reactor.setPropertyUpdateFrequency(updateEvery);
    if (label.startsWith("isothermal")) {
      reactor.setEnergyMode(neqsim.process.equipment.reactor.PlugFlowReactor.EnergyMode.ISOTHERMAL);
    }
    if (withBed) {
      // **Three millimetre pellets at a void fraction of 0.40 are the class's own defaults**,
      // and a bed is what makes the Ergun equation the pressure row rather than the empty-tube
      // Darcy branch - measured, the two differ by four orders of magnitude on this gas.
      neqsim.process.equipment.reactor.CatalystBed bed =
          new neqsim.process.equipment.reactor.CatalystBed(3.0, 0.40, 800.0);
      bed.setActivityFactor(activity);
      reactor.setCatalystBed(bed);
    }

    neqsim.process.equipment.reactor.KineticReaction reaction =
        new neqsim.process.equipment.reactor.KineticReaction("methanecombustion");
    reaction.addReactant("methane", 1.0, 1.0);
    reaction.addReactant("oxygen", 2.0, 1.0);
    reaction.addProduct("CO2", 1.0);
    reaction.addProduct("water", 2.0);
    reaction.setPreExponentialFactor(1.0e4);
    reaction.setActivationEnergy(80000.0);
    reaction.setHeatOfReaction(-802000.0);
    reactor.addReaction(reaction);

    reactor.run();
    System.out.println(label);
    print("feed", inlet);
    printOrEmpty("product", reactor.getOutletStream());
    System.out.println("conversion=" + reactor.getConversion());
    System.out.println("pressure_drop_bar=" + reactor.getPressureDrop());
    System.out.println("outlet_temperature_K=" + reactor.getOutletTemperature());
    System.out.println("residence_time_s=" + reactor.getResidenceTime());
    System.out.println("heat_duty_W=" + reactor.getHeatDuty());
    System.out.println("profile=");
    System.out.println(reactor.getAxialProfile().toCSV());
    System.out.println();
  }

  /// **Does a molar density carry the volume translation?** One state, two ways to read it.
  ///
  /// `KineticReaction.getConcentration` is `x_i * getPhase(0).getDensity("mol/m3")`, and this
  /// probe printed once already that a *mass* density asked with a unit is the Peneloux-corrected
  /// one. Whether the molar density behaves the same way is a measurement, not a reading.
  static void concentrationProbe() {
    SystemInterface fluid = new SystemPrEos(600.0, 5.0);
    fluid.addComponent("methane", 0.05);
    fluid.addComponent("oxygen", 0.10);
    fluid.addComponent("nitrogen", 0.85);
    fluid.setMixingRule(2);
    Stream stream = new Stream("probe", fluid);
    stream.setFlowRate(1.0, "mol/sec");
    stream.run();
    SystemInterface flashed = stream.getThermoSystem();
    double molarMass = flashed.getMolarMass("kg/mol");
    System.out.println("concentration_probe");
    System.out.println("molar_density_with_unit=" + flashed.getPhase(0).getDensity("mol/m3"));
    System.out.println("mass_density_with_unit=" + flashed.getPhase(0).getDensity("kg/m3"));
    System.out.println("mass_density_bare=" + flashed.getPhase(0).getDensity());
    System.out.println("molar_mass_kg_per_mol=" + molarMass);
    System.out.println("mass_over_molar=" + (flashed.getPhase(0).getDensity("kg/m3") / molarMass));

    // **The superficial velocity is the one intermediate a reader cannot guess.** The class
    // computes `getVolume("m3") / getNumberOfMoles() * getTotalNumberOfMoles() / totalArea`, and
    // whether the first call is extensive or per-mole is exactly what this prints.
    double area = Math.PI * 0.1 * 0.1 / 4.0;
    double volume = flashed.getVolume("m3");
    double numberOfMoles = flashed.getNumberOfMoles();
    double totalMoles = flashed.getTotalNumberOfMoles();
    System.out.println("velocity_probe");
    System.out.println("volume_m3=" + volume);
    System.out.println("number_of_moles=" + numberOfMoles);
    System.out.println("total_number_of_moles=" + totalMoles);
    System.out.println("tube_area_m2=" + area);
    System.out.println("superficial_velocity_m_s=" + (volume / numberOfMoles * totalMoles / area));
    System.out.println("gas_density_kg_m3=" + flashed.getDensity("kg/m3"));
    System.out.println("gas_phase_density_kg_m3=" + flashed.getPhase(0).getDensity("kg/m3"));
    System.out.println("gas_viscosity_pa_s=" + flashed.getPhase(0).getViscosity("kg/msec"));

    // **The energy balance's heat capacity flow is `getCp("J/molK") * totalMolFlow`**, and
    // azoth has no real-mixture molar Cp at all - only `eos.ideal_gas_cp` and a liquid
    // correlation. This prints the number the port would have to reproduce, beside the two
    // things a port could reach for instead: the ideal-gas Cp and a finite difference of the
    // molar enthalpy, which is what `Cp = (dH/dT)_P` says it is.
    System.out.println("cp_probe");
    System.out.println("system_cp_J_molK=" + flashed.getCp("J/molK"));
    System.out.println("phase_cp_J_molK=" + flashed.getPhase(0).getCp("J/molK"));
    System.out.println("phase_cp_molar_bare=" + flashed.getPhase(0).getCp());
    double h1 = enthalpyAt(flashAt(600.0, 5.0));
    double h2 = enthalpyAt(flashAt(600.001, 5.0));
    System.out.println("molar_enthalpy_600K=" + h1);
    System.out.println("finite_difference_dH_dT=" + ((h2 - h1) / 0.001));
    System.out.println();
  }

  /// A flashed single-point state, for the finite-difference Cp.
  static SystemInterface flashAt(double temperatureK, double pressureBara) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBara);
    fluid.addComponent("methane", 0.05);
    fluid.addComponent("oxygen", 0.10);
    fluid.addComponent("nitrogen", 0.85);
    fluid.setMixingRule(2);
    Stream stream = new Stream("fd", fluid);
    stream.setFlowRate(1.0, "mol/sec");
    stream.run();
    return stream.getThermoSystem();
  }

  /// The molar enthalpy of a flashed state, J/mol.
  static double enthalpyAt(SystemInterface flashed) {
    return flashed.getEnthalpy() / flashed.getTotalNumberOfMoles();
  }


  /// computes.
  static void flareRows() {
    String[] binary = new String[] { "methane", "n-butane" };
    double[] binaryZ = new double[] { 0.9, 0.1 };
    String[] sales = new String[] { "methane", "ethane", "propane", "n-butane", "nitrogen" };
    double[] salesZ = new double[] { 0.85, 0.07, 0.03, 0.02, 0.03 };
    String[] sour = new String[] { "methane", "CO2", "nitrogen" };
    double[] sourZ = new double[] { 0.8, 0.15, 0.05 };
    String[] heavy = new String[] { "methane", "n-butane", "n-heptane" };
    double[] heavyZ = new double[] { 0.6, 0.3, 0.1 };

    runFlare("binary", binary, binaryZ);
    runFlare("sales_gas", sales, salesZ);
    runFlare("sour_gas", sour, sourZ);
    runFlare("heavy", heavy, heavyZ);
    flarePerformanceRows(binary, binaryZ, sales, salesZ);
  }

  /// The radiation model, the dispersion surrogate and the capacity check.
  ///
  /// **The radiation model is a point source at the flame height**:
  /// `radiantFraction * duty / (4 pi (d^2 + h^2))`, and its inverse for a threshold. The
  /// surrogate is the tip's exit velocity and momentum flux from the inlet's own density;
  /// the capacity check divides the operating point by three optional design limits. **Each
  /// limit is stated as a multiple of the state it is checked against**, so the printed
  /// design value is the oracle the port is held to rather than a typed number.
  static void flarePerformanceRows(String[] binary, double[] binaryZ, String[] sales,
      double[] salesZ) {
    runFlarePerformance("performance_defaults", binary, binaryZ, 1.0, null, null, null, 30.0,
        4000.0, null, null, null);
    runFlarePerformance("performance_design_capacities", binary, binaryZ, 1.0, null, null, null,
        50.0, 1500.0, 0.8, 0.5, 0.5);
    runFlarePerformance("performance_hot_flame", sales, salesZ, 1.0, 60.0, 0.3, 0.6, 10.0, 2000.0,
        null, null, null);
    // **A threshold the horizontal distance can satisfy**: at 1 mol/s the 30 m flame height
    // puts every contour inside its own sphere, so the inverse answers zero. A hundredfold
    // flow and a 1 kW/m2 threshold put the contour outside it.
    runFlarePerformance("performance_high_duty", binary, binaryZ, 100.0, null, null, null, 60.0,
        1000.0, null, null, null);
    // The class's setters clamp rather than refuse: a negative height is zero, a fraction
    // above one is one, and a zero tip diameter is 1e-4 m.
    runFlarePerformance("performance_clamped_fields", binary, binaryZ, 1.0, -5.0, 5.0, 0.0, 30.0,
        4000.0, null, null, null);
  }

  static void runFlarePerformance(String label, String[] names, double[] z, double molPerSecond,
      Double flameHeight, Double radiantFraction, Double tipDiameter, double groundDistance,
      double fluxThreshold, Double heatFactor, Double massFactor, Double molarFactor) {
    Stream inlet = feed(names, z, 288.15, 1.01325, molPerSecond);
    neqsim.process.equipment.flare.Flare flare =
        new neqsim.process.equipment.flare.Flare("flare", inlet);
    if (flameHeight != null) {
      flare.setFlameHeight(flameHeight);
    }
    if (radiantFraction != null) {
      flare.setRadiantFraction(radiantFraction);
    }
    if (tipDiameter != null) {
      flare.setTipDiameter(tipDiameter);
    }
    flare.run();
    if (heatFactor != null) {
      flare.setDesignHeatDutyCapacity(heatFactor * flare.getHeatDuty(), "W");
    }
    if (massFactor != null) {
      flare.setDesignMassFlowCapacity(massFactor * inlet.getFlowRate("kg/sec"), "kg/sec");
    }
    if (molarFactor != null) {
      flare.setDesignMolarFlowCapacity(molarFactor * inlet.getFlowRate("mole/sec"), "mole/sec");
    }
    neqsim.process.equipment.flare.Flare.CapacityCheckResult capacity = flare.evaluateCapacity();
    neqsim.process.equipment.flare.dto.FlareDispersionSurrogateDTO dispersion =
        flare.getDispersionSurrogate();

    System.out.println(label);
    print("inlet", inlet);
    System.out.println("heat_duty_W=" + flare.getHeatDuty());
    System.out.println("co2_emission_kg_per_s=" + flare.getCO2Emission());
    System.out.println("ground_distance_m=" + groundDistance);
    System.out.println("radiation_flux_W_per_m2=" + flare.estimateRadiationHeatFlux(groundDistance));
    System.out.println("radiation_flux_threshold_W_per_m2=" + fluxThreshold);
    System.out.println("radiation_distance_m=" + flare.radiationDistanceForFlux(fluxThreshold));
    System.out.println("dispersion_exit_velocity_m_per_s=" + dispersion.getExitVelocityMs());
    System.out.println("dispersion_momentum_flux_Pa=" + dispersion.getMomentumFlux());
    System.out.println("dispersion_momentum_per_mass=" + dispersion.getMomentumPerMass());
    System.out.println(
        "dispersion_standard_volume_m3_per_s=" + dispersion.getStandardVolumeSm3PerSec());
    System.out.println("capacity_heat_duty_W=" + capacity.getHeatDutyW());
    System.out.println("capacity_design_heat_duty_W=" + capacity.getDesignHeatDutyW());
    System.out.println("capacity_heat_utilization=" + capacity.getHeatUtilization());
    System.out.println("capacity_mass_rate_kg_s=" + capacity.getMassRateKgS());
    System.out.println("capacity_design_mass_rate_kg_s=" + capacity.getDesignMassRateKgS());
    System.out.println("capacity_mass_utilization=" + capacity.getMassUtilization());
    System.out.println("capacity_molar_rate_mole_s=" + capacity.getMolarRateMoleS());
    System.out.println("capacity_design_molar_rate_mole_s=" + capacity.getDesignMolarRateMoleS());
    System.out.println("capacity_molar_utilization=" + capacity.getMolarUtilization());
    System.out.println("capacity_overloaded=" + capacity.isOverloaded());
    System.out.println();
  }

  static void runFlare(String label, String[] names, double[] z) {
    Stream inlet = feed(names, z, 288.15, 1.01325, 1.0);
    neqsim.process.equipment.flare.Flare flare =
        new neqsim.process.equipment.flare.Flare("flare", inlet);
    flare.run();
    System.out.println(label);
    print("inlet", inlet);
    printOrEmpty("product", flare.getOutletStream());
    System.out.println("heat_duty_W=" + flare.getHeatDuty());
    System.out.println("co2_emission_kg_per_s=" + flare.getCO2Emission());
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
    neqsim.process.mechanicaldesign.ejector.EjectorMechanicalDesign design = ejector.getMechanicalDesign();
    System.out.println("design_motive_nozzle_throat_area_m2=" + design.getMotiveNozzleThroatArea());
    System.out.println("design_motive_nozzle_diameter_m=" + design.getMotiveNozzleDiameter());
    System.out.println("design_suction_inlet_area_m2=" + design.getSuctionInletArea());
    System.out.println("design_suction_inlet_diameter_m=" + design.getSuctionInletDiameter());
    System.out.println("design_mixing_chamber_area_m2=" + design.getMixingChamberArea());
    System.out.println("design_mixing_chamber_diameter_m=" + design.getMixingChamberDiameter());
    System.out.println("design_diffuser_outlet_area_m2=" + design.getDiffuserOutletArea());
    System.out.println("design_diffuser_outlet_diameter_m=" + design.getDiffuserOutletDiameter());
    System.out.println("design_motive_nozzle_effective_length_m=" + design.getMotiveNozzleEffectiveLength());
    System.out.println("design_suction_inlet_length_m=" + design.getSuctionInletLength());
    System.out.println("design_mixing_chamber_length_m=" + design.getMixingChamberLength());
    System.out.println("design_diffuser_outlet_length_m=" + design.getDiffuserOutletLength());
    System.out.println("design_body_volume_m3=" + design.getBodyVolume());
    System.out.println("design_connected_piping_volume_m3=" + design.getConnectedPipingVolume());
    System.out.println("design_total_volume_m3=" + design.getTotalVolume());
    System.out.println("design_suction_connection_length_m=" + design.getSuctionConnectionLength());
    System.out.println("design_discharge_connection_length_m=" + design.getDischargeConnectionLength());
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
  /// The Gibbs reactor, and the reason this one prints a *trace* rather than an answer.
  ///
  /// **A Gibbs solve has a fixed point instead of a formula**, so the outlet composition
  /// alone cannot say whether a port reproduced the machine: two different solvers can agree
  /// to three figures on the answer and differ on every step that reached it. So each row
  /// prints the convergence trace - the iteration count, the final error, the element
  /// balance on the class's own seven element names, the Lagrange multipliers, the Gibbs
  /// energy history and the objective vector - and the outlet beside it.
  ///
  /// **Driven on PR, not SRK.** The class carries no cubic: it works on whatever fluid it is
  /// handed, and azoth's `Stream::mixture()` resolves PR with the databank's kij, so the
  /// oracle has to be the same fluid the port runs. This is the call the column probe makes
  /// for the same reason. The class's own JUnit fluids are SRK, which is their fluid and not
  /// this port's - `GibbsReactorTest.testComponentNotInDatabaseMolesUnchanged` is the
  /// methane/oxygen case below at `damping 0.01`, carried here on PR.
  ///
  /// **The phase and its root are printed because the port has to choose one.** The class
  /// reads `getPhase(0).getComponent(i).getFugacityCoefficient()` and never states which root
  /// that is, so the probe records the phase type, the compressibility and the fugacity
  /// coefficients it actually used instead of leaving a reader to infer it.
  static void gibbsReactorRows() {
    gibbsReactorRow("methane_oxygen_adiabatic",
        new String[] { "methane", "oxygen", "CO2", "CO", "water" },
        new double[] { 0.05, 0.5, 0.0, 0.0, 0.0 }, 298.15, 100.0, true, 0.01, 10000, 1e-3);
    gibbsReactorRow("methane_oxygen_isothermal",
        new String[] { "methane", "oxygen", "CO2", "CO", "water" },
        new double[] { 0.05, 0.5, 0.0, 0.0, 0.0 }, 298.15, 100.0, false, 0.05, 5000, 1e-3);
    gibbsReactorRow("ammonia_isothermal", new String[] { "hydrogen", "nitrogen", "ammonia" },
        new double[] { 1.5, 0.5, 0.0 }, 450.0, 300.0, false, 0.05, 5000, 1e-3);
    gibbsReactorRow("ammonia_adiabatic", new String[] { "hydrogen", "nitrogen", "ammonia" },
        new double[] { 1.5, 0.5, 0.0 }, 450.0, 300.0, true, 0.05, 5000, 1e-3);
    gibbsReactorRow("steam_methane_hot", new String[] { "methane", "water" },
        new double[] { 1.0, 1.0 }, 1000.0, 1.0, false, 0.05, 5000, 1e-3);
    // A feed whose species set is larger than the ones that react, so the variable and the
    // feed-excluded sets are both non-trivial and the argon row is exercised - `argon` is
    // the species whose atom count sits in the column the class calls `Ar` and the file
    // calls `Na`.
    // **A feed whose every species stays significant**, which is what a case needs: the
    // comparison between the two implementations is relative, so a species that settles at
    // `1e-15` in one and `1e-7` in the other cannot be asserted at all. The reverse water-gas
    // shift carries four species against three active elements - so nothing is over-determined -
    // and at these temperatures all four are tens of per cent.
    gibbsReactorRow("rwgs_1200k",
        new String[] { "CO2", "hydrogen", "CO", "water" },
        new double[] { 1.0, 1.0, 0.0, 0.0 }, 1200.0, 1.0, false, 0.05, 5000, 1e-3);
    gibbsReactorRow("rwgs_1000k",
        new String[] { "CO2", "hydrogen", "CO", "water" },
        new double[] { 1.0, 1.0, 0.0, 0.0 }, 1000.0, 1.0, false, 0.05, 5000, 1e-3);
    gibbsReactorRow("argon_in_feed", new String[] { "hydrogen", "oxygen", "water", "argon" },
        new double[] { 0.1, 1.0, 0.0, 0.05 }, 298.15, 50.0, true, 0.05, 5000, 1e-3);
  }

  /// **The step controls, on the states `GibbsReactorAlgorithmTest` solves.**
  ///
  /// `useArmijoLineSearch` and `useAdaptiveStepSize` are both `false` by default, so no row in
  /// `process_gibbs_reactor.tsv` reaches either. The class's own test drives them, and this
  /// capture prints the whole path rather than the assertion: the iteration count, the final
  /// error, the outlet and both diagnostic histories, so a port is held to the trajectory and
  /// not only to "it converged".
  ///
  /// **The fixed-damping row is the baseline the class's own test compares against**: the same
  /// feed, damping and gate as `armijo_methane_adiabatic`, with the line search off. The two
  /// must reach the same equilibrium, and printing both is what makes that checkable.
  static void gibbsAlgorithmRows() {
    gibbsAlgorithmRow("armijo_methane_adiabatic",
        new String[] { "methane", "oxygen", "CO2", "water" }, new double[] { 0.05, 0.5, 0.0, 0.0 },
        298.15, 100.0, true, 0.01, 10000, 1.0e-3, true, false);
    gibbsAlgorithmRow("fixed_methane_adiabatic",
        new String[] { "methane", "oxygen", "CO2", "water" }, new double[] { 0.05, 0.5, 0.0, 0.0 },
        298.15, 100.0, true, 0.01, 10000, 1.0e-3, false, false);
    gibbsAlgorithmRow("adaptive_hydrogen",
        new String[] { "hydrogen", "oxygen", "water" }, new double[] { 0.1, 1.0, 0.0 }, 298.15,
        1.0, false, 0.05, 5000, 1.0e-3, false, true);
    gibbsAlgorithmRow("armijo_adaptive_hydrogen",
        new String[] { "hydrogen", "oxygen", "water" }, new double[] { 0.1, 1.0, 0.0 }, 298.15,
        1.0, false, 0.05, 5000, 1.0e-3, true, true);
    // **An undamped step under the line search**, which is where the backtracks are: the class's
    // own Armijo row caps alpha at its `0.01` damping, so it accepts that alpha every pass.
    gibbsAlgorithmRow("armijo_ammonia_undamped",
        new String[] { "hydrogen", "nitrogen", "ammonia" }, new double[] { 1.5, 0.5, 0.0 }, 450.0,
        300.0, false, 1.0, 5000, 1.0e-3, true, false);
    // **Two rows on a state every species stays significant in**, which is what a case needs:
    // the reverse water-gas shift at 1000 K carries four species against three active elements,
    // so nothing sits on the `1e-6` floor and the comparison can be relative.
    gibbsAlgorithmRow("armijo_rwgs_undamped",
        new String[] { "CO2", "hydrogen", "CO", "water" }, new double[] { 1.0, 1.0, 0.0, 0.0 },
        1000.0, 1.0, false, 1.0, 5000, 1.0e-3, true, false);
    gibbsAlgorithmRow("adaptive_rwgs",
        new String[] { "CO2", "hydrogen", "CO", "water" }, new double[] { 1.0, 1.0, 0.0, 0.0 },
        1000.0, 1.0, false, 0.05, 5000, 1.0e-3, false, true);
    // **The two states `GibbsReactorAlgorithmTest` uses for the regularisation and the
    // diagnostic.** The flag is on in the first and off in the second; the condition history is
    // recorded on the first five passes either way.
    gibbsAlgorithmRow("regularized_methane",
        new String[] { "methane", "oxygen", "CO2", "CO", "water" },
        new double[] { 0.05, 0.5, 0.0, 0.0, 0.0 }, 298.15, 100.0, true, 0.01, 10000, 1.0e-3,
        false, false, true, 1.0e10, 1.0e-6);
    gibbsAlgorithmRow("diagnostic_ammonia",
        new String[] { "hydrogen", "nitrogen", "ammonia" }, new double[] { 1.5, 0.5, 0.0 }, 450.0,
        300.0, false, 0.05, 5000, 1.0e-3, false, false, false, 1.0e10, 1.0e-6);
  }

  static void gibbsAlgorithmRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, boolean adiabatic, double damping, int maxIterations, double tolerance,
      boolean armijo, boolean adaptive) {
    gibbsAlgorithmRow(label, names, z, temperatureK, pressureBara, adiabatic, damping,
        maxIterations, tolerance, armijo, adaptive, false, 1.0e10, 1.0e-6);
  }

  static void gibbsAlgorithmRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, boolean adiabatic, double damping, int maxIterations, double tolerance,
      boolean armijo, boolean adaptive, boolean regularization, double threshold, double tau) {
    Stream inlet = feed(names, z, temperatureK, pressureBara, 1.0);

    neqsim.process.equipment.reactor.GibbsReactor reactor =
        new neqsim.process.equipment.reactor.GibbsReactor("gibbs " + label, inlet);
    reactor.setUseAllDatabaseSpecies(false);
    reactor.setDampingComposition(damping);
    reactor.setMaxIterations(maxIterations);
    reactor.setConvergenceTolerance(tolerance);
    reactor.setUseArmijoLineSearch(armijo);
    reactor.setUseAdaptiveStepSize(adaptive);
    reactor.setUseRegularization(regularization);
    reactor.setRegularizationThreshold(threshold);
    reactor.setRegularizationTau(tau);
    reactor.setEnergyMode(adiabatic
        ? neqsim.process.equipment.reactor.GibbsReactor.EnergyMode.ADIABATIC
        : neqsim.process.equipment.reactor.GibbsReactor.EnergyMode.ISOTHERMAL);

    reactor.run();
    StreamInterface outlet = reactor.getOutletStream();

    System.out.println(label);
    print("feed", inlet);
    print("product", outlet);
    // **What the solve carried, beside what the outlet stream published.** The two are not the
    // same object in the class, and a row where they disagree says which one moved.
    StringBuilder moles = new StringBuilder("outlet_moles=");
    for (int i = 0; i < reactor.getOutletMole().size(); i++) {
      moles.append(inlet.getThermoSystem().getComponent(i).getComponentName()).append(":")
          .append(reactor.getOutletMole().get(i));
      if (i + 1 < reactor.getOutletMole().size()) {
        moles.append(" ");
      }
    }
    System.out.println(moles);
    System.out.println("energy_mode=" + reactor.getEnergyMode() + "\tdamping=" + damping
        + "\tmax_iterations=" + maxIterations + "\ttolerance=" + tolerance
        + "\tarmijo=" + armijo + "\tadaptive=" + adaptive);
    System.out.println("converged=" + reactor.hasConverged() + "\titerations="
        + reactor.getActualIterations() + "\tfinal_error=" + reactor.getFinalConvergenceError()
        + "\tmass_balance_converged=" + reactor.getMassBalanceConverged());
    System.out.println("gibbs_history=" + joined(reactor.getGibbsEnergyHistory()));
    System.out.println("step_size_history=" + joined(reactor.getStepSizeHistory()));
    System.out.println("regularization=" + regularization + "\tthreshold=" + threshold
        + "\ttau=" + tau);
    System.out.println("condition_number_history="
        + joined(reactor.getConditionNumberHistory()));
    System.out.println();
  }

  /// One history as the space-separated list the capture's other rows use.
  static String joined(java.util.List<Double> values) {
    StringBuilder text = new StringBuilder();
    for (int i = 0; i < values.size(); i++) {
      if (i > 0) {
        text.append(' ');
      }
      text.append(values.get(i));
    }
    return text.toString();
  }

  static void gibbsReactorRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, boolean adiabatic, double damping, int maxIterations, double tolerance) {
    Stream inlet = feed(names, z, temperatureK, pressureBara, 1.0);

    neqsim.process.equipment.reactor.GibbsReactor reactor =
        new neqsim.process.equipment.reactor.GibbsReactor("gibbs", inlet);
    reactor.setUseAllDatabaseSpecies(false);
    reactor.setDampingComposition(damping);
    reactor.setMaxIterations(maxIterations);
    reactor.setConvergenceTolerance(tolerance);
    reactor.setEnergyMode(adiabatic ? neqsim.process.equipment.reactor.GibbsReactor.EnergyMode.ADIABATIC
        : neqsim.process.equipment.reactor.GibbsReactor.EnergyMode.ISOTHERMAL);

    reactor.run();
    StreamInterface outlet = reactor.getOutletStream();

    System.out.println(label);
    print("feed", inlet);
    // **What the class actually took as its input.** `feed_n` above is the *stream's* flow
    // rate; `inlet_mole` is the list the solve reads, and the two are not the same number:
    // `GibbsReactor.run` fills it from `getComponent(i).getNumberOfMolesInPhase()`, which is
    // per phase and not the total. Printing it is the only way to know what the row was given.
    System.out.println("feed_fluid_total_moles=" + inlet.getThermoSystem().getTotalNumberOfMoles()
        + "\tfeed_phases=" + inlet.getThermoSystem().getNumberOfPhases());
    StringBuilder classInlet = new StringBuilder("class_inlet_mole=");
    for (int i = 0; i < reactor.getInletMole().size(); i++) {
      classInlet.append(inlet.getThermoSystem().getComponent(i).getComponentName()).append(":")
          .append(reactor.getInletMole().get(i));
      if (i + 1 < reactor.getInletMole().size()) {
        classInlet.append(" ");
      }
    }
    System.out.println(classInlet);
    // **The inlet's phase structure, which is the one the solve inherits.** `run` ends with
    // `getOutletStream().run(id)`, which flashes, so `outlet_phases` below describes the state
    // *after* the reactor - not the state the iterations saw. The solve's fluid is a clone of
    // the inlet's, and `init(3)` refreshes a phase's derivatives without ever adding one, so
    // the phase count that matters is this one.
    SystemInterface feedFluid = inlet.getThermoSystem();
    System.out.println("feed_phases=" + feedFluid.getNumberOfPhases() + "\tfeed_phase0_type="
        + feedFluid.getPhase(0).getPhaseTypeName() + "\tfeed_phase0_z=" + feedFluid.getPhase(0).getZ()
        + "\tfeed_phase0_moles=" + feedFluid.getPhase(0).getNumberOfMolesInPhase());
    printOrEmpty("product", outlet);
    System.out.println("energy_mode=" + reactor.getEnergyMode() + "\tdamping=" + damping
        + "\tmax_iterations=" + maxIterations + "\ttolerance=" + tolerance);
    System.out.println("converged=" + reactor.hasConverged() + "\titerations="
        + reactor.getActualIterations() + "\tfinal_error=" + reactor.getFinalConvergenceError()
        + "\tmass_balance_converged=" + reactor.getMassBalanceConverged());

    SystemInterface fluid = outlet.getThermoSystem();
    System.out.println("outlet_phases=" + fluid.getNumberOfPhases() + "\tphase0_type="
        + fluid.getPhase(0).getPhaseTypeName() + "\tphase0_z=" + fluid.getPhase(0).getZ()
        + "\tphase0_moles=" + fluid.getPhase(0).getNumberOfMolesInPhase());

    String[] elements = reactor.getElementNames();
    double[] balanceIn = reactor.getElementMoleBalanceIn();
    double[] balanceOut = reactor.getElementMoleBalanceOut();
    double[] balanceDiff = reactor.getElementMoleBalanceDiff();
    for (int i = 0; i < elements.length; i++) {
      System.out.println("element_" + elements[i] + "=" + balanceIn[i] + "\t" + balanceOut[i] + "\t"
          + balanceDiff[i]);
    }

    double[] lambda = reactor.getLagrangianMultipliers();
    StringBuilder multipliers = new StringBuilder("lambda=");
    for (int i = 0; i < lambda.length; i++) {
      multipliers.append(elements[i]).append(":").append(lambda[i]);
      if (i + 1 < lambda.length) {
        multipliers.append(" ");
      }
    }
    System.out.println(multipliers);

    // **The fugacity coefficients the objective actually read.** `getFugacityCoefficient(0)`
    // returns phase 0's array, and the class calls it once per objective evaluation without
    // ever stating which root or which state - so it is recorded rather than inferred.
    double[] phi = reactor.getFugacityCoefficient(0);
    StringBuilder fugacities = new StringBuilder("phi=");
    for (int i = 0; i < phi.length; i++) {
      fugacities.append(fluid.getPhase(0).getComponent(i).getName()).append(":").append(phi[i]);
      if (i + 1 < phi.length) {
        fugacities.append(" ");
      }
    }
    System.out.println(fugacities);

    StringBuilder outletMoles = new StringBuilder("outlet_moles=");
    for (int i = 0; i < reactor.getOutletMole().size(); i++) {
      outletMoles.append(fluid.getComponent(i).getComponentName()).append(":")
          .append(reactor.getOutletMole().get(i));
      if (i + 1 < reactor.getOutletMole().size()) {
        outletMoles.append(" ");
      }
    }
    System.out.println(outletMoles);

    StringBuilder objective = new StringBuilder("objective=");
    for (int i = 0; i < fluid.getNumberOfComponents(); i++) {
      String name = fluid.getComponent(i).getComponentName();
      objective.append(name).append(":").append(reactor.getObjectiveFunctionValues().get(name));
      if (i + 1 < fluid.getNumberOfComponents()) {
        objective.append(" ");
      }
    }
    System.out.println(objective);

    // The history is per iteration and is the only record of the path, so it is printed
    // whole rather than summarised - a hundred and one values is what a port has to match.
    java.util.List<Double> history = reactor.getGibbsEnergyHistory();
    StringBuilder gibbs = new StringBuilder("gibbs_history=");
    for (int i = 0; i < history.size(); i++) {
      gibbs.append(history.get(i));
      if (i + 1 < history.size()) {
        gibbs.append(" ");
      }
    }
    System.out.println(gibbs);
    System.out.println();
  }

  /// **One iteration of the Gibbs solve, so the map can be compared rather than the trajectory.**
  ///
  /// The six rows above show that the port tracks NeqSim to `1e-14` kJ/mol absolute and then
  /// separates, which says the disagreement is in the iteration rather than in the equations -
  /// and says nothing about *where*. `maxIterations = n` is the instrument: the loop's
  /// convergence branch is `(deltaXNorm < tol && iteration >= minIterations) || iteration ==
  /// maxIterations`, so at `n = 2` exactly one update is applied and the second pass stops. The
  /// state, the multipliers, the objective and the Jacobian that come out are one step from the
  /// same place the port is one step from.
  static void gibbsReactorSteps() {
    for (int iterations : new int[] { 2, 3 }) {
      gibbsReactorStepRow("ammonia_steps_" + iterations, iterations);
    }
  }

  static void gibbsReactorStepRow(String label, int maxIterations) {
    String[] names = new String[] { "hydrogen", "nitrogen", "ammonia" };
    Stream inlet = feed(names, new double[] { 1.5, 0.5, 0.0 }, 450.0, 300.0, 1.0);

    neqsim.process.equipment.reactor.GibbsReactor reactor =
        new neqsim.process.equipment.reactor.GibbsReactor("gibbs", inlet);
    reactor.setUseAllDatabaseSpecies(false);
    reactor.setMaxIterations(maxIterations);
    reactor.setConvergenceTolerance(1e-3);
    reactor.setEnergyMode(
        neqsim.process.equipment.reactor.GibbsReactor.EnergyMode.ISOTHERMAL);
    reactor.run();

    System.out.println(label);
    System.out.println("iterations=" + reactor.getActualIterations() + "\tconverged="
        + reactor.hasConverged() + "\tfinal_error=" + reactor.getFinalConvergenceError());

    SystemInterface fluid = reactor.getOutletStream().getThermoSystem();
    StringBuilder moles = new StringBuilder("outlet_moles=");
    for (int i = 0; i < reactor.getOutletMole().size(); i++) {
      moles.append(fluid.getComponent(i).getComponentName()).append(":")
          .append(reactor.getOutletMole().get(i));
      if (i + 1 < reactor.getOutletMole().size()) {
        moles.append(" ");
      }
    }
    System.out.println(moles);

    StringBuilder objective = new StringBuilder("objective=");
    for (int i = 0; i < fluid.getNumberOfComponents(); i++) {
      String name = fluid.getComponent(i).getComponentName();
      objective.append(name).append(":").append(reactor.getObjectiveFunctionValues().get(name));
      if (i + 1 < fluid.getNumberOfComponents()) {
        objective.append(" ");
      }
    }
    System.out.println(objective);

    String[] elements = reactor.getElementNames();
    double[] lambda = reactor.getLagrangianMultipliers();
    StringBuilder multipliers = new StringBuilder("lambda=");
    for (int i = 0; i < lambda.length; i++) {
      multipliers.append(elements[i]).append(":").append(lambda[i]);
      if (i + 1 < lambda.length) {
        multipliers.append(" ");
      }
    }
    System.out.println(multipliers);

    // **The Jacobian as it stands, which after `maxIterations` passes is the one the *last*
    // pass computed.** Its rows are the species variables then the active element balances, and
    // its column labels say which, so a reader can line it up without guessing.
    double[][] jacobian = reactor.getJacobianMatrix();
    java.util.List<String> rows = reactor.getJacobianRowLabels();
    java.util.List<String> cols = reactor.getJacobianColLabels();
    System.out.println("jacobian_rows=" + rows);
    System.out.println("jacobian_cols=" + cols);
    for (int i = 0; i < jacobian.length; i++) {
      StringBuilder row = new StringBuilder("jacobian_" + i + "=");
      for (int j = 0; j < jacobian[i].length; j++) {
        row.append(jacobian[i][j]);
        if (j + 1 < jacobian[i].length) {
          row.append(" ");
        }
      }
      System.out.println(row);
    }
    System.out.println();
  }

  /// **The one thing `ProcessSystem` is required for.** Every other row in this file drives a
  /// single unit or a stream, because `SimulationInterface.run()` defaults to
  /// `run(UUID.randomUUID())` and a unit runs standalone. A **recycle** needs the sequential loop
  /// around it - the tear has to be evaluated against the previous pass - so this is the first
  /// and only row that builds a `ProcessSystem`.
  ///
  /// **It is `specs/flowsheets/demo.toml`'s graph**: feed, mixer, pump, heater, separator, with the
  /// separator's liquid returned to the mixer. The wiring is NeqSim's own idiom - the recycle's
  /// inlet is the produced stream and the mixer takes the recycle's *outlet* - and the recycle is
  /// added last, which is where an insertion-order run puts it.
  ///
  /// **The converged state is the answer and the residuals are what make the convergence
  /// checkable.** A tear that lands in the same place by a different number of passes is a
  /// different machine, so the row prints both: every unit's outlet record, and the recycle's
  /// iteration count with its four residuals.
  static void flowsheetRows() {
    flowsheetRow("demo", 320.0);
    // **And the same graph with the heater cold enough to condense.** At 320 K and 20 bar this
    // feed is all vapour, so the separator's liquid is nil and the tear carries nothing - the
    // recycle converges on the zero-flow floor rather than on a comparison. At 250 K the tear
    // carries material and the loop has something to converge, which is the row worth having.
    flowsheetRow("demo_condensing", 250.0);
  }

  /**
   * The same graph with each acceleration method, at a tolerance the loop does not close at.
   *
   * **This row exists because the acceleration changes the answer, and that had to be measured
   * rather than assumed.** At 250 K the heated feed condenses and the tear carries material, and
   * with `flowTolerance` at `1e-8` direct substitution never closes - the recycled liquid has no
   * exit and grows by a fixed amount every pass. Broyden *does* close it, at a state with a
   * hundredfold smaller liquid, which is either the class's answer or a defect in it.
   */
  static void flowsheetAcceleratedRows() {
    flowsheetRow("demo_condensing_direct", 250.0, null, 0.0);
    flowsheetRow("demo_condensing_wegstein", 250.0, AccelerationMethod.WEGSTEIN, 1e-8);
    flowsheetRow("demo_condensing_broyden", 250.0, AccelerationMethod.BROYDEN, 1e-8);
  }

  static void flowsheetRow(String label, double heaterK) {
    flowsheetRow(label, heaterK, null, 0.0);
  }

  static void flowsheetRow(String label, double heaterK, AccelerationMethod acceleration,
      double flowTolerance) {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.9, 0.1 };

    SystemInterface fluid = new SystemPrEos(300.0, 5.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream feed = new Stream("feed_1", fluid);
    feed.setFlowRate(1.0, "mol/sec");
    feed.run();

    neqsim.process.equipment.mixer.Mixer mix1 = new neqsim.process.equipment.mixer.Mixer("mix1");
    mix1.addStream(feed);

    neqsim.process.equipment.pump.Pump p1 = new neqsim.process.equipment.pump.Pump("p1", mix1.getOutletStream());
    p1.setOutletPressure(20.0, "bara");
    p1.setIsentropicEfficiency(0.75);

    neqsim.process.equipment.heatexchanger.Heater hx1 = new neqsim.process.equipment.heatexchanger.Heater("hx1", p1.getOutletStream());
    hx1.setOutletTemperature(heaterK, "K");

    neqsim.process.equipment.separator.Separator sep1 = new neqsim.process.equipment.separator.Separator("sep1", hx1.getOutletStream());
    // The class takes the drop in bara and the entrainment as a fraction of the vapour's
    // moles carried into the liquid, on the `feed` basis the palette entry declares.
    sep1.setPressureDrop(0.5);
    sep1.setEntrainment(0.0, "mole", "feed", "gas", "oil");

    // The tear: the separator's liquid returns to the mixer.
    //
    // **The outlet stream is made here, and that is not decoration.** `Recycle.outletStream`
    // starts as `null` (`Recycle.java:51`) and is only assigned inside `run` (`:1246`), so the
    // idiom NeqSim's own example uses - `mixer.addStream(recycle.getOutletStream())` at
    // `simpleTopSideProcess2.java:47` - wires a **null**, and `Mixer.getCanonicalStreamIndices`
    // then dereferences it. Making the stream first is what `ProcessRecipesDocumentationTest:121`
    // does (`recycle.setOutletStream(recycleTear)`).
    neqsim.process.equipment.util.Recycle recycle1 = new neqsim.process.equipment.util.Recycle("recycle_1");
    recycle1.addStream(sep1.getLiquidOutStream());
    neqsim.process.equipment.stream.Stream tear =
        new neqsim.process.equipment.stream.Stream("recycle_1_tear", fluid.clone());
    tear.setFlowRate(0.0, "mol/sec");
    tear.run();
    recycle1.setOutletStream(tear);
    mix1.addStream(tear);

    neqsim.process.processmodel.ProcessSystem operations = new neqsim.process.processmodel.ProcessSystem();
    operations.add(feed);
    operations.add(mix1);
    operations.add(p1);
    operations.add(hx1);
    operations.add(sep1);
    operations.add(recycle1);
    // **`runSequential` is not what `run()` does by default.** `useOptimizedExecution` is `true`
    // (`ProcessSystem.java:273`), so `run(id)` takes `runOptimized` and the legacy sequential path
    // is reached only with this flag off - which is the path azoth's executor is a port of, and so
    // the one an oracle for it has to drive.
    if (acceleration != null) {
      recycle1.setAccelerationMethod(acceleration);
    }
    if (flowTolerance > 0.0) {
      recycle1.setFlowTolerance(flowTolerance);
    }
    operations.setUseOptimizedExecution(false);
    operations.run();

    System.out.println(label);
    System.out.println("heater_k=" + heaterK);
    System.out.println("execution=runSequential");
    print("feed_1", feed);
    print("mix1_product", mix1.getOutletStream());
    print("p1_outlet", p1.getOutletStream());
    print("hx1_outlet", hx1.getOutletStream());
    print("sep1_liquid", sep1.getLiquidOutStream());
    printOrEmpty("sep1_vapour", sep1.getGasOutStream());
    // **`isActive` is what makes the low-flow cutoff measurable rather than inferred.** A tear
    // whose inlet carries less than `minimumFlow` kg/hr is switched off by `deactivateOnLowFlow`,
    // and `runUnitProfiled` then skips it - so the outer loop runs a second pass that does not
    // advance `getIterations()`. Without this line the `demo` row's `recycle_iterations=1` reads
    // as a loop that stopped early rather than one that was switched off.
    System.out.println("recycle_iterations=" + recycle1.getIterations());
    System.out.println("recycle_active=" + recycle1.isActive());
    System.out.println("recycle_active=" + recycle1.isActive());
    System.out.println("recycle_solved=" + recycle1.solved());
    System.out.println("recycle_error_flow=" + recycle1.getErrorFlow());
    System.out.println("recycle_error_composition=" + recycle1.getErrorComposition());
    System.out.println("recycle_error_temperature=" + recycle1.getErrorTemperature());
    System.out.println("recycle_error_pressure=" + recycle1.getErrorPressure());
    System.out.println();
  }

  /// **The rate-based column, on two cubics.** `RateBasedPackedColumnTest` runs SRK, and this
  /// library's process layer resolves PR - so the class's own state is run on both, and the
  /// pair measures what the cubic moved. The PR rows are the ones the port is held to. The two
  /// TEG rows are `SystemSrkCPAstatoil` with mixing rule 10, which PR cannot represent, so they
  /// are evidence for a refusal rather than an oracle for a state.
  static void rateBasedRows() {
    rateBasedCo2WaterRow("co2_water_absorber", 0.10, 0.0, 6.0, "Pall-Ring-50", true, false);
    rateBasedCo2WaterRow("co2_water_absorber_srk", 0.10, 0.0, 6.0, "Pall-Ring-50", false, false);
    rateBasedCo2WaterRow("co2_water_stripper", 0.0, 0.04, 6.0, "Pall-Ring-50", true, false);
    rateBasedCo2WaterRow("co2_water_short", 0.10, 0.0, 1.0, "Pall-Ring-50", true, false);
    rateBasedCo2WaterRow("co2_water_tall", 0.10, 0.0, 8.0, "Pall-Ring-50", true, false);
    rateBasedCo2WaterRow("co2_water_mellapak", 0.10, 0.0, 6.0, "Mellapak-250Y", true, false);
    rateBasedCo2WaterRow("heat_transfer_disabled", 0.10, 0.0, 6.0, "Pall-Ring-50", true, true);
    rateBasedCo2WaterRow("zero_packed_height", 0.10, 0.0, 0.0, "Pall-Ring-50", true, false);
    rateBasedTegRow("teg_dehydration", 2000.0, 8.0, 5, 5.0);
    rateBasedTegRow("teg_circulation", 50.0, 12.0, 6, 20.0);
  }

  /// **`BILLET_SCHULTES_1999` in one capture, beside its own `ONDA_1968` baseline.** The class's
  /// second "correlation" is not one: it multiplies the two film coefficients by constants taken
  /// off the packing row, so one state run twice is the whole measurement. **Two rows of packings
  /// are needed since upstream PR #4405**: on the built-in `Pall-Ring-50` the constants are `1.0`
  /// and `0.40`, so the two correlations are the identity, while `Pall-Ring-50-plastic` carries
  /// `cp = 0.698` and `ch = 2.725` - the pair is `1.1354166` on `kGa` and `0.698` on `kLa`,
  /// neither floored at the class's `0.1`. The plastic pair is the state the pin's `Pall-Ring-50`
  /// resolved to.
  static void rateBasedBilletRows() {
    rateBasedCo2WaterRow("billet_onda", 0.10, 0.0, 6.0, "Pall-Ring-50", true, false,
        neqsim.process.equipment.distillation.RateBasedPackedColumn.MassTransferCorrelation.ONDA_1968);
    rateBasedCo2WaterRow("billet_schultes", 0.10, 0.0, 6.0, "Pall-Ring-50", true, false,
        neqsim.process.equipment.distillation.RateBasedPackedColumn.MassTransferCorrelation.BILLET_SCHULTES_1999);
    rateBasedCo2WaterRow("plastic_onda", 0.10, 0.0, 6.0, "Pall-Ring-50-plastic", true, false,
        neqsim.process.equipment.distillation.RateBasedPackedColumn.MassTransferCorrelation.ONDA_1968);
    rateBasedCo2WaterRow("plastic_schultes", 0.10, 0.0, 6.0, "Pall-Ring-50-plastic", true, false,
        neqsim.process.equipment.distillation.RateBasedPackedColumn.MassTransferCorrelation.BILLET_SCHULTES_1999);
  }

  /// **The two solvers this port does not carry, run on the state it does.**
  ///
  /// `billet_onda`'s own state - the PR row `process_rate_based_packed_column.tsv` holds and the
  /// port reproduces to `1e-4` - is run four ways: this port's pair as the control, then each of
  /// `SegmentSolver.SIMULTANEOUS_RESIDUAL` and `ColumnSolver.EQUATION_ORIENTED` alone, then both.
  /// **What decides whether either is a port or a non-port is not the class's source but this
  /// capture**: `EQUATION_ORIENTED` seeds from the fixed-point profile and then calls
  /// `acceptSolution` unconditionally, so it never reports a failure - the readable thing is
  /// `lastColumnResidualNorm` against `columnResidualTolerance`. A branch that lands above its own
  /// tolerance is a number with no oracle, and `SIMULTANEOUS_RESIDUAL` carries `@Disabled` on its
  /// own test upstream.
  static void rateBasedSolverRows() {
    neqsim.process.equipment.distillation.RateBasedPackedColumn.SegmentSolver sequential =
        neqsim.process.equipment.distillation.RateBasedPackedColumn.SegmentSolver.SEQUENTIAL_EXPLICIT;
    neqsim.process.equipment.distillation.RateBasedPackedColumn.SegmentSolver simultaneous =
        neqsim.process.equipment.distillation.RateBasedPackedColumn.SegmentSolver.SIMULTANEOUS_RESIDUAL;
    neqsim.process.equipment.distillation.RateBasedPackedColumn.ColumnSolver fixedPoint =
        neqsim.process.equipment.distillation.RateBasedPackedColumn.ColumnSolver.FIXED_POINT_PROFILE;
    neqsim.process.equipment.distillation.RateBasedPackedColumn.ColumnSolver equationOriented =
        neqsim.process.equipment.distillation.RateBasedPackedColumn.ColumnSolver.EQUATION_ORIENTED;
    rateBasedSolverRow("sequential_fixed_point", sequential, fixedPoint);
    rateBasedSolverRow("simultaneous_fixed_point", simultaneous, fixedPoint);
    rateBasedSolverRow("sequential_equation_oriented", sequential, equationOriented);
    rateBasedSolverRow("simultaneous_equation_oriented", simultaneous, equationOriented);
    // **The class's own equation-oriented test state**, which is `configuredColumn` at three
    // segments with the Newton capped and the homotopy removed - `testEquationOriented-
    // ColumnSolverSolvesColumnWideResiduals`' own four settings. It is here because that test
    // asserts a direction, a component balance and *finiteness*, and never that the norm met the
    // tolerance: if it does not converge at its own settings either, the branch is one no upstream
    // test holds to a residual.
    rateBasedSolverRow("equation_oriented_test_settings", sequential, equationOriented, 3, 2, 1,
        1.0e-5);
  }

  /// `billet_onda`'s state - `configuredColumn`'s, re-cased on PR - under one solver pair.
  static void rateBasedSolverRow(String label,
      neqsim.process.equipment.distillation.RateBasedPackedColumn.SegmentSolver segmentSolver,
      neqsim.process.equipment.distillation.RateBasedPackedColumn.ColumnSolver columnSolver) {
    rateBasedSolverRow(label, segmentSolver, columnSolver, 4, 8, 3, 1.0e-6);
  }

  /// The same row with the column-wide Newton's own knobs stated rather than left at the class's
  /// defaults - `maxColumnResidualIterations`, `columnHomotopySteps` and `columnResidualTolerance`.
  static void rateBasedSolverRow(String label,
      neqsim.process.equipment.distillation.RateBasedPackedColumn.SegmentSolver segmentSolver,
      neqsim.process.equipment.distillation.RateBasedPackedColumn.ColumnSolver columnSolver,
      int segments, int maxColumnResidualIterations, int homotopySteps,
      double columnResidualTolerance) {
    Stream gas = rateFeed("gas in", true, 313.15, 50.0, 1000.0, new String[] { "methane", "CO2" },
        new double[] { 0.9, 0.1 });
    Stream liquid = rateFeed("lean liquid", true, 303.15, 50.0, 2000.0,
        new String[] { "water", "CO2" }, new double[] { 1.0, 0.0 });
    neqsim.process.equipment.distillation.RateBasedPackedColumn column =
        new neqsim.process.equipment.distillation.RateBasedPackedColumn("rate based column", gas,
            liquid);
    column.setColumnDiameter(1.0);
    column.setPackedHeight(6.0);
    column.setNumberOfSegments(segments);
    column.setMaxIterations(20);
    column.setPackingType("Pall-Ring-50");
    column.setTransferComponents("CO2");
    column.setMassTransferCorrectionFactor(3.0);
    column.setConvergenceTolerance(1.0e-9);
    column.setMaxColumnResidualIterations(maxColumnResidualIterations);
    column.setColumnHomotopySteps(homotopySteps);
    column.setColumnResidualTolerance(columnResidualTolerance);
    column.setSegmentSolver(segmentSolver);
    column.setColumnSolver(columnSolver);
    String exception = "";
    try {
      column.run();
    } catch (RuntimeException ex) {
      exception = ex.getClass().getName() + ": " + ex.getMessage();
    }
    System.out.println(label);
    System.out.println("segment_solver=" + segmentSolver.name());
    System.out.println("column_solver=" + columnSolver.name());
    System.out.println("max_iterations=" + column.getMaxIterations());
    System.out.println("convergence_tolerance_mol_per_s=" + column.getConvergenceTolerance());
    System.out.println("max_column_residual_iterations=" + column.getMaxColumnResidualIterations());
    System.out.println("column_homotopy_steps=" + column.getColumnHomotopySteps());
    System.out.println("column_residual_tolerance=" + column.getColumnResidualTolerance());
    System.out.println("exception=" + exception);
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("convergence_residual=" + column.getLastConvergenceResidual());
    System.out.println("column_residual_norm=" + column.getLastColumnResidualNorm());
    System.out.println("column_residual_iterations=" + column.getLastColumnResidualIterations());
    System.out.println("gas_component_balance_residual="
        + column.getLastGasComponentBalanceResidual());
    System.out.println("liquid_component_balance_residual="
        + column.getLastLiquidComponentBalanceResidual());
    System.out.println("column_energy_balance_residual="
        + column.getLastColumnEnergyBalanceResidual());
    System.out.println("total_absolute_molar_transfer_mol_per_s="
        + column.getTotalAbsoluteMolarTransfer());
    print("gas_out", column.getGasOutStream());
    print("liquid_out", column.getLiquidOutStream());
    System.out.println("segment_count=" + column.getSegmentResults().size());
    for (neqsim.process.equipment.distillation.RateBasedPackedColumn.SegmentResult segment : column
        .getSegmentResults()) {
      String prefix = "segment" + segment.getSegmentNumber() + "_";
      System.out.println(prefix + "segment_solver=" + segment.getSegmentSolver());
      System.out.println(prefix + "residual_iterations=" + segment.getResidualIterations());
      System.out.println(prefix + "max_flux_residual_mol_per_s="
          + segment.getMaxFluxResidualMolPerSec());
      System.out.println(
          prefix + "heat_balance_residual_W=" + segment.getHeatBalanceResidualW());
      System.out.println(prefix + "net_molar_transfer_mol_per_s=" + segment.getNetMolarTransfer());
    }
    System.out.println();
  }

  /// `configuredColumn`'s own state and settings - the four segments, the 3.0 correction, the
  /// 1e-9 tolerance and `CO2` as the only transfer component - on the cubic the row names.
  static void rateBasedCo2WaterRow(String label, double gasCo2, double liquidCo2, double height,
      String packing, boolean pr, boolean heatOff) {
    rateBasedCo2WaterRow(label, gasCo2, liquidCo2, height, packing, pr, heatOff,
        neqsim.process.equipment.distillation.RateBasedPackedColumn.MassTransferCorrelation.ONDA_1968);
  }

  /// The same state with the mass-transfer correlation stated rather than left at the class's
  /// `ONDA_1968` field initialiser, which is why the rows above are unchanged by this overload.
  static void rateBasedCo2WaterRow(String label, double gasCo2, double liquidCo2, double height,
      String packing, boolean pr, boolean heatOff,
      neqsim.process.equipment.distillation.RateBasedPackedColumn.MassTransferCorrelation
          correlation) {
    Stream gas = rateFeed("gas in", pr, 313.15, 50.0, 1000.0, new String[] { "methane", "CO2" },
        new double[] { 1.0 - gasCo2, gasCo2 });
    Stream liquid = rateFeed("lean liquid", pr, 303.15, 50.0, 2000.0,
        new String[] { "water", "CO2" }, new double[] { 1.0 - liquidCo2, liquidCo2 });
    neqsim.process.equipment.distillation.RateBasedPackedColumn column =
        new neqsim.process.equipment.distillation.RateBasedPackedColumn("rate based column", gas,
            liquid);
    column.setColumnDiameter(1.0);
    column.setPackedHeight(height);
    column.setNumberOfSegments(4);
    // 60 rather than the class's 20: the packing fix moves the stripper to 21 passes.
    column.setMaxIterations(60);
    column.setPackingType(packing);
    column.setTransferComponents("CO2");
    column.setMassTransferCorrectionFactor(3.0);
    column.setConvergenceTolerance(1.0e-9);
    column.setMassTransferCorrelation(correlation);
    if (heatOff) {
      column.setHeatTransferModel(
          neqsim.process.equipment.distillation.RateBasedPackedColumn.HeatTransferModel.NONE);
    }
    column.run();
    printRateBased(label, column, gas, liquid, packing);
  }

  /// The class's own TEG contactor states, which are CPA-Statoil and are the refusal's evidence.
  static void rateBasedTegRow(String label, double tegKgPerHour, double height, int segments,
      double correction) {
    SystemInterface wet = new neqsim.thermo.system.SystemSrkCPAstatoil(313.15, 70.0);
    wet.addComponent("methane", 0.94);
    wet.addComponent("ethane", 0.04);
    wet.addComponent("propane", 0.018);
    wet.addComponent("water", 0.002);
    wet.createDatabase(true);
    wet.setMixingRule(10);
    Stream wetGas = new Stream("wet natural gas", wet);
    wetGas.setFlowRate(1000.0, "kg/hr");
    wetGas.run();

    SystemInterface teg = new neqsim.thermo.system.SystemSrkCPAstatoil(303.15, 70.0);
    teg.addComponent("TEG", 0.9999);
    teg.addComponent("water", 0.0001);
    teg.createDatabase(true);
    teg.setMixingRule(10);
    Stream leanTeg = new Stream("lean TEG", teg);
    leanTeg.setFlowRate(tegKgPerHour, "kg/hr");
    leanTeg.run();

    neqsim.process.equipment.distillation.RateBasedPackedColumn column =
        new neqsim.process.equipment.distillation.RateBasedPackedColumn("TEG contactor", wetGas,
            leanTeg);
    column.setColumnDiameter(1.2);
    column.setPackedHeight(height);
    column.setNumberOfSegments(segments);
    // 60 rather than the class's 25: the packing fix moves this profile to 33 passes.
    column.setMaxIterations(60);
    column.setPackingType("Mellapak-250Y");
    column.setTransferComponents("water");
    column.setMassTransferCorrectionFactor(correction);
    column.setConvergenceTolerance(1.0e-9);
    column.run();
    printRateBased(label, column, wetGas, leanTeg, "Mellapak-250Y");
  }

  /// **What `averageDiffusivity` and `binaryDiffusivity` actually read**, on the state the class
  /// flashes: the phase's effective diffusion vector and its binary matrix. The segment model
  /// falls back to `DEFAULT_GAS_DIFFUSIVITY`/`DEFAULT_LIQUID_DIFFUSIVITY` when the effective
  /// vector carries nothing positive, so whether the class computes a real one decides whether
  /// the whole transport runs on the correlations or on the fallbacks.
  static void printDiffusion(String port, StreamInterface stream) {
    SystemInterface fluid = stream.getThermoSystem().clone();
    try {
      new neqsim.thermodynamicoperations.ThermodynamicOperations(fluid).TPflash();
    } catch (RuntimeException ex) {
      fluid.init(3);
    }
    fluid.initProperties();
    for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
      StringBuilder effective = new StringBuilder(port + "_phase" + p + "_effective_diffusivity=");
      for (int i = 0; i < fluid.getPhase(p).getNumberOfComponents(); i++) {
        effective.append(fluid.getPhase(p).getPhysicalProperties().getEffectiveDiffusionCoefficient(i))
            .append(" ");
      }
      System.out.println(effective.toString().trim());
      for (int i = 0; i < fluid.getPhase(p).getNumberOfComponents(); i++) {
        for (int j = 0; j < fluid.getPhase(p).getNumberOfComponents(); j++) {
          if (i != j) {
            System.out.println(port + "_phase" + p + "_diffusivity_" + i + "_" + j + "="
                + fluid.getPhase(p).getPhysicalProperties().getDiffusionCoefficient(i, j));
          }
        }
      }
    }
  }

  /// A feed on the cubic a rate-based row names: SRK, which the class's own tests use, or the PR
  /// the probe's other feeds use.
  static Stream rateFeed(String name, boolean pr, double temperatureK, double pressureBara,
      double kgPerHour, String[] names, double[] z) {
    SystemInterface fluid = pr ? new SystemPrEos(temperatureK, pressureBara)
        : new neqsim.thermo.system.SystemSrkEos(temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream stream = new Stream(name, fluid);
    stream.setFlowRate(kgPerHour, "kg/hr");
    stream.run();
    return stream;
  }

  /// One rate-based row: the configuration, the four ports, the convergence scalar and the
  /// per-segment profile the port's layer diff reads.
  static void printRateBased(String label,
      neqsim.process.equipment.distillation.RateBasedPackedColumn column, Stream gasIn,
      Stream liquidIn, String packing) {
    System.out.println(label);
    System.out.println("packing=" + packing);
    System.out.println("column_diameter_m=" + column.getColumnDiameter());
    System.out.println("packed_height_m=" + column.getPackedHeight());
    System.out.println("number_of_segments=" + column.getNumberOfSegments());
    System.out.println("max_iterations=" + column.getMaxIterations());
    System.out.println("convergence_tolerance_mol_per_s=" + column.getConvergenceTolerance());
    System.out.println("mass_transfer_correction_factor=" + column.getMassTransferCorrectionFactor());
    System.out.println("mass_transfer_correlation=" + column.getMassTransferCorrelation().name());
    System.out.println("film_model=" + column.getFilmModel().name());
    System.out.println("heat_transfer_model=" + column.getHeatTransferModel().name());
    System.out.println("segment_solver=" + column.getSegmentSolver().name());
    System.out.println("column_solver=" + column.getColumnSolver().name());
    print("gas_in", gasIn);
    print("liquid_in", liquidIn);
    printDiffusion("gas_in", gasIn);
    printDiffusion("liquid_in", liquidIn);
    print("gas_out", column.getGasOutStream());
    print("liquid_out", column.getLiquidOutStream());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("convergence_residual_mol_per_s=" + column.getLastConvergenceResidual());
    System.out.println("total_absolute_molar_transfer_mol_per_s="
        + column.getTotalAbsoluteMolarTransfer());
    for (java.util.Map.Entry<String, Double> entry : column.getComponentTransferTotals().entrySet()) {
      System.out.println("transfer_total_" + entry.getKey() + "_mol_per_s=" + entry.getValue());
    }
    System.out.println("segment_count=" + column.getSegmentResults().size());
    for (neqsim.process.equipment.distillation.RateBasedPackedColumn.SegmentResult segment : column
        .getSegmentResults()) {
      String prefix = "segment" + segment.getSegmentNumber() + "_";
      System.out.println(prefix + "height_from_bottom_m=" + segment.getHeightFromBottom());
      System.out.println(prefix + "gas_temperature_K=" + segment.getGasTemperatureK());
      System.out.println(prefix + "liquid_temperature_K=" + segment.getLiquidTemperatureK());
      System.out.println(prefix + "gas_pressure_bara=" + segment.getGasPressureBar());
      System.out.println(prefix + "liquid_pressure_bara=" + segment.getLiquidPressureBar());
      System.out.println(prefix + "gas_molar_flow_mol_per_s=" + segment.getGasMolarFlow());
      System.out.println(prefix + "liquid_molar_flow_mol_per_s=" + segment.getLiquidMolarFlow());
      System.out.println(prefix + "gas_density_kg_per_m3=" + segment.getGasDensity());
      System.out.println(prefix + "liquid_density_kg_per_m3=" + segment.getLiquidDensity());
      System.out.println(prefix + "gas_viscosity_Pa_s=" + segment.getGasViscosity());
      System.out.println(prefix + "liquid_viscosity_Pa_s=" + segment.getLiquidViscosity());
      System.out.println(prefix + "gas_diffusivity_m2_per_s=" + segment.getGasDiffusivity());
      System.out.println(prefix + "liquid_diffusivity_m2_per_s=" + segment.getLiquidDiffusivity());
      System.out.println(prefix + "wetted_area_m2_per_m3=" + segment.getWettedArea());
      System.out.println(prefix + "kga=" + segment.getKGa());
      System.out.println(prefix + "kla=" + segment.getKLa());
      System.out.println(prefix + "gas_heat_transfer_coefficient_W_per_m3_K="
          + segment.getGasHeatTransferCoefficient());
      System.out.println(prefix + "liquid_heat_transfer_coefficient_W_per_m3_K="
          + segment.getLiquidHeatTransferCoefficient());
      System.out.println(prefix + "overall_heat_transfer_coefficient_W_per_m3_K="
          + segment.getOverallHeatTransferCoefficient());
      System.out.println(prefix + "interface_temperature_K=" + segment.getInterfaceTemperatureK());
      System.out.println(prefix + "heat_transfer_rate_W=" + segment.getHeatTransferRateW());
      System.out.println(prefix + "pressure_drop_per_meter_Pa=" + segment.getPressureDropPerMeter());
      System.out.println(prefix + "percent_flood=" + segment.getPercentFlood());
      System.out.println(prefix + "net_molar_transfer_mol_per_s=" + segment.getNetMolarTransfer());
      System.out.println(
          prefix + "enthalpy_balance_residual_W=" + segment.getEnthalpyBalanceResidualW());
      for (java.util.Map.Entry<String, Double> entry : segment.getInterfaceEquilibriumRatios()
          .entrySet()) {
        System.out.println(
            prefix + "interface_ratio_" + entry.getKey() + "=" + entry.getValue());
      }
      for (java.util.Map.Entry<String, Double> entry : segment.getComponentMoleTransfer()
          .entrySet()) {
        System.out.println(prefix + "transfer_" + entry.getKey() + "_mol_per_s=" + entry.getValue());
      }
    }
    System.out.println();
  }

  /// **The two capacity-limit families the columns carry, on solved states.**
  ///
  /// `DistillationColumn` computes `Fs = u·sqrt(rho_g)` against `maxAllowableFsFactor` (its own
  /// `2.5`); `AbsorptionColumn` computes the Souders-Brown `Ks = u·sqrt(rho_g/(rho_l - rho_g))`
  /// against `maxAllowableGasLoadFactor` (`0.15`) and overrides the Fs limit to `3.0` in its
  /// constructor. Neither family was measured before this probe: `absorption_column`,
  /// `stripping_column` and `distillation_column` all declare `max_allowable_gas_load_factor` and
  /// refuse it rather than reading it, and `getFsFactor()` is the quantity that would answer it.
  ///
  /// **Both families read the outlet streams, not a tray**, and both velocities are superficial
  /// over the **total** cross-section `pi·D^2/4` - which is what separates them from
  /// `TrayHydraulicsCalculator.getFsFactor()`, whose area is the net one less a downcomer. The
  /// volumetric flow is the *gas outlet's*, `getFlowRate("m3/sec")`; the Fs family takes the
  /// system's density and the K family `getPhase(0)`'s, so both are printed beside the lines they
  /// come from, along with `n·M/rho` as the port's own spelling of the same flow.
  ///
  /// **Rows that move one thing at a time**: the diameter (both velocities and the areas with it,
  /// the minimum diameters against), each stated limit (the utilizations and the inverses, not the
  /// factors, and one limit each side of the factor so both verdicts are in the capture), and
  /// rows that reflash an outlet, because the outlet is what both families read.
  static void capacityRows() {
    // The distillation column the port is already held to, at the class's own diameter of 1.0 m.
    capacityBinaryRow("binary_methane_butane_default_diameter", 1.0, null);
    // The same solved state at another diameter: both factors fall with the area, both minimum
    // diameters are invariant (they are what the factors are held to, not functions of D).
    capacityBinaryRow("binary_methane_butane_diameter_2", 2.0, null);
    capacityBinaryRow("binary_methane_butane_diameter_half", 0.5, null);
    // A stated limit either side of the factor the default diameter produces, so the capture
    // carries both verdicts.
    capacityBinaryRow("binary_methane_butane_fs_limit_5", 1.0, 5.0);
    capacityBinaryRow("binary_methane_butane_fs_limit_0_01", 1.0, 0.01);
    // The absorber, which carries both families and the subclass's own two defaults.
    capacityAbsorberRow("lean_oil_absorber_default_diameter", 1.0, null, null, false, "");
    capacityAbsorberRow("lean_oil_absorber_diameter_2", 2.0, null, null, false, "");
    capacityAbsorberRow("lean_oil_absorber_k_limit_0_30", 1.0, 0.30, null, false, "");
    capacityAbsorberRow("lean_oil_absorber_k_limit_0_005", 1.0, 0.005, null, false, "");
    capacityAbsorberRow("lean_oil_absorber_fs_limit_2_5", 1.0, null, 2.5, false, "");
    // **The same solved state with its liquid outlet reflashed to a single liquid**, which is the
    // two-sided control for the K family's density: on the row above, `getPhase(0)` of that outlet
    // is the *gas* it carries, `rho_l - rho_g` is `0.3` and the `10.0` fallback substitutes
    // `1000.0` - so closing the outlet to one liquid phase is what makes the real density visible.
    capacityAbsorberRow("lean_oil_absorber_liquid_out_liquid", 1.0, null, null, false, "liquid");
    // And the same outlet reflashed the other way, to a gas, which trips the fallback for the
    // opposite reason.
    capacityAbsorberRow("lean_oil_absorber_liquid_out_gas", 1.0, null, null, false, "gas");
    // **A two-phase gas outlet**, which separates the two families' density conventions: the Fs
    // family takes the system's, the K family `getPhase(0)`'s.
    capacityAbsorberRow("lean_oil_absorber_gas_out_two_phase", 1.0, null, null, false,
        "gas_two_phase");
    // The stripper, which is the same class under another name and inherits both defaults.
    capacityAbsorberRow("hydrocarbon_stripper_default_diameter", 1.0, null, null, true, "");
    // **The packed column, whose two diameters are not one.** `setColumnDiameter` writes the
    // class's own `columnDiameter` and `calcPackingHydraulics` *ends* by writing the inherited
    // `internalDiameter` from it - or from `getRequiredDiameter` when it is not stated - so the
    // Fs family reads a diameter the caller may never have named.
    capacityPackedRow("packed_distillation_binary_unstated", 2.0, "Pall-Ring-50", null, null);
    capacityPackedRow("packed_distillation_binary_column_diameter_0_5", 2.0, "Pall-Ring-50", 0.5,
        null);
    capacityPackedRow("packed_distillation_binary_fs_limit_5", 2.0, "Pall-Ring-50", 0.5, 5.0);
  }

  static void capacityPackedRow(String label, double packedHeightM, String packingType,
      Double columnDiameterM, Double fsLimit) {
    String[] names = new String[] { "methane", "n-butane" };
    double[] z = new double[] { 0.5, 0.5 };
    SystemInterface fluid = new SystemPrEos(300.0, 20.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    Stream inlet = new Stream("feed", fluid);
    inlet.setFlowRate(1000.0, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.PackedColumn column =
        new neqsim.process.equipment.distillation.PackedColumn(label, packedHeightM, packingType,
            true, true);
    column.addFeedStream(inlet, 2);
    column.setCondenserTemperature(-20.0, "C");
    column.setReboilerTemperature(100.0, "C");
    column.setTopPressure(19.0);
    column.setBottomPressure(20.0);
    column.setTemperatureTolerance(1.0e-6);
    column.setMaxNumberOfIterations(200, true);
    if (columnDiameterM != null) {
      column.setColumnDiameter(columnDiameterM);
    }
    if (fsLimit != null) {
      column.setMaxAllowableFsFactor(fsLimit);
    }
    neqsim.process.processmodel.ProcessSystem process = new neqsim.process.processmodel.ProcessSystem();
    process.add(inlet);
    process.add(column);
    process.run();
    // The inherited diameter is passed as the row's own, because that is what the Fs family
    // reads; `column_diameter_stated_m` is the class's other field, which does not reach it.
    printCapacity(label, "packed_column", column, column.getInternalDiameter(), columnDiameterM,
        fsLimit, null, "");
  }

  /// **`ColumnInternalsDesigner`, the trayed half, on solved columns of this probe's own.**
  ///
  /// **This is the class `TrayHydraulicsCalculator` is driven *by* - and it is not the tray
  /// calculator.** It walks the column's trays, picks the **controlling** one by the largest
  /// vapour mass flow, sizes a diameter off that tray, then runs one calculator per tray at the
  /// diameter it resolved and sums them: the total pressure drop, the two percent-flood extremes,
  /// the average efficiency and the per-tray results all come from that second loop, and the
  /// first calculator - the one the diameter was sized from - is discarded.
  ///
  /// **Three of its own lines are the port's, not the capture's, and each is measured here.**
  /// The pre-loop calculator is built **without a relative volatility** while the loop's carries
  /// one per tray, so the sizing branch runs at the class's default `2.0`. `getTrayProperties`
  /// takes the phase-0 composition for the gas side of the spread - the array is gas-first - and
  /// falls back to `600.0` for the liquid density when the fluid's phase 0 is a gas on a
  /// single-phase tray. And the spread itself is capped at `20.0`.
  ///
  /// **The ingredients are printed per tray beside every answer**, because the class's two
  /// extraction helpers are private: what they read off a tray is its two outlet streams and its
  /// own fluid's phases, and both are printed here as the class reads them.
  static void designerRows() {
    // The port's own oracle column, at the designer's own defaults.
    designerColumnRow("binary_designer_defaults", "sieve", 0.6, 0.05, 12.7, 0.1, 0.1, 0.8, -1.0);
    // The other two types the tray calculator names, at that same state.
    designerColumnRow("binary_designer_valve", "valve", 0.6, 0.05, 12.7, 0.1, 0.1, 0.8, -1.0);
    designerColumnRow("binary_designer_bubble_cap", "bubble-cap", 0.6, 0.05, 12.7, 0.1, 0.1, 0.8,
        -1.0);
    // A stated override, which replaces the sizing branch entirely.
    designerColumnRow("binary_designer_diameter_override_0_5", "sieve", 0.6, 0.05, 12.7, 0.1, 0.1,
        0.8, 0.5);
    // A looser flood fraction, which moves the sized diameter and every per-tray load with it.
    designerColumnRow("binary_designer_flood_0_7", "sieve", 0.6, 0.05, 12.7, 0.1, 0.1, 0.7, -1.0);
    // A different geometry, on the same state.
    designerColumnRow("binary_designer_wider_trays", "sieve", 0.75, 0.075, 6.0, 0.12, 0.15, 0.75,
        -1.0);
    // **The absorber, which has no ends**: its five trays are all middle trays, so the designer
    // walks five rather than the column's six.
    designerAbsorberRow("absorber_designer_defaults", "sieve", 0.6, 0.05, 12.7, 0.1, 0.1, 0.8,
        -1.0);
  }

  /// The binary column `process_column_capacity.tsv` already holds, driven by the designer.
  static void designerColumnRow(String label, String internalsType, double traySpacingM,
      double weirHeightM, double holeDiameterMm, double holeAreaFraction,
      double downcommerAreaFraction, double designFloodFraction, double diameterOverrideM) {
    SystemInterface fluid = new SystemPrEos(300.0, 20.0);
    fluid.addComponent("methane", 0.5);
    fluid.addComponent("n-butane", 0.5);
    fluid.setMixingRule(2);
    Stream inlet = new Stream("column feed", fluid);
    inlet.setFlowRate(1000.0, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn("col1", 4, true, true);
    column.addFeedStream(inlet, 2);
    column.setCondenserTemperature(-20.0, "C");
    column.setReboilerTemperature(100.0, "C");
    column.setTopPressure(19.0);
    column.setBottomPressure(20.0);
    column.setTemperatureTolerance(1.0e-6);
    column.setMaxNumberOfIterations(200, true);
    column.run();
    designerRow(label, "distillation_column", column, internalsType, traySpacingM, weirHeightM,
        holeDiameterMm, holeAreaFraction, downcommerAreaFraction, designFloodFraction,
        diameterOverrideM);
  }

  /// The lean-oil absorber, which is a column with no ends at all.
  static void designerAbsorberRow(String label, String internalsType, double traySpacingM,
      double weirHeightM, double holeDiameterMm, double holeAreaFraction,
      double downcommerAreaFraction, double designFloodFraction, double diameterOverrideM) {
    String[] names = new String[] { "methane", "ethane", "propane", "n-butane", "n-pentane",
        "n-heptane" };
    double[] gasZ = new double[] { 0.920, 0.040, 0.025, 0.010, 0.005, 0.0 };
    double[] pureHeptane = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
    double bara = 15.0;
    Stream gas = feedMass("gas", names, gasZ, 303.15, bara, 2000.0);
    Stream solvent = feedMass("solvent", names, pureHeptane, 293.15, bara, 600.0);

    neqsim.process.equipment.absorber.AbsorptionColumn column =
        new neqsim.process.equipment.absorber.AbsorptionColumn(label, 5);
    column.addGasInStream(gas);
    column.addSolventInStream(solvent);
    column.setTopPressure(bara);
    column.setBottomPressure(bara);
    column.setTemperatureTolerance(1.0e-4);
    column.setMassBalanceTolerance(5.0e-2);
    column.setEnthalpyBalanceTolerance(5.0e-2);
    column.setMaxNumberOfIterations(80, true);
    column.setSolverType(
        neqsim.process.equipment.distillation.DistillationColumn.SolverType.DIRECT_SUBSTITUTION);
    neqsim.process.processmodel.ProcessSystem process = new neqsim.process.processmodel.ProcessSystem();
    process.add(gas);
    process.add(solvent);
    process.add(column);
    process.run();
    designerRow(label, "absorption_column", column, internalsType, traySpacingM, weirHeightM,
        holeDiameterMm, holeAreaFraction, downcommerAreaFraction, designFloodFraction,
        diameterOverrideM);
  }

  static void designerRow(String label, String machine,
      neqsim.process.equipment.distillation.DistillationColumn column, String internalsType,
      double traySpacingM, double weirHeightM, double holeDiameterMm, double holeAreaFraction,
      double downcommerAreaFraction, double designFloodFraction, double diameterOverrideM) {
    neqsim.process.equipment.distillation.internals.ColumnInternalsDesigner designer =
        new neqsim.process.equipment.distillation.internals.ColumnInternalsDesigner(column);
    designer.setInternalsType(internalsType);
    designer.setTraySpacing(traySpacingM);
    designer.setWeirHeight(weirHeightM);
    designer.setHoleDiameter(holeDiameterMm);
    designer.setHoleAreaFraction(holeAreaFraction);
    designer.setDowncommerAreaFraction(downcommerAreaFraction);
    designer.setDesignFloodFraction(designFloodFraction);
    designer.setColumnDiameterOverride(diameterOverrideM);
    designer.calculate();

    System.out.println(label);
    System.out.println("machine=" + machine);
    System.out.println("internals_type=" + internalsType);
    System.out.println("tray_spacing_m=" + traySpacingM);
    System.out.println("weir_height_m=" + weirHeightM);
    System.out.println("hole_diameter_mm=" + holeDiameterMm);
    System.out.println("hole_area_fraction=" + holeAreaFraction);
    System.out.println("downcommer_area_fraction=" + downcommerAreaFraction);
    System.out.println("design_flood_fraction=" + designFloodFraction);
    System.out.println("column_diameter_override_m=" + diameterOverrideM);
    System.out.println("column_tray_count=" + column.getTrays().size());
    System.out.println("column_stage_count=" + column.getNumberOfTrays());
    // **What the designer reads off each tray, and what it resolves from it.** The two flows are
    // `getTrayFlows`'s own `kg/hr` over 3600; the five properties are `getTrayProperties`'s.
    for (int i = 0; i < column.getTrays().size(); i++) {
      neqsim.process.equipment.distillation.SimpleTray tray = column.getTrays().get(i);
      SystemInterface trayFluid = tray.getFluid();
      System.out.println("tray" + i + "_temperature_K=" + tray.getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + tray.getPressure());
      System.out.println("tray" + i + "_vapor_mass_flow_kg_per_s="
          + tray.getGasOutStream().getFlowRate("kg/hr") / 3600.0);
      System.out.println("tray" + i + "_liquid_mass_flow_kg_per_s="
          + tray.getLiquidOutStream().getFlowRate("kg/hr") / 3600.0);
      System.out.println("tray" + i + "_phases=" + trayFluid.getNumberOfPhases());
      trayFluid.initProperties();
      System.out.println("tray" + i + "_vapor_density_kg_per_m3="
          + (trayFluid.hasPhaseType("gas") ? trayFluid.getPhase("gas").getDensity("kg/m3")
              : trayFluid.getPhase(0).getDensity("kg/m3")));
      int liquidIndex = trayFluid.hasPhaseType("oil") ? trayFluid.getPhaseNumberOfPhase("oil")
          : (trayFluid.hasPhaseType("aqueous") ? trayFluid.getPhaseNumberOfPhase("aqueous") : 1);
      System.out.println("tray" + i + "_liquid_density_kg_per_m3="
          + trayFluid.getPhase(liquidIndex).getDensity("kg/m3"));
      System.out.println("tray" + i + "_liquid_viscosity_kg_per_msec="
          + trayFluid.getPhase(liquidIndex).getViscosity("kg/msec"));
      System.out.println("tray" + i + "_interphase_surface_tension_N_per_m="
          + trayFluid.getInterphaseProperties().getSurfaceTension(0, 1));
      System.out.println(
          "tray" + i + "_relative_volatility=" + relativeVolatility(trayFluid, liquidIndex));
    }
    // ---- The class's own answers.
    System.out.println("required_diameter_m=" + designer.getRequiredDiameter());
    System.out.println("controlling_tray_index=" + designer.getControllingTrayIndex());
    System.out.println("design_ok=" + designer.isDesignOk());
    System.out.println("max_percent_flood=" + designer.getMaxPercentFlood());
    System.out.println("min_percent_flood=" + designer.getMinPercentFlood());
    System.out.println("average_tray_efficiency=" + designer.getAverageTrayEfficiency());
    System.out.println("total_pressure_drop_Pa=" + designer.getTotalPressureDrop());
    System.out.println("total_pressure_drop_mbar=" + designer.getTotalPressureDropMbar());
    for (int i = 0; i < designer.getTrayResults().size(); i++) {
      neqsim.process.equipment.distillation.internals.TrayHydraulicsCalculator result =
          designer.getTrayResults().get(i);
      System.out.println("result" + i + "_percent_flood=" + result.getPercentFlood());
      System.out
          .println("result" + i + "_total_pressure_drop_Pa=" + result.getTotalTrayPressureDrop());
      System.out.println("result" + i + "_tray_efficiency=" + result.getTrayEfficiency());
      System.out.println("result" + i + "_design_ok=" + result.isDesignOk());
    }
    System.out.println();
  }

  /// `getTrayProperties`'s own relative-volatility spread, at a tray's fluid.
  ///
  /// The gas side is **phase 0**'s mole fraction - the class reads `getPhase(0)` and not the gas
  /// phase's - so on a gas-first system the two agree and on an all-liquid one they do not. The
  /// ratio is the largest `y_i/x_i` over the smallest, capped at `20.0`, and the default `2.0`
  /// stands wherever no component clears both `1e-10` thresholds.
  static double relativeVolatility(SystemInterface fluid, int liquidIndex) {
    double maxRatio = 0.0;
    double minRatio = 1.0e10;
    for (int i = 0; i < fluid.getPhase(0).getNumberOfComponents(); i++) {
      double xGas = fluid.getPhase(0).getComponent(i).getx();
      double xLiquid = fluid.getPhase(liquidIndex).getComponent(i).getx();
      if (xLiquid > 1.0e-10 && xGas > 1.0e-10) {
        double ratio = xGas / xLiquid;
        if (ratio > maxRatio) {
          maxRatio = ratio;
        }
        if (ratio < minRatio) {
          minRatio = ratio;
        }
      }
    }
    double relative = minRatio > 0.0 ? maxRatio / minRatio : 2.0;
    return Math.min(relative, 20.0);
  }

  /// **`hydraulicPressureDropCouplingEnabled`, which is a *solver* change and not a report.**
  ///
  /// `updatePressureProfileFromHydraulics` runs inside `updateColumnTearVariables`, the tear's
  /// own update: it builds a designer at the stated internals type, takes its
  /// **`getTotalPressureDrop`**, and calls `applyHydraulicPressureDrop`, which rewrites one end
  /// so the *difference* between the two equals that drop - `bottomTrayPressure =
  /// topTrayPressure + drop/1e5` where the top is positive, and `topTrayPressure =
  /// max(1e-6, bottomTrayPressure - drop/1e5)` otherwise - then rebuilds the profile. The
  /// residual it returns is the relative change it made to the end it rewrote.
  ///
  /// **The flag does not gate the coupling alone.** `hasActiveColumnTearVariables` is
  /// `!sideDrawSpecifications.isEmpty() || !pumparounds.isEmpty() ||
  /// hydraulicPressureDropCouplingEnabled` - so **turning the coupling on is what puts the
  /// column on `solveWithColumnTearVariables`' coordinated loop**, the same loop
  /// `process_column_tear.tsv` measured as not converging on a side-draw flow specification and
  /// a pumparound. The rows here are that measurement's twin, on states this port already
  /// reproduces.
  static void couplingRows() {
    // The control: the captured binary column with the coupling off, at the class's own tear
    // settings.
    couplingColumnRow("binary_coupling_off", false, "sieve", 1.0e-4, 30);
    // **And on**, which is the coordinated loop.
    couplingColumnRow("binary_coupling_on", true, "sieve", 1.0e-4, 30);
    // On at a looser tear tolerance, to see whether the loop converges when asked for less.
    couplingColumnRow("binary_coupling_on_tear_1e2", true, "sieve", 1.0e-2, 30);
    // On under the other internals type, whose per-tray drop is a different number.
    couplingColumnRow("binary_coupling_on_valve", true, "valve", 1.0e-4, 30);
    // **The absorber**, which is the same machine with no ends: its trays are all middle ones
    // and its two ends are not reboiler and condenser.
    couplingAbsorberRow("absorber_coupling_off", false, "sieve", 1.0e-4, 30);
    couplingAbsorberRow("absorber_coupling_on", true, "sieve", 1.0e-4, 30);
  }

  /// The captured binary column, driven with the coupling on or off.
  static void couplingColumnRow(String label, boolean enabled, String internalsType,
      double tearTolerance, int tearIterations) {
    SystemInterface fluid = new SystemPrEos(300.0, 20.0);
    fluid.addComponent("methane", 0.5);
    fluid.addComponent("n-butane", 0.5);
    fluid.setMixingRule(2);
    Stream inlet = new Stream("column feed", fluid);
    inlet.setFlowRate(1000.0, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn("col1", 4, true, true);
    column.addFeedStream(inlet, 2);
    column.setCondenserTemperature(-20.0, "C");
    column.setReboilerTemperature(100.0, "C");
    column.setTopPressure(19.0);
    column.setBottomPressure(20.0);
    column.setTemperatureTolerance(1.0e-6);
    column.setMaxNumberOfIterations(200, true);
    neqsim.process.processmodel.ProcessSystem process = new neqsim.process.processmodel.ProcessSystem();
    process.add(inlet);
    process.add(column);
    couplingRow(label, "distillation_column", column, process, enabled, internalsType,
        tearTolerance, tearIterations);
  }

  /// The lean-oil absorber, which is the same class with two inlets and no ends.
  static void couplingAbsorberRow(String label, boolean enabled, String internalsType,
      double tearTolerance, int tearIterations) {
    String[] names = new String[] { "methane", "ethane", "propane", "n-butane", "n-pentane",
        "n-heptane" };
    double[] gasZ = new double[] { 0.920, 0.040, 0.025, 0.010, 0.005, 0.0 };
    double[] pureHeptane = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
    double bara = 15.0;
    Stream gas = feedMass("gas", names, gasZ, 303.15, bara, 2000.0);
    Stream solvent = feedMass("solvent", names, pureHeptane, 293.15, bara, 600.0);

    neqsim.process.equipment.absorber.AbsorptionColumn column =
        new neqsim.process.equipment.absorber.AbsorptionColumn(label, 5);
    column.addGasInStream(gas);
    column.addSolventInStream(solvent);
    column.setTopPressure(bara);
    column.setBottomPressure(bara);
    column.setTemperatureTolerance(1.0e-4);
    column.setMassBalanceTolerance(5.0e-2);
    column.setEnthalpyBalanceTolerance(5.0e-2);
    column.setMaxNumberOfIterations(80, true);
    column.setSolverType(
        neqsim.process.equipment.distillation.DistillationColumn.SolverType.DIRECT_SUBSTITUTION);
    neqsim.process.processmodel.ProcessSystem process = new neqsim.process.processmodel.ProcessSystem();
    process.add(gas);
    process.add(solvent);
    process.add(column);
    couplingRow(label, "absorption_column", column, process, enabled, internalsType, tearTolerance,
        tearIterations);
  }

  static void couplingRow(String label, String machine,
      neqsim.process.equipment.distillation.DistillationColumn column,
      neqsim.process.processmodel.ProcessSystem process, boolean enabled, String internalsType,
      double tearTolerance, int tearIterations) {
    double statedTop = column.getTopPressure();
    double statedBottom = column.getBottomPressure();
    if (enabled) {
      column.enableHydraulicPressureDropCoupling(internalsType);
    }
    column.setColumnTearTolerance(tearTolerance);
    column.setMaxColumnTearIterations(tearIterations);
    process.run();

    System.out.println(label);
    System.out.println("machine=" + machine);
    System.out.println("coupling_enabled=" + enabled);
    System.out.println("internals_type=" + internalsType);
    System.out.println("tear_tolerance=" + tearTolerance);
    System.out.println("max_tear_iterations=" + tearIterations);
    System.out.println("top_pressure_stated_bara=" + statedTop);
    System.out.println("bottom_pressure_stated_bara=" + statedBottom);
    System.out.println("top_pressure_bara=" + column.getTopPressure());
    System.out.println("bottom_pressure_bara=" + column.getBottomPressure());
    System.out.println("solved=" + column.solved());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    System.out.println("mass_residual=" + column.getLastMassResidual());
    System.out.println("energy_residual=" + column.getLastEnergyResidual());
    // ---- The tear's own reporting, which is where a coupled solve says what it did.
    System.out.println("tear_converged=" + column.isLastColumnTearConverged());
    System.out.println("tear_iterations=" + column.getLastColumnTearIterationCount());
    System.out.println("tear_residual=" + column.getLastColumnTearResidual());
    System.out.println("tear_rejected_candidates=" + column.getLastColumnTearRejectedCandidateCount());
    System.out.println("tear_rollbacks=" + column.getLastColumnTearRollbackCount());
    System.out.println("tear_inner_iterations=" + column.getLastColumnTearInnerIterationCount());
    System.out.println("tear_candidate_history=" + column.getLastColumnTearCandidateHistory());
    System.out.println("last_hydraulic_pressure_drop_Pa=" + column.getLastHydraulicPressureDropPa());
    System.out.println("last_hydraulic_pressure_drop_residual="
        + column.getLastHydraulicPressureDropResidual());
    for (int i = 0; i < column.getNumberOfTrays(); i++) {
      System.out.println("tray" + i + "_temperature_K=" + column.getTray(i).getTemperature());
      System.out.println("tray" + i + "_pressure_bara=" + column.getTray(i).getPressure());
    }
    print("gas_out", column.getGasOutStream());
    print("liquid_out", column.getLiquidOutStream());
    System.out.println();
  }

  /// **`DistillationColumnMechanicalDesign`, which is the trayed column's own sizing.**
  ///
  /// `calcDesign` reads **tray 0** - the bottom end - and nothing else: the vapour outlet's molar
  /// flow and its fluid's density and molar mass, then the liquid outlet's. From those it sizes
  /// the diameter and the height, takes the wall thickness, and reports the tray efficiency, the
  /// two flooding factors, the weir loading and the two duties. **It is the class the designer's
  /// report feeds**: it reads the same internals tree the coupling does, so this capture and
  /// `process_internals_designer.tsv` are the two halves of one machine.
  ///
  /// **No capture in this directory held it**, and the class is 1,348 lines with a public surface
  /// of forty-odd getters, so what is printed here is the surface rather than a sample of it.
  static void mechanicalRows() {
    mechanicalBinaryRow("binary_mechanical_defaults", null, null);
    // **The stated-diameter and the valve rows exist to make two of the class's own branches
    // measurable.** A diameter override replaces the rating diameter the designer is driven at,
    // and that is the *only* thing it replaces: `floodingFactor`, `weirLength` and `weirLoading`
    // were computed one statement earlier, at the Souders-Brown diameter the override throws
    // away - so this row's three quantities are the defaults row's, to the last digit, beside a
    // `getColumnDiameter()` of `2.0`. A tray type of `valve` moves the Souders-Brown `kFactor`
    // from `0.1` to `0.12` and the designer's own internals type with it.
    mechanicalBinaryRow("binary_mechanical_override", 2.0, null);
    mechanicalBinaryRow("binary_mechanical_valve", null, "valve");
    mechanicalAbsorberRow("absorber_mechanical");
    // **`AbsorberMechanicalDesign` has no row, and the reason is measured**: its only constructor
    // call site in `src/main` is `SimpleAbsorber.getMechanicalDesign()` (`:383`), and
    // `SimpleAbsorber extends Separator` - so driving it on an `AbsorptionColumn`, which extends
    // `DistillationColumn`, throws `ClassCastException` from `PressureVesselDesignStandard:36`.
    // The class is the separator family's sizing for a machine this port does not carry, which is
    // why `process.absorption_column`'s mechanical design is the *column's* design and not this.
    // **And the stripper's own row**, appended last so the four rows above keep their indices. It
    // is the same state as the product capture's `hydrocarbon_stripper` row - the same subclass,
    // the same two feeds, the same gate - and it exists because the stripper inherits every
    // getter and overrides none: one absorber row cannot settle a second reading.
    mechanicalStripperRow("hydrocarbon_stripper_mechanical");
  }

  /// **The stripper's mechanical row**, on the state `absorberRows()` captures as
  /// `hydrocarbon_stripper`.
  static void mechanicalStripperRow(String label) {
    AbsorberCase unit = runAbsorber(label, STRIPPER_COMPONENTS, STRIPPER_GAS_Z, 150.0, 70.0,
        STRIPPER_COMPONENTS, STRIPPER_LIQUID_Z, 900.0, 70.0, 5, 12.0, 0.0, 1.0e-4,
        "DIRECT_SUBSTITUTION", false, null, null, 0.0);
    mechanicalRow(label, "stripping_column",
        (neqsim.process.equipment.distillation.DistillationColumn) unit.column, null, null);
  }

  /// The port's own binary column, at the class's own constructor defaults.
  static void mechanicalBinaryRow(String label, Double diameterOverride, String trayType) {
    SystemInterface fluid = new SystemPrEos(300.0, 20.0);
    fluid.addComponent("methane", 0.5);
    fluid.addComponent("n-butane", 0.5);
    fluid.setMixingRule(2);
    Stream inlet = new Stream("column feed", fluid);
    inlet.setFlowRate(1000.0, "kg/hr");
    inlet.run();

    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn("col1", 4, true, true);
    column.addFeedStream(inlet, 2);
    column.setCondenserTemperature(-20.0, "C");
    column.setReboilerTemperature(100.0, "C");
    column.setTopPressure(19.0);
    column.setBottomPressure(20.0);
    column.setTemperatureTolerance(1.0e-6);
    column.setMaxNumberOfIterations(200, true);
    neqsim.process.processmodel.ProcessSystem process = new neqsim.process.processmodel.ProcessSystem();
    process.add(inlet);
    process.add(column);
    process.run();
    mechanicalRow(label, "distillation_column", column, diameterOverride, trayType);
  }

  /// The lean-oil absorber, which is the same class on a column with no ends.
  static void mechanicalAbsorberRow(String label) {
    String[] names = new String[] { "methane", "ethane", "propane", "n-butane", "n-pentane",
        "n-heptane" };
    double[] gasZ = new double[] { 0.920, 0.040, 0.025, 0.010, 0.005, 0.0 };
    double[] pureHeptane = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
    double bara = 15.0;
    Stream gas = feedMass("gas", names, gasZ, 303.15, bara, 2000.0);
    Stream solvent = feedMass("solvent", names, pureHeptane, 293.15, bara, 600.0);

    neqsim.process.equipment.absorber.AbsorptionColumn column =
        new neqsim.process.equipment.absorber.AbsorptionColumn(label, 5);
    column.addGasInStream(gas);
    column.addSolventInStream(solvent);
    column.setTopPressure(bara);
    column.setBottomPressure(bara);
    column.setTemperatureTolerance(1.0e-4);
    column.setMassBalanceTolerance(5.0e-2);
    column.setEnthalpyBalanceTolerance(5.0e-2);
    column.setMaxNumberOfIterations(80, true);
    column.setSolverType(
        neqsim.process.equipment.distillation.DistillationColumn.SolverType.DIRECT_SUBSTITUTION);
    neqsim.process.processmodel.ProcessSystem process = new neqsim.process.processmodel.ProcessSystem();
    process.add(gas);
    process.add(solvent);
    process.add(column);
    process.run();
    mechanicalRow(label, "absorption_column", column, null, null);
  }

  static void mechanicalRow(String label, String machine,
      neqsim.process.equipment.distillation.DistillationColumn column, Double diameterOverride,
      String trayType) {
    neqsim.process.mechanicaldesign.distillation.DistillationColumnMechanicalDesign design =
        new neqsim.process.mechanicaldesign.distillation.DistillationColumnMechanicalDesign(column);
    if (diameterOverride != null) {
      design.setColumnDiameterOverride(diameterOverride);
    }
    if (trayType != null) {
      design.setTrayType(trayType);
    }
    design.calcDesign();

    System.out.println(label);
    System.out.println("machine=" + machine);
    // ---- What it reads: tray 0's two outlets, which is the whole of its input surface.
    neqsim.process.equipment.distillation.SimpleTray bottom = column.getTray(0);
    System.out.println("tray0_vapor_mol_per_hr="
        + bottom.getGasOutStream().getFlowRate("mol/hr"));
    System.out.println("tray0_vapor_density_kg_per_m3="
        + bottom.getGasOutStream().getFluid().getDensity("kg/m3"));
    System.out.println("tray0_vapor_molar_mass_kg_per_mol="
        + bottom.getGasOutStream().getFluid().getMolarMass() );
    System.out.println("tray0_liquid_mol_per_hr="
        + bottom.getLiquidOutStream().getFlowRate("mol/hr"));
    System.out.println("tray0_liquid_density_kg_per_m3="
        + bottom.getLiquidOutStream().getFluid().getDensity("kg/m3"));
    System.out.println("tray0_liquid_molar_mass_kg_per_mol="
        + bottom.getLiquidOutStream().getFluid().getMolarMass());
    // **And the same read after the fluid is asked to initialise**, which is the mechanism the
    // zero above is. Nothing else changes: the same stream, the same getter, one call apart.
    bottom.getLiquidOutStream().getFluid().initProperties();
    System.out.println("tray0_liquid_density_after_init_kg_per_m3="
        + bottom.getLiquidOutStream().getFluid().getDensity("kg/m3"));
    neqsim.process.mechanicaldesign.distillation.DistillationColumnMechanicalDesign again =
        new neqsim.process.mechanicaldesign.distillation.DistillationColumnMechanicalDesign(column);
    again.calcDesign();
    System.out.println("weir_loading_after_init=" + again.getWeirLoading());
    System.out.println("tray_pressure_drop_after_init_Pa=" + again.getTrayPressureDrop());
    System.out.println("column_diameter_after_init_m=" + again.getColumnDiameter());
    // ---- And what it answers.
    System.out.println("number_of_trays=" + design.getNumberOfTrays());
    System.out.println("actual_trays=" + design.getActualTrays());
    System.out.println("tray_efficiency=" + design.getTrayEfficiency());
    System.out.println("tray_spacing_m=" + design.getTraySpacing());
    System.out.println("tray_type=" + design.getTrayType());
    System.out.println("column_diameter_m=" + design.getColumnDiameter());
    System.out.println("column_height_m=" + design.getColumnHeight());
    System.out.println("column_wall_thickness_m=" + design.getColumnWallThickness());
    // **The `_m`, `_Pa` and `_W` suffixes on this block are the probe's, not the class's**, and
    // the three keys beside them state the class's own documented unit for the same number:
    // `getColumnWallThickness` is mm, `getTrayPressureDrop` mbar/tray, `getTotalPressureDrop` bar
    // and both duties kW. Reading the suffix as the class's unit is what made an earlier version
    // of the README entry claim four unit slips that were this file's own labels.
    System.out.println("column_wall_thickness_mm=" + design.getColumnWallThickness());
    System.out.println("flooding_factor=" + design.getFloodingFactor());
    System.out.println("max_flooding_factor=" + design.getMaxFloodingFactor());
    System.out.println("weir_loading=" + design.getWeirLoading());
    System.out.println("tray_pressure_drop_Pa=" + design.getTrayPressureDrop());
    System.out.println("tray_pressure_drop_mbar_per_tray=" + design.getTrayPressureDrop());
    System.out.println("total_pressure_drop_Pa=" + design.getTotalPressureDrop());
    System.out.println("total_pressure_drop_bar=" + design.getTotalPressureDrop());
    System.out.println("reboiler_duty_W=" + design.getReboilerDuty());
    System.out.println("reboiler_duty_kw=" + design.getReboilerDuty());
    System.out.println("condenser_duty_W=" + design.getCondenserDuty());
    System.out.println("condenser_duty_kw=" + design.getCondenserDuty());
    System.out.println("material_grade=" + design.getMaterialGrade());
    // ---- The four quantities the wall-thickness formula is built from, so its own arithmetic
    // is checkable from this capture rather than from the source alone.
    System.out.println("max_operation_pressure=" + design.getMaxOperationPressure());
    System.out.println("tensile_strength=" + design.getTensileStrength());
    System.out.println("joint_efficiency=" + design.getJointEfficiency());
    System.out.println("corrosion_allowance=" + design.getCorrosionAllowance());
    System.out.println("design_pressure_times_1_1=" + design.getMaxOperationPressure() * 1.1);
    System.out.println("contactor_internals_type=" + design.getContactorInternalsType());
    // ---- And the designer it built, so the two halves of `calcDesign` are separated here: the
    // Souders-Brown pass above and the internals tree below, which is what replaces its diameter.
    neqsim.process.equipment.distillation.internals.ColumnInternalsDesigner internals =
        design.getContactorInternalsDesigner();
    System.out.println("internals_required_diameter_m=" + internals.getRequiredDiameter());
    System.out.println("internals_max_percent_flood=" + internals.getMaxPercentFlood());
    System.out.println("internals_total_pressure_drop_pa=" + internals.getTotalPressureDrop());
    System.out.println("column_diameter_override_m=" + design.getColumnDiameterOverride());
    System.out.println("max_contactor_pressure_drop_bar=" + design.getMaxContactorPressureDropBar());
    System.out.println();
  }

  /// **`ReactorMechanicalDesign`, which is the reactor vessel's own sizing.**
  ///
  /// `calcDesign` reads the **inlet** stream of a `TwoPortInterface` and nothing else: its
  /// pressure, temperature, volumetric flow and the fluid's density, then sizes a catalyst volume
  /// from a space velocity, a vessel from a bed fill fraction and a target L/D, an ASME VIII Div 1
  /// shell and head, an Ergun bed drop and a weight. **No capture in this directory held the
  /// class.**
  ///
  /// **Three rows, and between them every branch that changes an answer.** The defaults row is
  /// the ammonia feed at 300 bara, where the GHSV division floors the catalyst volume at `0.01`
  /// m3 and the shell is thick. The CSTR row takes the LHSV branch instead, which is the only
  /// place `reactorType` is read, and its `1.5` 1/hr is chosen so the vessel is **above** the
  /// `0.5` m diameter floor: at these flows every default-GHSV row floors it, and a capture that
  /// only ever reads `0.5` does not measure the cube root. The high-pressure row feeds 1 bara with a stated
  /// `maxOperationPressure` of `150.0`, so `max(inlet, maxOperation)` takes the max branch, and
  /// its 299 m3/hr over the default `5000` keeps the GHSV division **above** the floor.
  ///
  /// **Five configured parameters have no getter.** `ghsv`, `lhsv`, `catalystBulkDensity`,
  /// `catalystParticleDiameter` and `bedVoidFraction` are setters only, so each row prints what
  /// it asked for under `requested_` - the class's own initialisers on the defaults row.
  static void reactorMechanicalRows() {
    reactorMechanicalRow("gibbs_fixed_bed_defaults",
        new String[] { "hydrogen", "nitrogen", "ammonia" }, new double[] { 1.5, 0.5, 0.0 },
        450.0, 300.0, 1.0, null, null, null);
    reactorMechanicalRow("gibbs_cstr",
        new String[] { "hydrogen", "nitrogen", "ammonia" }, new double[] { 1.5, 0.5, 0.0 },
        450.0, 300.0, 1.0, "CSTR", 1.5, null);
    reactorMechanicalRow("gibbs_high_pressure",
        new String[] { "CO2", "hydrogen", "CO", "water" }, new double[] { 1.0, 1.0, 0.0, 0.0 },
        1000.0, 1.0, 1.0, null, null, 150.0);
  }

  /// One reactor mechanical row: the feed the class's own Gibbs rows state, the reactor run,
  /// then the design on the inlet it read.
  static void reactorMechanicalRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, double molPerSecond, String reactorType, Double lhsv,
      Double maxOperationPressure) {
    Stream inlet = feed(names, z, temperatureK, pressureBara, molPerSecond);
    neqsim.process.equipment.reactor.GibbsReactor reactor =
        new neqsim.process.equipment.reactor.GibbsReactor("gibbs", inlet);
    reactor.setUseAllDatabaseSpecies(false);
    reactor.run();

    neqsim.process.mechanicaldesign.reactor.ReactorMechanicalDesign design =
        reactor.getMechanicalDesign();
    if (reactorType != null) {
      design.setReactorType(reactorType);
    }
    if (lhsv != null) {
      design.setLHSV(lhsv);
    }
    if (maxOperationPressure != null) {
      design.setMaxOperationPressure(maxOperationPressure);
    }
    design.calcDesign();

    System.out.println(label);
    // ---- What it reads: the inlet's own state, which is the whole of its input surface.
    System.out.println("inlet_p_bara=" + inlet.getPressure("bara"));
    System.out.println("inlet_t_c=" + inlet.getTemperature("C"));
    System.out.println("inlet_flow_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("inlet_density_kg_per_m3=" + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println(
        "inlet_molar_mass_kg_per_mol=" + inlet.getThermoSystem().getMolarMass());
    System.out.println("inlet_moles_per_sec=" + inlet.getFlowRate("mol/sec"));
    // ---- The configuration, from the getters where the class has one and from the probe where
    // it does not.
    System.out.println("requested_reactor_type="
        + (reactorType == null ? "FIXED_BED" : reactorType));
    System.out.println("requested_ghsv_per_hr=5000.0");
    System.out.println("requested_lhsv_per_hr=" + (lhsv == null ? 1.0 : lhsv));
    System.out.println("requested_catalyst_bulk_density_kg_per_m3=800.0");
    System.out.println("requested_catalyst_particle_diameter_mm=3.0");
    System.out.println("requested_bed_void_fraction=0.4");
    System.out.println("class_reactor_type=" + design.getReactorType());
    System.out.println("class_number_of_beds=" + design.getNumberOfBeds());
    System.out.println("class_max_operation_pressure=" + design.getMaxOperationPressure());
    System.out.println("class_max_operation_temperature=" + design.getMaxOperationTemperature());
    System.out.println("class_joint_efficiency=" + design.getJointEfficiency());
    System.out.println("class_corrosion_allowance=" + design.getCorrosionAllowance());
    // ---- And what it answers.
    System.out.println("vessel_diameter_m=" + design.getVesselDiameter());
    System.out.println("vessel_length_m=" + design.getVesselLength());
    System.out.println("shell_thickness_mm=" + design.getShellThickness());
    System.out.println("head_thickness_mm=" + design.getHeadThickness());
    System.out.println("design_pressure_bara=" + design.getDesignPressureBara());
    System.out.println("design_temperature_c=" + design.getDesignTemperatureC());
    System.out.println("catalyst_volume_m3=" + design.getCatalystVolume());
    System.out.println("catalyst_mass_kg=" + design.getCatalystMass());
    System.out.println("bed_pressure_drop_bar=" + design.getBedPressureDrop());
    System.out.println("empty_vessel_weight_kg=" + design.getEmptyVesselWeight());
    System.out.println("total_equipped_weight_kg=" + design.getTotalEquippedWeight());
    System.out.println("reactor_type=" + design.getReactorType());
    System.out.println("number_of_beds=" + design.getNumberOfBeds());
    System.out.println();
  }

  /// **`ControlValveSizing`, the `ValveMechanicalDesign` default, and the two solve modes.**
  ///
  /// The class's own sizing tests re-cased on PR: the rich-gas regression state
  /// (`ThrottlingValveSizingStateRegressionTest`), the methane state and a water liquid. Each row
  /// runs the valve unsized - the first `run` calls `calcKv` through `getMechanicalDesign()` -
  /// then prints what the sizing reads and what it answers, and then drives the two modes with the
  /// Kv it produced: `isCalcPressure` solves the outlet pressure from the trim, and the flow mode
  /// solves the mole flow from the drop.
  static void valveSizingRows() {
    String[] rich = new String[] { "methane", "ethane", "propane", "n-butane" };
    valveSizingRow("sizing_gas_rich", rich, new double[] { 0.1, 0.2, 0.4, 0.3 }, 296.25, 1.82925,
        14942.5, 1.41325, 100.0, false, false);
    valveSizingRow("sizing_gas_rich_isothermal", rich, new double[] { 0.1, 0.2, 0.4, 0.3 },
        296.25, 1.82925, 14942.5, 1.41325, 100.0, true, false);
    valveSizingRow("sizing_gas_rich_choked_cap", rich, new double[] { 0.1, 0.2, 0.4, 0.3 },
        296.25, 1.82925, 14942.5, 1.41325, 100.0, false, true);
    valveSizingRow("sizing_methane", new String[] { "methane" }, new double[] { 1.0 }, 293.15,
        10.0, 4749.1, 9.0, 100.0, false, false);
    valveSizingRow("sizing_water_liquid", new String[] { "water" }, new double[] { 1.0 }, 289.15,
        14.8, 1.69272e6, 4.46, 100.0, false, false);
  }

  static void valveSizingRow(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, double kgPerHour, double outletBara, double openingPercent,
      boolean isothermal, boolean allowChoked) {
    Stream inlet = feedMass("inlet", names, z, temperatureK, pressureBara, kgPerHour);

    neqsim.process.equipment.valve.ThrottlingValve valve =
        new neqsim.process.equipment.valve.ThrottlingValve(label, inlet);
    valve.setOutletPressure(outletBara);
    valve.setPercentValveOpening(openingPercent);
    valve.setIsoThermal(isothermal);
    valve.getMechanicalDesign().getValveSizingMethod().setAllowChoked(allowChoked);

    System.out.println(label);
    // ---- What the sizing reads.
    System.out.println("inlet_p_bara=" + inlet.getPressure("bara"));
    System.out.println("inlet_t_k=" + inlet.getTemperature("K"));
    System.out.println("inlet_mol_per_sec=" + inlet.getFlowRate("mol/sec"));
    System.out.println("inlet_kg_per_hr=" + inlet.getFlowRate("kg/hr"));
    System.out.println("inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("inlet_density_kg_per_m3="
        + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("inlet_z=" + inlet.getThermoSystem().getZ());
    System.out.println("inlet_gamma2=" + inlet.getThermoSystem().getGamma2());
    System.out.println("inlet_molar_mass_gr_per_mol="
        + inlet.getThermoSystem().getMolarMass("gr/mol"));
    StringBuilder feedZ = new StringBuilder("feed_z=");
    for (int i = 0; i < inlet.getThermoSystem().getNumberOfComponents(); i++) {
      feedZ.append(inlet.getThermoSystem().getComponent(i).getComponentName()).append(":")
          .append(inlet.getThermoSystem().getComponent(i).getz());
      if (i + 1 < inlet.getThermoSystem().getNumberOfComponents()) {
        feedZ.append(" ");
      }
    }
    System.out.println(feedZ);
    System.out.println("outlet_p_stated_bara=" + outletBara);
    System.out.println("opening_percent=" + openingPercent);
    System.out.println("is_isothermal=" + isothermal);
    System.out.println("allow_choked=" + allowChoked);

    // ---- The first run sizes the valve, because no Kv was set.
    valve.run();
    System.out.println("gas_valve=" + valve.isGasValve());
    System.out.println("kv=" + valve.getKv());
    System.out.println("cv=" + valve.getCv());
    System.out.println("cg=" + valve.getCg());
    java.util.Map<String, Object> sized =
        ((neqsim.process.mechanicaldesign.valve.ValveMechanicalDesign) valve.getMechanicalDesign())
            .calcValveSize();
    System.out.println("choked=" + sized.get("choked"));
    System.out.println("outlet_stream_p_bara=" + valve.getOutletStream().getPressure("bara"));
    System.out.println("outlet_stream_t_k=" + valve.getOutletStream().getTemperature("K"));
    System.out.println("calculated_molar_flow=" + valve.calculateMolarFlow());

    // ---- The pressure mode: the drop from the trim, at the same opening.
    valve.setIsCalcOutPressure(true);
    valve.run();
    System.out.println("solved_outlet_p_bara=" + valve.getOutletPressure());
    System.out.println("solved_outlet_stream_p_bara="
        + valve.getOutletStream().getPressure("bara"));
    System.out.println("solved_outlet_stream_t_k=" + valve.getOutletStream().getTemperature("K"));

    // ---- The flow mode, with the Kv the sizing produced stated explicitly.
    double sizedKv = valve.getKv();
    valve.setIsCalcOutPressure(false);
    valve.setKv(sizedKv);
    valve.setOutletPressure(outletBara);
    valve.run();
    System.out.println("stated_kv=" + valve.getKv());
    System.out.println("flow_mode_molar_flow=" + valve.calculateMolarFlow());
    System.out.println();
  }

  /// ``PumpChart`` and ``Pump``'s chart branch.
  ///
  /// **The chart's own arithmetic is fluid-free**, so those rows compare a reduced-polynomial
  /// fit and nothing else: the eight-speed table ``PumpChartTest`` states, read at the flow
  /// and speed its own ``Pump.run`` reaches, at both edges of the measured reduced-flow range,
  /// and past both. The density, NPSH and viscosity rows are ``PumpDensityCorrectionTest``,
  /// ``PumpNPSHCurveTest`` and ``PumpViscosityCorrectionTest``'s states.
  ///
  /// **The run rows use PR water, not the SRK water the class tests build.** The chart branch
  /// multiplies the head by the inlet density, and the port's fluid is PR; a row on SRK water
  /// would price the equation of state rather than the chart.
  static void pumpCurveRows() {
    neqsim.process.equipment.pump.PumpChart test =
        chartTest(new double[] {0.3, 1.0, 1.0, 1.0});

    chartFitRow("chart_test_run_point", test, 3974.8455, 12913.0);
    chartFitRow("chart_test_second_speed", test, 3599.6333, 10453.0);
    chartFitRow("chart_test_lowest_speed", test, 1636.5807, 8200.0);
    chartFitRow("chart_test_above_reduced_range", test, 6000.0, 12913.0);
    chartFitRow("chart_test_below_reduced_range", test, 100.0, 12913.0);
    System.out.println("chart_bep_flow_at_reference=" + test.getBestEfficiencyFlowRate());
    System.out.println("chart_bep_flow_at_1500=" + test.getBestEfficiencyFlowRate(1500.0));
    System.out.println();

    // The class's two head units are one specific energy, so the same three-point curve in
    // metre and in kJ/kg has to give one pressure rise. `testHeadUnitKJkg`'s own factor.
    neqsim.process.equipment.pump.PumpChart metre = threePointChart("meter");
    neqsim.process.equipment.pump.PumpChart kJkg = threePointChart("kJ/kg");
    System.out.println("head_unit_meter_at_40=" + metre.getHead(40.0, 500.0));
    System.out.println("head_unit_kjkg_at_40=" + kJkg.getHead(40.0, 500.0));
    System.out.println("head_unit_meter_delta_p_pa="
        + metre.getHead(40.0, 500.0) * 1000.0 * 9.80665);
    System.out.println("head_unit_kjkg_delta_p_pa="
        + kJkg.getHead(40.0, 500.0) * 1000.0 * 1000.0);
    System.out.println();

    neqsim.process.equipment.pump.PumpChart density =
        simpleChart(new double[] {18.0, 298.15, 1.0, 1.0, 998.0});
    System.out.println("density_reference=" + density.getReferenceDensity());
    System.out.println("density_has_correction=" + density.hasDensityCorrection());
    System.out.println("density_chart_head_at_40=" + density.getHead(40.0, 1000.0));
    System.out.println("density_corrected_at_998="
        + density.getCorrectedHead(40.0, 1000.0, 998.0));
    System.out.println("density_corrected_at_800="
        + density.getCorrectedHead(40.0, 1000.0, 800.0));
    System.out.println("density_corrected_at_1200="
        + density.getCorrectedHead(40.0, 1000.0, 1200.0));
    neqsim.process.equipment.pump.PumpChart noDensity =
        simpleChart(new double[] {18.0, 298.15, 1.0, 1.0});
    System.out.println("density_absent_reference=" + noDensity.getReferenceDensity());
    System.out.println("density_absent_has_correction=" + noDensity.hasDensityCorrection());
    System.out.println("density_absent_corrected_at_800="
        + noDensity.getCorrectedHead(40.0, 1000.0, 800.0));
    System.out.println();

    neqsim.process.equipment.pump.PumpChart npsh = npshChart();
    npsh.setNPSHCurve(new double[][] {{2.0, 2.2, 2.5, 3.0, 3.8, 4.8, 6.2, 8.0},
        {4.5, 4.95, 5.625, 6.75, 8.55, 10.8, 13.95, 18.0}});
    System.out.println("npsh_has_curve=" + npsh.hasNPSHCurve());
    System.out.println("npsh_at_40_1000=" + npsh.getNPSHRequired(40.0, 1000.0));
    System.out.println("npsh_at_60_1500=" + npsh.getNPSHRequired(60.0, 1500.0));
    System.out.println("npsh_at_35_1000=" + npsh.getNPSHRequired(35.0, 1000.0));
    System.out.println("npsh_at_150_1000=" + npsh.getNPSHRequired(150.0, 1000.0));
    System.out.println("npsh_at_1_1000=" + npsh.getNPSHRequired(1.0, 1000.0));
    System.out.println();

    neqsim.process.equipment.pump.PumpChart viscous =
        simpleChart(new double[] {18.0, 298.15, 1.0, 1.0, 998.0});
    viscous.setUseViscosityCorrection(true);
    for (double viscosity : new double[] {1.0, 50.0, 75.0, 100.0, 500.0}) {
      viscous.calculateViscosityCorrection(viscosity, 30.0, 88.0, 1000.0);
      System.out.println("viscosity=" + viscosity);
      System.out.println("viscosity_cq=" + viscous.getFlowCorrectionFactor());
      System.out.println("viscosity_ch=" + viscous.getHeadCorrectionFactor());
      System.out.println("viscosity_ceta=" + viscous.getEfficiencyCorrectionFactor());
      System.out.println("viscosity_base_head=" + viscous.getHead(30.0, 1000.0));
      System.out.println("viscosity_base_efficiency=" + viscous.getEfficiency(30.0, 1000.0));
      System.out.println("viscosity_corrected_head="
          + viscous.getViscosityCorrectedHead(viscous.getHead(30.0, 1000.0)));
      System.out.println("viscosity_corrected_efficiency="
          + viscous.getViscosityCorrectedEfficiency(viscous.getEfficiency(30.0, 1000.0)));
      System.out.println("viscosity_corrected_efficiency_method="
          + viscous.getCorrectedEfficiency(30.0, 1000.0, viscosity));
      // **A fresh chart, because the factors are cached.** `getFullyCorrectedHead` recalculates
      // them only when the viscosity moved by more than 0.1, so on the chart above it would
      // reuse the `flowBEP = 30, headBEP = 88` factors the loop seeded. A chart that has never
      // been corrected recomputes them at its own best-efficiency point, which is what
      // `Pump.run` meets - that is the row the port is held to.
      neqsim.process.equipment.pump.PumpChart fresh =
          simpleChart(new double[] {18.0, 298.15, 1.0, 1.0, 998.0});
      fresh.setUseViscosityCorrection(true);
      System.out.println("viscosity_run_fully_corrected_head="
          + fresh.getFullyCorrectedHead(30.0, 1000.0, 850.0, viscosity));
      System.out.println("viscosity_run_corrected_efficiency="
          + fresh.getCorrectedEfficiency(30.0, 1000.0, viscosity));
      System.out.println();
    }

    runChartRow("run_chart", new double[] {0.3, 1.0, 1.0, 1.0}, 12913.0);
    runChartRow("run_chart_reference_density", new double[] {18.0, 298.15, 1.0, 1.0, 998.0},
        12913.0);

    npshAvailableRow("npsh_available_water_25c_2bara", new String[] {"water"}, 298.15, 2.0);
    npshAvailableRow("npsh_available_water_60c_1p5bara", new String[] {"water"}, 333.15, 1.5);
    npshAvailableRow("npsh_available_water_80c_1p5bara", new String[] {"water"}, 353.15, 1.5);
    npshAvailableRow("npsh_available_hexane", new String[] {"n-hexane"}, 298.15, 2.0);

    Stream inlet = feedMass("inlet", new String[] {"water"}, new double[] {1.0}, 298.15, 2.0,
        39900.0);
    neqsim.process.equipment.pump.Pump pump =
        new neqsim.process.equipment.pump.Pump("p npsh", inlet);
    pump.getPumpChart().setCurves(new double[] {}, new double[] {1000.0, 1500.0},
        new double[][] {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0},
            {15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0}},
        new double[][] {{120.0, 118.0, 115.0, 110.0, 103.0, 94.0, 83.0, 70.0},
            {270.0, 265.5, 258.8, 247.5, 231.8, 211.5, 186.8, 157.5}},
        new double[][] {{60.0, 70.0, 78.0, 82.0, 81.0, 76.0, 68.0, 55.0},
            {62.0, 71.0, 79.0, 83.0, 82.0, 77.0, 69.0, 56.0}});
    pump.getPumpChart().setHeadUnit("meter");
    pump.getPumpChart().setNPSHCurve(new double[][] {{2.0, 2.2, 2.5, 3.0, 3.8, 4.8, 6.2, 8.0},
        {4.5, 4.95, 5.625, 6.75, 8.55, 10.8, 13.95, 18.0}});
    pump.setSpeed(1000.0);
    pump.setCheckNPSH(true);
    pump.setNPSHMargin(1.3);
    pump.run();
    System.out.println("npsh_run_label=npsh_curve_run");
    System.out.println("npsh_run_inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("npsh_run_inlet_density_kg_per_m3="
        + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("npsh_run_head_m=" + pump.getPumpChart().getCorrectedHead(
        inlet.getFlowRate("m3/hr"), 1000.0, inlet.getThermoSystem().getDensity("kg/m3")));
    System.out.println("npsh_run_efficiency_percent="
        + pump.getPumpChart().getEfficiency(inlet.getFlowRate("m3/hr"), 1000.0));
    print("inlet", inlet);
    print("outlet", pump.getOutletStream());
    System.out.println("npsh_run_power_kW=" + pump.getPower("kW"));
    System.out.println("npsh_run_required=" + pump.getNPSHRequired());
    System.out.println("npsh_run_available=" + pump.getNPSHAvailable());
    System.out.println("npsh_run_cavitating=" + pump.isCavitating());
    System.out.println("npsh_run_margin=" + pump.getNPSHMargin());
    System.out.println();
  }

  static neqsim.process.equipment.pump.PumpChart chartTest(double[] conditions) {
    neqsim.process.equipment.pump.PumpChart chart = new neqsim.process.equipment.pump.PumpChart();
    setChartTest(chart, conditions);
    return chart;
  }

  static void setChartTest(neqsim.process.equipment.pump.PumpChart chart, double[] conditions) {
    double[] speed = new double[] {12913, 12298, 11683, 11098, 10453, 9224, 8609, 8200};
    double[][] flow = new double[][] {
        {2789.1285, 3174.0375, 3689.2288, 4179.4503, 4570.2768, 4954.7728, 5246.0329, 5661.0331},
        {2571.1753, 2943.7254, 3440.2675, 3837.4448, 4253.0898, 4668.6643, 4997.1926, 5387.4952},
        {2415.3793, 2763.0706, 3141.7095, 3594.7436, 4047.6467, 4494.1889, 4853.7353, 5138.7858},
        {2247.2043, 2799.7342, 3178.3428, 3656.1551, 4102.778, 4394.1591, 4648.3224, 4840.4998},
        {2072.8397, 2463.9483, 2836.4078, 3202.5266, 3599.6333, 3978.0203, 4257.0022, 4517.345},
        {1835.9552, 2208.455, 2618.1322, 2940.8034, 3244.7852, 3530.1279, 3753.3738, 3895.9746},
        {1711.3386, 1965.8848, 2356.9431, 2685.9247, 3008.5154, 3337.2855, 3591.5092},
        {1636.5807, 2002.8708, 2338.0319, 2642.1245, 2896.4894, 3113.6264, 3274.8764, 3411.2977}};
    double[][] head = new double[][] {
        {80.0375, 78.8934, 76.2142, 71.8678, 67.0062, 60.6061, 53.0499, 39.728},
        {72.2122, 71.8369, 68.9009, 65.8341, 60.7167, 54.702, 47.2749, 35.7471},
        {65.1576, 64.5253, 62.6118, 59.1619, 54.0455, 47.0059, 39.195, 31.6387},
        {58.6154, 56.9627, 54.6647, 50.4462, 44.4322, 38.4144, 32.9084, 28.8109},
        {52.3295, 51.0573, 49.5283, 46.3326, 42.3685, 37.2502, 31.4884, 25.598},
        {40.6578, 39.6416, 37.6008, 34.6603, 30.9503, 27.1116, 23.2713, 20.4546},
        {35.2705, 34.6359, 32.7228, 31.0645, 27.0985, 22.7482, 18.0113},
        {32.192, 31.1756, 29.1329, 26.833, 23.8909, 21.3324, 18.7726, 16.3403}};
    double[][] efficiency = new double[][] {
        {77.2452238409573, 79.4154186459363, 80.737960012489, 80.5229826589649, 79.2210931638144,
            75.4719133864634, 69.6034181197298, 58.7322388482707},
        {77.0107837113504, 79.3069974136389, 80.8941189021135, 80.7190194665918, 79.5313242980328,
            75.5912622896367, 69.6846136362097, 60.0043057990909},
        {77.0043065299874, 79.1690958847856, 80.8038169975675, 80.6543975614197, 78.8532389102705,
            73.6664774270613, 66.2735600426727, 57.671664571658},
        {77.0716623789093, 80.4629750233093, 81.1390811169072, 79.6374242667478, 75.380928428817,
            69.5332969549779, 63.7997587622339, 58.8120614497758},
        {76.9705872525642, 79.8335492585324, 80.9468133671171, 80.5806471927835, 78.0462158225426,
            73.0403707523258, 66.5572286338589, 59.8624822515064},
        {77.5063036680357, 80.2056198362559, 81.0339108025933, 79.6085962687939, 76.3814534404405,
            70.8027503005902, 64.6437367160571, 60.5299349982342},
        {77.8175271586685, 80.065165942218, 81.0631362122632, 79.8955051771299, 76.1983240929369,
            69.289982774309, 60.8567149372229},
        {78.0924334304045, 80.9353551568667, 80.7904437766234, 78.8639325223295, 75.2170936751143,
            70.3105081673411, 65.5507568533569, 61.0391468300337}};
    chart.setCurves(conditions, speed, flow, head, efficiency);
    chart.setHeadUnit("meter");
  }

  static neqsim.process.equipment.pump.PumpChart simpleChart(double[] conditions) {
    neqsim.process.equipment.pump.PumpChart chart = new neqsim.process.equipment.pump.PumpChart();
    chart.setCurves(conditions, new double[] {1000.0},
        new double[][] {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0}},
        new double[][] {{120.0, 118.0, 115.0, 110.0, 103.0, 94.0, 83.0, 70.0}},
        new double[][] {{60.0, 70.0, 78.0, 82.0, 81.0, 76.0, 68.0, 55.0}});
    chart.setHeadUnit("meter");
    return chart;
  }

  static neqsim.process.equipment.pump.PumpChart threePointChart(String headUnit) {
    double factor = headUnit.equals("kJ/kg") ? 0.00981 : 1.0;
    neqsim.process.equipment.pump.PumpChart chart = new neqsim.process.equipment.pump.PumpChart();
    chart.setCurves(new double[] {}, new double[] {500.0}, new double[][] {{30.0, 40.0, 50.0}},
        new double[][] {{80.0 * factor, 70.0 * factor, 60.0 * factor}},
        new double[][] {{80.0, 80.0, 75.0}});
    chart.setHeadUnit(headUnit);
    return chart;
  }

  static neqsim.process.equipment.pump.PumpChart npshChart() {
    neqsim.process.equipment.pump.PumpChart chart = new neqsim.process.equipment.pump.PumpChart();
    chart.setCurves(new double[] {}, new double[] {1000.0, 1500.0},
        new double[][] {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0},
            {15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0}},
        new double[][] {{120.0, 118.0, 115.0, 110.0, 103.0, 94.0, 83.0, 70.0},
            {270.0, 265.5, 258.8, 247.5, 231.8, 211.5, 186.8, 157.5}},
        new double[][] {{60.0, 70.0, 78.0, 82.0, 81.0, 76.0, 68.0, 55.0},
            {62.0, 71.0, 79.0, 83.0, 82.0, 77.0, 69.0, 56.0}});
    chart.setHeadUnit("meter");
    return chart;
  }

  static void chartFitRow(String label, neqsim.process.equipment.pump.PumpChart chart,
      double flow, double speed) {
    System.out.println(label);
    System.out.println("flow_m3_per_hr=" + flow);
    System.out.println("speed_rpm=" + speed);
    System.out.println("head_m=" + chart.getHead(flow, speed));
    System.out.println("efficiency_percent=" + chart.getEfficiency(flow, speed));
    System.out.println();
  }

  static void runChartRow(String label, double[] conditions, double speed) {
    Stream inlet = feedMass("inlet", new String[] {"water"}, new double[] {1.0}, 293.15, 1.0,
        4.0e6);
    neqsim.process.equipment.pump.Pump pump =
        new neqsim.process.equipment.pump.Pump(label, inlet);
    setChartTest((neqsim.process.equipment.pump.PumpChart) pump.getPumpChart(), conditions);
    pump.setSpeed(speed);
    pump.run();

    System.out.println(label);
    System.out.println("inlet_m3_per_hr=" + inlet.getFlowRate("m3/hr"));
    System.out.println("inlet_density_kg_per_m3="
        + inlet.getThermoSystem().getDensity("kg/m3"));
    System.out.println("head_m=" + pump.getPumpChart().getCorrectedHead(
        inlet.getFlowRate("m3/hr"), speed, inlet.getThermoSystem().getDensity("kg/m3")));
    System.out.println("efficiency_percent=" + pump.getPumpChart().getEfficiency(
        inlet.getFlowRate("m3/hr"), speed));
    print("inlet", inlet);
    print("outlet", pump.getOutletStream());
    System.out.println("power_kW=" + pump.getPower("kW"));
    System.out.println();
  }

  static void npshAvailableRow(String label, String[] names, double temperatureK,
      double pressureBara) {
    Stream inlet = feedMass("inlet", names, new double[] {1.0}, temperatureK, pressureBara,
        1000.0);
    neqsim.process.equipment.pump.Pump pump =
        new neqsim.process.equipment.pump.Pump(label, inlet);
    pump.setOutletPressure(5.0);
    pump.setIsentropicEfficiency(0.75);
    pump.setCheckNPSH(true);
    pump.run();

    System.out.println(label);
    System.out.println("npsh_available_m=" + pump.getNPSHAvailable());
    System.out.println("npsh_required_m=" + pump.getNPSHRequired());
    System.out.println("npsh_cavitating=" + pump.isCavitating());
    System.out.println("npsh_margin=" + pump.getNPSHMargin());
    System.out.println();
  }

  static void capacityBinaryRow(String label, double diameterM, Double fsLimit) {
    SystemInterface fluid = new SystemPrEos(300.0, 20.0);
    fluid.addComponent("methane", 0.5);
    fluid.addComponent("n-butane", 0.5);
    fluid.setMixingRule(2);
    Stream inlet = new Stream("column feed", fluid);
    inlet.setFlowRate(1000.0, "kg/hr");
    inlet.run();

    // `binary_methane_butane_4_stages`, the row `process_column.tsv` carries and the port
    // reproduces exactly - so the capacity family is measured on a state that is already an
    // oracle rather than on a column of this probe's own.
    neqsim.process.equipment.distillation.DistillationColumn column =
        new neqsim.process.equipment.distillation.DistillationColumn("col1", 4, true, true);
    column.addFeedStream(inlet, 2);
    column.setCondenserTemperature(-20.0, "C");
    column.setReboilerTemperature(100.0, "C");
    column.setTopPressure(19.0);
    column.setBottomPressure(20.0);
    column.setTemperatureTolerance(1.0e-6);
    column.setMaxNumberOfIterations(200, true);
    column.setInternalDiameter(diameterM);
    if (fsLimit != null) {
      column.setMaxAllowableFsFactor(fsLimit);
    }
    column.run();
    printCapacity(label, "distillation_column", column, diameterM, null, fsLimit, null, "");
  }

  static void capacityAbsorberRow(String label, double diameterM, Double kLimit, Double fsLimit,
      boolean stripper, String perturb) {
    // The stripper's own five-component pair of streams, which is what its class test states.
    String[] names = stripper
        ? new String[] { "methane", "propane", "n-butane", "n-pentane", "n-heptane" }
        : new String[] { "methane", "ethane", "propane", "n-butane", "n-pentane", "n-heptane" };
    double[] gasZ = new double[] { 0.920, 0.040, 0.025, 0.010, 0.005, 0.0 };
    double[] pureHeptane = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
    double[] stripZ = new double[] { 0.9900, 0.0080, 0.0015, 0.0005, 0.0 };
    double[] richZ = new double[] { 0.02, 0.08, 0.12, 0.15, 0.63 };

    // `lean_oil_absorber` and `hydrocarbon_stripper`, the two unpinned rows of the class's own
    // tests, at the tolerances that capture solved them to.
    double bara = stripper ? 12.0 : 15.0;
    double tolerance = 1.0e-4;
    Stream gas = feedMass("gas", names, stripper ? stripZ : gasZ, stripper ? 343.15 : 303.15, bara,
        stripper ? 150.0 : 2000.0);
    Stream solvent = feedMass("solvent", names, stripper ? richZ : pureHeptane,
        stripper ? 343.15 : 293.15, bara, stripper ? 900.0 : 600.0);

    neqsim.process.equipment.absorber.AbsorptionColumn column = stripper
        ? new neqsim.process.equipment.absorber.StrippingColumn(label, 5)
        : new neqsim.process.equipment.absorber.AbsorptionColumn(label, 5);
    if (stripper) {
      ((neqsim.process.equipment.absorber.StrippingColumn) column).addStrippingGasStream(gas);
      ((neqsim.process.equipment.absorber.StrippingColumn) column).addRichLiquidStream(solvent);
    } else {
      column.addGasInStream(gas);
      column.addSolventInStream(solvent);
    }
    column.setTopPressure(bara);
    column.setBottomPressure(bara);
    column.setTemperatureTolerance(tolerance);
    column.setMassBalanceTolerance(5.0e-2);
    column.setEnthalpyBalanceTolerance(5.0e-2);
    column.setMaxNumberOfIterations(80, true);
    column.setSolverType(
        neqsim.process.equipment.distillation.DistillationColumn.SolverType.DIRECT_SUBSTITUTION);
    column.setInternalDiameter(diameterM);
    // **The K setter refuses a non-positive limit** rather than storing it, so a row that states
    // one is stating something the class accepted.
    if (kLimit != null) {
      column.setMaxAllowableGasLoadFactor(kLimit);
    }
    if (fsLimit != null) {
      column.setMaxAllowableFsFactor(fsLimit);
    }
    // Through a `ProcessSystem`, which is how the class's own tests drive it.
    neqsim.process.processmodel.ProcessSystem process = new neqsim.process.processmodel.ProcessSystem();
    process.add(gas);
    process.add(solvent);
    process.add(column);
    process.run();

    // **An outlet reflashed on its own**, which is the class's own state as it reads it: both
    // families take their ingredients from the outlet streams, so moving one is how the branch
    // that reads it is measured rather than read. `StreamInterface.run()` is the flash -
    // `SystemInterface` carries none in this revision.
    StreamInterface reflashed = null;
    if (perturb.equals("liquid")) {
      reflashed = column.getLiquidOutStream();
      reflashed.setTemperature(280.0);
      reflashed.setPressure(bara);
      reflashed.run();
    } else if (perturb.equals("gas")) {
      reflashed = column.getLiquidOutStream();
      reflashed.setTemperature(500.0);
      reflashed.setPressure(1.0);
      reflashed.run();
    } else if (perturb.equals("gas_two_phase")) {
      reflashed = column.getGasOutStream();
      reflashed.setTemperature(235.0);
      reflashed.setPressure(bara);
      reflashed.run();
    }
    printCapacity(label, stripper ? "stripping_column" : "absorption_column", column, diameterM,
        null, fsLimit, kLimit, reflashed == null ? "" : perturb);
  }

  static void printCapacity(String label, String machine,
      neqsim.process.equipment.distillation.DistillationColumn column, double statedDiameterM,
      Double statedColumnDiameterM, Double statedFsLimit, Double statedKLimit,
      String perturbedOutlet) {
    StreamInterface gasOut = column.getGasOutStream();
    StreamInterface liquidOut = column.getLiquidOutStream();
    SystemInterface gas = gasOut == null ? null : gasOut.getThermoSystem();
    SystemInterface liquid = liquidOut == null ? null : liquidOut.getThermoSystem();
    // The two property inits both K getters perform before they read, so the ingredients printed
    // below are the ones those getters see.
    if (gas != null) {
      gas.initPhysicalProperties(neqsim.physicalproperties.PhysicalPropertyType.MASS_DENSITY);
    }
    if (liquid != null) {
      liquid.initPhysicalProperties(neqsim.physicalproperties.PhysicalPropertyType.MASS_DENSITY);
    }

    boolean gasLoad = column instanceof neqsim.process.equipment.absorber.AbsorptionColumn;
    System.out.println(label);
    System.out.println("machine=" + machine);
    System.out.println("trays=" + column.getNumberOfTrays());
    System.out.println("internal_diameter_stated_m=" + statedDiameterM);
    System.out.println("internal_diameter_m=" + column.getInternalDiameter());
    System.out.println("column_diameter_stated_m=" + statedColumnDiameterM);
    System.out.println("max_allowable_fs_factor_stated=" + statedFsLimit);
    System.out.println("max_allowable_fs_factor=" + column.getMaxAllowableFsFactor());
    if (gasLoad) {
      neqsim.process.equipment.absorber.AbsorptionColumn absorber =
          (neqsim.process.equipment.absorber.AbsorptionColumn) column;
      System.out.println("max_allowable_gas_load_factor_stated=" + statedKLimit);
      System.out.println("max_allowable_gas_load_factor=" + absorber.getMaxAllowableGasLoadFactor());
    }
    System.out.println("reflashed_outlet=" + perturbedOutlet);
    System.out.println("solved=" + column.solved());
    System.out.println("iterations=" + column.getLastIterationCount());
    System.out.println("status=" + column.getLastSolveStatus());
    System.out.println("temperature_residual=" + column.getLastTemperatureResidual());
    printCapacitySide("gas_out", gas);
    printCapacitySide("liquid_out", liquid);
    System.out.println("gas_superficial_velocity_m_per_s=" + superficialVelocity(column, gas));
    // ---- The class's own answers, which the port is held to.
    System.out.println("fs_factor=" + column.getFsFactor());
    System.out.println("fs_factor_utilization=" + column.getFsFactorUtilization());
    System.out.println("fs_factor_within_design_limit=" + column.isFsFactorWithinDesignLimit());
    System.out.println("minimum_diameter_for_fs_limit_m=" + column.getMinimumDiameterForFsLimit());
    if (gasLoad) {
      neqsim.process.equipment.absorber.AbsorptionColumn absorber =
          (neqsim.process.equipment.absorber.AbsorptionColumn) column;
      System.out.println("gas_load_factor=" + absorber.getGasLoadFactor());
      System.out.println("gas_load_factor_utilization=" + absorber.getGasLoadFactorUtilization());
      System.out.println("gas_load_factor_within_design_limit="
          + absorber.isGasLoadFactorWithinDesignLimit());
      System.out.println("minimum_diameter_for_gas_load_limit_m="
          + absorber.getMinimumDiameterForGasLoadLimit());
    }
    System.out.println();
  }

  /// The two outlet states the capacity families read, as ingredients rather than as answers:
  /// the volumetric flow both velocities divide, and the two densities the Fs and K families each
  /// take their own copy of.
  static void printCapacitySide(String port, SystemInterface side) {
    if (side == null) {
      System.out.println(port + "=null");
      return;
    }
    double n = side.getFlowRate("mol/sec");
    double mass = n * side.getMolarMass();
    double density = side.getDensity("kg/m3");
    System.out.println(port + "_T_K=" + side.getTemperature());
    System.out.println(port + "_P_bara=" + side.getPressure());
    System.out.println(port + "_phases=" + side.getNumberOfPhases());
    System.out.println(port + "_phase0_type=" + side.getPhase(0).getPhaseTypeName());
    System.out.println(port + "_n_mol_per_s=" + n);
    System.out.println(port + "_molar_mass_kg_per_mol=" + side.getMolarMass());
    System.out.println(port + "_mass_flow_kg_per_s=" + mass);
    System.out.println(port + "_density_kg_per_m3=" + density);
    System.out.println(port + "_phase0_density_kg_per_m3="
        + side.getPhase(0).getPhysicalProperties().getDensity());
    System.out.println(port + "_volumetric_flow_m3_per_s=" + side.getFlowRate("m3/sec"));
    System.out.println(port + "_mass_over_density_m3_per_s=" + mass / density);
  }

  /// `pi·D^2/4` under the gas outlet's volumetric flow - the same quantity
  /// `getGasSuperficialVelocity` returns, printed from the ingredients so a row still carries a
  /// velocity if the getter is what moved.
  static double superficialVelocity(
      neqsim.process.equipment.distillation.DistillationColumn column, SystemInterface gas) {
    if (gas == null) {
      return 0.0;
    }
    return gas.getFlowRate("m3/sec") / (Math.PI * column.getInternalDiameter()
        * column.getInternalDiameter() / 4.0);
  }

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

  // === process.hem_release =====================================================

  static void hemReleaseRows() {
    // PR rows are the oracle: azoth's databank route resolves Peng-Robinson.
    hemRow("methane_50bara_choked", false, new String[] {"methane"}, new double[] {1.0}, 300.0,
        50.0, 0.01, 0.62, 101325.0);
    hemRow("methane_dilute_0_1bara", false, new String[] {"methane"}, new double[] {1.0}, 300.0,
        0.1, 0.01, 0.62, 500.0);
    for (double ratio : new double[] {0.05, 0.30, 0.70, 0.90, 0.98}) {
      hemRow("methane_dilute_back_" + ratio, false, new String[] {"methane"}, new double[] {1.0},
          300.0, 0.1, 0.01, 0.62, ratio * 10000.0);
    }
    hemRow("methane_5bara_unchoked", false, new String[] {"methane"}, new double[] {1.0}, 300.0,
        5.0, 0.01, 0.62, 490000.0);
    hemRow("methane_5bara_no_flow_at", false, new String[] {"methane"}, new double[] {1.0}, 300.0,
        5.0, 0.01, 0.62, 500000.0);
    hemRow("methane_5bara_no_flow_above", false, new String[] {"methane"}, new double[] {1.0},
        300.0, 5.0, 0.01, 0.62, 600000.0);
    hemRow("methane_ethane_5bara", false, new String[] {"methane", "ethane"},
        new double[] {1.0, 0.1}, 300.0, 5.0, 0.01, 0.62, 100000.0);
    hemRow("propane_20bara", false, new String[] {"propane"}, new double[] {1.0}, 300.0, 20.0,
        0.01, 0.62, 300000.0);
    hemRow("dense_co2_120bara", false, new String[] {"CO2"}, new double[] {1.0}, 310.0, 120.0,
        0.01, 0.62, 110.0e5);
    hemRow("cold_co2_unsupported", false, new String[] {"CO2"}, new double[] {1.0}, 210.0, 20.0,
        0.01, 0.62, 100000.0);
    hemRow("methane_scale_d01_cd050", false, new String[] {"methane"}, new double[] {1.0}, 300.0,
        5.0, 0.01, 0.5, 100000.0);
    hemRow("methane_scale_d02_cd100", false, new String[] {"methane"}, new double[] {1.0}, 300.0,
        5.0, 0.02, 1.0, 100000.0);
    hemRow("propane_butane_20bara", false, new String[] {"propane", "n-butane"},
        new double[] {0.8, 0.2}, 300.0, 20.0, 0.01, 0.62, 300000.0);

    // SRK rows are the pinned tests' own EOS, and evidence beside the PR oracle.
    hemRow("srk_methane_50bara_choked", true, new String[] {"methane"}, new double[] {1.0}, 300.0,
        50.0, 0.01, 0.62, 101325.0);
    hemRow("srk_propane_20bara", true, new String[] {"propane"}, new double[] {1.0}, 300.0, 20.0,
        0.01, 0.62, 300000.0);
    hemRow("srk_cold_co2_unsupported", true, new String[] {"CO2"}, new double[] {1.0}, 210.0,
        20.0, 0.01, 0.62, 100000.0);
    hemRow("srk_propane_butane_20bara", true, new String[] {"propane", "n-butane"},
        new double[] {0.8, 0.2}, 300.0, 20.0, 0.01, 0.62, 300000.0);

    // The two states no azoth input can express: an empty fluid and the hydrate flag.
    hemLifecycleRow("pr_empty_fluid", new SystemPrEos(300.0, 5.0), 0.01, 0.62, 100000.0);
    SystemInterface hydrate = hemFluid(false, new String[] {"methane"}, new double[] {1.0}, 300.0,
        5.0);
    hydrate.setHydrateCheck(true);
    hemLifecycleRow("pr_hydrate_flag", hydrate, 0.01, 0.62, 100000.0);
  }

  static SystemInterface hemFluid(boolean srk, String[] names, double[] amounts, double temperature,
      double pressureBar) {
    SystemInterface fluid =
        srk ? new SystemSrkEos(temperature, pressureBar) : new SystemPrEos(temperature, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], amounts[i]);
    }
    if (srk) {
      fluid.setMixingRule("classic");
    } else {
      fluid.setMixingRule(2);
    }
    new neqsim.thermodynamicoperations.ThermodynamicOperations(fluid).TPflash();
    fluid.init(3);
    return fluid;
  }

  static void hemRow(String label, boolean srk, String[] names, double[] amounts,
      double temperature, double pressureBar, double diameter, double coefficient, double back) {
    SystemInterface fluid = hemFluid(srk, names, amounts, temperature, pressureBar);
    ReleaseFlowRequest request = new ReleaseFlowRequest(fluid, diameter, coefficient, back);
    hemPrint(label, srk ? "srk" : "pr", fluid, request);
  }

  static void hemLifecycleRow(String label, SystemInterface fluid, double diameter,
      double coefficient, double back) {
    ReleaseFlowRequest request = new ReleaseFlowRequest(fluid, diameter, coefficient, back);
    hemPrint(label, "pr", fluid, request);
  }

  static void hemPrint(String label, String eos, SystemInterface fluid,
      ReleaseFlowRequest request) {
    ReleaseFlowResult result = new HomogeneousEquilibriumReleaseModel().calculate(request);
    System.out.println(label);
    System.out.println("eos=" + eos);
    System.out.println("input_t_k=" + fluid.getTemperature());
    System.out.println("input_p_bara=" + fluid.getPressure());
    double[] overall = fluid.getMolarComposition();
    StringBuilder z = new StringBuilder();
    for (int i = 0; i < overall.length; i++) {
      z.append(i == 0 ? "" : ",")
          .append(fluid.getPhase(0).getComponent(i).getName()).append(":").append(overall[i]);
    }
    System.out.println("input_z=" + z);
    System.out.println("diameter_m=" + request.getDiameterM());
    System.out.println("discharge_coefficient=" + request.getDischargeCoefficient());
    System.out.println("back_pressure_pa=" + request.getBackPressurePa());
    System.out.println("status=" + result.getStatus());
    System.out.println("model_id=" + result.getModelId());
    System.out.println("mass_flow_kg_s=" + (result.isUsable() ? result.getMassFlowRateKgS() : ""));
    System.out.println("choked=" + result.isChoked());
    System.out.println("throat_sound_speed_m_per_s="
        + (result.getThroatSoundSpeedMs() == null ? "" : result.getThroatSoundSpeedMs()));
    int d = 0;
    for (ReleaseFlowResult.Diagnostic diagnostic : result.getDiagnostics()) {
      System.out.println("diagnostic_" + d + "_code=" + diagnostic.getCode());
      System.out.println("diagnostic_" + d + "_message=" + diagnostic.getMessage());
      d++;
    }
    int i = 0;
    for (Map.Entry<ReleaseFlowResult.Station, ReleaseState> entry : result.getStations().entrySet()) {
      ReleaseState state = entry.getValue();
      String prefix = "station_" + i + "_";
      System.out.println(prefix + "name=" + entry.getKey().name());
      System.out.println(prefix + "p_pa=" + state.getPressurePa());
      System.out.println(prefix + "t_k=" + state.getTemperatureK());
      System.out.println(prefix + "density_kg_m3=" + state.getDensityKgM3());
      System.out.println(prefix + "h_j_per_kg=" + state.getEnthalpyJkg());
      System.out.println(prefix + "s_j_per_kgk=" + state.getEntropyJkgK());
      System.out.println(prefix + "velocity_m_per_s=" + state.getVelocityMs());
      System.out.println(prefix + "mass_flux_kg_m2_s=" + state.getMassFluxKgM2s());
      System.out.println(prefix + "gas_mass_fraction=" + state.getGasMassFraction());
      System.out.println(prefix + "component_mole_fractions="
          + hemFractions(state.getComponentMoleFractions()));
      System.out.println(prefix + "component_mass_fractions="
          + hemFractions(state.getComponentMassFractions()));
      System.out.println(prefix + "phase_mass_fractions="
          + hemFractions(state.getPhaseMassFractions()));
      i++;
    }
    System.out.println();
  }

  static String hemFractions(Map<String, Double> values) {
    StringBuilder out = new StringBuilder();
    for (Map.Entry<String, Double> entry : values.entrySet()) {
      if (out.length() > 0) {
        out.append(",");
      }
      out.append(entry.getKey()).append(":").append(entry.getValue());
    }
    return out.toString();
  }

  // === process.leak_blowdown ===================================================

  static void leakBlowdownRows() {
    leakBlowdownRow("pr_methane_50bara_10s_1s", false, new String[] {"methane"},
        new double[] {1.0}, 300.0, 50.0, 1.0, 0.01, 0.62, 101325.0, 10.0, 1.0);
    leakBlowdownRow("pr_methane_50bara_60s_1s", false, new String[] {"methane"},
        new double[] {1.0}, 300.0, 50.0, 1.0, 0.01, 0.62, 101325.0, 60.0, 1.0);
    leakBlowdownRow("pr_methane_export_30s_1s", false, new String[] {"methane"},
        new double[] {1.0}, 300.0, 50.0, 5.0, 0.02, 0.62, 101325.0, 30.0, 1.0);
    leakBlowdownRow("pr_propane_11bara_10s_1s", false, new String[] {"propane"},
        new double[] {1.0}, 300.0, 11.0, 0.01, 0.01, 0.62, 101325.0, 10.0, 1.0);
    leakBlowdownRow("pr_methane_no_forward_flow_3s", false, new String[] {"methane"},
        new double[] {1.0}, 300.0, 5.0, 1.0, 0.01, 0.62, 600000.0, 3.0, 1.0);
    leakBlowdownRow("pr_methane_2bara_subsonic_1s", false, new String[] {"methane"},
        new double[] {1.0}, 300.0, 2.0, 1.0, 0.01, 0.62, 190000.0, 1.0, 1.0);
    leakBlowdownRow("pr_methane_ethane_50bara_5s", false,
        new String[] {"methane", "ethane"}, new double[] {0.9, 0.1}, 300.0, 50.0, 1.0, 0.01, 0.62,
        101325.0, 5.0, 1.0);
    leakBlowdownRow("pr_single_step_depletion_100s", false, new String[] {"methane"},
        new double[] {1.0}, 300.0, 50.0, 1.0, 0.01, 0.62, 101325.0, 100.0, 100.0);
    // SRK rows are the pinned tests' own EOS, and evidence beside the PR oracle.
    leakBlowdownRow("srk_methane_50bara_60s_1s", true, new String[] {"methane"},
        new double[] {1.0}, 300.0, 50.0, 1.0, 0.01, 0.62, 101325.0, 60.0, 1.0);
    leakBlowdownRow("srk_methane_export_30s_1s", true, new String[] {"methane"},
        new double[] {1.0}, 300.0, 50.0, 5.0, 0.02, 0.62, 101325.0, 30.0, 1.0);
    leakBlowdownRow("srk_propane_11bara_10s_1s", true, new String[] {"propane"},
        new double[] {1.0}, 300.0, 11.0, 0.01, 0.01, 0.62, 101325.0, 10.0, 1.0);
  }

  static SystemInterface leakFluid(boolean srk, String[] names, double[] amounts,
      double temperature, double pressureBar) {
    SystemInterface fluid =
        srk ? new SystemSrkEos(temperature, pressureBar) : new SystemPrEos(temperature, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], amounts[i]);
    }
    if (srk) {
      fluid.setMixingRule("classic");
    } else {
      fluid.setMixingRule(2);
    }
    new neqsim.thermodynamicoperations.ThermodynamicOperations(fluid).TPflash();
    fluid.init(3);
    fluid.initPhysicalProperties();
    return fluid;
  }

  @SuppressWarnings({"checkstyle:MethodLength", "checkstyle:ParameterNumber"})
  static void leakBlowdownRow(String label, boolean srk, String[] names, double[] amounts,
      double temperature, double pressureBar, double volume, double diameter, double coefficient,
      double back, double duration, double step) {
    SystemInterface fluid = leakFluid(srk, names, amounts, temperature, pressureBar);
    LeakModel leak = LeakModel.builder().fluid(fluid).holeDiameter(diameter).vesselVolume(volume)
        .dischargeCoefficient(coefficient).backPressure(back).scenarioName(label).build();
    SourceTermResult result = leak.calculateSourceTerm(duration, step);
    System.out.println(label);
    System.out.println("eos=" + (srk ? "srk" : "pr"));
    System.out.println("input_t_k=" + fluid.getTemperature());
    System.out.println("input_p_bara=" + fluid.getPressure());
    double[] overall = fluid.getMolarComposition();
    StringBuilder z = new StringBuilder();
    for (int i = 0; i < overall.length; i++) {
      z.append(i == 0 ? "" : ",")
          .append(fluid.getPhase(0).getComponent(i).getName()).append(":").append(overall[i]);
    }
    System.out.println("input_z=" + z);
    System.out.println("vessel_volume_m3=" + volume);
    System.out.println("hole_diameter_m=" + diameter);
    System.out.println("discharge_coefficient=" + coefficient);
    System.out.println("back_pressure_pa=" + back);
    System.out.println("duration_s=" + duration);
    System.out.println("time_step_s=" + step);
    System.out.println("declared_points=" + result.getNumberOfPoints());
    System.out.println("initial_density_kg_m3=" + fluid.getDensity("kg/m3"));
    System.out.println("initial_gamma=" + fluid.getGamma());
    System.out.println("initial_beta=" + fluid.getBeta());
    System.out.println("initial_phase_count=" + fluid.getNumberOfPhases());
    System.out.println("initial_mass_flow_kg_s=" + leak.calculateMassFlowRate(fluid));
    System.out.println("initial_jet_velocity_m_per_s=" + leak.calculateJetVelocity(fluid));
    System.out.println("initial_jet_momentum_n=" + leak.calculateJetMomentum(fluid));
    System.out.println("initial_droplet_smd_m=" + leak.calculateDropletSMD(fluid));
    System.out.println("total_mass_released_kg=" + result.getTotalMassReleased());
    System.out.println("peak_mass_flow_kg_s=" + result.getPeakMassFlowRate());
    System.out.println("time_to_empty_s=" + result.getTimeToEmpty());
    double[] time = result.getTime();
    double[] rate = result.getMassFlowRate();
    double[] temp = result.getTemperature();
    double[] pressure = result.getPressure();
    double[] beta = result.getVaporMassFraction();
    double[] velocity = result.getJetVelocity();
    double[] momentum = result.getJetMomentum();
    double[] smd = result.getLiquidDropletSMD();
    int recorded = 0;
    for (int i = 0; i < time.length; i++) {
      if (temp[i] != 0.0) {
        recorded = i + 1;
      }
    }
    System.out.println("recorded_points=" + recorded);
    double area = Math.PI * diameter * diameter / 4.0;
    for (int i = 0; i < recorded; i++) {
      String prefix = "point_" + i + "_";
      System.out.println(prefix + "t_s=" + time[i]);
      System.out.println(prefix + "mass_flow_kg_s=" + rate[i]);
      System.out.println(prefix + "temperature_k=" + temp[i]);
      System.out.println(prefix + "pressure_pa=" + pressure[i]);
      System.out.println(prefix + "beta=" + beta[i]);
      System.out.println(prefix + "jet_velocity_m_per_s=" + velocity[i]);
      System.out.println(prefix + "jet_momentum_n=" + momentum[i]);
      System.out.println(prefix + "droplet_smd_m=" + smd[i]);
      SystemInterface probe = fluid.clone();
      probe.setTemperature(temp[i]);
      probe.setPressure(pressure[i] / 1e5);
      new neqsim.thermodynamicoperations.ThermodynamicOperations(probe).TPflash();
      probe.init(3);
      probe.initPhysicalProperties();
      System.out.println(prefix + "gamma=" + probe.getGamma());
      System.out.println(prefix + "phase_count=" + probe.getNumberOfPhases());
      StringBuilder types = new StringBuilder();
      StringBuilder densities = new StringBuilder();
      for (int j = 0; j < probe.getNumberOfPhases(); j++) {
        types.append(j == 0 ? "" : ",").append(probe.getPhase(j).getType().name());
        densities.append(j == 0 ? "" : ",").append(probe.getPhase(j).getDensity("kg/m3"));
      }
      System.out.println(prefix + "phase_types=" + types);
      System.out.println(prefix + "phase_densities_kg_m3=" + densities);
      System.out.println(prefix + "fresh_density_kg_m3=" + probe.getDensity("kg/m3"));
      double implied = velocity[i] > 0.0 ? rate[i] / (velocity[i] * area * coefficient) : Double.NaN;
      System.out.println(prefix + "implied_velocity_density_kg_m3=" + implied);
    }
    System.out.println();
  }

  // === process.gas_dispersion, process.gaussian_plume, process.probit =========

  static BoundaryConditions dispersionWeather(double ambientC, double wind, char stability,
      boolean offshore, double roughness) {
    return BoundaryConditions.builder().ambientTemperature(ambientC, "C").windSpeed(wind)
        .pasquillStabilityClass(stability).isOffshore(offshore).surfaceRoughness(roughness)
        .build();
  }

  static void gasDispersionRows() {
    // PR rows are the oracle: azoth's databank route resolves Peng-Robinson. The class's own
    // tests (`GasDispersionAnalyzerTest`), re-cased on PR.
    dispersionRow("pr_methane_gaussian", false, new String[] {"methane"}, new double[] {1.0},
        293.15, 20.0, 1.0, 0.0, null, null, 0.5, null, 0.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("pr_methane_sts0131_integral", false, new String[] {"methane"},
        new double[] {1.0}, 293.15, 20.0, 1.0, 0.0, null, null, 0.20, null, 0.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("pr_methane_sts0131_cfd", false, new String[] {"methane"}, new double[] {1.0},
        293.15, 20.0, 1.0, 0.0, null, null, 0.50, null, 0.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("pr_co2_dense", false, new String[] {"CO2"}, new double[] {1.0}, 293.15, 20.0,
        5.0, 0.0, null, null, 0.5, null, 0.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("pr_h2s_toxic_gaussian", false, new String[] {"methane", "H2S"},
        new double[] {0.98, 0.02}, 293.15, 20.0, 2.0, 0.0, null, null, 0.5, "H2S", 100.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("pr_methane_manual_lfl", false, new String[] {"methane"}, new double[] {1.0},
        293.15, 20.0, 1.0, 0.0, null, 0.05, 0.5, null, 0.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("pr_methane_forced_heavy", false, new String[] {"methane"},
        new double[] {1.0}, 293.15, 20.0, 1.0, 0.0,
        GasDispersionAnalyzer.ModelSelection.HEAVY_GAS, null, 0.5, null, 0.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("pr_co2_forced_gaussian", false, new String[] {"CO2"}, new double[] {1.0},
        293.15, 20.0, 5.0, 0.0, GasDispersionAnalyzer.ModelSelection.GAUSSIAN_PLUME, null, 0.5,
        null, 0.0, null, dispersionWeather(15.0, 5.0, 'D', false, 0.1));

    // `DispersionBenchmarkValidationTest`, re-cased on PR.
    dispersionRow("pr_benchmark_methane", false, new String[] {"methane"}, new double[] {1.0},
        288.15, 10.0, 1.0, 0.0, GasDispersionAnalyzer.ModelSelection.GAUSSIAN_PLUME, null, 0.5,
        null, 0.0, null, dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("pr_benchmark_co2_toxic", false, new String[] {"CO2"}, new double[] {1.0},
        288.15, 10.0, 10.0, 0.0, null, null, 0.5, "CO2", 50000.0, null,
        dispersionWeather(10.0, 2.0, 'F', false, 0.1));
    dispersionRow("pr_benchmark_propane", false, new String[] {"propane"}, new double[] {1.0},
        288.15, 10.0, 1.0, 0.0, null, null, 0.5, null, 0.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("pr_benchmark_h2s_forced_gaussian", false, new String[] {"methane", "H2S"},
        new double[] {0.90, 0.10}, 288.15, 10.0, 1.0, 0.0,
        GasDispersionAnalyzer.ModelSelection.GAUSSIAN_PLUME, null, 0.5, "H2S", 100.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));

    // The source-term row: the builder consumes the leak's peak rate.
    SystemInterface sourceFluid = hemFluid(false, new String[] {"methane"}, new double[] {1.0},
        293.15, 50.0);
    LeakModel leak = LeakModel.builder().fluid(sourceFluid).holeDiameter(0.01).vesselVolume(1.0)
        .dischargeCoefficient(0.62).backPressure(101325.0).scenarioName("source").build();
    SourceTermResult sourceTerm = leak.calculateSourceTerm(5.0, 1.0);
    dispersionRow("pr_source_term_rate", false, new String[] {"methane"}, new double[] {1.0},
        293.15, 50.0, Double.NaN, 0.0, null, null, 0.5, null, 0.0, sourceTerm,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));

    // SRK rows are the pinned tests' own EOS, and evidence beside the PR oracle.
    dispersionRow("srk_methane_gaussian", true, new String[] {"methane"}, new double[] {1.0},
        293.15, 20.0, 1.0, 0.0, null, null, 0.5, null, 0.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("srk_co2_dense", true, new String[] {"CO2"}, new double[] {1.0}, 293.15, 20.0,
        5.0, 0.0, null, null, 0.5, null, 0.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("srk_h2s_toxic_gaussian", true, new String[] {"methane", "H2S"},
        new double[] {0.98, 0.02}, 293.15, 20.0, 2.0, 0.0, null, null, 0.5, "H2S", 100.0, null,
        dispersionWeather(15.0, 5.0, 'D', false, 0.1));
    dispersionRow("srk_benchmark_co2_toxic", true, new String[] {"CO2"}, new double[] {1.0},
        288.15, 10.0, 10.0, 0.0, null, null, 0.5, "CO2", 50000.0, null,
        dispersionWeather(10.0, 2.0, 'F', false, 0.1));

    // The state the class refuses at construction, printed as evidence.
    System.out.println("pr_zero_rate_refused");
    System.out.println("eos=pr");
    System.out.println("rate_input_kg_s=0.0");
    try {
      GasDispersionAnalyzer.builder().scenarioName("zero").fluid(sourceFluid).massReleaseRate(0.0)
          .boundaryConditions(dispersionWeather(15.0, 5.0, 'D', false, 0.1)).build().analyze();
      System.out.println("refusal=none");
    } catch (RuntimeException ex) {
      System.out.println("refusal=" + ex.getClass().getSimpleName() + ": " + ex.getMessage());
    }
    System.out.println();
  }

  @SuppressWarnings({"checkstyle:MethodLength", "checkstyle:ParameterNumber"})
  static void dispersionRow(String label, boolean srk, String[] names, double[] amounts,
      double temperature, double pressureBar, double rate, double releaseHeight,
      GasDispersionAnalyzer.ModelSelection selection, Double manualLfl, double endpointFraction,
      String toxicComponent, double toxicThresholdPpm, SourceTermResult sourceTerm,
      BoundaryConditions conditions) {
    SystemInterface fluid = hemFluid(srk, names, amounts, temperature, pressureBar);
    GasDispersionAnalyzer.Builder builder = GasDispersionAnalyzer.builder().scenarioName(label)
        .fluid(fluid).boundaryConditions(conditions).releaseHeight(releaseHeight)
        .flammableEndpointFraction(endpointFraction);
    if (!Double.isNaN(rate)) {
      builder.massReleaseRate(rate);
    }
    if (selection != null) {
      builder.modelSelection(selection);
    }
    if (manualLfl != null) {
      builder.lowerFlammableLimit(manualLfl);
    }
    if (toxicComponent != null) {
      builder.toxicEndpoint(toxicComponent, toxicThresholdPpm);
    }
    if (sourceTerm != null) {
      builder.sourceTerm(sourceTerm);
    }
    GasDispersionResult result = builder.build().analyze();
    System.out.println(label);
    System.out.println("eos=" + (srk ? "srk" : "pr"));
    System.out.println("source_t_k=" + fluid.getTemperature());
    System.out.println("source_p_bara=" + fluid.getPressure());
    double[] overall = fluid.getMolarComposition();
    StringBuilder z = new StringBuilder();
    for (int i = 0; i < overall.length; i++) {
      z.append(i == 0 ? "" : ",")
          .append(fluid.getPhase(0).getComponent(i).getName()).append(":").append(overall[i]);
    }
    System.out.println("source_z=" + z);
    System.out.println("rate_input_kg_s=" + (Double.isNaN(rate) ? "" : rate));
    if (sourceTerm != null) {
      System.out.println("source_term_peak_kg_s=" + sourceTerm.getPeakMassFlowRate());
    }
    System.out.println("release_height_m=" + releaseHeight);
    System.out.println("model_selection=" + (selection == null ? "auto" : selection.name()));
    System.out.println("manual_lfl=" + (manualLfl == null ? "" : manualLfl));
    System.out.println("flammable_endpoint_fraction=" + endpointFraction);
    System.out.println("toxic_component=" + (toxicComponent == null ? "" : toxicComponent));
    System.out.println("toxic_threshold_ppm=" + (toxicComponent == null ? "" : toxicThresholdPpm));
    System.out.println("ambient_t_k=" + conditions.getAmbientTemperature());
    System.out.println("ambient_p_pa=" + conditions.getAtmosphericPressure());
    System.out.println("conditions_wind_m_per_s=" + conditions.getWindSpeed());
    System.out.println("conditions_stability=" + conditions.getPasquillStabilityClass());
    System.out.println("conditions_offshore=" + conditions.isOffshore());
    System.out.println("conditions_roughness_m=" + conditions.getSurfaceRoughness());
    System.out.println("selected_model=" + result.getSelectedModel());
    System.out.println("result_mass_release_rate_kg_s=" + result.getMassReleaseRateKgPerS());
    System.out.println(
        "flammable_mass_release_rate_kg_s=" + result.getFlammableMassReleaseRateKgPerS());
    System.out.println("source_density_kg_m3=" + result.getSourceDensityKgPerM3());
    System.out.println("air_density_kg_m3=" + result.getAirDensityKgPerM3());
    System.out.println("fuel_mole_fraction=" + result.getFuelMoleFraction());
    System.out.println("fuel_mass_fraction=" + result.getFuelMassFraction());
    System.out.println(
        "lfl_volume_fraction=" + result.getLowerFlammableLimitVolumeFraction());
    System.out.println(
        "endpoint_fraction_of_lfl=" + result.getFlammableEndpointFractionOfLfl());
    System.out.println(
        "distance_to_flammable_endpoint_m=" + result.getDistanceToFlammableEndpointM());
    System.out.println("distance_to_half_lfl_m=" + result.getDistanceToHalfLflM());
    System.out.println("distance_to_lfl_m=" + result.getDistanceToLflM());
    System.out.println("flammable_cloud_volume_m3=" + result.getFlammableCloudVolumeM3());
    System.out.println("toxic_component_name=" + result.getToxicComponentName());
    System.out.println("result_toxic_threshold_ppm=" + result.getToxicThresholdPpm());
    System.out.println("toxic_distance_m=" + result.getToxicDistanceM());
    System.out.println("result_wind_speed_m_per_s=" + result.getWindSpeedMPerS());
    System.out.println("result_stability_class=" + result.getStabilityClass());
    System.out.println("screening_basis=" + result.getScreeningBasis());
    System.out.println("has_flammable_cloud=" + result.hasFlammableCloud());
    System.out.println("has_toxic_endpoint=" + result.hasToxicEndpoint());
    System.out.println();
  }

  static void gaussianPlumeRows() {
    plumeRow("rural_d_100m", 1.0, 10.0, 5.0, GaussianPlume.Stability.D, GaussianPlume.Terrain.RURAL,
        100.0, 0.0, 0.0, 1.0e-3);
    plumeRow("rural_d_500m", 1.0, 10.0, 5.0, GaussianPlume.Stability.D, GaussianPlume.Terrain.RURAL,
        500.0, 0.0, 0.0, 1.0e-3);
    plumeRow("rural_d_1000m", 1.0, 10.0, 5.0, GaussianPlume.Stability.D,
        GaussianPlume.Terrain.RURAL, 1000.0, 0.0, 0.0, 1.0e-3);
    plumeRow("urban_b_500m", 1.0, 5.0, 3.0, GaussianPlume.Stability.B,
        GaussianPlume.Terrain.URBAN, 500.0, 0.0, 0.0, 1.0e-4);
    plumeRow("rural_d_offcentre", 1.0, 10.0, 5.0, GaussianPlume.Stability.D,
        GaussianPlume.Terrain.RURAL, 500.0, 25.0, 2.0, 1.0e-4);
    // The upstream factor-of-ten pair: ground-level source, 1e-3 threshold.
    plumeRow("rural_d_crossing", 10.0, 0.0, 4.0, GaussianPlume.Stability.D,
        GaussianPlume.Terrain.RURAL, 100.0, 0.0, 0.0, 1.0e-3);
    // The scan's two ends: never reached, and never fallen below (the 100000 m cap).
    plumeRow("rural_d_never_above", 1.0, 0.0, 5.0, GaussianPlume.Stability.D,
        GaussianPlume.Terrain.RURAL, 100.0, 0.0, 0.0, 1.0e6);
    plumeRow("rural_d_never_below", 10.0, 0.0, 4.0, GaussianPlume.Stability.D,
        GaussianPlume.Terrain.RURAL, 100.0, 0.0, 0.0, 1.0e-30);
    plumeRow("zero_emission", 0.0, 0.0, 5.0, GaussianPlume.Stability.D,
        GaussianPlume.Terrain.RURAL, 100.0, 0.0, 0.0, 1.0e-6);
  }

  static void plumeRow(String label, double emission, double stackHeight, double wind,
      GaussianPlume.Stability stability, GaussianPlume.Terrain terrain, double distance,
      double crosswind, double height, double threshold) {
    GaussianPlume plume =
        new GaussianPlume(emission, stackHeight, wind, stability, terrain);
    System.out.println(label);
    System.out.println("emission_rate_kg_s=" + emission);
    System.out.println("stack_height_m=" + stackHeight);
    System.out.println("wind_speed_m_per_s=" + wind);
    System.out.println("stability=" + stability.name());
    System.out.println("terrain=" + terrain.name());
    System.out.println("distance_m=" + distance);
    System.out.println("crosswind_m=" + crosswind);
    System.out.println("height_m=" + height);
    System.out.println("threshold_kg_m3=" + threshold);
    System.out.println("sigma_y_m=" + plume.sigmaY(distance));
    System.out.println("sigma_z_m=" + plume.sigmaZ(distance));
    System.out.println("concentration_kg_m3=" + plume.concentration(distance, crosswind, height));
    System.out.println(
        "centerline_ground_concentration_kg_m3=" + plume.centerlineGroundConcentration(distance));
    System.out.println("distance_to_threshold_m=" + plume.distanceToConcentration(threshold));
    System.out.println();
  }

  static void probitRows() {
    probitRow("custom_five_at_dose_one", 5.0, 0.0, 1.0, 1.0, Double.NaN);
    probitRow("custom_five_from_probit", 5.0, 0.0, 1.0, 1.0, 5.0);
    probitRow("h2s_low", -31.42, 3.008, 1.43, Math.pow(100.0, 1.43) * 10.0, Double.NaN);
    probitRow("h2s_high", -31.42, 3.008, 1.43, Math.pow(500.0, 1.43) * 30.0, Double.NaN);
    probitRow("thermal_low", -14.9, 2.56, 1.0, 60.0 * Math.pow(5000.0, 4.0 / 3.0) / 1.0e4,
        Double.NaN);
    probitRow("thermal_high", -14.9, 2.56, 1.0, 60.0 * Math.pow(35000.0, 4.0 / 3.0) / 1.0e4,
        Double.NaN);
    probitRow("chlorine_sample", -17.10, 1.69, 2.75, Math.pow(10.0, 2.75), Double.NaN);
    probitRow("ammonia_sample", -35.9, 1.85, 2.0, Math.pow(100.0, 2.0) * 5.0, Double.NaN);
    probitRow("blast_sample", -77.1, 6.91, 1.0, 20000.0, Double.NaN);
    // The class's non-positive-dose branch, printed as evidence.
    probitRow("zero_dose_evidence", 5.0, 0.0, 1.0, 0.0, Double.NaN);
  }

  static void probitRow(String label, double a, double b, double n, double dose, double y) {
    ProbitModel model = new ProbitModel(a, b, n);
    System.out.println(label);
    System.out.println("a=" + a);
    System.out.println("b=" + b);
    System.out.println("n=" + n);
    System.out.println("dose=" + (Double.isNaN(dose) ? "" : dose));
    System.out.println("probit=" + (Double.isNaN(dose) ? "" : model.probit(dose)));
    System.out.println(
        "probability_from_dose=" + (Double.isNaN(dose) ? "" : model.probabilityFromDose(dose)));
    System.out.println("probit_value=" + (Double.isNaN(y) ? "" : y));
    System.out.println("probability=" + (Double.isNaN(y) ? "" : model.probability(y)));
    System.out.println();
  }

  // === process.jet_fire, process.pool_fire, process.vce, process.bleve =======

  static void jetFireRows() {
    // The class's own test state (`FireModelsTest`): 2 kg/s, 50 MJ/kg, eta 0.20.
    jetFireRow("jet_fire_test_20m", 2.0, 50.0e6, 0.20, Double.NaN, 20.0, Double.NaN);
    jetFireRow("jet_fire_test_80m", 2.0, 50.0e6, 0.20, Double.NaN, 80.0, Double.NaN);
    jetFireRow("jet_fire_distance_12500", 2.0, 50.0e6, 0.20, Double.NaN, Double.NaN, 12500.0);
    // The transmissivity's humidity input, at one distance.
    jetFireRow("jet_fire_humidity_0", 2.0, 50.0e6, 0.20, 0.0, 50.0, Double.NaN);
    jetFireRow("jet_fire_humidity_05", 2.0, 50.0e6, 0.20, 0.5, 50.0, Double.NaN);
    jetFireRow("jet_fire_humidity_1", 2.0, 50.0e6, 0.20, 1.0, 50.0, Double.NaN);
    // Inside one metre the class evaluates tau at max(d, 1).
    jetFireRow("jet_fire_inside_one_metre", 2.0, 50.0e6, 0.20, 0.5, 0.5, Double.NaN);
    // A target small enough to test the 25-iteration cap.
    jetFireRow("jet_fire_small_target", 2.0, 50.0e6, 0.20, Double.NaN, Double.NaN, 10.0);
    // The class's zero-distance branch and its zero-target refusal, as evidence.
    jetFireRow("jet_fire_zero_distance_evidence", 2.0, 50.0e6, 0.20, Double.NaN, 0.0, Double.NaN);
    jetFireRow("jet_fire_zero_target_evidence", 2.0, 50.0e6, 0.20, Double.NaN, Double.NaN, 0.0);
  }

  static void jetFireRow(String label, double mass, double hoc, double eta, double humidity,
      double distance, double targetFlux) {
    JetFireModel model = new JetFireModel(mass, hoc, eta);
    if (!Double.isNaN(humidity)) {
      model.setHumidity(humidity);
    }
    System.out.println(label);
    System.out.println("mass_flow_kg_s=" + mass);
    System.out.println("heat_of_combustion_j_kg=" + hoc);
    System.out.println("radiative_fraction=" + eta);
    System.out.println("relative_humidity=" + (Double.isNaN(humidity) ? "default_0.5" : humidity));
    System.out.println("distance_m=" + (Double.isNaN(distance) ? "" : distance));
    System.out.println("target_flux_w_m2=" + (Double.isNaN(targetFlux) ? "" : targetFlux));
    System.out.println("total_radiative_power_w=" + model.totalRadiativePowerW());
    if (!Double.isNaN(distance)) {
      System.out.println("transmissivity=" + model.transmissivity(distance));
      System.out.println("incident_heat_flux_w_m2=" + model.incidentHeatFlux(distance));
    }
    if (!Double.isNaN(targetFlux)) {
      try {
        double solved = model.distanceToFlux(targetFlux);
        System.out.println("distance_to_flux_m=" + solved);
        System.out.println("flux_at_distance_w_m2=" + model.incidentHeatFlux(solved));
      } catch (RuntimeException ex) {
        System.out.println("refusal=" + ex.getClass().getSimpleName() + ": " + ex.getMessage());
      }
    }
    System.out.println();
  }

  static void poolFireRows() {
    // The class's own test state (`FireModelsTest`): 10 m, 0.06 kg/(m2 s), 45 MJ/kg,
    // 800 kg/m3, eta 0.30.
    poolFireRow("pool_fire_test_5m1", 10.0, 0.06, 45.0e6, 800.0, 0.30, 5.1, Double.NaN);
    poolFireRow("pool_fire_test_10m", 10.0, 0.06, 45.0e6, 800.0, 0.30, 10.0, Double.NaN);
    poolFireRow("pool_fire_test_20m", 10.0, 0.06, 45.0e6, 800.0, 0.30, 20.0, Double.NaN);
    poolFireRow("pool_fire_test_50m", 10.0, 0.06, 45.0e6, 800.0, 0.30, 50.0, Double.NaN);
    // 12.5 kW/m2 is above this pool's near-field flux; 10 kW/m2 is reachable.
    poolFireRow("pool_fire_distance_12500", 10.0, 0.06, 45.0e6, 800.0, 0.30, Double.NaN,
        12500.0);
    poolFireRow("pool_fire_distance_10000", 10.0, 0.06, 45.0e6, 800.0, 0.30, Double.NaN,
        10000.0);
    poolFireRow("pool_fire_zero_distance", 10.0, 0.06, 45.0e6, 800.0, 0.30, 0.0, Double.NaN);
    // A target above the model's near-field flux: the class answers its near bound.
    poolFireRow("pool_fire_unreachable_1mw", 10.0, 0.06, 45.0e6, 800.0, 0.30, Double.NaN,
        1.0e6);
    poolFireRow("pool_fire_zero_target_evidence", 10.0, 0.06, 45.0e6, 800.0, 0.30, Double.NaN,
        0.0);
    // The same geometry with heat of combustion, liquid density and radiative fraction
    // changed: the class reads none of them.
    poolFireRow("pool_fire_dead_parameters", 10.0, 0.06, 10.0e6, 500.0, 0.90, 20.0, Double.NaN);
  }

  static void poolFireRow(String label, double diameter, double burnRate, double hoc,
      double liquidDensity, double eta, double distance, double targetFlux) {
    PoolFireModel model = new PoolFireModel(diameter, burnRate, hoc, liquidDensity, eta);
    System.out.println(label);
    System.out.println("pool_diameter_m=" + diameter);
    System.out.println("burn_rate_kg_m2_s=" + burnRate);
    System.out.println("heat_of_combustion_j_kg=" + hoc);
    System.out.println("liquid_density_kg_m3=" + liquidDensity);
    System.out.println("radiative_fraction=" + eta);
    System.out.println("distance_m=" + (Double.isNaN(distance) ? "" : distance));
    System.out.println("target_flux_w_m2=" + (Double.isNaN(targetFlux) ? "" : targetFlux));
    System.out.println("total_burn_rate_kg_s=" + model.totalBurnRateKgPerS());
    System.out.println("flame_height_m=" + model.flameHeightM());
    System.out.println("surface_emissive_power_w_m2=" + model.surfaceEmissivePowerWPerM2());
    if (!Double.isNaN(distance)) {
      System.out.println("view_factor=" + model.viewFactorVertical(distance));
      System.out.println("incident_heat_flux_w_m2=" + model.incidentHeatFlux(distance));
    }
    if (!Double.isNaN(targetFlux)) {
      double solved = model.distanceToFlux(targetFlux);
      System.out.println("distance_to_flux_m=" + solved);
      System.out.println("flux_at_distance_w_m2=" + model.incidentHeatFlux(solved));
    }
    System.out.println();
  }

  static void vceRows() {
    // The class's own test state (`FireModelsTest`): 1000 kg, 46 MJ/kg, class 4.
    vceRow("vce_test_20m", 1000.0, 46.0e6, 4, 20.0, Double.NaN);
    vceRow("vce_test_200m", 1000.0, 46.0e6, 4, 200.0, Double.NaN);
    // A class-4 plateau is 10.1 kPa, so 20.7 kPa is unreachable and 5 kPa is reachable.
    vceRow("vce_distance_20700", 1000.0, 46.0e6, 4, Double.NaN, 20700.0);
    vceRow("vce_distance_5000", 1000.0, 46.0e6, 4, Double.NaN, 5000.0);
    for (int strength = 1; strength <= 10; strength++) {
      vceRow("vce_class_" + strength, 1000.0, 46.0e6, strength, 20.0, Double.NaN);
    }
    // A class-1 cloud cannot reach 1 MPa; the class answers its near bound.
    vceRow("vce_unreachable_1mpa_class1", 1000.0, 46.0e6, 1, Double.NaN, 1.0e6);
    vceRow("vce_zero_target_evidence", 1000.0, 46.0e6, 4, Double.NaN, 0.0);
  }

  static void vceRow(String label, double mass, double hoc, int strength, double distance,
      double targetOverpressure) {
    VCEModel model = new VCEModel(mass, hoc, strength);
    System.out.println(label);
    System.out.println("cloud_fuel_mass_kg=" + mass);
    System.out.println("heat_of_combustion_j_kg=" + hoc);
    System.out.println("strength_class=" + strength);
    System.out.println("distance_m=" + (Double.isNaN(distance) ? "" : distance));
    System.out.println(
        "target_overpressure_pa=" + (Double.isNaN(targetOverpressure) ? "" : targetOverpressure));
    System.out.println("combustion_energy_j=" + model.combustionEnergyJ());
    if (!Double.isNaN(distance)) {
      System.out.println("scaled_distance=" + model.scaledDistance(distance));
      System.out.println(
          "dimensionless_overpressure=" + model.dimensionlessOverpressure(distance));
      System.out.println("overpressure_pa=" + model.overpressurePa(distance));
    }
    if (!Double.isNaN(targetOverpressure)) {
      try {
        double solved = model.distanceToOverpressure(targetOverpressure);
        System.out.println("distance_to_overpressure_m=" + solved);
        System.out.println("overpressure_at_distance_pa=" + model.overpressurePa(solved));
      } catch (RuntimeException ex) {
        System.out.println("refusal=" + ex.getClass().getSimpleName() + ": " + ex.getMessage());
      }
    }
    System.out.println();
  }

  static void bleveRows() {
    // The class's own test state (`FireModelsTest`): 1000 kg, 46 MJ/kg, eta 0.30.
    bleveRow("bleve_test_geometry", 1000.0, 46.0e6, 0.30, Double.NaN, Double.NaN, Double.NaN);
    bleveRow("bleve_flux_100m", 1000.0, 46.0e6, 0.30, Double.NaN, 100.0, Double.NaN);
    // The same distance and humidity as the jet fire's tau rows: the helper is shared.
    bleveRow("bleve_humidity_05_50m", 1000.0, 46.0e6, 0.30, 0.5, 50.0, Double.NaN);
    bleveRow("bleve_humidity_0", 1000.0, 46.0e6, 0.30, 0.0, 50.0, Double.NaN);
    bleveRow("bleve_distance_12500", 1000.0, 46.0e6, 0.30, Double.NaN, Double.NaN, 12500.0);
    bleveRow("bleve_small_target", 1000.0, 46.0e6, 0.30, Double.NaN, Double.NaN, 10.0);
    // A target above the fireball's near-field flux: the class answers its near bound.
    bleveRow("bleve_unreachable_1mw", 1000.0, 46.0e6, 0.30, Double.NaN, Double.NaN, 1.0e6);
    // The class accepts a negative radiative fraction and answers a negative flux.
    bleveRow("bleve_negative_fraction", 1000.0, 46.0e6, -0.30, Double.NaN, 100.0, 12500.0);
    bleveRow("bleve_zero_distance_evidence", 1000.0, 46.0e6, 0.30, Double.NaN, 0.0, Double.NaN);
    bleveRow("bleve_zero_target_evidence", 1000.0, 46.0e6, 0.30, Double.NaN, Double.NaN, 0.0);
  }

  static void bleveRow(String label, double mass, double hoc, double eta, double humidity,
      double distance, double targetFlux) {
    BLEVECalculator model = new BLEVECalculator(mass, hoc, eta);
    if (!Double.isNaN(humidity)) {
      model.setHumidity(humidity);
    }
    System.out.println(label);
    System.out.println("fuel_mass_kg=" + mass);
    System.out.println("heat_of_combustion_j_kg=" + hoc);
    System.out.println("radiative_fraction=" + eta);
    System.out.println("relative_humidity=" + (Double.isNaN(humidity) ? "default_0.5" : humidity));
    System.out.println("distance_m=" + (Double.isNaN(distance) ? "" : distance));
    System.out.println("target_flux_w_m2=" + (Double.isNaN(targetFlux) ? "" : targetFlux));
    System.out.println("fireball_diameter_m=" + model.fireballDiameterM());
    System.out.println("fireball_duration_s=" + model.fireballDurationS());
    System.out.println("surface_emissive_power_w_m2=" + model.surfaceEmissivePowerWPerM2());
    if (!Double.isNaN(distance)) {
      System.out.println("transmissivity=" + model.transmissivity(distance));
      System.out.println("incident_heat_flux_w_m2=" + model.incidentHeatFlux(distance));
    }
    if (!Double.isNaN(targetFlux)) {
      double solved = model.distanceToFlux(targetFlux);
      System.out.println("distance_to_flux_m=" + solved);
      System.out.println("flux_at_distance_w_m2=" + model.incidentHeatFlux(solved));
    }
    System.out.println();
  }

  // === process.event_tree, process.consequence_analysis =====================

  static void eventTreeRows() {
    eventTreeRow("event_tree_two_branch", 1.0e-3, new String[] {"Ignition"},
        new double[] {0.1});
    eventTreeRow("event_tree_three_branch", 1.0e-2,
        new String[] {"Immediate ignition", "Delayed ignition", "ESD success"},
        new double[] {0.05, 0.10, 0.99});
    eventTreeRow("event_tree_no_branches", 1.0e-3, new String[] {}, new double[] {});
    eventTreeRow("event_tree_zero_and_one", 1.0e-3, new String[] {"Certain", "Impossible"},
        new double[] {1.0, 0.0});

    System.out.println("event_tree_negative_frequency_evidence");
    System.out.println("initiating_frequency_per_year=-1.0E-3");
    try {
      new EventTreeAnalyzer("bad", -1.0e-3);
      System.out.println("refusal=none");
    } catch (RuntimeException ex) {
      System.out.println("refusal=" + ex.getClass().getSimpleName() + ": " + ex.getMessage());
    }
    System.out.println();

    System.out.println("event_tree_probability_outside_evidence");
    System.out.println("probability_of_yes=1.5");
    try {
      EventTreeAnalyzer tree = new EventTreeAnalyzer("bad", 1.0e-3);
      tree.addBranch("Bad", 1.5);
      System.out.println("refusal=none");
    } catch (RuntimeException ex) {
      System.out.println("refusal=" + ex.getClass().getSimpleName() + ": " + ex.getMessage());
    }
    System.out.println();

    // The class's `1 << n` has no cap; the count is printed, never evaluated.
    System.out.println("event_tree_shift_evidence");
    System.out.println("branch_count=21");
    System.out.println("left_shift_outcome_count=" + (1 << 21));
    EventTreeAnalyzer shifted = new EventTreeAnalyzer("shift", 1.0e-3);
    for (int i = 0; i < 21; i++) {
      shifted.addBranch("B" + i, 0.5);
    }
    System.out.println("branches_added=21");
    System.out.println();
  }

  static void eventTreeRow(String label, double frequency, String[] names,
      double[] probabilities) {
    System.out.println(label);
    System.out.println("initiating_frequency_per_year=" + frequency);
    System.out.println("branch_count=" + names.length);
    for (int i = 0; i < names.length; i++) {
      System.out.println("branch" + i + "_name=" + names[i]);
      System.out.println("branch" + i + "_probability_of_yes=" + probabilities[i]);
    }
    EventTreeAnalyzer tree = new EventTreeAnalyzer(label, frequency);
    for (int i = 0; i < names.length; i++) {
      tree.addBranch(names[i], probabilities[i]);
    }
    java.util.List<EventTreeAnalyzer.Outcome> outcomes = tree.evaluate();
    System.out.println("outcome_count=" + outcomes.size());
    for (int i = 0; i < outcomes.size(); i++) {
      EventTreeAnalyzer.Outcome outcome = outcomes.get(i);
      System.out.println("outcome" + i + "_path=" + outcome.path);
      System.out.println("outcome" + i + "_probability=" + outcome.probability);
      System.out.println("outcome" + i + "_frequency_per_year=" + outcome.frequencyPerYear);
    }
    System.out.println();
  }

  static void consequenceAnalysisRows() {
    JetFireModel jet = new JetFireModel(30.0, 50.0e6, 0.25);
    ConsequenceAnalysisEngine jetEngine = new ConsequenceAnalysisEngine("Leak", 1.0e-4);
    jetEngine.addJetFire(0.05, jet, ProbitModel.thermalFatality(), 60.0);
    consequenceRow("consequence_analysis_jet_fire", 1.0e-4, 10.0, jetEngine,
        new String[] {"jet_fire"}, new double[] {jet.incidentHeatFlux(10.0)},
        new String[] {"outcome0_mass_flow_kg_s=30.0",
            "outcome0_heat_of_combustion_j_kg=5.0E7",
            "outcome0_radiative_fraction=0.25", "outcome0_probit=thermal",
            "outcome0_exposure_seconds=60.0", "outcome0_conditional_probability=0.05"});

    PoolFireModel pool = new PoolFireModel(10.0, 0.06, 45.0e6, 800.0, 0.30);
    ConsequenceAnalysisEngine poolEngine = new ConsequenceAnalysisEngine("Pool", 1.0e-4);
    poolEngine.addPoolFire(0.10, pool, ProbitModel.thermalFatality(), 60.0);
    consequenceRow("consequence_analysis_pool_fire", 1.0e-4, 20.0, poolEngine,
        new String[] {"pool_fire"}, new double[] {pool.incidentHeatFlux(20.0)},
        new String[] {"outcome0_pool_diameter_m=10.0", "outcome0_burn_rate_kg_m2_s=0.06",
            "outcome0_probit=thermal", "outcome0_exposure_seconds=60.0",
            "outcome0_conditional_probability=0.10"});

    VCEModel vce = new VCEModel(1000.0, 46.0e6, 4);
    ConsequenceAnalysisEngine vceEngine = new ConsequenceAnalysisEngine("VCE", 1.0e-4);
    vceEngine.addVCE(0.30, vce, ProbitModel.blastLungFatality());
    consequenceRow("consequence_analysis_vce", 1.0e-4, 20.0, vceEngine,
        new String[] {"vce"}, new double[] {vce.overpressurePa(20.0)},
        new String[] {"outcome0_cloud_fuel_mass_kg=1000.0",
            "outcome0_heat_of_combustion_j_kg=4.6E7", "outcome0_strength_class=4",
            "outcome0_probit=blast", "outcome0_conditional_probability=0.30"});

    GaussianPlume plume = new GaussianPlume(1.0, 10.0, 5.0, GaussianPlume.Stability.D,
        GaussianPlume.Terrain.RURAL);
    ConsequenceAnalysisEngine toxicEngine = new ConsequenceAnalysisEngine("Toxic", 1.0e-4);
    toxicEngine.addToxicDispersion(0.20, plume, ProbitModel.h2sFatality(), 0.0340819, 288.15,
        1.01325, 10.0);
    consequenceRow("consequence_analysis_toxic_dispersion", 1.0e-4, 50.0, toxicEngine,
        new String[] {"toxic_dispersion"},
        new double[] {plume.centerlineGroundConcentration(50.0)},
        new String[] {"outcome0_emission_rate_kg_s=1.0",
            "outcome0_effective_stack_height_m=10.0", "outcome0_wind_speed_m_per_s=5.0",
            "outcome0_stability=D", "outcome0_terrain=rural",
            "outcome0_molar_mass_kg_mol=0.0340819", "outcome0_temperature_k=288.15",
            "outcome0_pressure_bara=1.01325", "outcome0_exposure_minutes=10.0",
            "outcome0_exposure_seconds=600.0", "outcome0_probit=h2s",
            "outcome0_conditional_probability=0.20"});

    // All four branches in one engine, in caller order, at one receiver.
    PoolFireModel mixedPool = new PoolFireModel(10.0, 0.06, 45.0e6, 800.0, 0.30);
    VCEModel mixedVce = new VCEModel(1000.0, 46.0e6, 4);
    GaussianPlume mixedPlume = new GaussianPlume(1.0, 10.0, 5.0, GaussianPlume.Stability.D,
        GaussianPlume.Terrain.RURAL);
    ConsequenceAnalysisEngine mixed = new ConsequenceAnalysisEngine("Mixed", 1.0e-3);
    mixed.addJetFire(0.05, new JetFireModel(30.0, 50.0e6, 0.25),
        ProbitModel.thermalFatality(), 60.0);
    mixed.addPoolFire(0.10, mixedPool, ProbitModel.thermalFatality(), 60.0);
    mixed.addVCE(0.30, mixedVce, ProbitModel.blastLungFatality());
    mixed.addToxicDispersion(0.20, mixedPlume, ProbitModel.h2sFatality(), 0.0340819, 288.15,
        1.01325, 10.0);
    consequenceRow("consequence_analysis_mixed", 1.0e-3, 20.0, mixed,
        new String[] {"jet_fire", "pool_fire", "vce", "toxic_dispersion"},
        new double[] {new JetFireModel(30.0, 50.0e6, 0.25).incidentHeatFlux(20.0),
            mixedPool.incidentHeatFlux(20.0), mixedVce.overpressurePa(20.0),
            mixedPlume.centerlineGroundConcentration(20.0)},
        new String[] {"outcome0_conditional_probability=0.05",
            "outcome1_conditional_probability=0.10",
            "outcome2_conditional_probability=0.30",
            "outcome3_conditional_probability=0.20"});

    // A zero exposure makes the thermal dose zero, which the class answers as 0.0.
    ConsequenceAnalysisEngine zeroExposure = new ConsequenceAnalysisEngine("Zero exposure", 1.0e-4);
    zeroExposure.addJetFire(0.05, new JetFireModel(30.0, 50.0e6, 0.25),
        ProbitModel.thermalFatality(), 0.0);
    consequenceRow("consequence_analysis_zero_exposure", 1.0e-4, 10.0, zeroExposure,
        new String[] {"jet_fire"},
        new double[] {new JetFireModel(30.0, 50.0e6, 0.25).incidentHeatFlux(10.0)},
        new String[] {"outcome0_exposure_seconds=0.0"});

    // A zero emission makes the plume concentration zero.
    GaussianPlume silent = new GaussianPlume(0.0, 10.0, 5.0, GaussianPlume.Stability.D,
        GaussianPlume.Terrain.RURAL);
    ConsequenceAnalysisEngine zeroEmission = new ConsequenceAnalysisEngine("Zero emission", 1.0e-4);
    zeroEmission.addToxicDispersion(0.20, silent, ProbitModel.h2sFatality(), 0.0340819,
        288.15, 1.01325, 10.0);
    consequenceRow("consequence_analysis_zero_emission", 1.0e-4, 50.0, zeroEmission,
        new String[] {"toxic_dispersion"},
        new double[] {silent.centerlineGroundConcentration(50.0)},
        new String[] {"outcome0_emission_rate_kg_s=0.0", "outcome0_exposure_seconds=600.0"});

    // The class answers an infinite flux at a zero receiver distance.
    ConsequenceAnalysisEngine zeroDistance = new ConsequenceAnalysisEngine("Zero distance", 1.0e-4);
    JetFireModel zeroDistanceJet = new JetFireModel(30.0, 50.0e6, 0.25);
    zeroDistance.addJetFire(0.05, zeroDistanceJet, ProbitModel.thermalFatality(), 60.0);
    consequenceRow("consequence_analysis_zero_distance_evidence", 1.0e-4, 0.0, zeroDistance,
        new String[] {"jet_fire"},
        new double[] {zeroDistanceJet.incidentHeatFlux(0.0)},
        new String[] {"outcome0_receiver_distance_m=0.0"});

    System.out.println("consequence_analysis_negative_probability_evidence");
    System.out.println("conditional_probability=-0.1");
    try {
      ConsequenceAnalysisEngine bad = new ConsequenceAnalysisEngine("Bad", 1.0e-4);
      bad.addJetFire(-0.1, new JetFireModel(30.0, 50.0e6, 0.25),
          ProbitModel.thermalFatality(), 60.0);
      System.out.println("refusal=none");
    } catch (RuntimeException ex) {
      System.out.println("refusal=" + ex.getClass().getSimpleName() + ": " + ex.getMessage());
    }
    System.out.println();
  }

  static void consequenceRow(String label, double releaseFrequency, double receiver,
      ConsequenceAnalysisEngine engine, String[] types, double[] intermediate, String[] inputs) {
    System.out.println(label);
    System.out.println("release_frequency_per_year=" + releaseFrequency);
    System.out.println("receiver_distance_m=" + receiver);
    for (String input : inputs) {
      System.out.println(input);
    }
    for (int i = 0; i < types.length; i++) {
      System.out.println("outcome" + i + "_outcome_type=" + types[i]);
      System.out.println("outcome" + i + "_intermediate=" + intermediate[i]);
    }
    java.util.List<ConsequenceAnalysisEngine.OutcomeResult> rows = engine.evaluate(receiver);
    System.out.println("outcome_count=" + rows.size());
    for (int i = 0; i < rows.size(); i++) {
      ConsequenceAnalysisEngine.OutcomeResult row = rows.get(i);
      System.out.println("outcome" + i + "_outcome=" + row.outcomeName);
      System.out.println("outcome" + i + "_frequency_per_year=" + row.outcomeFrequencyPerYear);
      System.out.println("outcome" + i + "_fatality_probability=" + row.fatalityProbability);
      System.out.println(
          "outcome" + i + "_fatality_frequency_per_year=" + row.fatalityFrequencyPerYear);
    }
    System.out.println(
        "individual_fatality_risk_per_year=" + engine.individualFatalityRiskPerYear(receiver));
    System.out.println();
  }
}
