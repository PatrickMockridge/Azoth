// FlowAssuranceRemainderProbe - the reachable remainder of pvtsimulation/flowassurance at the pin.
//
// Six calculators in one JVM run: the de Waard-Milliams CO2/H2S/S8 corrosion correlation, the
// API RP 14E / DNV RP O501 erosion screen, the emulsion viscosity models, the lumped pipeline
// cooldown march, the hydrate risk map over a profile, and the surface cooldown no-touch-time
// analyser. Each section prints the class's own units; the cases convert once.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar FlowAssuranceRemainderProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar FlowAssuranceRemainderProbe > captures/flow_assurance_remainder_probe.tsv

import java.util.List;
import neqsim.pvtsimulation.flowassurance.DeWaardMilliamsCorrosion;
import neqsim.pvtsimulation.flowassurance.EmulsionViscosityCalculator;
import neqsim.pvtsimulation.flowassurance.ErosionPredictionCalculator;
import neqsim.pvtsimulation.flowassurance.HydrateRiskMapper;
import neqsim.pvtsimulation.flowassurance.PipelineCooldownCalculator;
import neqsim.pvtsimulation.flowassurance.SurfCooldownAnalyzer;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class FlowAssuranceRemainderProbe {

  // ─── de Waard-Milliams ────────────────────────────────────────────────────

  static void dew(String label, double tC, double pco2, double ph, double ph2s, double v,
      double ie, double glycol, double sulfur, double life) {
    DeWaardMilliamsCorrosion m = new DeWaardMilliamsCorrosion();
    m.setTemperatureCelsius(tC);
    m.setCO2PartialPressure(pco2);
    m.setPH(ph);
    m.setH2SPartialPressure(ph2s);
    m.setFlowVelocity(v);
    m.setInhibitorEfficiency(ie);
    m.setGlycolFraction(glycol);
    m.setSulfurDepositionRate(sulfur);
    double corrected = m.calculateCorrosionRate();
    System.out.println("state=" + label
        + " baseline_mmyr=" + m.calculateBaselineRate()
        + " ph_factor=" + m.getPHCorrectionFactor()
        + " scale_factor=" + m.getScaleCorrectionFactor()
        + " glycol_factor=" + m.getGlycolCorrectionFactor()
        + " flow_factor=" + m.getFlowCorrectionFactor()
        + " corrected_mmyr=" + corrected
        + " sulfur_rate_mmyr=" + m.getSulfurCorrosionRate()
        + " fes_rate_mmyr=" + m.getFeSCorrosionRate()
        + " feco3_si=" + m.getFeCO3SaturationIndex()
        + " feco3_scale_forming=" + (m.getFeCO3SaturationIndex() > 1.0)
        + " sour=" + m.isSourService()
        + " severity=" + m.getCorrosionSeverity()
        + " morphology=" + m.getFeSScaleMorphology()
        + " allowance_mm=" + m.estimateCorrosionAllowance(life));
  }

  // ─── erosion ──────────────────────────────────────────────────────────────

  static void erosion(String label, double rho, double v, double d, double wallMm, double sandKgDay,
      double dpMm, double rhoP, double mu, String geometry, double bend, String material,
      double apiC, double caMm, double life) {
    ErosionPredictionCalculator e = new ErosionPredictionCalculator();
    e.setMixtureDensity(rho);
    e.setMixtureVelocity(v);
    e.setPipeDiameter(d);
    e.setWallThickness(wallMm);
    e.setSandRate(sandKgDay);
    e.setSandParticleDiameter(dpMm);
    e.setSandParticleDensity(rhoP);
    e.setMixtureViscosity(mu);
    e.setGeometry(geometry);
    e.setBendRadiusRatio(bend);
    e.setPipeMaterial(material);
    e.setApiCFactor(apiC);
    e.setCorrosionAllowance(caMm);
    e.setDesignLife(life);
    e.calculate();
    System.out.println("state=" + label
        + " erosional_velocity_ms=" + e.getErosionalVelocity()
        + " velocity_ratio=" + e.getVelocityRatio()
        + " erosion_rate_mmyr=" + e.getErosionRate()
        + " cumulative_mm=" + e.getCumulativeErosion()
        + " remaining_mm=" + e.getRemainingWallThickness()
        + " within_api=" + e.isWithinApiLimits()
        + " within_erosion=" + e.isWithinErosionLimits()
        + " risk=" + e.getRiskLevel());
  }

  // ─── emulsion ─────────────────────────────────────────────────────────────

  static void emulsion(String label, double muO, double muW, double rhoO, double rhoW, double wc,
      double phiMax, String model, double tightness, boolean demulsifier, double efficiency) {
    EmulsionViscosityCalculator e = new EmulsionViscosityCalculator();
    e.setOilViscosity(muO);
    e.setWaterViscosity(muW);
    e.setOilDensity(rhoO);
    e.setWaterDensity(rhoW);
    e.setWaterCut(wc);
    e.setMaxPackingFraction(phiMax);
    e.setModel(model);
    e.setTightnessFactor(tightness);
    e.setDemulsifierPresent(demulsifier);
    e.setDemulsifierEfficiency(efficiency);
    e.calculate();
    System.out.println("state=" + label
        + " effective_viscosity_cp=" + e.getEffectiveViscosity()
        + " relative_viscosity=" + e.getRelativeViscosity()
        + " viscosity_ratio=" + e.getViscosityRatio()
        + " inversion_water_cut=" + e.getInversionWaterCut()
        + " emulsion_type=" + e.getEmulsionType()
        + " inverted=" + e.isInverted());
  }

  // ─── pipeline cooldown ────────────────────────────────────────────────────

  static void cooldown(String label, double diameter, double wall, double insulation, double coating,
      double insulationK, Double uValue, Double htc, double t0, double ambient, double rho, double cp,
      double dtMinutes, double totalHours, Double targetK, Double steelDensity, Double steelCp) {
    PipelineCooldownCalculator c = new PipelineCooldownCalculator();
    c.setInternalDiameter(diameter);
    c.setWallThickness(wall);
    c.setInsulationThickness(insulation);
    c.setCoatingThickness(coating);
    c.setInsulationConductivity(insulationK);
    if (steelDensity != null) {
      c.setSteelDensity(steelDensity);
    }
    if (steelCp != null) {
      c.setSteelSpecificHeat(steelCp);
    }
    if (uValue != null) {
      c.setOverallUValue(uValue);
    } else {
      c.useLayerCalculation();
      c.setExternalHTC(htc);
    }
    c.setInitialFluidTemperature(t0);
    c.setAmbientTemperature(ambient);
    c.setFluidDensity(rho);
    c.setFluidSpecificHeat(cp);
    c.setTimeStepMinutes(dtMinutes);
    c.setTotalTimeHours(totalHours);
    c.calculate();
    double[] times = c.getTimeHours();
    double[] temperatures = c.getFluidTemperature();
    System.out.println("state=" + label
        + " steps=" + (times.length - 1)
        + " time_constant_hours=" + c.getTimeConstantHours()
        + " calculated_u=" + c.getCalculatedUValue()
        + " target_hours=" + (targetK == null ? "absent" : c.getTimeToReachTemperature(targetK))
        + " temperature_at_1h=" + c.getTemperatureAtTime(1.0)
        + " temperature_at_10h=" + c.getTemperatureAtTime(10.0));
    if (times.length <= 25) {
      System.out.println("  time_hours=" + join(times) + " fluid_temperature_k=" + join(temperatures));
    }
  }

  static String join(double[] values) {
    StringBuilder text = new StringBuilder("[");
    for (int i = 0; i < values.length; i++) {
      if (i > 0) {
        text.append(",");
      }
      text.append(values[i]);
    }
    return text.append("]").toString();
  }

  // ─── hydrate risk map ─────────────────────────────────────────────────────

  static SystemInterface wetGas(double water) {
    return wetGasAt(water, 120.0);
  }

  static SystemInterface wetGasAt(double water, double pressureBara) {
    SystemSrkEos fluid = new SystemSrkEos(338.15, pressureBara);
    fluid.addComponent("methane", 0.82);
    fluid.addComponent("ethane", 0.06);
    fluid.addComponent("propane", 0.03);
    fluid.addComponent("water", water);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    return fluid;
  }

  static SystemInterface dryGas() {
    SystemSrkEos fluid = new SystemSrkEos(338.15, 120.0);
    fluid.addComponent("methane", 0.90);
    fluid.addComponent("ethane", 0.07);
    fluid.addComponent("propane", 0.03);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    return fluid;
  }

  static void hydrateMap(String label, SystemInterface fluid, double[][] points,
      Double highRisk, Double mediumRisk) {
    HydrateRiskMapper mapper = new HydrateRiskMapper(fluid);
    for (double[] point : points) {
      mapper.addProfilePoint(point[0], point[1], point[2]);
    }
    if (highRisk != null) {
      mapper.setRiskThresholds(highRisk, mediumRisk);
    }
    HydrateRiskMapper.RiskProfile profile = mapper.calculate();
    System.out.println("state=" + label
        + " overall_risk=" + profile.getOverallRisk()
        + " minimum_subcooling_c=" + profile.getMinimumSubcoolingC()
        + " critical_count=" + profile.getCriticalPointCount()
        + " unknown_count=" + profile.getUnknownPointCount()
        + " total_points=" + profile.getPoints().size()
        + " failure_reasons=" + profile.getFailureReasons().size());
    List<HydrateRiskMapper.RiskPoint> rows = profile.getPoints();
    for (int i = 0; i < rows.size(); i++) {
      HydrateRiskMapper.RiskPoint row = rows.get(i);
      System.out.println("  point=" + i
          + " distance_km=" + row.distanceKm
          + " pressure_bara=" + row.pressureBara
          + " actual_temperature_c=" + row.actualTemperatureC
          + " hydrate_temperature_c=" + row.hydrateTemperatureC
          + " subcooling_c=" + row.subcoolingC
          + " risk=" + row.riskLevel);
    }
  }

  // ─── surface cooldown ─────────────────────────────────────────────────────

  static void surf(String label, SystemInterface fluid, Double uValue, Double htc,
      double insulationK, double insulation, double seabedC, Double operatingC, Double marginK,
      Double requiredHours, double totalHours, double dtMinutes, boolean requireWater) {
    SurfCooldownAnalyzer analyzer = new SurfCooldownAnalyzer(fluid);
    analyzer.setInternalDiameter(0.254);
    analyzer.setWallThickness(0.0159);
    analyzer.setInsulationThickness(insulation);
    analyzer.setCoatingThickness(0.006);
    analyzer.setInsulationConductivity(insulationK);
    if (uValue != null) {
      analyzer.setOverallUValue(uValue);
    } else {
      analyzer.setExternalHTC(htc);
    }
    analyzer.setSeabedTemperature(seabedC);
    if (operatingC != null) {
      analyzer.setOperatingTemperature(operatingC);
    }
    if (marginK != null) {
      analyzer.setHydrateMargin(marginK);
    }
    if (requiredHours != null) {
      analyzer.setRequiredNoTouchTimeHours(requiredHours);
    }
    analyzer.setTotalTimeHours(totalHours);
    analyzer.setTimeStepMinutes(dtMinutes);
    analyzer.setRequireWater(requireWater);
    try {
      analyzer.calculate();
    } catch (RuntimeException failure) {
      System.out.println("state=" + label + " status=refused reason="
          + failure.getClass().getSimpleName() + ":" + failure.getMessage());
      return;
    }
    PipelineCooldownCalculator inner = analyzer.getCooldownCalculator();
    System.out.println("state=" + label
        + " verdict=" + analyzer.getVerdict()
        + " water_present=" + analyzer.isWaterPresent()
        + " water_mole_fraction=" + analyzer.getWaterMoleFraction()
        + " no_touch_hours=" + analyzer.getNoTouchTimeHours()
        + " initial_temperature_k=" + analyzer.getInitialFluidTemperatureK()
        + " hydrate_temperature_k=" + analyzer.getHydrateEquilibriumTemperatureK()
        + " fluid_density=" + analyzer.getFluidDensity()
        + " fluid_cp=" + analyzer.getFluidSpecificHeat()
        + " time_constant_hours=" + analyzer.getTimeConstantHours()
        + " calculated_u=" + (inner == null ? "absent" : inner.getCalculatedUValue()));
  }

  public static void main(String[] args) {
    System.out.println("pvtsimulation/flowassurance remainder at f0c7436c69");

    System.out.println("# --- DeWaardMilliamsCorrosion ---");
    dew("default", 60, 2, 4.5, 0, 2, 0, 0, 0, 25);
    dew("co2_zero", 25, 0, 3, 0, 1, 0, 0, 0, 25);
    dew("severity_low", 25, 0.001, 3, 0, 1, 0, 0, 0, 25);
    dew("severity_medium", 25, 0.05, 3, 0, 1, 0, 0, 0, 25);
    dew("severity_high", 25, 0.5, 3, 0, 1, 0, 0, 0, 25);
    dew("severity_very_high", 25, 1.0, 3, 0, 1, 0, 0, 0, 25);
    dew("ph_below_sat", 60, 2, 3.0, 0, 2, 0, 0, 0, 25);
    dew("ph_above_sat", 60, 2, 6.0, 0, 2, 0, 0, 0, 25);
    dew("scale_50", 50, 2, 6, 0, 2, 0, 0, 0, 25);
    dew("scale_60", 60, 2, 6, 0, 2, 0, 0, 0, 25);
    dew("scale_70", 70, 2, 6, 0, 2, 0, 0, 0, 25);
    dew("scale_100_ph6", 100, 2, 6, 0, 2, 0, 0, 0, 25);
    dew("scale_160_ph7", 160, 2, 7, 0, 2, 0, 0, 0, 25);
    dew("glycol_0_3", 60, 2, 4.5, 0, 2, 0, 0.3, 0, 25);
    dew("glycol_0_8", 60, 2, 4.5, 0, 2, 0, 0.8, 0, 25);
    dew("glycol_0_81", 60, 2, 4.5, 0, 2, 0, 0.81, 0, 25);
    dew("velocity_1", 60, 2, 4.5, 0, 1, 0, 0, 0, 25);
    dew("velocity_2", 60, 2, 4.5, 0, 2, 0, 0, 0, 25);
    dew("velocity_10", 60, 2, 4.5, 0, 10, 0, 0, 0, 25);
    dew("inhibitor_0_8", 60, 2, 4.5, 0, 2, 0.8, 0, 0, 25);
    dew("inhibitor_1_5", 60, 2, 4.5, 0, 2, 1.5, 0, 0, 25);
    dew("inhibitor_negative", 60, 2, 4.5, 0, 2, -0.5, 0, 0, 25);
    dew("sulfur", 60, 2, 4.5, 0, 2, 0, 0, 0.5, 25);
    dew("sour_threshold_below", 60, 2, 4.5, 0.003, 2, 0, 0, 0, 25);
    dew("sour_threshold_above", 60, 2, 4.5, 0.0031, 2, 0, 0, 0, 25);
    dew("fes_cap", 200, 2, 4.5, 100, 2, 0, 0, 0, 25);
    dew("morph_none", 150, 2, 4.5, 0, 2, 0, 0, 0, 25);
    dew("morph_50", 50, 2, 4.5, 0.1, 2, 0, 0, 0, 25);
    dew("morph_80", 80, 2, 4.5, 0.1, 2, 0, 0, 0, 25);
    dew("morph_150", 150, 2, 4.5, 0.1, 2, 0, 0, 0, 25);
    dew("feco3_high_ph", 90, 1, 7, 0, 2, 0, 0, 0, 25);

    System.out.println("# --- ErosionPredictionCalculator ---");
    erosion("api_base", 100, 10, 0.1524, 10, 0, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);
    erosion("api_c150", 100, 10, 0.1524, 10, 0, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 150, 3, 25);
    erosion("api_v5", 100, 5, 0.1524, 10, 0, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);
    erosion("api_v20", 100, 20, 0.1524, 10, 0, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);
    erosion("density_zero", 0, 10, 0.1524, 10, 0, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);
    erosion("elbow_rd1", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 1.0, "carbon_steel", 100, 3, 25);
    erosion("elbow_rd1_5", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);
    erosion("elbow_rd2", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 2.0, "carbon_steel", 100, 3, 25);
    erosion("elbow_rd5", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 5.0, "carbon_steel", 100, 3, 25);
    erosion("tee", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "tee", 1.5, "carbon_steel", 100, 3, 25);
    erosion("blind_tee", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "blind_tee", 1.5, "carbon_steel", 100, 3, 25);
    erosion("reducer", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "reducer", 1.5, "carbon_steel", 100, 3, 25);
    erosion("choke", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "choke", 1.5, "carbon_steel", 100, 3, 25);
    erosion("weld", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "weld", 1.5, "carbon_steel", 100, 3, 25);
    erosion("straight", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "straight", 1.5, "carbon_steel", 100, 3, 25);
    erosion("dp_0_5", 100, 10, 0.1524, 10, 50, 0.5, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);
    erosion("dp_1_5", 100, 10, 0.1524, 10, 50, 1.5, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);
    erosion("duplex", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 1.5, "duplex_steel", 100, 3, 25);
    erosion("super_duplex", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 1.5, "super_duplex", 100, 3, 25);
    erosion("13cr", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 1.5, "13cr", 100, 3, 25);
    erosion("inconel", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 1.5, "inconel", 100, 3, 25);
    erosion("titanium", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 1.5, "titanium", 100, 3, 25);
    erosion("tungsten_carbide", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 1.5, "tungsten_carbide", 100, 3, 25);
    erosion("eroded", 100, 10, 0.1524, 1, 500, 0.25, 2650, 0.001, "choke", 1.5, "carbon_steel", 100, 3, 25);
    erosion("life_50", 100, 10, 0.1524, 10, 50, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 50);
    erosion("rate_medium", 100, 5, 0.1524, 10, 500, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);
    erosion("rate_high", 100, 5, 0.1524, 10, 2000, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);
    erosion("rate_critical", 100, 5, 0.1524, 10, 4000, 0.25, 2650, 0.001, "elbow", 1.5, "carbon_steel", 100, 3, 25);

    System.out.println("# --- EmulsionViscosityCalculator ---");
    emulsion("default_pal", 10, 0.5, 850, 1020, 0.30, 0.74, "pal_rhodes", 0.5, false, 0);
    emulsion("einstein", 10, 0.5, 850, 1020, 0.01, 0.74, "einstein", 0.5, false, 0);
    emulsion("taylor", 10, 0.5, 850, 1020, 0.01, 0.74, "taylor", 0.5, false, 0);
    emulsion("brinkman", 10, 0.5, 850, 1020, 0.30, 0.74, "brinkman", 0.5, false, 0);
    emulsion("pal_jammed", 0.001, 1.0, 1500, 500, 0.80, 0.74, "pal_rhodes", 0.0, false, 0);
    emulsion("woelflin_loose", 10, 0.5, 850, 1020, 0.30, 0.74, "woelflin", 0.2, false, 0);
    emulsion("woelflin_medium", 10, 0.5, 850, 1020, 0.30, 0.74, "woelflin", 0.5, false, 0);
    emulsion("woelflin_tight", 10, 0.5, 850, 1020, 0.30, 0.74, "woelflin", 0.8, false, 0);
    emulsion("richardson_low", 10, 0.5, 850, 1020, 0.30, 0.74, "richardson", 0.2, false, 0);
    emulsion("richardson_high", 10, 0.5, 850, 1020, 0.30, 0.74, "richardson", 1.0, false, 0);
    emulsion("inverted", 10, 0.5, 850, 1020, 0.5, 0.74, "pal_rhodes", 0.5, false, 0);
    emulsion("low_clamp", 1000, 0.5, 500, 1000, 0.30, 0.74, "pal_rhodes", 0.0, false, 0);
    emulsion("high_clamp", 0.001, 1.0, 1500, 500, 0.30, 0.74, "pal_rhodes", 1.0, false, 0);
    emulsion("wc_zero", 10, 0.5, 850, 1020, 0.0, 0.74, "pal_rhodes", 0.5, false, 0);
    emulsion("wc_one", 10, 0.5, 850, 1020, 1.0, 0.74, "pal_rhodes", 0.5, false, 0);
    emulsion("tight_0_75", 10, 0.5, 850, 1020, 0.30, 0.74, "pal_rhodes", 0.75, false, 0);
    emulsion("tight_1", 10, 0.5, 850, 1020, 0.30, 0.74, "pal_rhodes", 1.0, false, 0);
    emulsion("demulsifier_0_6", 10, 0.5, 850, 1020, 0.30, 0.74, "pal_rhodes", 0.5, true, 0.6);
    emulsion("demulsifier_zero_eff", 10, 0.5, 850, 1020, 0.30, 0.74, "pal_rhodes", 0.5, true, 0.0);

    System.out.println("# --- PipelineCooldownCalculator ---");
    cooldown("u3_default", 0.254, 0.0127, 0.050, 0.006, 0.17, 3.0, null, 353.15, 277.15, 750, 2200, 5, 48, null, null, null);
    cooldown("u1_5", 0.254, 0.0127, 0.050, 0.006, 0.17, 1.5, null, 353.15, 277.15, 750, 2200, 5, 48, null, null, null);
    cooldown("u6", 0.254, 0.0127, 0.050, 0.006, 0.17, 6.0, null, 353.15, 277.15, 750, 2200, 5, 48, null, null, null);
    cooldown("layers_ins50", 0.254, 0.0127, 0.050, 0.006, 0.17, null, 300.0, 353.15, 277.15, 750, 2200, 5, 48, null, null, null);
    cooldown("layers_bare", 0.254, 0.0127, 0.0, 0.006, 0.17, null, 300.0, 353.15, 277.15, 750, 2200, 5, 48, null, null, null);
    cooldown("layers_htc500", 0.254, 0.0127, 0.050, 0.006, 0.17, null, 500.0, 353.15, 277.15, 750, 2200, 5, 48, null, null, null);
    cooldown("ceil_steps", 0.254, 0.0127, 0.050, 0.006, 0.17, 3.0, null, 353.15, 277.15, 750, 2200, 180, 10, null, null, null);
    cooldown("target_exact", 0.254, 0.0127, 0.050, 0.006, 0.17, 3.0, null, 353.15, 277.15, 750, 2200, 60, 10, 353.15, null, null);
    cooldown("target_unreached", 0.254, 0.0127, 0.050, 0.006, 0.17, 3.0, null, 353.15, 277.15, 750, 2200, 60, 10, 270.0, null, null);
    cooldown("heat_up", 0.254, 0.0127, 0.050, 0.006, 0.17, 3.0, null, 277.15, 300.0, 750, 2200, 60, 10, 300.0, null, null);
    cooldown("steel_8000_450", 0.254, 0.0127, 0.050, 0.006, 0.17, 3.0, null, 353.15, 277.15, 750, 2200, 60, 10, null, 8000.0, 450.0);
    cooldown("coating_20mm", 0.254, 0.0127, 0.050, 0.020, 0.17, 3.0, null, 353.15, 277.15, 750, 2200, 60, 10, null, null, null);
    cooldown("coarse_1h", 0.254, 0.0127, 0.050, 0.006, 0.17, 3.0, null, 353.15, 277.15, 750, 2200, 60, 10, null, null, null);

    System.out.println("# --- HydrateRiskMapper ---");
    double[][] fourPoints = { {0, 100, 30}, {10, 95, 25}, {25, 80, 15}, {50, 60, 8} };
    double[][] bandPoints = { {0, 120, 25}, {5, 110, 20}, {10, 100, 12}, {20, 90, 20} };
    hydrateMap("four_points", wetGas(0.09), fourPoints, null, null);
    hydrateMap("bands", wetGas(0.09), bandPoints, null, null);
    hydrateMap("custom_thresholds", wetGas(0.09), fourPoints, 1.0, 2.0);
    hydrateMap("one_point", wetGas(0.09), new double[][] { {0, 100, 20} }, null, null);
    hydrateMap("dry", dryGas(), fourPoints, null, null);
    hydrateMap("water_zero", wetGas(0.0), fourPoints, null, null);

    System.out.println("# --- SurfCooldownAnalyzer ---");
    surf("wet_u2_5", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, 8.0, 72, 5, false);
    surf("wet_u1_5", wetGas(0.09), 1.5, null, 0.22, 0.060, 4.0, null, 3.0, 8.0, 72, 5, false);
    surf("wet_u6", wetGas(0.09), 6.0, null, 0.22, 0.060, 4.0, null, 3.0, 8.0, 72, 5, false);
    surf("layers", wetGas(0.09), null, 500.0, 0.22, 0.060, 4.0, null, 3.0, 8.0, 72, 5, false);
    surf("no_required", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, null, 72, 5, false);
    surf("req_half", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, 4.0, 72, 5, false);
    surf("req_1_2t", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, 9.6, 72, 5, false);
    surf("req_2t", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, 16.0, 72, 5, false);
    surf("operating_override", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, 35.0, 3.0, 8.0, 72, 5, false);
    surf("initial_below_target", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, 20.0, 3.0, 8.0, 72, 5, false);
    surf("seabed_300k", wetGas(0.09), 2.5, null, 0.22, 0.060, 26.85, null, 3.0, 8.0, 72, 5, false);
    surf("dry", dryGas(), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, 8.0, 72, 5, false);
    surf("dry_require", dryGas(), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, 8.0, 72, 5, true);
    surf("wet_require", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, 8.0, 72, 5, true);
    surf("cold_operating", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, 15.0, 3.0, 8.0, 72, 5, false);
    surf("low_pressure", wetGasAt(0.09, 10.0), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, 8.0, 72, 5, false);
    surf("margin_1", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, null, 1.0, 8.0, 72, 5, false);
    surf("margin_5", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, null, 5.0, 8.0, 72, 5, false);
    surf("dt_600", wetGas(0.09), 2.5, null, 0.22, 0.060, 4.0, null, 3.0, 8.0, 72, 10, false);
  }
}
