// `ScalePredictionCalculator` on the states the port will hold: the upstream tests'
// brines, the Ksp pins, the auto-pH path and the pressure-correction boundary. Every
// intermediate the class exposes is printed, so the port can be diffed key by key.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ScalePredictionProbe.java
//     java -cp .:neqsim-f0c7436.jar ScalePredictionProbe > captures/scale_prediction_probe.tsv

import java.lang.reflect.Field;
import neqsim.pvtsimulation.flowassurance.ScalePredictionCalculator;

public class ScalePredictionProbe {

  public static void main(String[] args) {
    row("upstream_defaults", new ScalePredictionCalculator());

    ScalePredictionCalculator ksp = new ScalePredictionCalculator();
    ksp.setTemperatureCelsius(25.0);
    ksp.setPressureBara(1.013);
    row("ksp_at_25c_1atm", ksp);

    row("north_sea_manual_ph_80c_100bara",
        northSea(false, 6.5));
    row("north_sea_auto_ph", northSea(true, 6.5));

    ScalePredictionCalculator barite = bariteBrine(false, 6.5, 1.013);
    row("barite_supersaturated_manual_ph", barite);
    row("barite_supersaturated_auto_ph", bariteBrine(true, 6.5, 1.013));
    row("barite_just_above_atmosphere", bariteBrine(false, 6.5, 1.014));

    row("ion_pairing_plain", ionPairing(false));
    row("ion_pairing_mg_na", ionPairing(true));

    ScalePredictionCalculator secondOrder = new ScalePredictionCalculator();
    secondOrder.setCalciumConcentration(400.0);
    secondOrder.setBicarbonateConcentration(150.0);
    secondOrder.setSulphateConcentration(10.0);
    secondOrder.setTotalDissolvedSolids(35000.0);
    secondOrder.setTemperatureCelsius(80.0);
    secondOrder.setPressureBara(400.0);
    secondOrder.setCO2PartialPressure(2.0);
    secondOrder.setPH(6.5);
    secondOrder.setSecondOrderPressureCorrection(false);
    row("second_order_off", secondOrder);
    secondOrder.setSecondOrderPressureCorrection(true);
    row("second_order_on", secondOrder);

    ScalePredictionCalculator absent = new ScalePredictionCalculator();
    absent.setCalciumConcentration(400.0);
    absent.setBicarbonateConcentration(150.0);
    absent.setSulphateConcentration(10.0);
    absent.setTotalDissolvedSolids(35000.0);
    absent.setTemperatureCelsius(80.0);
    absent.setPressureBara(100.0);
    absent.setCO2PartialPressure(2.0);
    absent.setPH(6.5);
    row("absent_barium_strontium_iron", absent);

    ScalePredictionCalculator noBicarbonate = new ScalePredictionCalculator();
    noBicarbonate.setCalciumConcentration(400.0);
    noBicarbonate.setBicarbonateConcentration(0.0);
    noBicarbonate.setSulphateConcentration(10.0);
    noBicarbonate.setTotalDissolvedSolids(35000.0);
    noBicarbonate.setTemperatureCelsius(80.0);
    noBicarbonate.setPressureBara(100.0);
    noBicarbonate.setCO2PartialPressure(2.0);
    noBicarbonate.enableAutoPH();
    row("zero_bicarbonate_auto_ph", noBicarbonate);

    ScalePredictionCalculator srfe = new ScalePredictionCalculator();
    srfe.setCalciumConcentration(400.0);
    srfe.setStrontiumConcentration(300.0);
    srfe.setIronConcentration(50.0);
    srfe.setBicarbonateConcentration(150.0);
    srfe.setSulphateConcentration(200.0);
    srfe.setTotalDissolvedSolids(35000.0);
    srfe.setTemperatureCelsius(80.0);
    srfe.setPressureBara(100.0);
    srfe.setCO2PartialPressure(2.0);
    srfe.setPH(6.5);
    row("strontium_and_iron", srfe);
  }

  static ScalePredictionCalculator northSea(boolean autoPh, double ph) {
    ScalePredictionCalculator calc = new ScalePredictionCalculator();
    calc.setCalciumConcentration(400.0);
    calc.setBariumConcentration(0.0);
    calc.setStrontiumConcentration(0.0);
    calc.setIronConcentration(0.0);
    calc.setBicarbonateConcentration(150.0);
    calc.setSulphateConcentration(10.0);
    calc.setTotalDissolvedSolids(35000.0);
    calc.setTemperatureCelsius(80.0);
    calc.setPressureBara(100.0);
    calc.setCO2PartialPressure(2.0);
    calc.setPH(ph);
    if (autoPh) {
      calc.enableAutoPH();
    }
    return calc;
  }

  static ScalePredictionCalculator bariteBrine(boolean autoPh, double ph, double pressureBara) {
    ScalePredictionCalculator calc = new ScalePredictionCalculator();
    calc.setCalciumConcentration(100.0);
    calc.setBariumConcentration(200.0);
    calc.setStrontiumConcentration(0.0);
    calc.setIronConcentration(0.0);
    calc.setBicarbonateConcentration(100.0);
    calc.setSulphateConcentration(2000.0);
    calc.setTotalDissolvedSolids(40000.0);
    calc.setTemperatureCelsius(60.0);
    calc.setPressureBara(pressureBara);
    calc.setCO2PartialPressure(0.5);
    calc.setPH(ph);
    if (autoPh) {
      calc.enableAutoPH();
    }
    return calc;
  }

  static ScalePredictionCalculator ionPairing(boolean withMagnesiumAndSodium) {
    ScalePredictionCalculator calc = new ScalePredictionCalculator();
    calc.setCalciumConcentration(100.0);
    calc.setBariumConcentration(50.0);
    calc.setStrontiumConcentration(0.0);
    calc.setIronConcentration(0.0);
    calc.setBicarbonateConcentration(100.0);
    calc.setSulphateConcentration(2000.0);
    calc.setTotalDissolvedSolids(50000.0);
    calc.setMagnesiumConcentration(withMagnesiumAndSodium ? 1300.0 : 0.0);
    calc.setSodiumConcentration(withMagnesiumAndSodium ? 11000.0 : 0.0);
    calc.setTemperatureCelsius(60.0);
    calc.setPressureBara(1.013);
    calc.setCO2PartialPressure(0.5);
    calc.setPH(6.5);
    return calc;
  }

  static void row(String label, ScalePredictionCalculator calc) {
    calc.calculate();
    StringBuilder out = new StringBuilder();
    out.append(label).append('\n');
    out.append("temperature_K=").append(calc.getTemperatureKelvin()).append('\n');
    out.append("ionic_strength_mol_per_l=").append(calc.getIonicStrengthMolPerL()).append('\n');
    out.append("debye_huckel_a=").append(calc.getDebyeHuckelAParameter()).append('\n');
    out.append("divalent_gamma=").append(calc.getDivalentActivityCoefficient()).append('\n');
    out.append("ph_used=").append(ph(calc)).append('\n');
    out.append("free_calcium_mol_per_l=").append(calc.getFreeCalciumMolPerL()).append('\n');
    out.append("free_sulphate_mol_per_l=").append(calc.getFreeSulphateMolPerL()).append('\n');
    out.append("free_carbonate_mol_per_l=").append(calc.getFreeCarbonateMolPerL()).append('\n');
    out.append("total_barium_mol_per_l=").append(calc.getTotalBariumMolPerL()).append('\n');
    out.append("total_strontium_mol_per_l=").append(calc.getTotalStrontiumMolPerL()).append('\n');
    out.append("total_iron_mol_per_l=").append(calc.getTotalIronMolPerL()).append('\n');
    out.append("total_sodium_mol_per_l=").append(calc.getTotalSodiumMolPerL()).append('\n');
    out.append("tds_mg_per_l=").append(calc.getTotalDissolvedSolidsMgPerL()).append('\n');
    out.append("ksp_calcite=").append(calc.getKspCalcite()).append('\n');
    out.append("ksp_barite=").append(calc.getKspBarite()).append('\n');
    out.append("ksp_celestite=").append(calc.getKspCelestite()).append('\n');
    out.append("ksp_anhydrite=").append(calc.getKspAnhydrite()).append('\n');
    out.append("ksp_siderite=").append(calc.getKspSiderite()).append('\n');
    out.append("log10_ksp_calcite=").append(Math.log10(calc.getKspCalcite())).append('\n');
    out.append("log10_ksp_barite=").append(Math.log10(calc.getKspBarite())).append('\n');
    out.append("log10_ksp_celestite=").append(Math.log10(calc.getKspCelestite())).append('\n');
    out.append("log10_ksp_anhydrite=").append(Math.log10(calc.getKspAnhydrite())).append('\n');
    out.append("log10_ksp_siderite=").append(Math.log10(calc.getKspSiderite())).append('\n');
    out.append("si_caco3=").append(calc.getCaCO3SaturationIndex()).append('\n');
    out.append("si_baso4=").append(calc.getBaSO4SaturationIndex()).append('\n');
    out.append("si_srso4=").append(calc.getSrSO4SaturationIndex()).append('\n');
    out.append("si_caso4=").append(calc.getCaSO4SaturationIndex()).append('\n');
    out.append("si_feco3=").append(calc.getFeCO3SaturationIndex()).append('\n');
    out.append("scale_risks=").append(String.join("|", calc.getScaleRisks())).append('\n');
    out.append("has_scaling_risk=").append(calc.hasScalingRisk()).append('\n');
    System.out.print(out);
    System.out.println();
  }

  static double ph(ScalePredictionCalculator calc) {
    try {
      Field field = ScalePredictionCalculator.class.getDeclaredField("pH");
      field.setAccessible(true);
      return (Double) field.get(calc);
    } catch (ReflectiveOperationException ex) {
      throw new AssertionError("pH is not readable", ex);
    }
  }
}
