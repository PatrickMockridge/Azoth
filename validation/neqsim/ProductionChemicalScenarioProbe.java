// ProductionChemicalScenarioProbe - the treatment-aware scale scenario at the pinned commit.
//
// `ProductionChemicalScaleScenario` applies pH-adjuster equivalents to a closed carbonate
// balance, consumes dissolved H2S by a scavenger's capacity, and evaluates the coupled
// mineral equilibrium before and after. The chemicals are `ProductionChemical` records; the
// probe builds them through the class's factory methods and one bespoke record, and prints
// the treated pH and bicarbonate, the dose equivalents, the scavenger capacity and residual
// H2S, the baseline and treated indices, the change, and the warning map.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar ProductionChemicalScenarioProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProductionChemicalScenarioProbe > captures/production_chemical_scenario_probe.tsv

import java.util.Map;
import neqsim.process.chemistry.ProductionChemical;
import neqsim.process.chemistry.scale.ProductionChemicalScaleScenario;
import neqsim.pvtsimulation.flowassurance.MultiMineralScaleEquilibrium;

public class ProductionChemicalScenarioProbe {

  static ProductionChemicalScaleScenario scenario() {
    return new ProductionChemicalScaleScenario()
        .setTemperatureCelsius(60.0).setPressureBara(100.0).setPH(6.5)
        .setCalciumMgL(2000.0).setBariumMgL(50.0).setStrontiumMgL(100.0)
        .setIronMgL(5.0).setMagnesiumMgL(500.0).setSodiumMgL(20000.0)
        .setSulphateMgL(500.0).setBicarbonateMgL(800.0)
        .setTotalDissolvedSolidsMgL(80000.0).setDissolvedH2SMgL(20.0)
        .setCO2PartialPressureBar(2.0)
        .setActivityModel(MultiMineralScaleEquilibrium.ActivityModel.DAVIES);
  }

  static void row(String label, ProductionChemicalScaleScenario scenario) {
    scenario.evaluate();
    System.out.println("  " + label + " treated_ph=" + scenario.getTreatedPH()
        + " residual_h2s_mg_l=" + scenario.getResidualDissolvedH2SMgL()
        + " delta_caco3=" + scenario.getSaturationIndexChange("CaCO3")
        + " delta_baso4=" + scenario.getSaturationIndexChange("BaSO4"));
    Map<String, Double> baseline =
        (Map<String, Double>) scenario.toMap().get("baselineSaturationIndex");
    for (String mineral : baseline.keySet()) {
      System.out.println("    " + mineral + " baseline_si=" + baseline.get(mineral)
          + " treated_si=" + scenario.getTreatedSaturationIndex(mineral) + " change="
          + scenario.getSaturationIndexChange(mineral));
    }
    for (Map.Entry<String, String> warning : scenario.getWarnings().entrySet()) {
      System.out.println("    warning=" + warning.getKey() + " message=" + warning.getValue());
    }
    Map<String, Object> chemistry =
        (Map<String, Object>) scenario.toMap().get("treatedChemistry");
    System.out.println("    chemistry treated_ph=" + chemistry.get("treatedPH")
        + " treated_bicarbonate_mg_l=" + chemistry.get("treatedBicarbonateMgL")
        + " dose_equivalents_mol_l=" + chemistry.get("baseDoseEquivalentsMolL")
        + " scavenger_capacity_mg_l=" + chemistry.get("scavengerCapacityH2SMgL")
        + " residual_h2s_mg_l=" + chemistry.get("residualDissolvedH2SMgL"));
  }

  public static void main(String[] args) {
    System.out.println("production_chemical_scenario");

    System.out.println(" ph_adjusters");
    row("caustic_100", scenario().addChemical(ProductionChemical.causticPHAdjuster("caustic", 100.0)));
    row("mdea_200", scenario().addChemical(ProductionChemical.phStabilizer("mdea", 200.0)));
    ProductionChemical hydrochloric = ProductionChemical.acid("acid", 15.0);
    hydrochloric.setDosagePpm(100.0);
    row("hcl_100", scenario().addChemical(hydrochloric));

    System.out.println(" scavengers");
    row("triazine_100", scenario().addChemical(ProductionChemical.h2sScavenger("scavenger", 100.0)));
    row("triazine_10", scenario().addChemical(ProductionChemical.h2sScavenger("scavenger", 10.0)));

    System.out.println(" inhibitors");
    row("scale_inhibitor_10", scenario().addChemical(ProductionChemical.scaleInhibitor("si", 10.0)));
    row("corrosion_inhibitor_10",
        scenario().addChemical(ProductionChemical.corrosionInhibitor("ci", 10.0)));

    System.out.println(" corners");
    row("no_bicarbonate", scenario().setBicarbonateMgL(0.0));
    row("huge_caustic",
        scenario().addChemical(ProductionChemical.causticPHAdjuster("caustic", 5000.0)));
    row("bdot", scenario().addChemical(ProductionChemical.phStabilizer("mdea", 50.0))
        .setActivityModel(MultiMineralScaleEquilibrium.ActivityModel.BDOT));
    row("pitzer",
        scenario().addChemical(ProductionChemical.phStabilizer("mdea", 50.0))
            .setActivityModel(MultiMineralScaleEquilibrium.ActivityModel.PITZER_BINARY)
            .setPitzerIonicStrengthMolal(1.2));

    ProductionChemical bespoke = new ProductionChemical("bespoke", ProductionChemical.ChemicalType.PH_ADJUSTER);
    bespoke.setActiveIngredient("KOH");
    bespoke.setActiveWtPct(45.0);
    bespoke.setDosagePpm(100.0);
    row("unsupported_adjuster", scenario().addChemical(bespoke));
  }
}
