// `PackingHydraulicsCalculator` driven directly, on states of its own.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar PackingProbe.java
//     java -cp .:neqsim-f0c7436.jar PackingProbe > captures/packing_probe.tsv
//
// **This is the rate-based column's hydraulics, and the equilibrium packed column's capture already
// holds its output** - `process_packed_column.tsv` prints `kga`, `kla`, the wetted area, the
// percent flood, the flooding velocity and the pressure drop of a solved column. But those rows do
// not print the *inputs* the calculator reads: the two mass flows, the four transport properties,
// the two diffusivities and the packing's geometry. So this probe states them and prints both
// sides, which is what makes the port's case an oracle rather than a self-consistency check.
//
// **Three of the four correlations are Onda's, and the fourth is Leva's.** The wetted area, `kGa`
// and `kLa` are Onda (1968); the bed pressure drop is Leva's; the flooding velocity is the Eckert
// GPDC fit; and the HETP is a two-resistance combination of the two film coefficients with a
// stripping factor of one.

import neqsim.process.equipment.distillation.internals.PackingHydraulicsCalculator;
import neqsim.process.equipment.distillation.internals.PackingSpecificationLibrary;

public class PackingProbe {

  public static void main(String[] args) {
    // `size` is the second capture this probe writes - the sizing path, which no row below takes.
    if (args.length > 0 && args[0].equals("size")) {
      sizeRows();
      return;
    }
    // The class's own default packing, at the states the rate-based column's tests run: a
    // 1 m column, the absorber's flows and transport properties.
    row("pall_ring_50_absorber", "Pall-Ring-50", 1.0, 5.0, 0.35, 3.5, 45.0, 990.0, 1.8e-5, 6.5e-4,
        0.072, 3.3e-7, 1.9e-9);
    // A structured packing, which is the other category and a different geometry.
    row("mellapak_250y_absorber", "Mellapak-250Y", 1.2, 8.0, 0.28, 4.0, 30.0, 1000.0, 1.7e-5,
        5.5e-4, 0.065, 4.0e-7, 2.2e-9);
    // A higher liquid load, which is what moves the wetted area and the flood.
    row("pall_ring_50_high_liquid", "Pall-Ring-50", 1.0, 5.0, 0.35, 12.0, 45.0, 990.0, 1.8e-5,
        6.5e-4, 0.072, 3.3e-7, 1.9e-9);
    // A taller bed at the same state, which the HETP's theoretical stages read.
    row("pall_ring_50_tall", "Pall-Ring-50", 1.0, 20.0, 0.35, 3.5, 45.0, 990.0, 1.8e-5, 6.5e-4,
        0.072, 3.3e-7, 1.9e-9);
  }

  /// **The sizing path, which no row above takes.** `sizeColumnDiameter` floods a **1.0 m** trial,
  /// applies `designFloodFraction` to get a design velocity, sizes the area from the vapour's
  /// volumetric flow, and rounds up to the class's own table of standard diameters. The four rows
  /// above *set* a diameter and the packed column's capture reaches this only through a solved
  /// middle tray, so nothing in this directory measures the arithmetic itself. Every state here
  /// prints the trial flooding velocity the class leaves in its own field and the diameter it
  /// answers with; the design velocity, the area and the unrounded diameter are derived from
  /// those two and are labelled as derived.
  static void sizeRows() {
    sizeRow("pall_ring_50_absorber_sized", "Pall-Ring-50", 5.0, 0.35, 3.5, 45.0, 990.0, 1.8e-5,
        6.5e-4, 0.072, 3.3e-7, 1.9e-9, 0.70);
    sizeRow("mellapak_250y_absorber_sized", "Mellapak-250Y", 8.0, 0.28, 4.0, 30.0, 1000.0, 1.7e-5,
        5.5e-4, 0.065, 4.0e-7, 2.2e-9, 0.70);
    sizeRow("pall_ring_50_high_liquid_sized", "Pall-Ring-50", 5.0, 0.35, 12.0, 45.0, 990.0,
        1.8e-5, 6.5e-4, 0.072, 3.3e-7, 1.9e-9, 0.70);
    // A heavy vapour at a low flow, which sizes small enough to reach the table's floor.
    sizeRow("pall_ring_50_light_load_sized", "Pall-Ring-50", 5.0, 0.02, 3.5, 45.0, 990.0, 1.8e-5,
        6.5e-4, 0.072, 3.3e-7, 1.9e-9, 0.70);
    // **The flood fraction is the one knob the sizing reads and the calculation does not**, so a
    // second value of it on the same state is the whole measurement of that claim.
    sizeRow("pall_ring_50_absorber_sized_half_flood", "Pall-Ring-50", 5.0, 0.35, 3.5, 45.0, 990.0,
        1.8e-5, 6.5e-4, 0.072, 3.3e-7, 1.9e-9, 0.50);
  }

  static void sizeRow(String label, String packing, double packedHeightM,
      double vaporMassFlowKgPerS, double liquidMassFlowKgPerS, double vaporDensity,
      double liquidDensity, double vaporViscosity, double liquidViscosity, double surfaceTension,
      double vaporDiffusivity, double liquidDiffusivity, double designFloodFraction) {
    PackingHydraulicsCalculator hydraulics = new PackingHydraulicsCalculator();
    hydraulics.setPackingSpecification(PackingSpecificationLibrary.getOrDefault(packing));
    hydraulics.setPackedHeight(packedHeightM);
    hydraulics.setVaporMassFlow(vaporMassFlowKgPerS);
    hydraulics.setLiquidMassFlow(liquidMassFlowKgPerS);
    hydraulics.setVaporDensity(vaporDensity);
    hydraulics.setLiquidDensity(liquidDensity);
    hydraulics.setVaporViscosity(vaporViscosity);
    hydraulics.setLiquidViscosity(liquidViscosity);
    hydraulics.setSurfaceTension(surfaceTension);
    hydraulics.setVaporDiffusivity(vaporDiffusivity);
    hydraulics.setLiquidDiffusivity(liquidDiffusivity);
    hydraulics.setDesignFloodFraction(designFloodFraction);
    double sized = hydraulics.sizeColumnDiameter();

    System.out.println(label);
    System.out.println("packing=" + packing);
    System.out.println("packed_height_m=" + packedHeightM);
    System.out.println("vapor_mass_flow_kg_per_s=" + vaporMassFlowKgPerS);
    System.out.println("liquid_mass_flow_kg_per_s=" + liquidMassFlowKgPerS);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("vapor_viscosity_Pa_s=" + vaporViscosity);
    System.out.println("liquid_viscosity_Pa_s=" + liquidViscosity);
    System.out.println("surface_tension_N_per_m=" + surfaceTension);
    System.out.println("design_flood_fraction=" + designFloodFraction);
    System.out.println("trial_column_diameter_m=" + 1.0);
    // **The class's own two fields after the call**: the flood it computed at the trial and the
    // diameter it answered with. Everything below them is derived from the pair and says so.
    System.out.println("trial_flooding_velocity_m_per_s=" + hydraulics.getFloodingVelocity());
    System.out.println("sized_diameter_m=" + sized);
    double designVelocity = hydraulics.getFloodingVelocity() * designFloodFraction;
    double vaporVolumetricFlow = vaporMassFlowKgPerS / Math.max(vaporDensity, 0.01);
    System.out.println("derived_design_velocity_m_per_s=" + designVelocity);
    System.out.println("derived_vapor_volumetric_flow_m3_per_s=" + vaporVolumetricFlow);
    System.out.println(
        "derived_required_area_m2=" + (designVelocity <= 0 ? 1.0 : vaporVolumetricFlow / designVelocity));
    System.out.println("derived_required_diameter_m="
        + (designVelocity <= 0 ? 1.0
            : Math.sqrt(4.0 * (vaporVolumetricFlow / designVelocity) / Math.PI)));
    System.out.println();
  }

  static void row(String label, String packing, double diameterM, double packedHeightM,
      double vaporMassFlowKgPerS, double liquidMassFlowKgPerS, double vaporDensity,
      double liquidDensity, double vaporViscosity, double liquidViscosity, double surfaceTension,
      double vaporDiffusivity, double liquidDiffusivity) {
    PackingHydraulicsCalculator hydraulics = new PackingHydraulicsCalculator();
    hydraulics.setPackingSpecification(PackingSpecificationLibrary.getOrDefault(packing));
    hydraulics.setColumnDiameter(diameterM);
    hydraulics.setPackedHeight(packedHeightM);
    hydraulics.setVaporMassFlow(vaporMassFlowKgPerS);
    hydraulics.setLiquidMassFlow(liquidMassFlowKgPerS);
    hydraulics.setVaporDensity(vaporDensity);
    hydraulics.setLiquidDensity(liquidDensity);
    hydraulics.setVaporViscosity(vaporViscosity);
    hydraulics.setLiquidViscosity(liquidViscosity);
    hydraulics.setSurfaceTension(surfaceTension);
    hydraulics.setVaporDiffusivity(vaporDiffusivity);
    hydraulics.setLiquidDiffusivity(liquidDiffusivity);
    hydraulics.calculate();

    System.out.println(label);
    System.out.println("packing=" + packing);
    System.out.println("column_diameter_m=" + diameterM);
    System.out.println("packed_height_m=" + packedHeightM);
    System.out.println("vapor_mass_flow_kg_per_s=" + vaporMassFlowKgPerS);
    System.out.println("liquid_mass_flow_kg_per_s=" + liquidMassFlowKgPerS);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("vapor_viscosity_Pa_s=" + vaporViscosity);
    System.out.println("liquid_viscosity_Pa_s=" + liquidViscosity);
    System.out.println("surface_tension_N_per_m=" + surfaceTension);
    System.out.println("vapor_diffusivity_m2_per_s=" + vaporDiffusivity);
    System.out.println("liquid_diffusivity_m2_per_s=" + liquidDiffusivity);
    // The packing the library resolved, so the port's own table is held to the same numbers.
    System.out.println("packing_name=" + hydraulics.getPackingName());
    System.out.println("category=" + hydraulics.getPackingCategory());
    System.out.println("specific_surface_area=" + hydraulics.getSpecificSurfaceArea());
    System.out.println("void_fraction=" + hydraulics.getVoidFraction());
    System.out.println("packing_factor=" + hydraulics.getPackingFactor());
    System.out.println("fs_factor=" + hydraulics.getFsFactor());
    System.out.println("flooding_velocity_m_per_s=" + hydraulics.getFloodingVelocity());
    System.out.println("percent_flood=" + hydraulics.getPercentFlood());
    System.out.println("packing_pressure_drop_Pa=" + hydraulics.getPressureDropPerMeter());
    System.out.println("wetted_area_m2_per_m3=" + hydraulics.getWettedArea());
    System.out.println("kga=" + hydraulics.getKGa());
    System.out.println("kla=" + hydraulics.getKLa());
    System.out.println("hetp_m=" + hydraulics.getHETP());
    System.out.println("theoretical_stages=" + hydraulics.getNumberOfTheoreticalStages());
    System.out.println("htu_g_m=" + hydraulics.getHtuG());
    System.out.println("htu_l_m=" + hydraulics.getHtuL());
    System.out.println("htu_og_m=" + hydraulics.getHtuOG());
    System.out.println("hydraulics_ok=" + hydraulics.isDesignOk());
    System.out.println("wetting_ok=" + hydraulics.isWettingOk());
    System.out.println();
  }
}
