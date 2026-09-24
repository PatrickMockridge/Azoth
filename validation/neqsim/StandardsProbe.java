// ISO 6976 gas quality, for the `standards.iso6976` id.
//
//     javac -proc:none -cp neqsim-f0c7436.jar StandardsProbe.java
//     java -cp .:neqsim-f0c7436.jar StandardsProbe > captures/process_gas_quality.tsv
//
// The only file in the repository that is not part of either implementation, and here on
// purpose: it is the ground truth the port is measured against.
//
// # What it drives
//
// `Standard_ISO6976` - **the original class, not `Standard_ISO6976_2016`** - because that is
// the one `Stream.LCV()` builds and the one `unit_ops.flare`'s duty comes from. It queries
// NeqSim's `ISO6976constants` table, which is the file
// `databank/sources/neqsim/ISO6976constants.csv`, compiled here to
// `data/standards/iso6976.csv`.
//
// The reference type is `molar`, so the calorific values come back per mole rather than per
// cubic metre: the port reports molar quantities and composes the volumetric ones from them,
// and the two are comparable only if the probe states which one it asked for. The
// compression factor, the relative density and the two densities are returned before that
// scaling, so the reference type does not move them.
//
// The reference state is `real`, which is what `Stream.LCV()` sets.
//
// # The rows
//
// Four gases at the 15 C metering and 25 C combustion references the standard's own examples
// use, and the flare's pair - 0 C and 15.55 C - on one of them, because that is the call
// `Stream.LCV()` makes and the one the port's consumer actually reads.

import neqsim.standards.gasquality.Standard_ISO6976;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

public class StandardsProbe {
  public static void main(String[] args) {
    // **The ISO 6976 constants are a database table, not a component property.** NeqSim's
    // tables are built from its bundled CSVs on demand, and `createTemporaryTables` is what
    // asks for them; without it the lookup returns no rows and the standard answers zero.
    neqsim.util.database.NeqSimDataBase.createTemporaryTables();
    String[] methane = new String[] { "methane" };
    double[] methaneZ = new double[] { 1.0 };

    String[] binary = new String[] { "methane", "n-butane" };
    double[] binaryZ = new double[] { 0.9, 0.1 };

    String[] sales = new String[] { "methane", "ethane", "propane", "n-butane", "nitrogen" };
    double[] salesZ = new double[] { 0.85, 0.07, 0.03, 0.02, 0.03 };

    String[] sour = new String[] { "methane", "CO2", "nitrogen" };
    double[] sourZ = new double[] { 0.8, 0.15, 0.05 };

    row("pure_methane_15_25", methane, methaneZ, 15.0, 25.0);
    row("methane_butane_15_25", binary, binaryZ, 15.0, 25.0);
    row("sales_gas_15_25", sales, salesZ, 15.0, 25.0);
    row("sour_gas_15_25", sour, sourZ, 15.0, 25.0);
    // **The flare's own pair**, on the gas a flare would carry: `Stream.LCV()` builds
    // `new Standard_ISO6976(fluid, 0, 15.55, "volume")`.
    row("methane_butane_0_15_55", binary, binaryZ, 0.0, 15.55);
    row("sales_gas_0_15_55", sales, salesZ, 0.0, 15.55);
    row("pure_methane_0_15_55", methane, methaneZ, 0.0, 15.55);
  }

  static void row(String label, String[] names, double[] z, double volumetricC, double energyC) {
    SystemInterface fluid = new SystemPrEos(288.15, 1.01325);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], z[i]);
    }
    fluid.setMixingRule(2);
    // **The fractions the standard reads are the phase's, and a fluid that has not been
    // initialised has none**: `getComponent(i).getz()` is zero until `init(0)` fills it, and
    // the standard then sums zeros and answers a compression factor of one. NeqSim's own
    // test does not call this because its fluid has been through a flash first.
    fluid.init(0);

    Standard_ISO6976 standard = new Standard_ISO6976(fluid.clone(), volumetricC, energyC, "molar");
    standard.setReferenceState("real");
    standard.calculate();

    System.out.println(label);
    System.out.println("volumetric_reference_C=" + volumetricC);
    System.out.println("energy_reference_C=" + energyC);
    // **The class's own units**: molar mass in g/mol and the calorific values in kJ/mol,
    // because `getValue` returns the table's numbers unscaled. The compiled table is in SI
    // and the port's outputs are too, so the cases divide the first by 1000 and multiply the
    // second by 1000 - stated here because a capture read beside a case would otherwise look
    // off by three orders of magnitude.
    System.out.println("molar_mass_g_per_mol=" + standard.getValue("MolarMass"));
    System.out.println("compression_factor=" + standard.getValue("CompressionFactor"));
    System.out.println("relative_density=" + standard.getValue("RelativeDensity"));
    System.out.println("density_ideal=" + standard.getValue("DensityIdeal"));
    System.out.println("density_real=" + standard.getValue("DensityReal"));
    System.out.println("superior_calorific_value_kj_per_mol=" + standard.getValue("SuperiorCalorificValue"));
    System.out.println("inferior_calorific_value_kj_per_mol=" + standard.getValue("InferiorCalorificValue"));
    System.out.println();
  }
}
