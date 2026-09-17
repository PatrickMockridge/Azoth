import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Drives NeqSim's `VUflashSingleComp` on the states `eos.vu_flash_single_comp` carries cases for.
 *
 * <p>A pure component in the two-phase region: the pressure fixes the saturation line, and an
 * internal energy between the saturated liquid's and the vapour's fixes the split. The probe
 * finds the saturation temperature itself, prints the two saturated internal energies, and then
 * hands the halfway value back to `ThermodynamicOperations.VUflash` - which routes a pure
 * component to `VUflashSingleComp` - and prints the temperature and vapour fraction it returns.
 *
 * <p>The `Vspec` NeqSim takes is never read by that class; the volume follows from the pressure
 * and the internal energy. The probe passes one anyway, because the signature requires it.
 */
public class VuSingleCompProbe {

  private static void one(String name, double pressureBar) {
    SystemInterface sat = new SystemPrEos(300.0, pressureBar);
    sat.addComponent(name, 1.0);
    sat.setMixingRule("classic");
    sat.setAttractiveTerm(1);
    try {
      new ThermodynamicOperations(sat).bubblePointTemperatureFlash();
    } catch (Exception e) {
      System.out.printf("%-10s %6.1f bar  bubblePointTemperatureFlash THREW %s%n", name, pressureBar,
          e.getClass().getSimpleName());
      return;
    }
    double tsat = sat.getTemperature();
    sat.init(2);
    double uGas = sat.getPhase(0).getInternalEnergy("J/mol");
    double uLiq = sat.getPhase(1).getInternalEnergy("J/mol");
    double uSpec = 0.5 * (uLiq + uGas);
    System.out.printf("%-10s %6.1f bar  Tsat=%.12f  u_liq=%.9f  u_gas=%.9f  u_spec=%.9f%n", name,
        pressureBar, tsat, uLiq, uGas, uSpec);

    SystemInterface flash = new SystemPrEos(298.0, pressureBar);
    flash.addComponent(name, 1.0);
    flash.setMixingRule("classic");
    flash.setAttractiveTerm(1);
    try {
      new ThermodynamicOperations(flash).VUflash(1.0, uSpec);
      System.out.printf("            -> T=%.12f beta=%.12f P=%.6f bar%n", flash.getTemperature(),
          flash.getBeta(), flash.getPressure());
    } catch (Exception e) {
      System.out.printf("            -> VUflash THREW %s: %s%n", e.getClass().getSimpleName(),
          e.getMessage());
    }
  }

  public static void main(String[] args) {
    one("propane", 10.0);
    one("propane", 20.0);
    one("n-butane", 2.0);
  }
}
