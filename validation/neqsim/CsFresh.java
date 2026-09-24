import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class CsFresh {
  public static void main(String[] args) {
    double[][] rows = {
      { 0.962671905697446, 0.029469548133595282, 0.007858546168958742 },
      { 0.02036659877800409, 0.5804480651731161, 0.39918533604887985 },
    };
    String[] labels = { "overhead", "bottoms" };
    for (int r = 0; r < 2; r++) {
      SystemInterface f = new SystemPrEos(300.0, 20.0);
      f.addComponent("methane", rows[r][0]);
      f.addComponent("n-butane", rows[r][1]);
      f.addComponent("n-pentane", rows[r][2]);
      f.setMixingRule(2);
      f.setTotalFlowRate(1.0, "mol/sec");
      f.init(0);
      new ThermodynamicOperations(f).TPflash();
      f.initProperties();
      System.out.println(labels[r] + " phases=" + f.getNumberOfPhases() + " h="
          + f.getEnthalpy() / f.getTotalNumberOfMoles() + " z=" + f.getPhase(0).getZ());
    }
  }
}
