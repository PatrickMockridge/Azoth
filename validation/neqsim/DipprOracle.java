import neqsim.thermo.system.SystemSrkEos;

public class DipprOracle {
  public static void main(String[] args) {
    double[] temps = {248.15, 298.15, 350.0};
    String[] names = {"i-pentane", "propanePVTsim", "nbutanePVTsim", "Piperazine", "COS"};
    SystemSrkEos s = new SystemSrkEos(300.0, 1.0);
    for (String n : names) s.addComponent(n, 1.0);
    s.createDatabase(true);
    for (double t : temps) {
      for (int i = 0; i < s.getPhase(0).getNumberOfComponents(); i++) {
        var c = s.getPhase(0).getComponent(i);
        System.out.printf("T=%.2f %-16s %.17g bar%n", t, c.getName(), c.getAntoineVaporPressure(t));
      }
    }
  }
}
