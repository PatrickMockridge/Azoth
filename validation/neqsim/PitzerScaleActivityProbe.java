// `PitzerScaleActivityModel`: the Hamer-Wu 1972 NaCl table, the eight salts' parameter sets
// (read by reflection, so the transcription surface is visible), the nine ion coefficients at
// a strength grid, and the two products the equilibrium solver reads.
//
//     javac -proc:none -cp neqsim-f0c7436.jar PitzerScaleActivityProbe.java
//     java -cp .:neqsim-f0c7436.jar PitzerScaleActivityProbe > captures/pitzer_scale_activity_probe.tsv

import java.lang.reflect.Field;
import java.util.Arrays;
import neqsim.pvtsimulation.flowassurance.PitzerScaleActivityModel;
import neqsim.pvtsimulation.flowassurance.ScalePredictionCalculator;

public class PitzerScaleActivityProbe {

  /** Hamer and Wu (1972) Table 16, the 29 points NeqSim's own test reads. */
  static final double[][] HAMER_WU = {
      {0.001, 0.965}, {0.002, 0.952}, {0.005, 0.928}, {0.010, 0.903}, {0.020, 0.872},
      {0.050, 0.822}, {0.100, 0.779}, {0.200, 0.734}, {0.300, 0.709}, {0.400, 0.693},
      {0.500, 0.681}, {0.600, 0.673}, {0.700, 0.667}, {0.800, 0.662}, {0.900, 0.659},
      {1.000, 0.657}, {1.200, 0.655}, {1.400, 0.656}, {1.600, 0.658}, {1.800, 0.662},
      {2.000, 0.668}, {2.500, 0.688}, {3.000, 0.714}, {3.500, 0.746}, {4.000, 0.783},
      {4.500, 0.826}, {5.000, 0.874}, {5.500, 0.928}, {6.000, 0.986},
  };

  public static void main(String[] args) throws Exception {
    ScalePredictionCalculator predictor = new ScalePredictionCalculator();
    predictor.setTemperatureCelsius(25.0);
    double a = predictor.getDebyeHuckelAParameter();
    System.out.println("debye_huckel_a_25c=" + a);

    hamerWu(a);
    parameterTable();
    ionGrid(a);
    meanAtMolalities(a);
    products(a);
  }

  static void hamerWu(double a) {
    double maximumAbsoluteError = 0.0;
    double sumAbsoluteError = 0.0;
    double pitzerSquared = 0.0;
    double daviesSquared = 0.0;
    double bdotSquared = 0.0;
    int highSalinityCount = 0;
    System.out.println("hamer_wu_table");
    for (double[] point : HAMER_WU) {
      double molality = point[0];
      double reference = point[1];
      PitzerScaleActivityModel model = new PitzerScaleActivityModel(298.15, molality, a);
      double calculated = model.getMeanActivityCoefficient(PitzerScaleActivityModel.Salt.NACL, molality);
      double error = Math.abs(calculated - reference);
      maximumAbsoluteError = Math.max(maximumAbsoluteError, error);
      sumAbsoluteError += error;
      if (molality >= 0.5) {
        double sqrtI = Math.sqrt(molality);
        double davies = Math.pow(10.0, -a * (sqrtI / (1.0 + sqrtI) - 0.3 * molality));
        double bdot = Math.pow(10.0, -a * sqrtI / (1.0 + 5.0 * 0.3283 * sqrtI) + 0.041 * molality);
        pitzerSquared += error * error;
        daviesSquared += Math.pow(davies - reference, 2.0);
        bdotSquared += Math.pow(bdot - reference, 2.0);
        highSalinityCount++;
      }
      System.out.println("  molality=" + molality + " reference=" + reference + " calculated="
          + calculated + " error=" + error);
    }
    System.out.println("hamer_wu_count=" + HAMER_WU.length);
    System.out.println("hamer_wu_max_abs_error=" + maximumAbsoluteError);
    System.out.println("hamer_wu_mean_abs_error=" + (sumAbsoluteError / HAMER_WU.length));
    System.out.println("hamer_wu_pitzer_rmse=" + Math.sqrt(pitzerSquared / highSalinityCount));
    System.out.println("hamer_wu_davies_rmse=" + Math.sqrt(daviesSquared / highSalinityCount));
    System.out.println("hamer_wu_bdot_rmse=" + Math.sqrt(bdotSquared / highSalinityCount));
    System.out.println();
  }

  static void parameterTable() throws Exception {
    System.out.println("salt_parameters");
    for (PitzerScaleActivityModel.Salt salt : PitzerScaleActivityModel.Salt.values()) {
      StringBuilder out = new StringBuilder("  salt=").append(salt.name());
      for (String field : new String[] {"cationStoichiometry", "anionStoichiometry", "cationCharge",
          "anionCharge", "alpha1", "alpha2"}) {
        out.append(' ').append(field).append('=').append(value(salt, field));
      }
      for (String field : new String[] {"beta0", "beta1", "beta2", "cPhi"}) {
        out.append(' ').append(field).append('=').append(Arrays.toString((double[]) value(salt, field)));
      }
      System.out.println(out);
    }
    System.out.println();
  }

  static void ionGrid(double a) {
    double[] strengths = {1.0e-12, 0.001, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0};
    System.out.println("ion_coefficients");
    for (double strength : strengths) {
      PitzerScaleActivityModel model = new PitzerScaleActivityModel(298.15, strength, a);
      StringBuilder out = new StringBuilder("  ionic_strength=").append(strength);
      for (PitzerScaleActivityModel.Ion ion : PitzerScaleActivityModel.Ion.values()) {
        out.append(' ').append(ion.name()).append('=')
            .append(model.getIonActivityCoefficient(ion));
      }
      System.out.println(out);
    }
    System.out.println();
  }

  static void meanAtMolalities(double a) {
    double[] strengths = {0.0, 1.0, 3.0};
    System.out.println("salt_mean_coefficients");
    for (double strength : strengths) {
      PitzerScaleActivityModel model = new PitzerScaleActivityModel(298.15, strength, a);
      StringBuilder out = new StringBuilder("  ionic_strength=").append(strength);
      for (PitzerScaleActivityModel.Salt salt : PitzerScaleActivityModel.Salt.values()) {
        out.append(' ').append(salt.name()).append('=')
            .append(model.getMeanActivityCoefficient(salt, strength == 0.0 ? 1.0 : strength));
      }
      System.out.println(out);
    }
    System.out.println("  zero_molality=" + new PitzerScaleActivityModel(298.15, 0.0, a)
        .getMeanActivityCoefficient(PitzerScaleActivityModel.Salt.NACL, 0.0));
    System.out.println();
  }

  static void products(double a) {
    System.out.println("activity_products");
    for (double strength : new double[] {1.0, 2.0}) {
      PitzerScaleActivityModel model = new PitzerScaleActivityModel(298.15, strength, a);
      System.out.println("  ionic_strength=" + strength + " baso4="
          + model.getActivityCoefficientProduct(PitzerScaleActivityModel.Ion.BARIUM,
              PitzerScaleActivityModel.Ion.SULPHATE)
          + " caco3="
          + model.getActivityCoefficientProduct(PitzerScaleActivityModel.Ion.CALCIUM,
              PitzerScaleActivityModel.Ion.CARBONATE));
    }
  }

  static Object value(Object target, String name) throws Exception {
    Field field = PitzerScaleActivityModel.Salt.class.getDeclaredField(name);
    field.setAccessible(true);
    return field.get(target);
  }
}
