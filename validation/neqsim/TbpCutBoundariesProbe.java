// `OilAssayCharacterisation.addTBPCutBoundariesKelvin` - a cumulative volume curve cut into an
// assay, one row per interval.
//
// The curve is stated as **boundary** yields rather than interval yields: N cumulative percentages
// and N boiling points give N-1 intervals, and the density vector is one shorter than either. Each
// interval's volume fraction is `(cumulative[i+1] - cumulative[i]) / 100`, its gravity is the one
// stated beside it, and the boiling point it is characterised by is the **midpoint** of the pair
// that bounds it - `withBoilingRangeKelvin` keeps all three, so the range is not lost to the
// average.
//
//   javac -proc:none -cp /path/to/neqsim-f0c7436.jar TbpCutBoundariesProbe.java
//   java -cp .:/path/to/neqsim-f0c7436.jar TbpCutBoundariesProbe > captures/tbp_cut_boundaries_probe.tsv

import java.util.List;
import neqsim.thermo.characterization.OilAssayCharacterisation;
import neqsim.thermo.characterization.OilAssayCharacterisation.AssayCut;
import neqsim.thermo.system.SystemSrkEos;

public class TbpCutBoundariesProbe {

  /** One curve, and the cuts the class makes of it. */
  private static void sweep(
      String label,
      String prefix,
      double[] cumulativeVolumePercent,
      double[] boilingPointKelvin,
      double[] specificGravity) {
    OilAssayCharacterisation assay =
        new OilAssayCharacterisation(new SystemSrkEos(298.15, 1.0));

    System.out.printf("# curve = %s%n", label);
    System.out.printf("prefix = %s%n", prefix);
    System.out.printf("cumulative_volume_percent = [%s]%n", vector(cumulativeVolumePercent));
    System.out.printf("boiling_point = [%s]%n", vector(boilingPointKelvin));
    System.out.printf("specific_gravity = [%s]%n", vector(specificGravity));

    try {
      assay.addTBPCutBoundariesKelvin(prefix, cumulativeVolumePercent, boilingPointKelvin, specificGravity);
    } catch (Exception ex) {
      System.out.printf("# refused: %s%n", ex.getMessage());
      return;
    }

    List<AssayCut> cuts = assay.getCuts();
    System.out.printf("cut_count = %d%n", cuts.size());
    StringBuilder names = new StringBuilder();
    StringBuilder fractions = new StringBuilder();
    StringBuilder gravities = new StringBuilder();
    StringBuilder lower = new StringBuilder();
    StringBuilder upper = new StringBuilder();
    StringBuilder average = new StringBuilder();
    for (int i = 0; i < cuts.size(); i++) {
      AssayCut cut = cuts.get(i);
      if (i > 0) {
        names.append(", ");
        fractions.append(", ");
        gravities.append(", ");
        lower.append(", ");
        upper.append(", ");
        average.append(", ");
      }
      names.append(cut.getName());
      fractions.append(String.format("%.15g", cut.getVolumeFraction()));
      gravities.append(String.format("%.15g", cut.resolveDensity()));
      lower.append(String.format("%.15g", cut.getLowerBoilingPointKelvin()));
      upper.append(String.format("%.15g", cut.getUpperBoilingPointKelvin()));
      average.append(String.format("%.15g", cut.resolveAverageBoilingPoint()));
    }
    System.out.printf("cut_names = [%s]%n", names);
    System.out.printf("volume_fraction = [%s]%n", fractions);
    System.out.printf("specific_gravity = [%s]%n", gravities);
    System.out.printf("lower_boiling_point = [%s]%n", lower);
    System.out.printf("upper_boiling_point = [%s]%n", upper);
    System.out.printf("average_boiling_point = [%s]%n", average);
  }

  private static String vector(double[] values) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < values.length; i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(String.format("%.15g", values[i]));
    }
    return out.toString();
  }

  /** The class's own refusals, each printed rather than thrown. */
  private static void refused(String label, double[] cumulative, double[] boiling, double[] gravity) {
    System.out.printf("# curve = %s%n", label);
    try {
      new OilAssayCharacterisation(new SystemSrkEos(298.15, 1.0))
          .addTBPCutBoundariesKelvin("TBP", cumulative, boiling, gravity);
      System.out.printf("# refused: nothing was thrown%n");
    } catch (Exception ex) {
      System.out.printf("# refused: %s%n", ex.getMessage());
    }
  }

  public static void main(String[] args) {
    // A five-boundary TBP curve: four cuts, four gravities, and the yields the intervals imply.
    sweep(
        "four_cuts",
        "TBP",
        new double[] {0.0, 20.0, 45.0, 75.0, 100.0},
        new double[] {320.0, 380.0, 440.0, 520.0, 620.0},
        new double[] {0.68, 0.76, 0.82, 0.88});

    // The fewest boundaries the class accepts: two yields, two temperatures, one gravity.
    sweep(
        "one_cut",
        "C",
        new double[] {0.0, 100.0},
        new double[] {350.0, 600.0},
        new double[] {0.80});

    // A prefix and a yield that are not whole numbers, so the naming and the arithmetic are both
    // visible rather than rounded away.
    sweep(
        "an_uneven_curve",
        "NAPHTHA_",
        new double[] {0.0, 12.5, 37.5, 100.0},
        new double[] {300.0, 350.0, 410.0, 480.0},
        new double[] {0.66, 0.74, 0.83});

    // **The refusals, one label each.** An N-length density vector where N-1 is required; a curve
    // that does not start at zero; yields that do not increase; a temperature that does not; and a
    // gravity that is not positive.
    refused(
        "too_many_gravities",
        new double[] {0.0, 50.0, 100.0},
        new double[] {320.0, 420.0, 520.0},
        new double[] {0.70, 0.80, 0.90});
    refused(
        "a_curve_that_does_not_start_at_zero",
        new double[] {5.0, 50.0, 100.0},
        new double[] {320.0, 420.0, 520.0},
        new double[] {0.70, 0.80});
    refused(
        "yields_that_do_not_increase",
        new double[] {0.0, 60.0, 50.0, 100.0},
        new double[] {320.0, 380.0, 440.0, 520.0},
        new double[] {0.70, 0.80, 0.90});
    refused(
        "a_temperature_that_does_not_increase",
        new double[] {0.0, 50.0, 100.0},
        new double[] {320.0, 300.0, 520.0},
        new double[] {0.70, 0.80});
    refused(
        "a_gravity_that_is_not_positive",
        new double[] {0.0, 50.0, 100.0},
        new double[] {320.0, 420.0, 520.0},
        new double[] {0.70, 0.0});
  }
}
