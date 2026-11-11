// `TubeInsertModel` driven directly, on the class's own factories and Reynolds boundaries.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar TubeInsertProbe.java
//     java -cp .:neqsim-f0c7436.jar TubeInsertProbe > captures/tube_insert_probe.tsv
//
// **The twisted tape's blockage factor cancels the tape thickness** (`t / (t / 0.05)` is
// `0.05` for any non-zero thickness), and the coiled wire's pitch cancels the same way, so
// those two parameters are carried and unread. The Re = 2300 and Re = 500 rows pin the
// branch boundaries, and the laminar rows pin that their correlations read neither Re nor Pr.

import neqsim.process.mechanicaldesign.heatexchanger.TubeInsertModel;

public class TubeInsertProbe {

  public static void main(String[] args) {
    // 1. The plain tube.
    row("none_default", new TubeInsertModel(), 10000.0, 5.0);

    // 2-7. Twisted tape across the regime boundary, the twist-ratio clamp and thickness.
    row("twisted_tape_turbulent", TubeInsertModel.createTwistedTape(5.0), 10000.0, 5.0);
    row("twisted_tape_laminar", TubeInsertModel.createTwistedTape(5.0), 1000.0, 5.0);
    row("twisted_tape_at_2300", TubeInsertModel.createTwistedTape(5.0), 2300.0, 5.0);
    row("twisted_tape_y2", TubeInsertModel.createTwistedTape(2.0), 10000.0, 5.0);
    row("twisted_tape_y1_clamped", TubeInsertModel.createTwistedTape(1.0), 10000.0, 5.0);
    TubeInsertModel thickTape = TubeInsertModel.createTwistedTape(5.0);
    thickTape.setTapeThickness(0.003);
    row("twisted_tape_thick", thickTape, 10000.0, 5.0);

    // 8-12. Wire matrix across its two boundaries and at two densities.
    row("wire_matrix_turbulent", TubeInsertModel.createWireMatrix(0.5), 10000.0, 5.0);
    row("wire_matrix_at_2300", TubeInsertModel.createWireMatrix(0.5), 2300.0, 5.0);
    row("wire_matrix_transition", TubeInsertModel.createWireMatrix(0.5), 1000.0, 5.0);
    row("wire_matrix_at_500", TubeInsertModel.createWireMatrix(0.5), 500.0, 5.0);
    row("wire_matrix_laminar", TubeInsertModel.createWireMatrix(0.5), 400.0, 5.0);
    row("wire_matrix_dense", TubeInsertModel.createWireMatrix(0.7), 10000.0, 5.0);

    // 13-16. Coiled wire across the boundary, two roughnesses, and a zero Reynolds number.
    row("coiled_wire_turbulent", TubeInsertModel.createCoiledWire(0.02, 45.0), 10000.0, 5.0);
    row("coiled_wire_laminar", TubeInsertModel.createCoiledWire(0.02, 45.0), 1000.0, 5.0);
    row("coiled_wire_rough", TubeInsertModel.createCoiledWire(0.05, 60.0), 10000.0, 5.0);
    row("coiled_wire_pr_two", TubeInsertModel.createCoiledWire(0.02, 45.0), 10000.0, 2.0);
    row("twisted_tape_zero_re", TubeInsertModel.createTwistedTape(5.0), 0.0, 5.0);

    System.out.flush();
  }

  static void row(String label, TubeInsertModel model, double re, double pr) {
    double hRatio = model.getHeatTransferEnhancementRatio(re, pr);
    double fRatio = model.getPressureDropPenaltyRatio(re);
    double pec = model.getPerformanceEvaluationCriteria(re, pr);
    double[] applied = model.applyEnhancement(1000.0, 5000.0, re, pr);
    System.out.println(label);
    System.out.println("insert_type=" + model.toMap(re, pr).get("insertType"));
    System.out.println("reynolds_number=" + re);
    System.out.println("prandtl_number=" + pr);
    System.out.println("heat_transfer_enhancement_ratio=" + hRatio);
    System.out.println("pressure_drop_penalty_ratio=" + fRatio);
    System.out.println("performance_evaluation_criteria=" + pec);
    System.out.println("enhanced_heat_transfer_coefficient=" + applied[0]);
    System.out.println("enhanced_pressure_drop=" + applied[1]);
    System.out.println();
  }
}
