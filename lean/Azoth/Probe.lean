import Azoth.Dim

open Units (Dimension)

example : Azoth.Dim.ofExponents [0, 0, 0, 0, 0, 0, 0] = (0 : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [1, 0, 0, 0, 0, 0, 0] = (Dimension.Length : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [1, 0, 0, 0, 0, 0, 0] = (Dimension.Length : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [2, 0, 0, 0, 0, 0, 0] = (Dimension.Area : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [3, 0, -1, 0, 0, 0, 0] = (Dimension.Volume / Dimension.Time : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [0, 1, -1, 0, 0, 0, 0] = (Dimension.Mass / Dimension.Time : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [0, 0, -1, 0, 0, 1, 0] = (Dimension.AmountOfSubstance / Dimension.Time : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [-3, 1, 0, 0, 0, 0, 0] = (Dimension.Mass / Dimension.Volume : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [1, 0, -1, 0, 0, 0, 0] = (Dimension.Speed : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [-1, 1, -2, 0, 0, 0, 0] = (Dimension.Pressure : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [-1, 1, -1, 0, 0, 0, 0] = (Dimension.Pressure * Dimension.Time : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [0, 0, 0, 0, 1, 0, 0] = (Dimension.Temperature : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [2, 1, -3, 0, 0, 0, 0] = (Dimension.Power : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [2, 0, -2, 0, -1, 0, 0] = (Dimension.Energy / (Dimension.Mass * Dimension.Temperature) : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [1, 1, -3, 0, -1, 0, 0] = (Dimension.Power / (Dimension.Length * Dimension.Temperature) : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [0, 1, -3, 0, -1, 0, 0] = (Dimension.Power / (Dimension.Area * Dimension.Temperature) : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [0, 1, 0, 0, 0, -1, 0] = (Dimension.Mass / Dimension.AmountOfSubstance : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [3, 0, 0, 0, 0, -1, 0] = (Dimension.Volume / Dimension.AmountOfSubstance : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [2, 1, -2, 0, 0, -1, 0] = (Dimension.Energy / Dimension.AmountOfSubstance : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [2, 1, -2, 0, -1, -1, 0] = (Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature) : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [2, 1, -2, 0, -2, -1, 0] = (Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature ^ 2) : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [2, 1, -2, 0, -3, -1, 0] = (Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature ^ 3) : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [2, 1, -2, 0, -4, -1, 0] = (Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature ^ 4) : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

example : Azoth.Dim.ofExponents [2, 1, -2, 0, -5, -1, 0] = (Dimension.Energy / (Dimension.AmountOfSubstance * Dimension.Temperature ^ 5) : Dimension) := by
  simp only [Azoth.Dim.ofExponents, Azoth.Dim.ofExponentsOn, Azoth.slots,
        Units.Dimension.Acceleration, Units.Dimension.AmountOfSubstance,
        Units.Dimension.Area, Units.Dimension.Energy, Units.Dimension.Force,
        Units.Dimension.Length, Units.Dimension.Mass, Units.Dimension.Power,
        Units.Dimension.Pressure, Units.Dimension.Speed, Units.Dimension.Temperature,
        Units.Dimension.Time, Units.Dimension.Volume, Units.Dimension.ofString,
        Units.Dimension.div_eq_sub, Units.Dimension.mul_eq_add,
        Units.Dimension.npow_eq_nsmul, sub_eq_add_neg,
        List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
        List.sum_cons, List.sum_nil] <;> simp <;> module

