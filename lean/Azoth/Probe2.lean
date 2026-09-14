import Azoth.Dim
import LeanUnits.Framework.Dimensions.Lemmas

/-- A theorem that leans on the vendored library's own lemma. -/
theorem via_vendor (d : Units.Dimension) : d * 1 = d := Units.Dimension.mul_one d
