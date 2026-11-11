//! `Component.createComponent`'s parachor selection, as upstream PR #4004 fixes it.

use azoth_eos::Cubic;
use azoth_eos::databank;
use azoth_eos::furst_electrolyte::electrolyte_cpa_statoil_mixture_of;

/// **A CPA component reads `PARACHOR_CPA`; a cubic one reads `PARACHOR`.**
///
/// `Component.createComponent` tests `this instanceof ComponentSrkCPA`, which before the fix
/// was two string comparisons and missed every subclass - so the electrolyte CPA components
/// read the ordinary column. Methanol is the substance the class's own test names: `73.0695`
/// on the cubic column and `87.543852818` on the CPA one, and the two mixtures below are the
/// two answers.
#[test]
fn cpa_components_read_the_cpa_parachor_column() {
    let names = ["methanol", "water"];
    let (cubic, _) = databank::mixture_of(&names, Cubic::Srk, None).expect("the cubic mixture");
    let (cpa, _) =
        databank::associating_mixture_of(&names, Cubic::Srk, None).expect("the CPA mixture");
    let (electrolyte, _) =
        electrolyte_cpa_statoil_mixture_of(&names, None).expect("the electrolyte CPA mixture");

    let parachors = |mixture: &azoth_eos::Mixture| -> Vec<f64> {
        mixture.components().iter().map(|c| c.parachor).collect()
    };
    let cubic = parachors(&cubic);
    assert!((cubic[0] - 73.069_489_067).abs() < 1.0e-9);
    assert!((cubic[1] - 52.878_622_6).abs() < 1.0e-9);

    // The fix's own number, and the water row beside it.
    let cpa = parachors(&cpa);
    assert!((cpa[0] - 87.543_852_818).abs() < 1.0e-9);
    assert!((cpa[1] - 52.419_867_266).abs() < 1.0e-9);

    // The electrolyte subclass reads the same column as the base CPA model, which is the
    // `instanceof` the fix restored.
    let electrolyte = parachors(&electrolyte);
    assert_eq!(electrolyte, cpa);
}
