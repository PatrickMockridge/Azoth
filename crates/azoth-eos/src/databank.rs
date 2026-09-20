//! The component databank: a mixture built from substance names.
//!
//! `data/components/components.csv` and `kij.csv` are generated from NeqSim's `COMP.csv`
//! and `INTER.csv` by `tools/gen_databank.py`, embedded with `include_str!` so a wheel
//! cannot find a different file or none. The Python side opens the same file from the
//! repository and `python/tests/test_data_agreement.py` compares the embedded bytes
//! against it.
//!
//! A keycard overrides the embedded tables by name, parameter by parameter: a card naming
//! only `omega` keeps the shipped `Tc` and `Pc`, and a name the databank does not have is
//! added - with no heat-capacity coefficients, because a card supplies the parameters a
//! cubic needs. An overlay is a value a caller passes and not a store; [`crate::card`]
//! reads a file and produces one.

use crate::Alpha;
use std::collections::HashMap;
use std::sync::OnceLock;

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, Result};

use crate::association::{AssociationCubic, AssociationRecord, SiteScheme, UmrCpaRecord};
use crate::bwrs::BwrsCoefficients;
use crate::cubic::Cubic;
use crate::mixing_rule::MixingRule;
use crate::mixture::{Component, Mixture};
use crate::molar_enthalpy_entropy::IdealGasModel;

/// The compiled component table, generated from NeqSim's `COMP.csv`.
const COMPONENTS_CSV: &str = include_str!("../../../data/components/components.csv");

/// The compiled interaction table, generated from NeqSim's `INTER.csv`.
const KIJ_CSV: &str = include_str!("../../../data/components/kij.csv");

/// The compiled UNIFAC tables, generated from NeqSim's `UNIFACcomp.csv`,
/// `UNIFACGroupParam.csv` and `UNIFACInterParam.csv`.
const UNIFAC_COMP_CSV: &str = include_str!("../../../data/components/UNIFACcomp.csv");
const UNIFAC_GROUP_CSV: &str = include_str!("../../../data/components/UNIFACGroupParam.csv");
const UNIFAC_INTER_CSV: &str = include_str!("../../../data/components/UNIFACInterParam.csv");

/// The compiled UNIFAC-PSRK interaction tables, generated from NeqSim's
/// `UNIFACInterParamB.csv` and `UNIFACInterParamC.csv`, the `b` and `c` of
/// `a_mn(T) = a + b T + c T^2`.
const UNIFAC_INTER_B_CSV: &str = include_str!("../../../data/components/UNIFACInterParamB.csv");
const UNIFAC_INTER_C_CSV: &str = include_str!("../../../data/components/UNIFACInterParamC.csv");

/// The compiled UNIFAC-UMR-PRU tables, generated from NeqSim's `UNIFACcompUMRPRU.csv`
/// and the six `UNIFACInterParam{A,B,C}_UMR{,MC}.csv` files.
const UNIFAC_COMP_UMRPRU_CSV: &str = include_str!("../../../data/components/UNIFACcompUMRPRU.csv");
const UNIFAC_A_UMR_CSV: &str = include_str!("../../../data/components/UNIFACInterParamA_UMR.csv");
const UNIFAC_A_UMRMC_CSV: &str =
    include_str!("../../../data/components/UNIFACInterParamA_UMRMC.csv");
const UNIFAC_B_UMR_CSV: &str = include_str!("../../../data/components/UNIFACInterParamB_UMR.csv");
const UNIFAC_B_UMRMC_CSV: &str =
    include_str!("../../../data/components/UNIFACInterParamB_UMRMC.csv");
const UNIFAC_C_UMR_CSV: &str = include_str!("../../../data/components/UNIFACInterParamC_UMR.csv");
const UNIFAC_C_UMRMC_CSV: &str =
    include_str!("../../../data/components/UNIFACInterParamC_UMRMC.csv");

/// The compiled MBWR-32 coefficient table, generated from NeqSim's `MBWR32param.csv`.
const MBWR32_CSV: &str = include_str!("../../../data/components/mbwr32.csv");

/// The compiled Pitzer pair parameters, generated from NeqSim's `PitzerParameters.csv`.
///
/// The electrolyte models' second table: `COMP.csv` says what an ion *is*, and this says
/// how a pair of them interacts. `PhasePitzer` is its reader upstream, and the tranche's
/// other activity models read their own tables beside it.
const PITZER_PARAMETERS_CSV: &str = include_str!("../../../data/components/PitzerParameters.csv");

/// The compiled salt table, generated from NeqSim's `COMPSALT.csv`.
const COMPSALT_CSV: &str = include_str!("../../../data/components/COMPSALT.csv");

/// Repo-relative path of the component table, which is how Python addresses the same
/// file. A constant rather than a string restated at the call site for the reason the
/// whole databank exists: two copies of a path can disagree.
pub const COMPONENTS_PATH: &str = "data/components/components.csv";

/// Repo-relative path of the interaction table.
pub const KIJ_PATH: &str = "data/components/kij.csv";

/// Repo-relative paths of the three UNIFAC tables.
pub const UNIFAC_COMP_PATH: &str = "data/components/UNIFACcomp.csv";
pub const UNIFAC_GROUP_PATH: &str = "data/components/UNIFACGroupParam.csv";
pub const UNIFAC_INTER_PATH: &str = "data/components/UNIFACInterParam.csv";
pub const UNIFAC_INTER_B_PATH: &str = "data/components/UNIFACInterParamB.csv";
pub const UNIFAC_INTER_C_PATH: &str = "data/components/UNIFACInterParamC.csv";
pub const UNIFAC_COMP_UMRPRU_PATH: &str = "data/components/UNIFACcompUMRPRU.csv";
pub const UNIFAC_A_UMR_PATH: &str = "data/components/UNIFACInterParamA_UMR.csv";
pub const UNIFAC_A_UMRMC_PATH: &str = "data/components/UNIFACInterParamA_UMRMC.csv";
pub const UNIFAC_B_UMR_PATH: &str = "data/components/UNIFACInterParamB_UMR.csv";
pub const UNIFAC_B_UMRMC_PATH: &str = "data/components/UNIFACInterParamB_UMRMC.csv";
pub const UNIFAC_C_UMR_PATH: &str = "data/components/UNIFACInterParamC_UMR.csv";
pub const UNIFAC_C_UMRMC_PATH: &str = "data/components/UNIFACInterParamC_UMRMC.csv";
pub const MBWR32_PATH: &str = "data/components/mbwr32.csv";

/// The exact bytes this build embedded for the component table.
///
/// Exposed so the Python side can compare bytes rather than parsed values - parsed
/// values agree across two *different* files, which is what a stale copy bundled into
/// a wheel looks like. `python/tests/test_data_agreement.py` is the caller.
#[must_use]
pub fn embedded_components() -> &'static str {
    COMPONENTS_CSV
}

/// The exact bytes this build embedded for the interaction table.
#[must_use]
pub fn embedded_kij() -> &'static str {
    KIJ_CSV
}

/// The exact bytes this build embedded for each UNIFAC table, so the Python side can
/// compare bytes rather than parsed values.
#[must_use]
pub fn embedded_unifac_comp() -> &'static str {
    UNIFAC_COMP_CSV
}

#[must_use]
pub fn embedded_unifac_group() -> &'static str {
    UNIFAC_GROUP_CSV
}

#[must_use]
pub fn embedded_unifac_inter() -> &'static str {
    UNIFAC_INTER_CSV
}

#[must_use]
pub fn embedded_unifac_inter_b() -> &'static str {
    UNIFAC_INTER_B_CSV
}

#[must_use]
pub fn embedded_unifac_inter_c() -> &'static str {
    UNIFAC_INTER_C_CSV
}

/// The exact bytes this build embedded for each UNIFAC-UMR-PRU table.
#[must_use]
pub fn embedded_unifac_comp_umrpru() -> &'static str {
    UNIFAC_COMP_UMRPRU_CSV
}

#[must_use]
pub fn embedded_unifac_umrpru_matrix(set: UmrpruSet, term: UmrpruTerm) -> &'static str {
    match (set, term) {
        (UmrpruSet::Umr, UmrpruTerm::A) => UNIFAC_A_UMR_CSV,
        (UmrpruSet::Umrmc, UmrpruTerm::A) => UNIFAC_A_UMRMC_CSV,
        (UmrpruSet::Umr, UmrpruTerm::B) => UNIFAC_B_UMR_CSV,
        (UmrpruSet::Umrmc, UmrpruTerm::B) => UNIFAC_B_UMRMC_CSV,
        (UmrpruSet::Umr, UmrpruTerm::C) => UNIFAC_C_UMR_CSV,
        (UmrpruSet::Umrmc, UmrpruTerm::C) => UNIFAC_C_UMRMC_CSV,
    }
}

/// The exact bytes this build embedded for the MBWR-32 table.
#[must_use]
pub fn embedded_mbwr32() -> &'static str {
    MBWR32_CSV
}

/// The `REFERENCESTATETYPE` value NeqSim treats as the symmetric (Raoult) reference.
///
/// `ComponentGE.fugcoef` compares the component's column against exactly this string and
/// falls to the Henry's-law branch for anything else - including the literal `0.0` two
/// of the table's rows carry.
pub const SOLVENT: &str = "solvent";

/// The `COMPTYPE` value whose rows carry no critical constants of their own.
///
/// `COMP.csv` files 62 rows under it, and for 27 of them the four critical columns hold
/// one default set - `Tc = 444.15` as the file states it, which the generator's 273.15
/// makes `717.3 K`, then `Pc = 290.89 bar`, `omega = 0.344` and `Vc = 99 cm3/mol`. The
/// remaining 35 inherit a neutral parent's numbers (`MDEA+` carries MDEA's) or carry
/// values with no stated source. `Entry::class` is where the tag reaches this side, and
/// [`mixture_of`] refuses rather than reading the filler.
pub const ION: &str = "ion";

/// One substance's constants, in the units the compiled table holds them in.
#[derive(Debug, Clone, PartialEq)]
pub struct Entry {
    /// The name it is looked up by, lower case.
    pub name: String,
    /// Critical temperature, in K.
    pub tc: f64,
    /// Critical pressure, in Pa.
    pub pc: f64,
    /// Acentric factor, dimensionless.
    pub omega: f64,
    /// The `Cp` polynomial's five coefficients, in J/(mol*K**n). `None` for a substance
    /// an overlay added: a card supplies the parameters a cubic needs. `mixture_of`
    /// refuses such a name rather than defaulting to zeros.
    pub cp: Option<[f64; 5]>,
    /// The Antoine coefficients `A`-`E` and the form they belong to, as the table
    /// states them.
    ///
    /// Read for `eos.ge_nrtl_phase`, whose liquid fugacity is `gamma_i P0_i / P`: `P0_i`
    /// is this correlation, and the phase needs it per component. `None` for a substance
    /// an overlay added, which carries the parameters a *cubic* needs and not these.
    pub antoine: Option<([f64; 5], String)>,
    /// Molar mass, in kg/mol.
    pub molar_mass: Option<f64>,
    /// Critical molar volume, in m³/mol.
    pub critical_volume: Option<f64>,
    /// Dipole moment, in debye.
    pub dipole: Option<f64>,
    /// Which reference state the component's activity model is written against, as
    /// NeqSim's `REFERENCESTATETYPE` column states it.
    ///
    /// Read by the activity-coefficient phases, whose fugacity coefficient is *not* one
    /// expression: `ComponentGE.fugcoef` returns `gamma_i P0_i / P` for a component
    /// tagged `solvent` and a Henry's-law coefficient for anything else. The values in
    /// the table are `solvent`, `solute`, and the literal `0.0` - which is neither, and
    /// is what NeqSim's own reader compares against `"solvent"` too, so a component
    /// carrying it takes the same branch NeqSim gives it.
    pub reference_state: String,
    /// The substance class, exactly as NeqSim's `COMPTYPE` column states it.
    ///
    /// Read for one thing, and it is a refusal rather than a term: [`mixture_of`] builds
    /// no cubic over an [`ION`], because the table's critical columns for those rows are
    /// a default rather than a measurement. Never computed with.
    pub class: String,
    /// The association parameters, or `None` for a component the table gives no scheme.
    ///
    /// 181 of the 348 compiled rows carry `0` in the scheme column and are non-associating;
    /// the rest name `1A`, `2A`, `2B` or `4C`.
    pub association: Option<AssociationRecord>,
    /// PC-SAFT's segment number `m`, dimensionless.
    ///
    /// **Zero means the table carries no PC-SAFT set**, which is how 106 of the 348 rows
    /// spell absence - the column is never blank. A model must refuse a zero rather than
    /// compute with it, the rule `ComponentSrkCPA`'s `|aCPA| > 1e-6` guard already follows.
    pub m_saft: f64,
    /// PC-SAFT's segment diameter `sigma`, in metres.
    ///
    /// NeqSim's `COMP.csv` carries it in ångström and divides by `1e10` as it reads it
    /// (`Component.java:548`); this table is compiled SI, so the conversion has happened
    /// once and is not repeated here.
    pub sigma_saft: f64,
    /// PC-SAFT's segment energy over Boltzmann's constant, `epsilon/k`, in K.
    pub epsik_saft: f64,
    /// `UMRCPA_MC1..5`, the five-parameter Mathias-Copeman coefficients.
    ///
    /// **Every component may carry these, associating or not.** NeqSim reads them in
    /// `Component.createComponent` (line 404) before any attractive term is chosen, and
    /// `ComponentUMRCPA.setAttractiveTerm` installs term 22 for a component whose set
    /// has them - a hydrocarbon as readily as a glycol. 23 of the 348 rows carry one.
    pub umrcpa_mc: [f64; 5],
    /// SAFT-VR-Mie's repulsive exponent `lambda_r`, dimensionless.
    ///
    /// **This one is an absence marker too, on exactly the rows `m_mie` is.** The table
    /// carries the standard `12` on 336 of its 348 rows and a fitted exponent on the 12
    /// that have a set, so `lambda_r_mie == 12.0` and `m_mie == 0.0` hold on the same 336
    /// rows. Two markers for one fact, which is why `m_mie` is the one read.
    pub lambda_r_mie: f64,
    /// SAFT-VR-Mie's attractive exponent `lambda_a`, dimensionless, `6` on 347 rows.
    ///
    /// Not a marker: the fitted rows carry `6` as well, apart from `co2` at `5.055`.
    pub lambda_a_mie: f64,
    /// SAFT-VR-Mie's segment number `m`, dimensionless.
    ///
    /// **Zero means the table carries no SAFT-VR-Mie set**, which is how 336 of the 348
    /// rows spell absence - 12 carry one, the light alkanes with `co2`, nitrogen and
    /// water. `sigma_mie` and `epsik_mie` are zero on exactly those rows.
    pub m_mie: f64,
    /// SAFT-VR-Mie's segment diameter `sigma`, in metres.
    pub sigma_mie: f64,
    /// SAFT-VR-Mie's segment energy over Boltzmann's constant, `epsilon/k`, in K.
    pub epsik_mie: f64,
    /// The ionic charge, in units of the elementary charge; zero for a neutral.
    ///
    /// **Zero does not mean "not an ion".** Four rows the table classes [`ION`] carry
    /// it: `h+pzcoo-` is a zwitterion and `caco3`, `nacl` and `cacl2` are neutral salts
    /// filed with the ions. Class is the marker; this is the value.
    pub ionic_charge: f64,
    /// The Deshmukh-Mather ion diameter, in ångström.
    ///
    /// **The table's absence marker is `0.0`, not a blank**, which is what 18 of the 62
    /// [`ION`] rows carry. A model that needs one must refuse a zero: a zero diameter is
    /// the infinite-dilution limit of the term it belongs to, so a plausible number would
    /// come out of it.
    pub deshmukh_mather_diameter: f64,
    /// The four coefficients of NeqSim's Henry correlation.
    ///
    /// **All four are zero on 296 of the 348 rows**, which is how the table spells "no
    /// correlation", so a caller must refuse rather than evaluate - see
    /// [`crate::henry`], which is the surface a model takes them through.
    pub henry: crate::henry::HenryRecord,
    /// `DIELECTRICPARAMETER1..5`, the dielectric mixing rule's coefficients.
    ///
    /// Carried exactly as NeqSim stores them. The unit is not established - the manifest
    /// records `neqsim-internal` - so no conversion is applied and none is guessed.
    pub dielectric: [f64; 5],
}

impl Entry {
    /// The cubic's record for this substance.
    ///
    /// Carries the association parameters when the table gives it a scheme, because a
    /// substance's description is not only its critical constants - see
    /// [`AssociationRecord`], which is why an associating equation of state cannot be
    /// built from `Tc` and `Pc` alone.
    pub fn component(&self) -> Result<Component> {
        Ok(
            Component::new(kelvins(self.tc), pascals(self.pc), self.omega)?
                .with_molar_mass(self.molar_mass)
                .with_association(self.association.clone()),
        )
    }
}

/// One substance as an overlay states it: each parameter named, or none.
///
/// Every field is optional, and that is the rule a keycard follows rather than a
/// convenience: a card naming only `omega` keeps the shipped `Tc` and `Pc`. A record
/// that replaced the shipped one whole would make a user correcting one value restate
/// the others, and lose them silently if they did not.
///
/// **A closed list.** It widens when the component model does, and
/// `specs/schema/component.schema.json` declares the cubic's parameters while
/// `specs/schema/association.schema.json` declares the association's. This is the second
/// place after `keycard.COMPONENT_PARAMETERS` that says which parameters a *cubic* reads,
/// and the two are held together by a test.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ComponentOverride {
    /// Critical temperature, in K.
    pub tc: Option<f64>,
    /// Critical pressure, in Pa.
    pub pc: Option<f64>,
    /// Acentric factor, dimensionless.
    pub omega: Option<f64>,
    /// The five ideal-gas heat-capacity coefficients, all or none.
    pub cp: Option<[f64; 5]>,
    /// The association this card states, or `None` to keep the table's.
    pub association: Option<AssociationOverride>,
    /// The ionic charge as a **charge number**, in units of the elementary charge.
    ///
    /// The card states this in `dimensionless`, which is the quantity every electrolyte
    /// model consumes: Pitzer's ionic strength is `1/2 sum m z^2` and NeqSim's own column
    /// holds the same number. **Zero is not "not an ion"** - `caco3` is a neutral salt
    /// filed with them - so this is a value and [`Self::ion`] is the marker.
    pub ionic_charge: Option<f64>,
    /// The Deshmukh-Mather ion diameter, in **metres**.
    ///
    /// The card states metres and the databank carries ångström, which is the one
    /// crossing between them: `ComponentDesmukhMather` multiplies its own column by
    /// `1e-10` at the point of use, so the two are the same quantity one conversion apart.
    pub deshmukh_mather_diameter: Option<f64>,
    /// The five dielectric-constant coefficients `d0`..`d4`, in the units their positions
    /// in `epsilon(T) = d0 + d1/T + d2 T + d3 T^2 + d4 T^3` give them.
    ///
    /// All five or none, for the reason [`Self::cp`] is: three of five is not a
    /// polynomial, and a partial set is not a value any model reads.
    pub dielectric: Option<[f64; 5]>,
    /// Whether the substance is an ion, stated only by a card that **adds** one.
    ///
    /// `None` is "the card says nothing", which leaves the table's class in force. A card
    /// cannot reclassify a substance the databank already carries: the class is what
    /// decides whether a cubic may be built at all, so a card able to clear it could hand
    /// a cubic the filler the refusal exists to keep out of one.
    ///
    /// A card-added ion needs no `Tc`, `Pc` or `omega` - an ion has no meaningful critical
    /// constants, and [`mixture_of`] refuses one - so this is also what exempts it from
    /// [`Self::is_complete`].
    pub ion: Option<bool>,
}

impl ComponentOverride {
    /// Whether this override names every parameter a cubic needs. Only asked of a
    /// substance the table does not have; one it does have is completed from the table.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.tc.is_some() && self.pc.is_some() && self.omega.is_some()
    }
}

/// One substance's association, as an overlay states it.
///
/// The card's values are in **SI** - the unit `specs/schema/association.schema.json`
/// declares for each, which is the unit a fitted CPA parameter set is published in - and
/// [`AssociationRecord`] holds the two fitted cubic constants in the source table's
/// internal scale. [`Self::applied_to`] is the one crossing between the two, so the factor
/// appears once rather than at every read.
///
/// Every parameter but the scheme is optional, for the reason [`ComponentOverride`] is: a
/// card regressing only `epsilon` and `kappa_AB` keeps the table's fitted set, and one
/// regressing the whole set states it. A card that states no fitted attraction leaves the
/// cubic's own in place, which is NeqSim's `|aCPA| > 1e-6` guard reached from the other
/// side rather than a different rule.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AssociationOverride {
    /// The site scheme, which a card must state: it is the one parameter that says
    /// whether the substance associates at all.
    pub scheme: SiteScheme,
    /// The association energy `eps`, in J/mol.
    pub energy: Option<f64>,
    /// `kappa_AB` for the SRK family.
    pub volume_srk: Option<f64>,
    /// The fitted attraction for SRK-CPA, in `Pa m**6/mol**2`.
    pub a_srk: Option<f64>,
    /// The fitted covolume for SRK-CPA, in `m**3/mol`.
    pub b_srk: Option<f64>,
    /// The SRK alpha correlation's `m`.
    pub m_srk: Option<f64>,
    /// `kappa_AB` for the PR family.
    pub volume_pr: Option<f64>,
    /// The fitted attraction for PR-CPA, in `Pa m**6/mol**2`.
    pub a_pr: Option<f64>,
    /// The fitted covolume for PR-CPA, in `m**3/mol`.
    pub b_pr: Option<f64>,
    /// The PR alpha correlation's `m`.
    pub m_pr: Option<f64>,
}

impl AssociationOverride {
    /// The record this override states, taking every parameter it does not name from
    /// `base` - which is the table's, or `None` for a substance that has none.
    ///
    /// **The site count is the table's where the scheme is, and the scheme's where the
    /// scheme is not.** A record carries a count *beside* its scheme because the two
    /// disagree upstream - the table's `1A` rows carry a count of zero - and a card that
    /// restated the scheme would otherwise resolve that disagreement silently, which is
    /// what the field exists to keep visible. A card naming a *different* scheme is
    /// stating a different molecule and has no count to inherit, so the new scheme's is
    /// the only one there is.
    #[must_use]
    pub fn applied_to(&self, base: Option<&AssociationRecord>) -> AssociationRecord {
        let kept = |pick: fn(&AssociationRecord) -> f64| base.map_or(0.0, pick);
        let internal = |si: Option<f64>, pick: fn(&AssociationRecord) -> f64| match si {
            Some(value) => value / crate::association::NEQSIM_INTERNAL_TO_SI,
            None => kept(pick),
        };
        let sites = match base {
            Some(record) if record.scheme == self.scheme => record.sites,
            _ => self.scheme.site_count() as u32,
        };
        AssociationRecord {
            scheme: self.scheme,
            sites,
            energy: self.energy.unwrap_or_else(|| kept(|r| r.energy)),
            volume_srk: self.volume_srk.unwrap_or_else(|| kept(|r| r.volume_srk)),
            a_srk: internal(self.a_srk, |r| r.a_srk),
            b_srk: internal(self.b_srk, |r| r.b_srk),
            m_srk: self.m_srk.unwrap_or_else(|| kept(|r| r.m_srk)),
            volume_pr: self.volume_pr.unwrap_or_else(|| kept(|r| r.volume_pr)),
            a_pr: internal(self.a_pr, |r| r.a_pr),
            b_pr: internal(self.b_pr, |r| r.b_pr),
            m_pr: self.m_pr.unwrap_or_else(|| kept(|r| r.m_pr)),
            racket_z: kept(|r| r.racket_z),
            volume_correction: kept(|r| r.volume_correction),
            umr_cpa: base.and_then(|r| r.umr_cpa.clone()),
        }
    }
}

/// A keycard's data, as a value a caller passes. [`crate::card::Card::overlay`] is where
/// a card file becomes one; a caller with the values in hand builds one with
/// [`Overlay::new`], so a carded lookup needs no file. Nothing holds one.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Overlay {
    /// Overrides by lower-cased name, and additions the table does not have.
    components: HashMap<String, ComponentOverride>,
    /// Overrides by lower-cased pair, stored both ways round.
    kij: HashMap<(String, String), f64>,
    /// Overrides of the *associating* interaction columns, keyed by pair and family.
    ///
    /// A separate map because it is a separate column: NeqSim's CPA rule reads
    /// `cpakij_SRK`/`cpakij_PR` and a classical mixture reads `KIJSRK` or `KIJPR`
    /// depending on its cubic, and on water/methanol they differ by a factor of two. A
    /// card stating one is not stating the other, so neither may stand in for the other -
    /// and the family is part of the key because the two `cpakij` columns are two fits
    /// rather than one converted.
    cpa_kij: HashMap<(String, String, AssociationCubic), f64>,
}

impl Overlay {
    /// An overlay that changes nothing.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Override one substance's parameters, or add a substance by name.
    pub fn set_component(&mut self, name: &str, parameters: ComponentOverride) -> &mut Self {
        self.components
            .insert(name.trim().to_lowercase(), parameters);
        self
    }

    /// Override one substance's association, or state one for a substance the table gives
    /// none by name.
    ///
    /// A read-modify-write, so a card stating a substance's critical constants and its
    /// association in two sections resolves to one override rather than the second
    /// replacing the first.
    pub fn set_association(&mut self, name: &str, association: AssociationOverride) -> &mut Self {
        let key = name.trim().to_lowercase();
        let mut parameters = self.components.get(&key).copied().unwrap_or_default();
        parameters.association = Some(association);
        self.components.insert(key, parameters);
        self
    }

    /// This overlay's statement about a substance's association, if it makes one.
    #[must_use]
    pub fn association(&self, name: &str) -> Option<&AssociationOverride> {
        self.components
            .get(&name.trim().to_lowercase())
            .and_then(|parameters| parameters.association.as_ref())
    }

    /// Override one pair's interaction parameter.
    ///
    /// **One number covers both cubic columns.** A card naming `kij` has said what the
    /// pair is; asking it to say so again per cubic would make the two halves separately
    /// overridable by accident, and a card that overrode only one would leave the other
    /// reading the table.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if both names are the same substance. `Mixture::new`
    ///   refuses a non-zero diagonal too, but only for a pair whose *both* names are in
    ///   the mixture being built - so a self-pair on a name nothing builds would be stored
    ///   and read by nothing, which is the failure this crate refuses everywhere.
    pub fn set_kij(&mut self, first: &str, second: &str, value: f64) -> Result<&mut Self> {
        let a = first.trim().to_lowercase();
        let b = second.trim().to_lowercase();
        if a == b {
            return Err(AzothError::invalid_input(
                "kij",
                format!("`{a}` does not interact with itself"),
            ));
        }
        // Stored both ways, exactly as the table is, so a caller need not know which
        // name came first.
        self.kij.insert((a.clone(), b.clone()), value);
        self.kij.insert((b, a), value);
        Ok(self)
    }

    /// Override one pair's *associating* interaction parameter, at one cubic family.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if both names are the same substance, for the reason
    ///   [`Self::set_kij`] refuses one.
    pub fn set_cpa_kij(
        &mut self,
        first: &str,
        second: &str,
        cubic: AssociationCubic,
        value: f64,
    ) -> Result<&mut Self> {
        let a = first.trim().to_lowercase();
        let b = second.trim().to_lowercase();
        if a == b {
            return Err(AzothError::invalid_input(
                "cpa_kij",
                format!("`{a}` does not interact with itself"),
            ));
        }
        self.cpa_kij.insert((a.clone(), b.clone(), cubic), value);
        self.cpa_kij.insert((b, a, cubic), value);
        Ok(self)
    }

    /// This overlay's statement about a pair's associating interaction at one family, if
    /// it makes one.
    #[must_use]
    pub fn cpa_kij_value(&self, first: &str, second: &str, cubic: AssociationCubic) -> Option<f64> {
        self.cpa_kij
            .get(&(
                first.trim().to_lowercase(),
                second.trim().to_lowercase(),
                cubic,
            ))
            .copied()
    }

    /// This overlay's statement about a substance, or `None` if it makes none.
    #[must_use]
    pub fn component(&self, name: &str) -> Option<&ComponentOverride> {
        self.components.get(&name.trim().to_lowercase())
    }

    /// This overlay's interaction parameter for a pair, or `None` if it states none.
    ///
    /// `None` and `Some(0.0)` are different answers, and the difference is a caller's
    /// deliberate reset to ideal mixing. Anything that filters a zero here silently
    /// undoes it.
    #[must_use]
    pub fn kij(&self, first: &str, second: &str) -> Option<f64> {
        self.kij
            .get(&(first.trim().to_lowercase(), second.trim().to_lowercase()))
            .copied()
    }

    /// Every pair this overlay states, each once, with the lower name first.
    ///
    /// The de-duplication is the opposite of the *storage*, which keeps both orderings
    /// so a caller need not know which name came first. This is the display direction:
    /// one row per pair.
    #[must_use]
    pub fn kij_pairs(&self) -> Vec<(String, String)> {
        let mut out: Vec<(String, String)> = self
            .kij
            .keys()
            .filter(|(first, second)| first < second)
            .cloned()
            .collect();
        out.sort();
        out
    }

    /// The names this overlay adds to the table, in no particular order.
    #[must_use]
    pub fn component_names(&self) -> Vec<&str> {
        self.components.keys().map(String::as_str).collect()
    }

    /// Whether this overlay changes nothing.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.components.is_empty() && self.kij.is_empty() && self.cpa_kij.is_empty()
    }
}

/// One interaction-pair record from `INTER.csv`, keyed by the ordered pair of names
/// exactly as they appear in the file.
///
/// Every directional field is stored **for this key**: the reversed key carries the
/// reversed value, done once where the table is parsed rather than at each read.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Interaction {
    /// NeqSim's `KIJPR`, the binary interaction parameter its **Peng-Robinson** phase
    /// reads, symmetric.
    ///
    /// **The cubic's column is chosen by cubic, and this one is only PR's.** NeqSim
    /// selects on the phase class - `phase.getClass().getName().equals(
    /// "neqsim.thermo.phase.PhasePrEos")` - and every other phase reads `KIJSRK`. The two
    /// are different fits rather than one converted into the other: water/methane is
    /// `0.651` here against `0.45`, benzene/methane is `0.0` against `0.0209`, and 194
    /// of the 957 in-scope pairs differ at all.
    kij_pr: f64,
    /// NeqSim's `KIJSRK`, the column its Soave-Redlich-Kwong phase and every other cubic
    /// reads. See [`Self::kij_pr`] for why the two are not interchangeable.
    kij_srk: f64,
    /// NeqSim's `cpakij_SRK`, the interaction parameter its **CPA** mixing rule reads for
    /// the Soave family.
    ///
    /// A *different column* from [`Self::kij`], not a refinement of it: water/methanol is
    /// `-0.153` here against `KIJPR`'s `-0.0789`, and reading the wrong one changes the
    /// mixture's attraction by 3.4% and its root by 0.077%.
    cpa_kij_srk: f64,
    /// NeqSim's `cpakij_PR`, the CPA rule's parameter for the Peng-Robinson family.
    cpa_kij_pr: f64,
    /// NeqSim's `KIJPCSAFT`, the interaction parameter its PC-SAFT phase reads.
    ///
    /// A third *different* fit, not a refinement of the other two: methane/n-butane is
    /// `0.022` here against the `0.01289789` its SRK and PR columns share, and 42 of the
    /// 957 in-scope pairs carry one at all.
    ///
    /// **NeqSim only reads it when it is driven the one way it can be driven.** Its own
    /// default mixing rule leaves the matrix null and a *mixture* throws a
    /// NullPointerException rather than using a zero - a pure fluid never asks for a pair
    /// and runs either way - so `setMixingRule("classic")` is what a caller must do, and
    /// the branch that serves it reads this column. There is no usable configuration in
    /// which NeqSim's PC-SAFT ignores it.
    pcsaft_kij: f64,
    /// The NRTL non-randomness parameter, symmetric.
    alpha: f64,
    /// The NRTL energy parameter for the ordered pair `(first, second)`: the `g_ij` in
    /// `tau_ij = g_ij / T`, in Kelvin. Directional - `g_ij` differs from `g_ji`.
    gij: f64,
    /// Whether this pair carries Huron-Vidal's *fitted* parameters, from `HVTYPE`.
    ///
    /// NeqSim reads the same column into `classicOrHV` and tests
    /// `mixRule[j][i].trim().equals("HV")`. `HV` is the fitted NRTL pair and `Classic`
    /// is the cubic's own excess energy, which is not a fallback but the other branch of
    /// the model.
    hv: bool,
    /// The Huron-Vidal non-randomness for this ordered pair, symmetric.
    hv_alpha: f64,
    /// The Huron-Vidal fitted energy `Dij` for this ordered pair, in Kelvin. Directional.
    hv_dij: f64,
    /// The Huron-Vidal temperature coefficient `DijT` for this ordered pair. Directional.
    hv_dij_t: f64,
    /// Whether this pair carries Wong-Sandler's fitted parameters, from `WSTYPE`.
    ws: bool,
    /// The Wong-Sandler temperature coefficient `DijT` for this ordered pair. Directional.
    ///
    /// A *different* column from `hv_dij_t`: NeqSim loads the HV one into `HVDijT` and
    /// the WS one into `NRTLDijT`, and the two rules are handed the corresponding array.
    ws_dij_t: f64,
    /// The Wong-Sandler interaction parameter, from `KIJWSunifac`.
    ///
    /// `KIJWSunifac` and not `KIJWS`: NeqSim assigns `WSintparam` from `kijWS` and
    /// overwrites it from `KIJWSunifac` on the next line, so the first read has never
    /// had an effect.
    kij_ws: f64,
}

/// The two parsed tables: substances by name, and interaction parameters by pair.
type Tables = (
    HashMap<String, Entry>,
    HashMap<(String, String), Interaction>,
);

/// The parsed tables, read once.
fn tables() -> &'static Tables {
    static TABLES: OnceLock<Tables> = OnceLock::new();
    TABLES.get_or_init(|| {
        (
            parse_components().expect("the embedded component table should parse"),
            parse_kij().expect("the embedded interaction table should parse"),
        )
    })
}

/// The column index of `name` in a header row, or a named failure.
fn column(header: &csv::StringRecord, name: &str) -> Result<usize> {
    header
        .iter()
        .position(|field| field == name)
        .ok_or_else(|| AzothError::InvalidInput {
            field: "databank".to_string(),
            reason: format!(
                "the compiled table has no `{name}` column; its header is {}",
                header.iter().collect::<Vec<_>>().join(", ")
            ),
        })
}

/// A field of one record, parsed as a number.
fn number(record: &csv::StringRecord, index: usize, column: &str, row: usize) -> Result<f64> {
    let raw = record.get(index).unwrap_or("").trim();
    raw.parse::<f64>().map_err(|_| AzothError::InvalidInput {
        field: "databank".to_string(),
        reason: format!("row {row}: `{column}` is {raw:?}, which is not a number"),
    })
}

/// A `csv` failure as this crate's error.
///
/// Unreachable in practice - the table is embedded, so a malformed one is a build
/// defect rather than a caller condition - and it is a `Result` anyway, because this
/// crate does not panic on data and because the table is a build artefact whose parse
/// is fallible in principle. It is not here for a keycard: a card is read by
/// [`crate::card`], which reports what is wrong with the document, and an overlay
/// arrives here already resolved.
fn csv_failure(error: csv::Error) -> AzothError {
    AzothError::InvalidInput {
        field: "databank".to_string(),
        reason: format!("the embedded table could not be read: {error}"),
    }
}

/// One parameter's value, where an empty cell means zero.
///
/// **Both upstream tables use a blank and a `0` interchangeably for "no parameter", so a
/// blank infers nothing the row does not already state.** Four cells in the compiled
/// tables are blank: `associationboundingvolume_pr` on `h2so4`, `hno3` and `c2h4-`, whose
/// every sibling association cell is already `0.0`; and `nrtlgij` on the `k+`/`h2s` pair,
/// which carries `KIJPR = 1` with no NRTL energy beside it.
///
/// A *malformed* value still refuses - the distinction is between an absent parameter and
/// a broken one, which [`number`] keeps.
fn optional_number(
    record: &csv::StringRecord,
    index: usize,
    name: &str,
    row: usize,
) -> Result<f64> {
    if record.get(index).unwrap_or("").trim().is_empty() {
        return Ok(0.0);
    }
    number(record, index, name, row)
}

/// One row's association parameters, or `None` where it carries no scheme.
///
/// The scheme column is text and `associationsites` a count. A row naming `"0"` - the
/// table's marker for a component with no scheme - is non-associating, and so is a name
/// this library does not carry; the two are distinguished because only the first is the
/// table's intent, but both mean the same thing to a model.
fn parse_association(
    record: &csv::StringRecord,
    index: &HashMap<&str, usize>,
    row: usize,
) -> Result<Option<AssociationRecord>> {
    let raw = record.get(index["associationscheme"]).unwrap_or("").trim();
    let Some(scheme) = SiteScheme::from_databank_name(raw) else {
        return Ok(None);
    };
    let fitted = |name: &str| -> Result<f64> { optional_number(record, index[name], name, row) };
    let umr_cpa = parse_umr_cpa(record, index, row)?;
    Ok(Some(AssociationRecord {
        scheme,
        sites: fitted("associationsites")? as u32,
        energy: fitted("associationenergy")?,
        volume_srk: fitted("associationboundingvolume_srk")?,
        a_srk: fitted("acpa_srk")?,
        b_srk: fitted("bcpa_srk")?,
        m_srk: fitted("mcpa_srk")?,
        volume_pr: fitted("associationboundingvolume_pr")?,
        a_pr: fitted("acpa_pr")?,
        b_pr: fitted("bcpa_pr")?,
        m_pr: fitted("mcpa_pr")?,
        racket_z: fitted("racketzcpa")?,
        volume_correction: fitted("volcorrcpa_t")?,
        umr_cpa,
    }))
}

/// One row's `UMRCPA_*` set, or `None` where the row does not carry one.
///
/// NeqSim's guard is two-part (`Component.java:534`): `UMRCPA_associating` is one **and**
/// `|UMRCPA_a0| > 1e-20`. Seven of the 348 rows pass it - water, methanol, ethanol and
/// the four glycols - and a row that fails it keeps the PR family's fitted values, which
/// is what `ComponentUMRCPA` reads first for being a `ComponentPR`.
fn parse_umr_cpa(
    record: &csv::StringRecord,
    index: &HashMap<&str, usize>,
    row: usize,
) -> Result<Option<UmrCpaRecord>> {
    let number = |name: &str| -> Result<f64> { optional_number(record, index[name], name, row) };
    let associating = number("umrcpa_associating")?;
    let a0 = number("umrcpa_a0")?;
    if associating != 1.0 || a0.abs() <= 1.0e-20 {
        return Ok(None);
    }
    Ok(Some(UmrCpaRecord {
        a0,
        b: number("umrcpa_b")?,
        energy: number("umrcpa_assocenergy")?,
        volume: number("umrcpa_assocvolume")?,
        racket_z: number("umrcpa_racketz")?,
    }))
}

fn parse_components() -> Result<HashMap<String, Entry>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(COMPONENTS_CSV.as_bytes());
    let header = reader.headers().map_err(csv_failure)?.clone();
    let mut index = HashMap::new();
    for name in [
        "name",
        "tc_k",
        "pc_pa",
        "acentric_factor",
        "cpa",
        "cpb",
        "cpc",
        "cpd",
        "cpe",
        "molar_mass_kg_per_mol",
        "critical_volume_m3_per_mol",
        "dipole_moment_debye",
        "antoine_type",
        "antoinea",
        "antoineb",
        "antoinec",
        "antoined",
        "antoinee",
        "referencestatetype",
        "comptype",
        "ioniccharge",
        "deshmationicdiameter",
        "henrycoef1",
        "henrycoef2",
        "henrycoef3",
        "henrycoef4",
        "dielectricparameter1",
        "dielectricparameter2",
        "dielectricparameter3",
        "dielectricparameter4",
        "dielectricparameter5",
        "associationscheme",
        "associationsites",
        "associationenergy",
        "associationboundingvolume_srk",
        "associationboundingvolume_pr",
        "msaft",
        "sigma_saft_m",
        "epsiksaft",
        "lambdarsaftvrmie",
        "lambdaasaftvrmie",
        "msaftvrmie",
        "sigma_saft_vr_mie_m",
        "epsiksaftvrmie",
        "acpa_srk",
        "bcpa_srk",
        "mcpa_srk",
        "acpa_pr",
        "bcpa_pr",
        "mcpa_pr",
        "racketzcpa",
        "volcorrcpa_t",
        "umrcpa_mc1",
        "umrcpa_mc2",
        "umrcpa_mc3",
        "umrcpa_mc4",
        "umrcpa_mc5",
        "umrcpa_a0",
        "umrcpa_b",
        "umrcpa_assocenergy",
        "umrcpa_assocvolume",
        "umrcpa_associating",
        "umrcpa_racketz",
    ] {
        index.insert(name, column(&header, name)?);
    }

    let mut out = HashMap::new();
    for (offset, record) in reader.records().enumerate() {
        let record = record.map_err(csv_failure)?;
        let name = record
            .get(index["name"])
            .unwrap_or("")
            .trim()
            .to_lowercase();
        if name.is_empty() {
            continue;
        }
        let row = offset + 2; // one for the header, one to count from one
        let mut cp = [0.0; 5];
        for (slot, key) in cp.iter_mut().zip(["cpa", "cpb", "cpc", "cpd", "cpe"]) {
            *slot = number(&record, index[key], key, row)?;
        }
        out.insert(
            name.clone(),
            Entry {
                name,
                tc: number(&record, index["tc_k"], "tc_k", row)?,
                pc: number(&record, index["pc_pa"], "pc_pa", row)?,
                omega: number(&record, index["acentric_factor"], "acentric_factor", row)?,
                // `Some` for everything the table carries: it holds the polynomial for
                // every row it has, and `mixture_of` refuses a name without one.
                cp: Some(cp),
                molar_mass: Some(number(
                    &record,
                    index["molar_mass_kg_per_mol"],
                    "molar_mass_kg_per_mol",
                    row,
                )?),
                critical_volume: Some(number(
                    &record,
                    index["critical_volume_m3_per_mol"],
                    "critical_volume_m3_per_mol",
                    row,
                )?),
                dipole: Some(number(
                    &record,
                    index["dipole_moment_debye"],
                    "dipole_moment_debye",
                    row,
                )?),
                antoine: Some((
                    [
                        number(&record, index["antoinea"], "antoinea", row)?,
                        number(&record, index["antoineb"], "antoineb", row)?,
                        number(&record, index["antoinec"], "antoinec", row)?,
                        number(&record, index["antoined"], "antoined", row)?,
                        number(&record, index["antoinee"], "antoinee", row)?,
                    ],
                    record
                        .get(index["antoine_type"])
                        .unwrap_or("")
                        .trim()
                        .to_string(),
                )),
                reference_state: record
                    .get(index["referencestatetype"])
                    .unwrap_or("")
                    .trim()
                    .to_string(),
                class: record
                    .get(index["comptype"])
                    .unwrap_or("")
                    .trim()
                    .to_lowercase(),
                association: parse_association(&record, &index, row)?,
                m_saft: number(&record, index["msaft"], "msaft", row)?,
                sigma_saft: number(&record, index["sigma_saft_m"], "sigma_saft_m", row)?,
                epsik_saft: number(&record, index["epsiksaft"], "epsiksaft", row)?,
                umrcpa_mc: {
                    let mut mc = [0.0; 5];
                    for (k, slot) in mc.iter_mut().enumerate() {
                        let column = format!("umrcpa_mc{}", k + 1);
                        *slot = number(&record, index[column.as_str()], &column, row)?;
                    }
                    mc
                },
                lambda_r_mie: number(&record, index["lambdarsaftvrmie"], "lambdarsaftvrmie", row)?,
                lambda_a_mie: number(&record, index["lambdaasaftvrmie"], "lambdaasaftvrmie", row)?,
                m_mie: number(&record, index["msaftvrmie"], "msaftvrmie", row)?,
                sigma_mie: number(
                    &record,
                    index["sigma_saft_vr_mie_m"],
                    "sigma_saft_vr_mie_m",
                    row,
                )?,
                epsik_mie: number(&record, index["epsiksaftvrmie"], "epsiksaftvrmie", row)?,
                ionic_charge: number(&record, index["ioniccharge"], "ioniccharge", row)?,
                deshmukh_mather_diameter: number(
                    &record,
                    index["deshmationicdiameter"],
                    "deshmationicdiameter",
                    row,
                )?,
                henry: crate::henry::HenryRecord {
                    h0: number(&record, index["henrycoef1"], "henrycoef1", row)?,
                    h1: number(&record, index["henrycoef2"], "henrycoef2", row)?,
                    h2: number(&record, index["henrycoef3"], "henrycoef3", row)?,
                    h3: number(&record, index["henrycoef4"], "henrycoef4", row)?,
                },
                dielectric: {
                    let mut dielectric = [0.0; 5];
                    for (k, slot) in dielectric.iter_mut().enumerate() {
                        let column = format!("dielectricparameter{}", k + 1);
                        *slot = number(&record, index[column.as_str()], &column, row)?;
                    }
                    dielectric
                },
            },
        );
    }
    Ok(out)
}

fn parse_kij() -> Result<HashMap<(String, String), Interaction>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(KIJ_CSV.as_bytes());
    let header = reader.headers().map_err(csv_failure)?.clone();
    let mut index = HashMap::new();
    for name in [
        "component_a",
        "component_b",
        "kij_pr",
        "kijsrk",
        "nrtlalpha",
        "nrtlgij",
        "nrtlgji",
        "hvtype",
        "hvalpha",
        "hvgij",
        "hvgji",
        "hvgijt",
        "hvgjit",
        "wstype",
        "wsgijt",
        "wsgjit",
        "kijwsunifac",
        "cpakij_srk",
        "cpakij_pr",
        "kijpcsaft",
    ] {
        index.insert(name, column(&header, name)?);
    }

    let mut out = HashMap::new();
    for (offset, record) in reader.records().enumerate() {
        let record = record.map_err(csv_failure)?;
        let a = record
            .get(index["component_a"])
            .unwrap_or("")
            .trim()
            .to_lowercase();
        let b = record
            .get(index["component_b"])
            .unwrap_or("")
            .trim()
            .to_lowercase();
        if a.is_empty() || b.is_empty() {
            continue;
        }
        let raw = record.get(index["kij_pr"]).unwrap_or("").trim();
        if raw.is_empty() {
            continue;
        }
        let value = raw.parse::<f64>().map_err(|_| AzothError::InvalidInput {
            field: "databank".to_string(),
            reason: format!("row {}: `kij_pr` is {raw:?}", offset + 2),
        })?;
        let row = offset + 2;
        // `optional_number` rather than `number`: a pair may carry a `KIJPR` and no NRTL
        // energy beside it, which is the one blank cell left in this table. See its own
        // comment for the four.
        let alpha = optional_number(&record, index["nrtlalpha"], "nrtlalpha", row)?;
        let gij = optional_number(&record, index["nrtlgij"], "nrtlgij", row)?;
        let gji = optional_number(&record, index["nrtlgji"], "nrtlgji", row)?;

        let hv = selector(&record, index["hvtype"], "hvtype", row)? == "HV";
        let ws = selector(&record, index["wstype"], "wstype", row)? == "WS";
        let hv_alpha = optional_number(&record, index["hvalpha"], "hvalpha", row)?;
        let hv_dij = optional_number(&record, index["hvgij"], "hvgij", row)?;
        let hv_dji = optional_number(&record, index["hvgji"], "hvgji", row)?;
        let hv_dij_t = optional_number(&record, index["hvgijt"], "hvgijt", row)?;
        let hv_dji_t = optional_number(&record, index["hvgjit"], "hvgjit", row)?;
        let ws_dij_t = optional_number(&record, index["wsgijt"], "wsgijt", row)?;
        let ws_dji_t = optional_number(&record, index["wsgjit"], "wsgjit", row)?;
        let kij_ws = optional_number(&record, index["kijwsunifac"], "kijwsunifac", row)?;
        let kij_srk = optional_number(&record, index["kijsrk"], "kijsrk", row)?;
        let cpa_kij_srk = optional_number(&record, index["cpakij_srk"], "cpakij_srk", row)?;
        let cpa_kij_pr = optional_number(&record, index["cpakij_pr"], "cpakij_pr", row)?;
        let pcsaft_kij = optional_number(&record, index["kijpcsaft"], "kijpcsaft", row)?;

        // `kij`, `alpha`, `hv_alpha` and the two selectors are symmetric, stored both
        // ways round so a caller need not know which name came first. `gij`, `hv_dij`,
        // `hv_dij_t` and `ws_dij_t` are directional, so the reversed key carries the
        // reversed value.
        out.insert(
            (a.clone(), b.clone()),
            Interaction {
                kij_pr: value,
                kij_srk,
                cpa_kij_srk,
                cpa_kij_pr,
                pcsaft_kij,
                alpha,
                gij,
                hv,
                hv_alpha,
                hv_dij,
                hv_dij_t,
                ws,
                ws_dij_t,
                kij_ws,
            },
        );
        out.insert(
            (b, a),
            Interaction {
                kij_pr: value,
                kij_srk,
                cpa_kij_srk,
                cpa_kij_pr,
                pcsaft_kij,
                alpha,
                gij: gji,
                hv,
                hv_alpha,
                hv_dij: hv_dji,
                hv_dij_t: hv_dji_t,
                ws,
                ws_dij_t: ws_dji_t,
                kij_ws,
            },
        );
    }
    Ok(out)
}

/// One text field of a record, trimmed.
///
/// NeqSim trims and nothing else - `mixRule[j][i].trim().equals("HV")` - so the
/// comparison this feeds is case-sensitive against the two spellings the file uses.
fn selector(record: &csv::StringRecord, index: usize, column: &str, row: usize) -> Result<String> {
    let raw = record.get(index).unwrap_or("").trim();
    if raw.is_empty() {
        return Err(AzothError::InvalidInput {
            field: "databank".to_string(),
            reason: format!("row {row}: `{column}` is empty; every row names its variant"),
        });
    }
    Ok(raw.to_string())
}

/// One substance's constants, or a failure naming it.
///
/// `overlay` is the card this call reads; `None` means the data this crate ships. An
/// override is applied here, so every caller resolves a name the same way. Owned rather
/// than `&'static`, because a substance an overlay adds is in no table.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if neither the table nor the overlay carries
///   the substance, or if the overlay adds one without every parameter a cubic reads.
pub fn entry(name: &str, overlay: Option<&Overlay>) -> Result<Entry> {
    let key = name.trim().to_lowercase();
    let base = tables().0.get(&key).cloned();
    let over = overlay.and_then(|o| o.component(&key));

    match (base, over) {
        (None, None) => Err(AzothError::property_unavailable(
            key,
            "critical constants".to_string(),
            "neither the databank nor the keycard has it. Names come from NeqSim's COMP.csv, \
             carried in data/components/components.csv; a keycard adds one by name"
                .to_string(),
        )),
        (Some(base), None) => Ok(base),
        (None, Some(over)) => {
            // **An ion is exempt from the cubic's parameters**, and that is the whole point
            // of the flag: an ion has no meaningful `Tc`, `Pc` or `omega`, so requiring them
            // would make a user invent the same filler NeqSim's own table carries - which
            // `mixture_of` refuses for exactly that reason.
            let ion = over.ion == Some(true);
            if !ion && !over.is_complete() {
                let mut missing = Vec::new();
                if over.tc.is_none() {
                    missing.push("Tc");
                }
                if over.pc.is_none() {
                    missing.push("Pc");
                }
                if over.omega.is_none() {
                    missing.push("omega");
                }
                return Err(AzothError::property_unavailable(
                    key,
                    "critical constants".to_string(),
                    format!(
                        "the keycard adds it but is missing {missing:?}. A substance the \
                         databank does not have needs every parameter a cubic reads, because \
                         completing it from a similar one would be inventing data"
                    ),
                ));
            }
            Ok(Entry {
                name: key,
                // Checked complete above for anything that is not an ion; for an ion the
                // defaults are deliberate - the zero states that no critical constant was
                // given, and `mixture_of` refuses one rather than reading the zero.
                tc: over.tc.unwrap_or_default(),
                pc: over.pc.unwrap_or_default(),
                omega: over.omega.unwrap_or_default(),
                // The card may also supply the polynomial, in which case the substance
                // has an enthalpy path; without it, it is a cubic only.
                cp: over.cp,
                // A card states the parameters a cubic reads; it carries no molar mass,
                // critical volume, dipole or Antoine coefficients, so a card-added
                // substance has none.
                molar_mass: None,
                critical_volume: None,
                dipole: None,
                antoine: None,
                association: over.association.as_ref().map(|a| a.applied_to(None)),
                // **A card carries no PC-SAFT set**, and zero is how this table spells
                // absence: a card states the parameters a cubic reads, and a model that
                // needed `m`, `sigma` and `epsilon/k` would refuse rather than invent them.
                m_saft: 0.0,
                sigma_saft: 0.0,
                epsik_saft: 0.0,
                umrcpa_mc: [0.0; 5],
                lambda_r_mie: 0.0,
                lambda_a_mie: 0.0,
                m_mie: 0.0,
                sigma_mie: 0.0,
                epsik_mie: 0.0,
                // Named rather than left blank: a card states a substance a cubic can
                // describe, and a cubic has no reference state. The activity-coefficient
                // phases read this, so a blank would have to mean something.
                reference_state: SOLVENT.to_string(),
                // `other` rather than blank, for the same reason, and it is the class that
                // permits a cubic: a card states the parameters a cubic reads, so the
                // substance it adds is one, and `mixture_of`'s [`ION`] refusal lets it
                // through. Everything below is the electrolyte data a card does not
                // state, and the zeros are its absence rather than a value.
                class: if ion {
                    ION.to_string()
                } else {
                    "other".to_string()
                },
                // A card carries no Henry correlation: `ComponentOverride` is a closed list
                // and an overlay states the parameters a cubic or an activity model reads.
                henry: crate::henry::HenryRecord {
                    h0: 0.0,
                    h1: 0.0,
                    h2: 0.0,
                    h3: 0.0,
                },
                ionic_charge: over.ionic_charge.unwrap_or_default(),
                // The card states metres and the table holds ångström; this is the one
                // crossing, so the factor appears once rather than at every read.
                deshmukh_mather_diameter: over
                    .deshmukh_mather_diameter
                    .map_or(0.0, |metres| metres * 1.0e10),
                dielectric: over.dielectric.unwrap_or([0.0; 5]),
            })
        }
        (Some(base), Some(over)) => Ok(Entry {
            // Parameter by parameter: what the overlay names, else what ships. A card
            // overriding one value does not restate, and does not lose, the others.
            tc: over.tc.unwrap_or(base.tc),
            pc: over.pc.unwrap_or(base.pc),
            omega: over.omega.unwrap_or(base.omega),
            cp: over.cp.or(base.cp),
            molar_mass: base.molar_mass,
            critical_volume: base.critical_volume,
            dipole: base.dipole,
            antoine: base.antoine,
            // `ComponentOverride` is a closed list and carries no PC-SAFT set, so the
            // table's survives a card untouched - the whole point of naming one parameter
            // rather than restating a record.
            m_saft: base.m_saft,
            sigma_saft: base.sigma_saft,
            epsik_saft: base.epsik_saft,
            umrcpa_mc: base.umrcpa_mc,
            lambda_r_mie: base.lambda_r_mie,
            lambda_a_mie: base.lambda_a_mie,
            m_mie: base.m_mie,
            sigma_mie: base.sigma_mie,
            epsik_mie: base.epsik_mie,
            reference_state: base.reference_state,
            // The table's, for the reason the PC-SAFT set below is: `ComponentOverride` is
            // a closed list, so a card correcting `Tc` does not turn an ion into a cubic.
            class: base.class,
            henry: base.henry,
            ionic_charge: over.ionic_charge.unwrap_or(base.ionic_charge),
            deshmukh_mather_diameter: over
                .deshmukh_mather_diameter
                .map_or(base.deshmukh_mather_diameter, |metres| metres * 1.0e10),
            dielectric: over.dielectric.unwrap_or(base.dielectric),
            // The card's scheme wins over the table's, and every parameter the card does
            // not name is the table's: `applied_to` is the one place the two are merged.
            association: over
                .association
                .as_ref()
                .map(|a| a.applied_to(base.association.as_ref()))
                .or(base.association),
            name: base.name,
        }),
    }
}

/// The binary interaction parameter for a pair, or zero.
///
/// Zero rather than a failure: an absent pair is the ideal-mixture default, which is
/// what NeqSim's own reader substitutes. An overlay's zero wins over a fitted value,
/// because overriding a pair back to ideal mixing is a caller stating something.
#[must_use]
pub fn kij(first: &str, second: &str, cubic: Cubic, overlay: Option<&Overlay>) -> f64 {
    if let Some(value) = overlay.and_then(|o| o.kij(first, second)) {
        // A card names one number and it wins for both columns: a caller writing `kij`
        // has said what the pair is, and asking them to say it twice would make the two
        // halves separately overridable by accident. The CPA columns beside it *are*
        // family-split, because the two families are separate fits of a different
        // quantity rather than two readings of the same one.
        return value;
    }
    tables()
        .1
        .get(&(first.trim().to_lowercase(), second.trim().to_lowercase()))
        .map_or(0.0, |interaction| match cubic {
            // NeqSim's own rule: exactly `PhasePrEos` reads `KIJPR` and every other
            // phase - SRK, RK, and even `PhaseUMRCPA`, which extends `PhasePrEos`
            // without being it - reads `KIJSRK`.
            Cubic::Pr => interaction.kij_pr,
            _ => interaction.kij_srk,
        })
}

/// The Wilke-Chang association parameter for a solvent, by name.
///
/// A name resolves through a lowercased lookup against NeqSim's table; a solvent
/// that is not listed - a non-associated one, most hydrocarbons - is 1.0. NeqSim's
/// table mixes cases: the six uppercase keys (`MEG`, `DEG`, `TEG`, `MDEA`, `MEA`,
/// `DEA`) are never matched by the lowercased lookup, so they fall through to 1.0
/// exactly as NeqSim's own `getAssociationParameter` does.
#[must_use]
pub fn wilke_chang_phi(name: &str) -> f64 {
    match name.trim().to_lowercase().as_str() {
        "water" | "h2o" | "d2o" => 2.26,
        "methanol" => 1.9,
        "ethanol" => 1.5,
        "1-propanol" | "2-propanol" => 1.2,
        "1-butanol" | "n-butanol" => 1.0,
        "acetic acid" => 1.3,
        "formic acid" => 1.6,
        _ => 1.0,
    }
}

/// The MBWR-32 coefficients of a substance by name, if the table carries it.
///
/// Only methane and ethane have MBWR-32 parameters in NeqSim's `mbwr32param`; every
/// other name is `None`, and a BWRS model refuses it rather than estimating a density.
#[must_use]
pub fn bwrs_coefficients(name: &str) -> Option<BwrsCoefficients> {
    fn parse() -> HashMap<String, BwrsCoefficients> {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(MBWR32_CSV.as_bytes());
        let mut out = HashMap::new();
        for record in reader.records().flatten() {
            // Columns: id, name, a0..a31, rhoc.
            let name = record.get(1).unwrap_or("").trim().to_lowercase();
            if name.is_empty() {
                continue;
            }
            let mut a = [0.0; 32];
            for (i, slot) in a.iter_mut().enumerate() {
                *slot = record
                    .get(2 + i)
                    .unwrap_or("")
                    .trim()
                    .parse()
                    .unwrap_or(0.0);
            }
            let rhoc = record.get(34).unwrap_or("").trim().parse().unwrap_or(0.0);
            out.insert(name, BwrsCoefficients { a, rhoc });
        }
        out
    }
    static TABLE: OnceLock<HashMap<String, BwrsCoefficients>> = OnceLock::new();
    TABLE
        .get_or_init(parse)
        .get(&name.trim().to_lowercase())
        .copied()
}

/// One ion pair's Pitzer parameters, as `PitzerParameters.csv` states them.
///
/// The three binary coefficients are **temperature-dependent**, and the form is NeqSim's
/// own (`PhasePitzer`, the comment above its `beta0T1` field):
///
/// ```text
/// beta0(T) = beta0_25 + t1 (1/T - 1/Tr) + t2 ln(T/Tr),   Tr = 298.15 K
/// ```
///
/// and likewise for `beta1` and `Cphi`. The `*_25` names are the values at `Tr` and the
/// `*_t` pairs are `(t1, t2)`, so a caller evaluates the correlation rather than reading
/// a number that is only right at one temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct PitzerRecord {
    /// The first ion, lower case, exactly as the file spells it: `na+`, `cl-`, `so4--`.
    pub ion1: String,
    /// The second ion. Which of the two is the cation is the caller's to know, and
    /// [`pitzer_pair`] answers either order.
    pub ion2: String,
    /// `beta0` at 298.15 K, dimensionless.
    pub beta0_25: f64,
    /// `beta1` at 298.15 K, dimensionless.
    pub beta1_25: f64,
    /// `Cphi` at 298.15 K, dimensionless.
    pub cphi_25: f64,
    /// `(t1, t2)` of `beta0(T)`; see the type's own comment for the form.
    pub beta0_t: [f64; 2],
    /// `(t1, t2)` of `beta1(T)`.
    pub beta1_t: [f64; 2],
    /// `(t1, t2)` of `Cphi(T)`.
    pub cphi_t: [f64; 2],
    /// `beta2`, the term a 2:2 pair carries and a 1:1 pair does not.
    ///
    /// **Nonzero on exactly four rows**: `ca++`, `mg++`, `sr++` and `fe++` against
    /// `so4--`, all negative, and zero on the other 26. So this column is not a
    /// parameter most pairs have - a caller evaluating it must know which pairs do.
    pub beta2_25: f64,
    /// `theta`, the same-sign ion mixing parameter (Harvie and Weare, 1984).
    ///
    /// **Zero on every row of the shipped table.** There are no same-sign pairs in it at
    /// all - all 30 are cation with anion - so Pitzer's same-sign term has no parameters
    /// here, and the zero is the absence rather than a fitted ideal solution. A model
    /// that needs it must refuse rather than evaluating a zero.
    pub theta: f64,
    /// `psi`, the ternary parameter of a cation-cation-anion or anion-anion-cation
    /// triplet, stated against the pair whose interaction it modifies.
    ///
    /// **Zero on every row too**, for the reason [`Self::theta`] gives: a ternary
    /// parameter needs a same-sign pair to hang off, and the table carries none.
    pub psi_common_ion: f64,
    /// The pair's fitted validity range, in K. **A caller outside it is extrapolating**,
    /// and the table states the range rather than leaving it to be discovered.
    pub t_min: f64,
    /// The upper end of the same range, in K.
    pub t_max: f64,
    /// The literature reference the file cites for this pair.
    pub reference: String,
}

/// Every row of the Pitzer pair table, in file order.
#[must_use]
pub fn pitzer_parameters() -> &'static [PitzerRecord] {
    fn parse() -> Vec<PitzerRecord> {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(PITZER_PARAMETERS_CSV.as_bytes());
        let index = column_index(PITZER_PARAMETERS_CSV);
        let mut out = Vec::new();
        for record in reader.records().flatten() {
            let field = |name: &str| -> String { field_of(&index, &record, name).to_owned() };
            let number = |name: &str| -> f64 {
                field_of(&index, &record, name)
                    .trim()
                    .parse()
                    .unwrap_or_default()
            };
            let ion1 = field("ion1").to_lowercase();
            if ion1.is_empty() {
                continue;
            }
            out.push(PitzerRecord {
                ion1,
                ion2: field("ion2").to_lowercase(),
                beta0_25: number("beta0_25"),
                beta1_25: number("beta1_25"),
                cphi_25: number("cphi_25"),
                beta0_t: [number("beta0_t1"), number("beta0_t2")],
                beta1_t: [number("beta1_t1"), number("beta1_t2")],
                cphi_t: [number("cphi_t1"), number("cphi_t2")],
                beta2_25: number("beta2_25"),
                theta: number("theta"),
                psi_common_ion: number("psi_common_ion"),
                t_min: number("tmin"),
                t_max: number("tmax"),
                reference: field("reference"),
            });
        }
        out
    }
    static TABLE: OnceLock<Vec<PitzerRecord>> = OnceLock::new();
    TABLE.get_or_init(parse)
}

/// The record for an ion pair, either order round.
///
/// Both orders because the table stores a pair once and the physical interaction is
/// symmetric: a caller holding `("cl-", "na+")` should not have to know the file wrote
/// sodium first. `None` for a pair the table does not carry, which a model must refuse
/// rather than evaluate at zero - see [`pitzer_parameters`].
#[must_use]
pub fn pitzer_pair(first: &str, second: &str) -> Option<&'static PitzerRecord> {
    let (a, b) = (first.trim().to_lowercase(), second.trim().to_lowercase());
    pitzer_parameters()
        .iter()
        .find(|r| (r.ion1 == a && r.ion2 == b) || (r.ion1 == b && r.ion2 == a))
}

/// One salt's record, as `COMPSALT.csv` states it.
#[derive(Debug, Clone, PartialEq)]
pub struct SaltRecord {
    /// The salt's name as the file spells it: `NaCl`, `CaCO3`.
    pub name: String,
    /// The cation, lower case, matching a [`crate::databank::entry`] name.
    pub cation: String,
    /// The anion, likewise.
    pub anion: String,
    /// How many cations one formula unit dissociates into.
    pub cation_stoichiometry: f64,
    /// How many anions.
    pub anion_stoichiometry: f64,
    /// The five solubility-product coefficients. **Their form is not established here**
    /// - the manifest records the column as `neqsim-internal` - so they are carried as
    ///   the file states them and a model that needs the correlation must read NeqSim's
    ///   `ChemicalReactionOperations` rather than assume a polynomial.
    pub ksp: [f64; 5],
    /// `Vdelta`, the molar volume change on dissolution. Unit unestablished, as above.
    pub volume_delta: f64,
    /// `waterstoc`, how many waters of hydration the dissolution carries. Zero on every
    /// row of the shipped table, which is the file's marker for a salt with none.
    pub water_stoichiometry: f64,
}

/// Every row of the salt table, in file order.
#[must_use]
pub fn salts() -> &'static [SaltRecord] {
    fn parse() -> Vec<SaltRecord> {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(COMPSALT_CSV.as_bytes());
        let index = column_index(COMPSALT_CSV);
        let mut out = Vec::new();
        for record in reader.records().flatten() {
            let field = |name: &str| -> String { field_of(&index, &record, name).to_owned() };
            let number = |name: &str| -> f64 {
                field_of(&index, &record, name)
                    .trim()
                    .parse()
                    .unwrap_or_default()
            };
            let name = field("saltname");
            if name.is_empty() {
                continue;
            }
            out.push(SaltRecord {
                name,
                cation: field("ion1").to_lowercase(),
                anion: field("ion2").to_lowercase(),
                cation_stoichiometry: number("stoc1"),
                anion_stoichiometry: number("stoc2"),
                ksp: [
                    number("kspwater"),
                    number("kspwater2"),
                    number("kspwater3"),
                    number("kspwater4"),
                    number("kspwater5"),
                ],
                volume_delta: number("vdelta"),
                water_stoichiometry: number("waterstoc"),
            });
        }
        out
    }
    static TABLE: OnceLock<Vec<SaltRecord>> = OnceLock::new();
    TABLE.get_or_init(parse)
}

/// The record for a salt by name, matched without regard to case or surrounding space.
#[must_use]
pub fn salt(name: &str) -> Option<&'static SaltRecord> {
    let key = name.trim().to_lowercase();
    salts().iter().find(|s| s.name.to_lowercase() == key)
}

/// A compiled column's position, by the name this project writes it under.
///
/// Read off the file's own first line rather than listed, so a column the manifest renames
/// resolves by the new name and one it stops carrying resolves to nothing.
fn column_index(text: &str) -> HashMap<&str, usize> {
    text.lines()
        .next()
        .unwrap_or("")
        .split(',')
        .enumerate()
        .map(|(position, name)| (name.trim(), position))
        .collect()
}

/// A column's value by the compiled header name it was written under.
///
/// **By name rather than by position**, because both tables' headers are the manifest's
/// `as` names and a column inserted in the middle should resolve to the wrong *name*
/// rather than shift every value one field left. A name a table does not carry is the
/// empty string, which the numeric readers turn into a zero - the table's own marker for
/// a parameter it does not state, and never a fabricated value.
fn field_of<'a>(
    index: &HashMap<&str, usize>,
    record: &'a csv::StringRecord,
    name: &str,
) -> &'a str {
    index
        .get(name)
        .and_then(|position| record.get(*position))
        .unwrap_or("")
}

/// Every substance name available, sorted: the table plus whatever an overlay adds.
#[must_use]
pub fn names(overlay: Option<&Overlay>) -> Vec<String> {
    let mut out: Vec<String> = tables().0.keys().cloned().collect();
    if let Some(overlay) = overlay {
        out.extend(overlay.components.keys().cloned());
    }
    out.sort();
    out.dedup();
    out
}

/// Every substance, ordered by name.
///
/// Ordered rather than in file order because the two implementations must be able to
/// compare them one for one, and a sort is the one ordering both can reproduce without
/// agreeing on how the file happens to be laid out.
#[must_use]
pub fn all_entries() -> Vec<&'static Entry> {
    let mut out: Vec<&Entry> = tables().0.values().collect();
    out.sort_by(|left, right| left.name.cmp(&right.name));
    out
}

/// Every interaction pair the databank carries, ordered, each pair once, as
/// `(first, second, KIJPR, KIJSRK)`.
///
/// **Both columns, because which one a mixture uses is decided by its cubic.** A caller
/// that collapses them has chosen a fluid for a cubic it does not know about yet.
///
/// The table is stored both ways round so a caller need not know which name came first;
/// this un-does that by keeping only the ordering where the first name sorts lower.
#[must_use]
pub fn all_kij() -> Vec<(String, String, f64, f64)> {
    let mut out: Vec<(String, String, f64, f64)> = tables()
        .1
        .iter()
        .filter(|((a, b), _)| a < b)
        .map(|((a, b), interaction)| {
            (
                a.clone(),
                b.clone(),
                interaction.kij_pr,
                interaction.kij_srk,
            )
        })
        .collect();
    out.sort_by(|left, right| (&left.0, &left.1).cmp(&(&right.0, &right.1)));
    out
}

/// The CPA interaction matrix for a list of names, flattened row-major.
///
/// **`cpakij_SRK` / `cpakij_PR`, not `KIJPR`.** NeqSim's `CPAMixingRuleHandler` reads the
/// `cpa` columns and `SystemSrkCPA` runs `setMixingRule(10)`, so an associating mixture's
/// attraction is mixed with these and a classical mixture's with `KIJPR`. On
/// water/methanol the two differ by a factor of two - `-0.153` against `-0.0789` - which
/// is a 3.4% difference in the mixture's `A` and 0.077% in its root.
///
/// An absent pair is zero, the ideal-mixture default NeqSim substitutes.
#[must_use]
pub fn cpa_kij(
    names: &[&str],
    cubic: crate::association::AssociationCubic,
    overlay: Option<&Overlay>,
) -> Vec<f64> {
    let n = names.len();
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            if let Some(value) = overlay.and_then(|o| o.cpa_kij_value(names[i], names[j], cubic)) {
                out[i * n + j] = value;
                continue;
            }
            let key = (
                names[i].trim().to_lowercase(),
                names[j].trim().to_lowercase(),
            );
            out[i * n + j] = tables().1.get(&key).map_or(0.0, |interaction| match cubic {
                crate::association::AssociationCubic::Srk => interaction.cpa_kij_srk,
                crate::association::AssociationCubic::Pr => interaction.cpa_kij_pr,
            });
        }
    }
    out
}

/// The PC-SAFT interaction matrix for a list of names, flattened row-major.
///
/// **`KIJPCSAFT`, a third fit.** Methane/n-butane is `0.022` here against the
/// `0.01289789` its SRK and PR columns share, and only 42 of the 516 in-scope pairs carry
/// one at all - so this is a sparse column of its own rather than a placeholder.
///
/// **NeqSim reads it whenever its PC-SAFT runs at all.** Its default mixing rule leaves
/// the matrix null and a mixture throws a NullPointerException before reaching a `k_ij`,
/// so the only configuration a caller can use is the classic one - whose branch is keyed
/// on the PC-SAFT phase class and reads this column.
///
/// An absent pair is zero, the ideal-mixture default NeqSim substitutes.
///
/// **A card cannot state one yet.** The keycard's pair record carries the cubic's `kij` and
/// the two CPA columns, and a third would be a schema change with no model needing it until
/// this one is read through a card - so the column is the table's alone for now rather than
/// a map nothing populates.
#[must_use]
pub fn pcsaft_kij(names: &[&str]) -> Vec<f64> {
    let n = names.len();
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let key = (
                names[i].trim().to_lowercase(),
                names[j].trim().to_lowercase(),
            );
            out[i * n + j] = tables()
                .1
                .get(&key)
                .map_or(0.0, |interaction| interaction.pcsaft_kij);
        }
    }
    out
}

/// The mixture an associating model runs on, from names.
///
/// **One call, because the four steps it replaces were four chances to be wrong.**
/// `SystemSrkCPA` mixes with `cpakij_SRK` and substitutes the fitted `a` and `b` in place
/// of the cubic's; a classical mixture reads the cubic's own column. The two are not
/// interchangeable and neither is a refinement of the other, so this resolves the whole
/// mixture rather than leaving a caller to assemble [`mixture_of`], [`cpa_kij`],
/// `with_cubic` and `with_association` in the right order with the right arguments.
///
/// The cubic is stated once and read three ways from here - the interaction column (see
/// [`mixture_of`]), the fitted association family and the root finder - which is the
/// point: water's `kappa_AB` is 0.0692 for SRK against 0.046473789 for PR, so a mixture
/// that took its cubic from one call and its association from another is a different
/// fluid.
///
/// # Errors
/// * As [`mixture_of`], plus the association's own refusals: a mixture whose association
///   is provably inert, or a component whose fitted set is negative.
pub fn associating_mixture_of(
    names: &[&str],
    cubic: Cubic,
    overlay: Option<&Overlay>,
) -> Result<(Mixture, IdealGasModel)> {
    let family = AssociationCubic::of(cubic);
    let (mixture, ideal_gas) = mixture_of(names, cubic, overlay)?;
    let mixture = mixture
        .with_mixing_rule(MixingRule::Classic {
            kij: cpa_kij(names, family, overlay),
        })
        .with_cubic(cubic)
        .with_association()?;
    Ok((mixture, ideal_gas))
}

/// The UMR-CPA fluid: Peng-Robinson, the UMR mixing rule and the association.
///
/// NeqSim's `SystemUMRCPAEoS`, and the only mixture in this library whose attraction is
/// mixed by a universal rule rather than an interaction matrix. Four choices are made
/// here and nowhere else:
///
/// * the cubic is Peng-Robinson, which `PhaseUMRCPA extends PhasePrEos` states;
/// * the mixing rule is [`MixingRule::Umr`] with the `_umrmc` UNIFAC tables, which is
///   what `SystemUMRCPAEoS`'s construction and `SystemThermo`'s `"UMR-CPA"` selector
///   both name;
/// * the interaction matrix is **zero**. NeqSim's UMR rule reads no `kij` column for
///   its `alpha_mix`, and the probe's own reading - `A = n B R T alpha_mix` with
///   `B = sum_i x_i b_i` - is reproduced with every pair at zero;
/// * the alpha is the five-parameter Mathias-Copeman term (`UMRCPA_MC1..5`), which is
///   `ComponentUMRCPA.setAttractiveTerm`'s term 22 and the only alpha under which the
///   substituted `a` and `b` are the ones NeqSim solves with.
///
/// # Errors
/// * As [`associating_mixture_of`], plus
///   [`AzothError::PropertyUnavailable`] if a name has no `UNIFACcompUMRPRU` group
///   decomposition, which the rule cannot mix without.
pub fn umr_cpa_mixture_of(names: &[&str]) -> Result<(Mixture, IdealGasModel)> {
    let (base, ideal_gas) = mixture_of(names, Cubic::Pr, None)?;
    let n = names.len();
    // **Every component needs its `UMRCPA_MC` set, and one without is refused.** NeqSim's
    // `ComponentUMRCPA.setAttractiveTerm` falls back to term 19 seeded with `MCPR1..3`
    // for a non-associating component with none, and to term 1 with `mCPA` for an
    // associating one - and it chooses **per component**, which this library's
    // mixture-level alpha cannot express: a mixture of a component that carries the set
    // and one that does not would need two attractive terms at once. So the model reads
    // the set it can express and refuses what it cannot, rather than giving a component
    // a different attractive term than NeqSim gives it.
    let mut components = base.components().to_vec();
    for (i, name) in names.iter().enumerate() {
        let mc = entry(name, None)?.umrcpa_mc;
        if mc.iter().all(|c| c.abs() <= 1.0e-20) {
            return Err(AzothError::property_unavailable(
                name.trim().to_lowercase(),
                "UMRCPA_MC1..5".to_string(),
                "carries no UMR-CPA Mathias-Copeman set; the UMR-CPA model's attractive \
                 term is that set, and NeqSim's per-component fallback to term 19 or term \
                 1 is not expressible as one mixture-level alpha"
                    .to_string(),
            ));
        }
        components[i] = components[i].clone().with_alpha_params(mc.to_vec());
    }
    let unifac = unifac_umrpru_parameters(names, UmrpruSet::Umrmc)?;
    let mixture = Mixture::new(components, vec![0.0; n * n])?
        .with_cubic(Cubic::Pr)
        .with_alpha(Alpha::MatCop5PrUmr)
        .with_mixing_rule(MixingRule::Umr {
            kij: vec![0.0; n * n],
            unifac,
        })
        .with_association()?;
    Ok((mixture, ideal_gas))
}

/// The NRTL non-randomness matrix for a list of names, flattened row-major.
///
/// `alpha[i][j] = alpha_ij`, symmetric; the diagonal and any pair the table does not
/// carry are `0.0`, the ideal-mixture default.
#[must_use]
pub fn nrtl_alpha(names: &[&str]) -> Vec<f64> {
    let table = &tables().1;
    let n = names.len();
    let mut out = vec![0.0; n * n];
    for (i, a) in names.iter().enumerate() {
        for (j, b) in names.iter().enumerate() {
            if i != j {
                out[i * n + j] = table
                    .get(&(a.trim().to_lowercase(), b.trim().to_lowercase()))
                    .map_or(0.0, |interaction| interaction.alpha);
            }
        }
    }
    out
}

/// The NRTL energy matrix `Dij` for a list of names, flattened row-major.
///
/// `dij[i][j] = g_ij`, the Kelvin energy in `tau_ij = g_ij / T`; directional, so
/// `dij[i][j]` and `dij[j][i]` differ in general. The diagonal and any absent pair are
/// `0.0`.
#[must_use]
pub fn nrtl_dij(names: &[&str]) -> Vec<f64> {
    let table = &tables().1;
    let n = names.len();
    let mut out = vec![0.0; n * n];
    for (i, a) in names.iter().enumerate() {
        for (j, b) in names.iter().enumerate() {
            if i != j {
                out[i * n + j] = table
                    .get(&(a.trim().to_lowercase(), b.trim().to_lowercase()))
                    .map_or(0.0, |interaction| interaction.gij);
            }
        }
    }
    out
}

/// The resolved Huron-Vidal parameters for a mixture, each matrix flattened row-major.
///
/// The `MixingRule::HuronVidal` variant's five matrices, resolved by name rather than
/// built by hand. The pairing is NeqSim's: a pair the interaction table marks `HV`
/// carries the fitted NRTL parameters, and every other pair carries the cubic's own
/// excess energy, which [`MixingRule::HuronVidal`] computes from `kij`, `a` and `b`.
#[derive(Debug, Clone, PartialEq)]
pub struct HuronVidalParameters {
    /// The cubic interaction matrix, as in [`MixingRule::Classic`].
    pub kij: Vec<f64>,
    /// The fitted NRTL energy `Dij` in kelvin, `N x N` row-major. Directional.
    pub hv_gij: Vec<f64>,
    /// The fitted temperature coefficient `DijT`, `N x N` row-major. Directional.
    pub hv_gij_t: Vec<f64>,
    /// The fitted non-randomness `alpha`, `N x N` row-major. Symmetric.
    pub hv_alpha: Vec<f64>,
    /// One flag per interaction: `true` where `HVTYPE` says `HV`.
    pub hv_pairs: Vec<bool>,
}

/// The resolved Wong-Sandler parameters for a mixture, each matrix flattened row-major.
///
/// The same shape as [`HuronVidalParameters`] with the rule's own `kij` and `DijT`: the
/// selectors come from `WSTYPE`, the interaction from `KIJWSunifac`, and the temperature
/// coefficient from `WSGIJT`/`WSGJIT`.
#[derive(Debug, Clone, PartialEq)]
pub struct WongSandlerParameters {
    /// The interaction matrix the rule's `b_mix` reads, from `KIJWSunifac`.
    pub kij: Vec<f64>,
    /// The fitted NRTL energy `Dij` in kelvin, `N x N` row-major. Directional.
    pub hv_gij: Vec<f64>,
    /// The fitted temperature coefficient `DijT`, `N x N` row-major. Directional.
    pub hv_gij_t: Vec<f64>,
    /// The fitted non-randomness `alpha`, `N x N` row-major. Symmetric.
    pub hv_alpha: Vec<f64>,
    /// One flag per interaction: `true` where `WSTYPE` says `WS`.
    pub hv_pairs: Vec<bool>,
}

/// The Huron-Vidal parameters for a list of names, from the interaction table.
///
/// The resolution the rule leaves to its caller otherwise, and the one the manifest
/// described as missing: NeqSim's `EosMixingRuleHandler` reads `HVTYPE` into
/// `classicOrHV` and tests `mixRule[j][i].trim().equals("HV")`, which is exactly
/// [`HuronVidalParameters::hv_pairs`].
///
/// Each component must be in the databank; a pair the interaction table does not carry
/// resolves to `Classic` with zero fitted parameters, which is the cubic's own excess
/// energy rather than an ideal mixture.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name is in neither the databank nor the
///   overlay.
pub fn huron_vidal_parameters(
    names: &[&str],
    cubic: Cubic,
    overlay: Option<&Overlay>,
) -> Result<HuronVidalParameters> {
    for name in names {
        entry(name, overlay)?;
    }
    let table = &tables().1;
    let n = names.len();
    let mut out = HuronVidalParameters {
        kij: vec![0.0; n * n],
        hv_gij: vec![0.0; n * n],
        hv_gij_t: vec![0.0; n * n],
        hv_alpha: vec![0.0; n * n],
        hv_pairs: vec![false; n * n],
    };
    for (i, a) in names.iter().enumerate() {
        for (j, b) in names.iter().enumerate() {
            if i == j {
                continue;
            }
            let key = (a.trim().to_lowercase(), b.trim().to_lowercase());
            let Some(interaction) = table.get(&key) else {
                continue;
            };
            out.kij[i * n + j] = match cubic {
                Cubic::Pr => interaction.kij_pr,
                _ => interaction.kij_srk,
            };
            out.hv_pairs[i * n + j] = interaction.hv;
            out.hv_alpha[i * n + j] = interaction.hv_alpha;
            out.hv_gij[i * n + j] = interaction.hv_dij;
            out.hv_gij_t[i * n + j] = interaction.hv_dij_t;
        }
    }
    Ok(out)
}

/// The Wong-Sandler parameters for a list of names, from the interaction table.
///
/// `WSTYPE` selects the pairs, `KIJWSunifac` is the interaction and `WSGIJT`/`WSGJIT`
/// the temperature coefficient. `KIJWS` is not read, and cannot be: NeqSim assigns
/// `WSintparam` from it and overwrites that from `KIJWSunifac` on the next line.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name is in neither the databank nor the
///   overlay.
pub fn wong_sandler_parameters(
    names: &[&str],
    overlay: Option<&Overlay>,
) -> Result<WongSandlerParameters> {
    for name in names {
        entry(name, overlay)?;
    }
    let table = &tables().1;
    let n = names.len();
    let mut out = WongSandlerParameters {
        kij: vec![0.0; n * n],
        hv_gij: vec![0.0; n * n],
        hv_gij_t: vec![0.0; n * n],
        hv_alpha: vec![0.0; n * n],
        hv_pairs: vec![false; n * n],
    };
    for (i, a) in names.iter().enumerate() {
        for (j, b) in names.iter().enumerate() {
            if i == j {
                continue;
            }
            let key = (a.trim().to_lowercase(), b.trim().to_lowercase());
            let Some(interaction) = table.get(&key) else {
                continue;
            };
            out.kij[i * n + j] = interaction.kij_ws;
            out.hv_pairs[i * n + j] = interaction.ws;
            out.hv_alpha[i * n + j] = interaction.hv_alpha;
            out.hv_gij[i * n + j] = interaction.hv_dij;
            out.hv_gij_t[i * n + j] = interaction.ws_dij_t;
        }
    }
    Ok(out)
}

/// One component's pure-liquid vapour-pressure correlation, as the table states it.
///
/// The coefficients, the label naming which of NeqSim's four Antoine forms they belong
/// to, and the critical constants the Wagner form reads. Carried whole rather than
/// resolved to a `p_sat`, because the phase evaluates it at its own temperature.
#[derive(Debug, Clone, PartialEq)]
pub struct AntoineRecord {
    /// NeqSim's own label, which `form_from_type` turns into the arithmetic.
    pub antoine_type: String,
    /// The five coefficients `A`-`E`, in NeqSim's internal scale.
    pub coefficients: [f64; 5],
    /// Critical temperature, in K.
    pub tc: f64,
    /// Critical pressure, in Pa.
    pub pc: f64,
}

/// The resolved parameters of an NRTL activity-coefficient *phase*.
///
/// The NRTL matrices [`NrtlParameters`] carries, beside what a phase needs and an
/// activity coefficient does not: each component's pure-liquid vapour pressure, because
/// the phase's fugacity coefficient is `gamma_i P0_i / P`.
#[derive(Debug, Clone, PartialEq)]
pub struct GeNrtlPhaseParameters {
    /// `alpha[i][j]`, `N x N` row-major. Symmetric with a zero diagonal.
    pub alpha: Vec<f64>,
    /// `dij[i][j] = g_ij` in kelvin, `N x N` row-major. Directional.
    pub dij: Vec<f64>,
    /// One vapour-pressure correlation per component, in order.
    pub antoine: Vec<AntoineRecord>,
}

/// The resolved parameters of a Van Laar acid activity-coefficient *phase*.
///
/// [`VanLaarAcidParameters`]' acid identity, beside each component's Antoine columns and
/// the critical constants they need. The columns are carried for *every* component, but
/// they are read only for the ones with `acid_index == 0`: the three modelled acids take
/// their `P0` from [`crate::nitric_sulfuric_acid_vapor_pressure`] instead, which is the
/// whole point of the phase - see its spec.
#[derive(Debug, Clone, PartialEq)]
pub struct GeVanLaarAcidPhaseParameters {
    /// Per component: `1` water, `2` nitric acid, `3` sulfuric acid, `0` for a species
    /// the model does not cover.
    pub acid_index: Vec<u8>,
    /// One vapour-pressure correlation per component, in order. Read only where
    /// `acid_index` is zero.
    pub antoine: Vec<AntoineRecord>,
}

/// The Van Laar acid phase parameters for a list of names.
///
/// **This resolver does not refuse a Henry's-law solute, and its four siblings do.** The
/// refusal there exists because `ComponentGE.fugcoef` branches on `referenceStateType`, so
/// computing the Raoult expression for a `solute` component would be the wrong branch with
/// no symptom. `ComponentGEVanLaarAcid` *overrides* `fugcoef` precisely to ignore the tag -
/// and it has to, because both acids are tagged `solute` and carry an all-zero Antoine row.
/// Refusing them here would refuse the model's own subject.
///
/// What it does refuse is a component the model does not cover (`acid_index == 0`) that
/// also has no Antoine correlation, because then its `P0` has no source at all.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name is in neither the databank nor the
///   overlay, or if an uncovered component carries no Antoine correlation.
pub fn ge_van_laar_acid_phase_parameters(
    names: &[&str],
    overlay: Option<&Overlay>,
) -> Result<GeVanLaarAcidPhaseParameters> {
    let acid = van_laar_acid_parameters(names, overlay)?;
    let mut antoine = Vec::with_capacity(names.len());
    for (name, &index) in names.iter().zip(&acid.acid_index) {
        let record = entry(name, overlay)?;
        let Some((coefficients, antoine_type)) = record.antoine.clone() else {
            return Err(AzothError::property_unavailable(
                record.name,
                "Antoine vapour-pressure coefficients".to_string(),
                "a species the acid model does not cover takes its `P0` from Antoine, \
                 since only water, nitric acid and sulfuric acid have a correlation \
                 here"
                    .to_string(),
            ));
        };
        // Only an *uncovered* component needs a real correlation; the three acids take
        // theirs from the calc. Checked here rather than at evaluation so the refusal
        // names the component.
        if index == 0 && coefficients.iter().all(|value| *value == 0.0) {
            return Err(AzothError::property_unavailable(
                record.name,
                "Antoine vapour-pressure coefficients".to_string(),
                "a species the acid model does not cover takes its `P0` from Antoine, \
                 and the databank carries none for this one"
                    .to_string(),
            ));
        }
        antoine.push(AntoineRecord {
            antoine_type,
            coefficients,
            tc: record.tc,
            pc: record.pc,
        });
    }
    Ok(GeVanLaarAcidPhaseParameters {
        acid_index: acid.acid_index,
        antoine,
    })
}

/// The resolved parameters of a UNIQUAC activity-coefficient *phase*.
///
/// The volume and surface parameters [`UniquacParameters`] carries, beside what a phase
/// needs and an activity coefficient does not: each component's pure-liquid vapour
/// pressure, because the phase's fugacity coefficient is `gamma_i P0_i / P`.
///
/// `aij` is *not* here. It is the caller's - no upstream table carries a UNIQUAC
/// interaction matrix - and it arrives as its own argument, the way it does for
/// [`crate::uniquac_activity_coefficients`].
#[derive(Debug, Clone, PartialEq)]
pub struct GeUniquacPhaseParameters {
    /// The van der Waals volume parameter `r_i` of each component.
    pub r: Vec<f64>,
    /// The van der Waals surface-area parameter `q_i` of each component.
    pub q: Vec<f64>,
    /// One vapour-pressure correlation per component, in order.
    pub antoine: Vec<AntoineRecord>,
}

/// The UNIQUAC phase parameters for a list of names.
///
/// The `r` and `q` resolve through [`uniquac_parameters`], so the phase and
/// `eos.uniquac_activity_coefficients` cannot disagree about them; this adds the
/// per-component vapour pressure beside them.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name has no UNIFAC group assignment or no
///   Antoine correlation.
/// * [`AzothError::InvalidInput`] if a component is tagged a Henry's-law solute, which
///   this phase does not implement.
pub fn ge_uniquac_phase_parameters(
    names: &[&str],
    overlay: Option<&Overlay>,
) -> Result<GeUniquacPhaseParameters> {
    let uniquac = uniquac_parameters(names)?;
    Ok(GeUniquacPhaseParameters {
        r: uniquac.r,
        q: uniquac.q,
        antoine: phase_antoine(names, overlay)?,
    })
}

/// The resolved parameters of a Wilson activity-coefficient *phase*.
///
/// Only the vapour-pressure columns: the Wilson correlation reads the component's molar
/// mass and critical temperature, which are what a [`Mixture`] carries and
/// [`databank::mixture_of`] already resolves, so the phase takes that mixture beside this
/// record rather than a second copy of the same two vectors.
#[derive(Debug, Clone, PartialEq)]
pub struct GeWilsonPhaseParameters {
    /// One vapour-pressure correlation per component, in order.
    pub antoine: Vec<AntoineRecord>,
}

/// The Wilson phase parameters for a list of names.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name is in neither the databank nor the
///   overlay, or if it carries no Antoine correlation.
/// * [`AzothError::InvalidInput`] if a component is tagged a Henry's-law solute, which
///   this phase does not implement.
pub fn ge_wilson_phase_parameters(
    names: &[&str],
    overlay: Option<&Overlay>,
) -> Result<GeWilsonPhaseParameters> {
    Ok(GeWilsonPhaseParameters {
        antoine: phase_antoine(names, overlay)?,
    })
}

/// The resolved parameters of a UNIFAC activity-coefficient *phase*.
///
/// The group tables [`UnifacParameters`] carries, beside what a phase needs and an
/// activity coefficient does not: each component's pure-liquid vapour pressure, because
/// the phase's fugacity coefficient is `gamma_i P0_i / P`.
#[derive(Debug, Clone, PartialEq)]
pub struct GeUnifacPhaseParameters {
    /// Per-component group counts, `N x G` row-major.
    pub groups: Vec<f64>,
    /// The volume `R` of each group, length `G`.
    pub group_r: Vec<f64>,
    /// The surface area `Q` of each group, length `G`.
    pub group_q: Vec<f64>,
    /// The main-group interaction matrix, `G x G` row-major, in Kelvin.
    pub aij: Vec<f64>,
    /// One vapour-pressure correlation per component, in order.
    pub antoine: Vec<AntoineRecord>,
}

/// The UNIFAC phase parameters for a list of names.
///
/// The group tables resolve through [`unifac_parameters`], so the phase and
/// `eos.unifac_activity_coefficients` cannot disagree about them; this adds the
/// per-component vapour pressure beside them.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name is in neither the databank nor the
///   overlay, if it has no UNIFAC group assignment, or if it carries no Antoine
///   correlation.
/// * [`AzothError::InvalidInput`] if a component is tagged a Henry's-law solute, which
///   this phase does not implement.
pub fn ge_unifac_phase_parameters(
    names: &[&str],
    overlay: Option<&Overlay>,
) -> Result<GeUnifacPhaseParameters> {
    let unifac = unifac_parameters(names)?;
    Ok(GeUnifacPhaseParameters {
        groups: unifac.groups,
        group_r: unifac.group_r,
        group_q: unifac.group_q,
        aij: unifac.aij,
        antoine: phase_antoine(names, overlay)?,
    })
}

/// The per-component vapour-pressure records a phase reads, or a refusal.
///
/// Shared by every activity-coefficient phase, because the two things it checks are the
/// same for all of them: a component tagged anything but `solvent` takes a Henry's-law
/// coefficient in `ComponentGE.fugcoef`, which this library does not implement, and a
/// component with no correlation has no `P0` to compose.
fn phase_antoine(names: &[&str], overlay: Option<&Overlay>) -> Result<Vec<AntoineRecord>> {
    let mut antoine = Vec::with_capacity(names.len());
    for name in names {
        let entry = entry(name, overlay)?;
        if entry.reference_state != SOLVENT {
            return Err(AzothError::invalid_input(
                "components",
                format!(
                    "`{}` is tagged `referenceStateType = {}` in NeqSim's component \
                     database, so `ComponentGE.fugcoef` gives it a Henry's-law fugacity \
                     coefficient rather than `gamma_i P0_i / P`. That branch is not \
                     ported, and this phase is the Raoult one. The substances the \
                     databank tags `solvent` - water, the alcohols, the glycols - are \
                     the ones it describes.",
                    entry.name, entry.reference_state
                ),
            ));
        }
        let Some((coefficients, antoine_type)) = entry.antoine else {
            return Err(AzothError::property_unavailable(
                entry.name,
                "Antoine vapour-pressure coefficients".to_string(),
                "the phase's fugacity coefficient is `gamma_i P0_i / P`, so every \
                 component needs a correlation; a keycard supplies the parameters a \
                 cubic reads and not these"
                    .to_string(),
            ));
        };
        antoine.push(AntoineRecord {
            antoine_type,
            coefficients,
            tc: entry.tc,
            pc: entry.pc,
        });
    }
    Ok(antoine)
}

/// The NRTL phase parameters for a list of names.
///
/// The activity-coefficient matrices resolve through [`nrtl_parameters`], so the phase
/// and `eos.nrtl_activity_coefficients` cannot disagree about them; this adds the
/// per-component vapour pressure beside them.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name is in neither the databank nor the
///   overlay, or if the databank carries it without an Antoine correlation - which is
///   what an overlay-added substance is, since a card states the parameters a cubic
///   reads.
pub fn ge_nrtl_phase_parameters(
    names: &[&str],
    overlay: Option<&Overlay>,
) -> Result<GeNrtlPhaseParameters> {
    let nrtl = nrtl_parameters(names, overlay)?;
    Ok(GeNrtlPhaseParameters {
        alpha: nrtl.alpha,
        dij: nrtl.dij,
        antoine: phase_antoine(names, overlay)?,
    })
}

/// The resolved NRTL parameters for a mixture, both matrices flattened row-major.
///
/// A caller-supplied record the way [`Component`] is, and carrying **no names** for the
/// same reason: [`nrtl_parameters`] does the lookup, and a model that looked up its own
/// inputs would answer from a file the caller never mentioned.
#[derive(Debug, Clone, PartialEq)]
pub struct NrtlParameters {
    /// `alpha[i][j]`, `N x N` row-major. Symmetric with a zero diagonal.
    pub alpha: Vec<f64>,
    /// `dij[i][j] = g_ij` in kelvin, `N x N` row-major. Directional.
    pub dij: Vec<f64>,
}

/// The NRTL `alpha` and `Dij` matrices for a list of names, flattened row-major.
///
/// The resolution `eos.nrtl_activity_coefficients` leaves to its caller, the way
/// [`bwrs_coefficients`] resolves the MBWR-32 set and [`mixture_of`] resolves a cubic's
/// constants. Both matrices come from the same `INTER.csv` row, so `alpha[i][j]` and
/// `alpha[j][i]` are the same number while `dij[i][j]` and `dij[j][i]` are not.
///
/// A pair the interaction table does not carry is `0.0`, which is an ideal interaction -
/// what NeqSim's NRTL does with an absent row, and why a pair the table has never seen
/// returns `gamma = 1`. That is a quiet answer for a real substance, and the spec says
/// so; a name that is in *neither* table is refused instead, because a typo is a mistake
/// rather than a mixture the table happens not to cover.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name is in neither the databank nor the
///   overlay.
pub fn nrtl_parameters(names: &[&str], overlay: Option<&Overlay>) -> Result<NrtlParameters> {
    for name in names {
        entry(name, overlay)?;
    }
    Ok(NrtlParameters {
        alpha: nrtl_alpha(names),
        dij: nrtl_dij(names),
    })
}

/// The resolved UNIFAC inputs for a mixture, each matrix flattened row-major.
///
/// `groups` is `N x G` (one row per component, one column per group), `group_r` and
/// `group_q` are length `G`, and `aij` is `G x G` (`a_{mn}` in Kelvin). `G` is the
/// union of the named components' subgroups, sorted by subgroup number, with absent
/// groups counted zero.
#[derive(Debug, Clone, PartialEq)]
pub struct UnifacParameters {
    /// Per-component group counts, `N x G` row-major.
    pub groups: Vec<f64>,
    /// The volume `R` of each group, length `G`.
    pub group_r: Vec<f64>,
    /// The surface area `Q` of each group, length `G`.
    pub group_q: Vec<f64>,
    /// The main-group interaction matrix, `G x G` row-major, in Kelvin.
    pub aij: Vec<f64>,
}

/// Which of UNIFAC-UMR-PRU's two parameter sets to read.
///
/// NeqSim decides this from `getComponent(0).getAttractiveTermNumber()` - the
/// `_umrmc` tables when it is 13, 19 or 22 - and that field does not exist in this
/// library, so the choice arrives as a named input instead of being inferred from the
/// mixture. It is a statement about which equation of state the caller is pairing the
/// activity model with, which is a thing the caller knows and a component record does
/// not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmrpruSet {
    /// The original UMR tables.
    Umr,
    /// The tables the Mathias-Copeman UMR-PRU variants share.
    Umrmc,
}

/// One of the three terms of UNIFAC-UMR-PRU's interaction,
/// `a_mn(T) = a_mn + b_mn (T - 298.15) + c_mn (T - 298.15)^2`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmrpruTerm {
    /// The constant term.
    A,
    /// The linear term.
    B,
    /// The quadratic term.
    C,
}

/// The resolved UNIFAC-UMR-PRU inputs for a mixture, each matrix flattened row-major.
///
/// The interaction is evaluated about 298.15 K rather than about zero, which is the
/// one place this differs from [`UnifacPsrkParameters`] beyond the tables it reads.
#[derive(Debug, Clone, PartialEq)]
pub struct UnifacUmrpruParameters {
    /// Per-component group counts, `N x G` row-major.
    pub groups: Vec<f64>,
    /// The volume `R` of each group, length `G`.
    pub group_r: Vec<f64>,
    /// The surface area `Q` of each group, length `G`.
    pub group_q: Vec<f64>,
    /// The constant term of the interaction, `G x G` row-major, in Kelvin.
    pub aij: Vec<f64>,
    /// The linear term, `G x G` row-major, in Kelvin per Kelvin.
    pub bij: Vec<f64>,
    /// The quadratic term, `G x G` row-major, in Kelvin per Kelvin squared.
    pub cij: Vec<f64>,
}

/// The UNIFAC-UMR-PRU inputs for a list of names, resolved from the vendored tables.
///
/// The group decomposition comes from `UNIFACcompUMRPRU.csv`, which carries 139
/// subgroups rather than the 133 of the plain UNIFAC table, and the interaction from
/// the `_umr` or `_umrmc` set `set` names. The group constants are the shared
/// `UNIFACGroupParam.csv`: only the decomposition and the interaction differ.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `names` is empty.
/// * [`AzothError::PropertyUnavailable`] if a name has no UMR-PRU group assignment.
pub fn unifac_umrpru_parameters(names: &[&str], set: UmrpruSet) -> Result<UnifacUmrpruParameters> {
    let tables = umrpru_tables();
    let basis = unifac_basis(names, &tables.group, &tables.members)?;
    let g = basis.union.len();

    let [a, b, c] = tables.matrices(set);
    let mut aij = vec![0.0; g * g];
    let mut bij = vec![0.0; g * g];
    let mut cij = vec![0.0; g * g];
    for (k, &s) in basis.union.iter().enumerate() {
        let (_, _, main) = *tables.group.get(&s).expect("a subgroup with a member row");
        for (m, &t) in basis.union.iter().enumerate() {
            let (_, _, main_t) = *tables.group.get(&t).expect("a subgroup with a member row");
            let key = (main, main_t);
            aij[k * g + m] = a.get(&key).copied().unwrap_or(0.0);
            bij[k * g + m] = b.get(&key).copied().unwrap_or(0.0);
            cij[k * g + m] = c.get(&key).copied().unwrap_or(0.0);
        }
    }

    Ok(UnifacUmrpruParameters {
        groups: basis.groups,
        group_r: basis.group_r,
        group_q: basis.group_q,
        aij,
        bij,
        cij,
    })
}

/// The parsed UNIFAC-UMR-PRU tables.
struct UmrpruTables {
    /// Subgroup number -> `(R, Q, main group)`, from the shared group table.
    group: HashMap<i64, (f64, f64, i64)>,
    /// Lower-cased name -> `[(subgroup, count)]`, from `UNIFACcompUMRPRU.csv`.
    members: HashMap<String, Vec<(i64, i64)>>,
    /// The six interaction tables, `[a, b, c]` for each of the two sets.
    umr: [HashMap<(i64, i64), f64>; 3],
    umrmc: [HashMap<(i64, i64), f64>; 3],
}

impl UmrpruTables {
    /// The `[a, b, c]` matrices of one set.
    fn matrices(&self, set: UmrpruSet) -> &[HashMap<(i64, i64), f64>; 3] {
        match set {
            UmrpruSet::Umr => &self.umr,
            UmrpruSet::Umrmc => &self.umrmc,
        }
    }
}

fn umrpru_tables() -> &'static UmrpruTables {
    static TABLES: OnceLock<UmrpruTables> = OnceLock::new();
    TABLES.get_or_init(|| parse_umrpru().expect("the embedded UMR-PRU tables should parse"))
}

/// The UMR-PRU tables: the shared group constants, the 139-subgroup decomposition, and
/// the six interaction matrices.
fn parse_umrpru() -> Result<UmrpruTables> {
    let group = parse_group_constants(UNIFAC_GROUP_CSV)?;
    let members = parse_members(UNIFAC_COMP_UMRPRU_CSV, 140)?;
    Ok(UmrpruTables {
        group,
        members,
        umr: [
            parse_interaction_table(UNIFAC_A_UMR_CSV)?,
            parse_interaction_table(UNIFAC_B_UMR_CSV)?,
            parse_interaction_table(UNIFAC_C_UMR_CSV)?,
        ],
        umrmc: [
            parse_interaction_table(UNIFAC_A_UMRMC_CSV)?,
            parse_interaction_table(UNIFAC_B_UMRMC_CSV)?,
            parse_interaction_table(UNIFAC_C_UMRMC_CSV)?,
        ],
    })
}

/// The resolved UNIQUAC volume and surface parameters for a mixture.
///
/// A caller-supplied record the way [`Component`] is, and carrying **no names** for the
/// same reason: [`uniquac_parameters`] does the lookup.
#[derive(Debug, Clone, PartialEq)]
pub struct UniquacParameters {
    /// The van der Waals volume parameter `r_i` of each component.
    pub r: Vec<f64>,
    /// The van der Waals surface-area parameter `q_i` of each component.
    pub q: Vec<f64>,
}

/// The UNIQUAC `r` and `q` for a list of names, one entry per component.
///
/// `r_i = sum_k nu_ik R_k` and `q_i = sum_k nu_ik Q_k`, the group sums
/// [`crate::unifac_activity_coefficients`] forms internally - which is what NeqSim's
/// `ComponentGEUnifac.getR`/`getQ` compute, and what a UNIQUAC `r`/`q` means when no
/// fitted value exists.
///
/// **Not** NeqSim's `rUNIQUAQ`/`qUNIQUAQ` columns, which `ComponentGEUniquac` reads.
/// Those are `0.0` for 109 of the 112 components the table carries - only water, acetic
/// acid and `H2S` have values - so a UNIQUAC built from them divides by zero for almost
/// every real mixture. The spec's assumptions state this rather than leaving the choice
/// looking arbitrary.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name has no group decomposition.
pub fn uniquac_parameters(names: &[&str]) -> Result<UniquacParameters> {
    let params = unifac_parameters(names)?;
    let g = params.group_r.len();
    let n = names.len();
    let mut r = vec![0.0; n];
    let mut q = vec![0.0; n];
    for i in 0..n {
        for k in 0..g {
            r[i] += params.groups[i * g + k] * params.group_r[k];
            q[i] += params.groups[i * g + k] * params.group_q[k];
        }
    }
    Ok(UniquacParameters { r, q })
}

/// The Taleb acid identity of each component of a mixture.
///
/// A caller-supplied record the way [`Component`] is, and carrying **no names** for the
/// same reason: [`van_laar_acid_parameters`] does the lookup.
#[derive(Debug, Clone, PartialEq)]
pub struct VanLaarAcidParameters {
    /// Per component: `1` water, `2` nitric acid, `3` sulfuric acid, `0` for a species
    /// the model does not cover.
    pub acid_index: Vec<u8>,
}

/// The Taleb (1996) acid identity of each name, in the components' order.
///
/// `1` water, `2` nitric acid, `3` sulfuric acid, `0` for anything else. NeqSim's
/// `ComponentGEVanLaarAcid.acidIndexOf` recognises several spellings of each - the
/// formulae as well as the names - and so does this, because a caller writing `HNO3`
/// means the same substance as one writing `nitric acid`. Both are in the component
/// databank, so both resolve to the same entry.
///
/// A name the *databank* does not carry is refused. A name it does carry but that is
/// not one of the three acids resolves to `0`, which is not an error:
/// `eos.van_laar_acid_activity_coefficients` gives such a component its penalty, which
/// is what NeqSim does with a dissolved carrier gas.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if a name is in neither the databank nor the
///   overlay.
pub fn van_laar_acid_parameters(
    names: &[&str],
    overlay: Option<&Overlay>,
) -> Result<VanLaarAcidParameters> {
    let mut acid_index = Vec::with_capacity(names.len());
    for name in names {
        entry(name, overlay)?;
        acid_index.push(acid_index_of(name));
    }
    Ok(VanLaarAcidParameters { acid_index })
}

/// The Taleb acid index of a name: `1` water, `2` nitric acid, `3` sulfuric acid, else
/// `0`. The spellings NeqSim's `ComponentGEVanLaarAcid.acidIndexOf` accepts.
fn acid_index_of(name: &str) -> u8 {
    match name.trim().to_lowercase().as_str() {
        "water" | "h2o" => 1,
        "nitric acid" | "hno3" => 2,
        "sulfuric acid" | "sulphuric acid" | "h2so4" => 3,
        _ => 0,
    }
}

/// The parsed UNIFAC tables: group constants by subgroup, main-group interactions,
/// and per-component group memberships.
struct UnifacTables {
    /// Subgroup number -> `(R, Q, main group)`.
    group: HashMap<i64, (f64, f64, i64)>,
    /// `(main group, main group)` -> `a_mn`.
    aij: HashMap<(i64, i64), f64>,
    /// `(main group, main group)` -> `b_mn`, the linear term of UNIFAC-PSRK's
    /// temperature-dependent interaction.
    bij: HashMap<(i64, i64), f64>,
    /// `(main group, main group)` -> `c_mn`, the quadratic term of the same.
    cij: HashMap<(i64, i64), f64>,
    /// Lower-cased name -> `[(subgroup, count)]`.
    members: HashMap<String, Vec<(i64, i64)>>,
}

fn unifac_tables() -> &'static UnifacTables {
    static TABLES: OnceLock<UnifacTables> = OnceLock::new();
    TABLES.get_or_init(|| parse_unifac().expect("the embedded UNIFAC tables should parse"))
}

/// One integer field of a record, parsed strictly.
fn integer(record: &csv::StringRecord, index: usize, column: &str, row: usize) -> Result<i64> {
    let raw = record.get(index).unwrap_or("").trim();
    raw.parse::<i64>().map_err(|_| AzothError::InvalidInput {
        field: "databank".to_string(),
        reason: format!("row {row}: `{column}` is {raw:?}, which is not an integer"),
    })
}

fn parse_unifac() -> Result<UnifacTables> {
    let group = parse_group_constants(UNIFAC_GROUP_CSV)?;

    let aij = parse_interaction_table(UNIFAC_INTER_CSV)?;
    let bij = parse_interaction_table(UNIFAC_INTER_B_CSV)?;
    let cij = parse_interaction_table(UNIFAC_INTER_C_CSV)?;

    let members = parse_members(UNIFAC_COMP_CSV, 140)?;

    Ok(UnifacTables {
        group,
        aij,
        bij,
        cij,
        members,
    })
}

/// One group-constant table: subgroup number -> `(R, Q, main group)`.
fn parse_group_constants(text: &str) -> Result<HashMap<i64, (f64, f64, i64)>> {
    let mut out = HashMap::new();
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(text.as_bytes());
    let header = reader.headers().map_err(csv_failure)?.clone();
    let mut idx = HashMap::new();
    for name in ["secondary", "volumer", "surfareaq", "main"] {
        idx.insert(name, column(&header, name)?);
    }
    for (offset, record) in reader.records().enumerate() {
        let record = record.map_err(csv_failure)?;
        let row = offset + 2;
        out.insert(
            integer(&record, idx["secondary"], "secondary", row)?,
            (
                number(&record, idx["volumer"], "volumer", row)?,
                number(&record, idx["surfareaq"], "surfareaq", row)?,
                integer(&record, idx["main"], "main", row)?,
            ),
        );
    }
    Ok(out)
}

/// One per-component group decomposition, keyed by lower-cased name.
///
/// `subgroups` is how many `subN` columns to read: 140 for both of NeqSim's
/// decompositions, whose rows carry `sub1`..`sub140` whether or not the model's own
/// subclass loop stops earlier.
fn parse_members(text: &str, subgroups: usize) -> Result<HashMap<String, Vec<(i64, i64)>>> {
    let mut out = HashMap::new();
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(text.as_bytes());
    let header = reader.headers().map_err(csv_failure)?.clone();
    let name_idx = column(&header, "name")?;
    let mut sub_idx: HashMap<i64, usize> = HashMap::new();
    for s in 1..=subgroups {
        let key = i64::try_from(s).map_err(|_| AzothError::InvalidInput {
            field: "databank".to_string(),
            reason: format!("`{s}` subgroups is more than the decomposition can number"),
        })?;
        sub_idx.insert(key, column(&header, &format!("sub{s}"))?);
    }
    for (offset, record) in reader.records().enumerate() {
        let record = record.map_err(csv_failure)?;
        let name = record.get(name_idx).unwrap_or("").trim().to_lowercase();
        if name.is_empty() {
            continue;
        }
        let mut subs = Vec::new();
        for (s, &idx) in &sub_idx {
            let count = integer(&record, idx, &format!("sub{s}"), offset + 2)?;
            if count > 0 {
                subs.push((*s, count));
            }
        }
        out.insert(name, subs);
    }
    Ok(out)
}

/// One main-group interaction table, keyed by `(main group, main group)`.
///
/// The three of them - `a`, and UNIFAC-PSRK's `b` and `c` - have the same shape, so
/// they are read the same way rather than three times over.
fn parse_interaction_table(text: &str) -> Result<HashMap<(i64, i64), f64>> {
    let mut out = HashMap::new();
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(text.as_bytes());
    let header = reader.headers().map_err(csv_failure)?.clone();
    let main_idx = column(&header, "maingroup")?;
    let mut n_idx = HashMap::new();
    for n in 1..=64 {
        n_idx.insert(n, column(&header, &format!("n{n}"))?);
    }
    for (offset, record) in reader.records().enumerate() {
        let record = record.map_err(csv_failure)?;
        let row = offset + 2;
        let main = integer(&record, main_idx, "maingroup", row)?;
        for (n, &idx) in &n_idx {
            let raw = record.get(idx).unwrap_or("").trim();
            let value = if raw.is_empty() {
                0.0
            } else {
                raw.parse::<f64>().map_err(|_| AzothError::InvalidInput {
                    field: "databank".to_string(),
                    reason: format!("row {row}: `n{n}` is {raw:?}, which is not a number"),
                })?
            };
            out.insert((main, *n), value);
        }
    }
    Ok(out)
}

/// The UNIFAC inputs for a list of names, resolved from the vendored group tables.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `names` is empty.
/// * [`AzothError::PropertyUnavailable`] if a name has no UNIFAC group assignment.
pub fn unifac_parameters(names: &[&str]) -> Result<UnifacParameters> {
    let tables = unifac_tables();
    let basis = unifac_basis(names, &tables.group, &tables.members)?;
    let g = basis.union.len();

    let mut aij = vec![0.0; g * g];
    for (k, &s) in basis.union.iter().enumerate() {
        let (_, _, main) = *tables.group.get(&s).expect("a subgroup with a member row");
        for (m, &t) in basis.union.iter().enumerate() {
            let (_, _, main_t) = *tables.group.get(&t).expect("a subgroup with a member row");
            aij[k * g + m] = tables.aij.get(&(main, main_t)).copied().unwrap_or(0.0);
        }
    }

    Ok(UnifacParameters {
        groups: basis.groups,
        group_r: basis.group_r,
        group_q: basis.group_q,
        aij,
    })
}

/// The resolved UNIFAC-PSRK inputs for a mixture, each matrix flattened row-major.
///
/// The same basis as [`UnifacParameters`], with the interaction split into the three
/// terms UNIFAC-PSRK fits separately: `a_mn(T) = a_mn + b_mn T + c_mn T^2`.
#[derive(Debug, Clone, PartialEq)]
pub struct UnifacPsrkParameters {
    /// Per-component group counts, `N x G` row-major.
    pub groups: Vec<f64>,
    /// The volume `R` of each group, length `G`.
    pub group_r: Vec<f64>,
    /// The surface area `Q` of each group, length `G`.
    pub group_q: Vec<f64>,
    /// The constant term of the interaction, `G x G` row-major, in Kelvin.
    pub aij: Vec<f64>,
    /// The linear term, `G x G` row-major, in Kelvin per Kelvin.
    pub bij: Vec<f64>,
    /// The quadratic term, `G x G` row-major, in Kelvin per Kelvin squared.
    pub cij: Vec<f64>,
}

/// The UNIFAC-PSRK inputs for a list of names, resolved from the vendored tables.
///
/// `NeqSim`'s `ComponentGEUnifacPSRK.calcaij` is
/// `aij + bij * T + cij * T^2`, reading `UNIFACInterParamB` and `UNIFACInterParamC`
/// beside the `a` table the plain UNIFAC model uses. The three are resolved over the
/// same group union, so the matrices line up column for column.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `names` is empty.
/// * [`AzothError::PropertyUnavailable`] if a name has no UNIFAC group assignment.
pub fn unifac_psrk_parameters(names: &[&str]) -> Result<UnifacPsrkParameters> {
    let tables = unifac_tables();
    let basis = unifac_basis(names, &tables.group, &tables.members)?;
    let g = basis.union.len();

    let mut aij = vec![0.0; g * g];
    let mut bij = vec![0.0; g * g];
    let mut cij = vec![0.0; g * g];
    for (k, &s) in basis.union.iter().enumerate() {
        let (_, _, main) = *tables.group.get(&s).expect("a subgroup with a member row");
        for (m, &t) in basis.union.iter().enumerate() {
            let (_, _, main_t) = *tables.group.get(&t).expect("a subgroup with a member row");
            let key = (main, main_t);
            aij[k * g + m] = tables.aij.get(&key).copied().unwrap_or(0.0);
            bij[k * g + m] = tables.bij.get(&key).copied().unwrap_or(0.0);
            cij[k * g + m] = tables.cij.get(&key).copied().unwrap_or(0.0);
        }
    }

    Ok(UnifacPsrkParameters {
        groups: basis.groups,
        group_r: basis.group_r,
        group_q: basis.group_q,
        aij,
        bij,
        cij,
    })
}

/// The group basis both UNIFAC resolutions share: the per-component counts, the
/// per-group constants, and the union of subgroups the columns are ordered by.
///
/// Built once and handed to both models rather than called twice, so the two cannot
/// order their columns differently - a difference that would permute every interaction
/// matrix without changing anything a case could see.
pub(crate) struct UnifacBasis {
    /// Per-component group counts, `N x G` row-major.
    pub groups: Vec<f64>,
    /// The volume `R` of each group, length `G`.
    pub group_r: Vec<f64>,
    /// The surface area `Q` of each group, length `G`.
    pub group_q: Vec<f64>,
    /// The subgroups the columns are ordered by, ascending.
    pub union: Vec<i64>,
}

fn unifac_basis(
    names: &[&str],
    group: &HashMap<i64, (f64, f64, i64)>,
    members: &HashMap<String, Vec<(i64, i64)>>,
) -> Result<UnifacBasis> {
    let n = names.len();
    if n == 0 {
        return Err(AzothError::invalid_input(
            "components",
            "a mixture needs at least one component",
        ));
    }
    let mut union: Vec<i64> = Vec::new();
    let mut resolved: Vec<&Vec<(i64, i64)>> = Vec::with_capacity(n);
    for name in names {
        let key = name.trim().to_lowercase();
        let subs = members.get(&key).ok_or_else(|| {
            AzothError::property_unavailable(
                key,
                "UNIFAC group assignment".to_string(),
                "not in UNIFACcomp.csv; a UNIFAC activity coefficient needs a group \
                 decomposition for every component"
                    .to_string(),
            )
        })?;
        for &(s, _) in subs {
            if !union.contains(&s) {
                union.push(s);
            }
        }
        resolved.push(subs);
    }
    union.sort_unstable();
    let g = union.len();

    let mut group_r = vec![0.0; g];
    let mut group_q = vec![0.0; g];
    for (k, &s) in union.iter().enumerate() {
        let (r, q, _) = *group.get(&s).expect("a subgroup with a member row");
        group_r[k] = r;
        group_q[k] = q;
    }

    let mut groups = vec![0.0; n * g];
    for (i, subs) in resolved.iter().enumerate() {
        for &(s, count) in subs.iter() {
            let k = union
                .iter()
                .position(|&x| x == s)
                .expect("a member subgroup in the union");
            groups[i * g + k] = count as f64;
        }
    }
    Ok(UnifacBasis {
        groups,
        group_r,
        group_q,
        union,
    })
}

/// A mixture and its ideal-gas model, built from substance names, which come back
/// together because a mixture without heat-capacity coefficients cannot produce an
/// enthalpy.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `names` is empty.
/// * [`AzothError::PropertyUnavailable`] if a name is in neither source, or if one has
///   no heat-capacity coefficients - which is what an overlay-added substance is.
/// * Propagates [`Component::new`]'s range checks.
pub fn mixture_of(
    names: &[&str],
    cubic: Cubic,
    overlay: Option<&Overlay>,
) -> Result<(Mixture, IdealGasModel)> {
    if names.is_empty() {
        return Err(AzothError::invalid_input(
            "components",
            "a mixture needs at least one component",
        ));
    }
    let entries: Vec<Entry> = names
        .iter()
        .map(|name| entry(name, overlay))
        .collect::<Result<_>>()?;

    let ions: Vec<&str> = entries
        .iter()
        .filter(|e| e.class == ION)
        .map(|e| e.name.as_str())
        .collect();
    if !ions.is_empty() {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "{} {} an ion, and a cubic has no notion of one. NeqSim's table fills its \
                 ion rows' critical columns with a default - one shared set on 27 of the 62, \
                 a neutral parent's numbers on most of the rest - so a cubic built from them \
                 returns plausible-looking wrong numbers, and a card that adds an ion has no \
                 better constants to give. An ion belongs to the electrolyte models, which \
                 read its charge and diameter instead.",
                ions.join(", "),
                if ions.len() == 1 { "is" } else { "are" },
            ),
        ));
    }

    let missing: Vec<&str> = entries
        .iter()
        .filter(|e| e.cp.is_none())
        .map(|e| e.name.as_str())
        .collect();
    if !missing.is_empty() {
        return Err(AzothError::property_unavailable(
            missing.join(", "),
            "heat-capacity coefficients".to_string(),
            "the databank carries them for every substance it ships; one a keycard adds \
             needs its own, because a cubic needs `Tc`, `Pc` and `omega` and an enthalpy \
             needs the polynomial as well"
                .to_string(),
        ));
    }

    let components = entries
        .iter()
        .map(Entry::component)
        .collect::<Result<Vec<_>>>()?;

    // Flattened row-major, symmetric with a zero diagonal - the shape `Mixture::new`
    // validates. The diagonal is zero because a component does not interact with
    // itself, and neither source has a self-pair to look up.
    let n = entries.len();
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let value = kij(&entries[i].name, &entries[j].name, cubic, overlay);
            matrix[i * n + j] = value;
            matrix[j * n + i] = value;
        }
    }

    // Every entry has a polynomial: the filter above refuses the mixture otherwise.
    let coefficient = |index: usize| -> Vec<f64> {
        entries
            .iter()
            .map(|e| e.cp.map_or(0.0, |cp| cp[index]))
            .collect()
    };
    let ideal_gas = IdealGasModel {
        cp_a: coefficient(0),
        cp_b: coefficient(1),
        cp_c: coefficient(2),
        cp_d: coefficient(3),
        cp_e: coefficient(4),
    };

    // **The cubic is set here, not left to the caller.** It is not decoration: the `kij`
    // above was read from *this* cubic's column, so a mixture that then changed cubic
    // would carry the other family's interaction matrix. Setting it here is what makes
    // the pair impossible to mismatched - the resolver has the cubic and uses it.
    Ok((
        Mixture::new(components, matrix)?.with_cubic(cubic),
        ideal_gas,
    ))
}
