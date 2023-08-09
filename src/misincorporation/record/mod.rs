use std::fmt::{self, Display, Formatter};
use crate::genome::{ChrName, Orientation, Strand, Position};

use serde::Deserialize;

/// A *partially* deserialized CSV record from a [mapDamage-v2](https://github.com/ginolhac/mapDamage)'s
/// [`misincorporation.txt`](https://ginolhac.github.io/mapDamage/#a4) output file.
/// This struct only keeps track of the most prevalent misincorporation patterns, i.e.:
/// - `C>T` transitions, at the [`Orientation::FivePrime`] end.
/// - `G>A` transitions, at the [`Orientation::ThreePrime`] end.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct MisincorporationRecord {
    #[serde(rename(deserialize = "Chr"))] pub chromosome: ChrName,
    #[serde(rename(deserialize = "End"))] pub end       : Orientation,
    #[serde(rename(deserialize = "Std"))] pub strand    : Strand,
    #[serde(rename(deserialize = "Pos"))] pub position  : Position, 
    #[serde(rename(deserialize = "C"))]   pub c_counts  : usize,
    #[serde(rename(deserialize = "G"))]   pub g_counts  : usize,
    #[serde(rename(deserialize = "C>T"))] pub c_to_t    : usize, 
    #[serde(rename(deserialize = "G>A"))] pub g_to_a    : usize,
}

impl MisincorporationRecord {
    /// Return the relative C>T frequency 
    /// This is computed as the number of observed `C>T`, divided by the number of observed `C`.
    pub(crate) fn c_to_t_freq(&self) -> f32 {
        self.c_to_t as f32 / self.c_counts as f32
    }

    /// Return the relative `G>A` frequency.
    /// This is computed as the number of observed `G>A`, divided by the number of observed `G`.
    pub(crate) fn g_to_a_freq(&self) -> f32 {
        self.g_to_a as f32 / self.g_counts as f32
    }

    /// Return the misincorporation frequency we're ***really*** interested in:  
    /// - If this entry is [`Orientation::FivePrime`]  -> return `C>T` relative frequency 
    ///   (see [`MisincorporationRecord::c_to_freq()`](MisincorporationRecord::c_to_t_freq))
    /// - If this entry is [`Orientation::ThreePrime`] -> return `G>A` relative frequency 
    ///   (see [`MisincorporationRecord::g_to_a_freq()`](MisincorporationRecord::g_t_a_freq))
    pub fn target_freq(&self) -> f32 {
        match self.end {
            Orientation::FivePrime  => self.c_to_t_freq(),
            Orientation::ThreePrime => self.g_to_a_freq(),
        }
    }
}

impl Display for MisincorporationRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        format!("Chr: {:<15} End: {: <3} Strand: {: <2} Pos: {: <4} C: {: <9} G: {: <9} C>T: {: <9} ({: <9.6}) G>A: {: <9} ({: <9.6})", 
            self.chromosome,
            self.end,
            self.strand,
            self.position,
            self.c_counts,
            self.g_counts,
            self.c_to_t,
            self.c_to_t_freq(),
            self.g_to_a,
            self.g_to_a_freq(),
        ).fmt(f)
    }
}

/// # Behavior:
/// Genate a dummy [`MisincorporationRecord`] from the provided parameters:
/// - `$chr`: A &str representation of a chromosome name 
/// - `$strand`: A [`Strand`] enum variant
/// - `$end`: An [`Orientation`] enum variant
/// - `$pos`: A u32 representing the base pair position
/// - `$counts`: The total number of nucleotide counts, found at `$pos`
/// - `$mis`: the misincorporation frequency found at `$pos` (either `C>T` or `G>A`, depending on the value of `$end`)
/// # Example
/// ```
/// use genome::{Strand, Orientation};
/// let record: MisincorporationRecord = mis_record!(
///     "chr1",
///     Strand::Reverse,
///     Orientation::FivePrime,
///     3,
///     10_000,
///     0.28
/// );
/// // Since we generated an Orientation::FivePrime record, we expect the g_to_a field to be set to 0, 
/// // and a c_to_t that is proportional to the provided $counts and $mis
/// assert_eq!(record.c_to_t, ((10_000.0/2.0) * 0.28).floor() as usize);
/// assert_eq!(record.g_to_a, 0);
/// ```
/// 
#[cfg(test)] 
#[macro_export]
macro_rules! mis_record {
        ($chr:expr, $strand:expr, $end:expr, $pos:expr, $counts:expr, $mis:expr) => {
            MisincorporationRecord {
                chromosome: $crate::genome::ChrName::new($chr),
                strand: $strand,
                end: $end,
                position: $crate::genome::Position::new($pos),
                c_counts: $counts / 2,
                c_to_t: if $end == $crate::genome::Orientation::FivePrime  { (($counts/2) as f64 * $mis).floor() as usize } else {0},
                g_counts: $counts / 2,
                g_to_a: if $end == $crate::genome::Orientation::ThreePrime { (($counts/2) as f64 * $mis).floor() as usize } else {0},

            }
        }
    }


#[cfg(test)]
mod test {
    use super::*;

    use serde_test::{Token, assert_de_tokens};

    /// # Behavior
    /// Construct an array of deserialization tokens ([`serde_test::Token`]) from the provided paramaters:
    /// - `$chr`: A &str representation of a chromosome name 
    /// - `$strand`: A [`Strand`] enum variant
    /// - `$end`: An [`Orientation`] enum variant
    /// - `$pos`: A u32 representing the base pair position
    /// - `$counts`: The total number of nucleotide counts, found at `$pos`
    /// - `$mis`: the misincorporation frequency found at `$pos` (either `C>T` or `G>A`, depending on the value of `$end`)
    /// 
    /// # Example
    /// ```
    /// use genome::{Strand, Orientation};
    /// let de_tokens: [Token; 28] = mis_tokens!("chr1", Strand::Reverse, Orientation::FivePrime,  1, 10_000, 0.5);
    /// assert_eq!(de_tokens.len(), 28);
    /// ```
    macro_rules! mis_tokens {
        ($chr:expr, $strand:expr, $end:expr, $pos:expr, $counts:expr, $mis:expr) => {[
            Token::Struct{name: "MisincorporationRecord", len: 8},
            Token::String("Chr"),
            Token::TupleStruct{name: "ChrName", len: 1},
            Token::String($chr),
            Token::TupleStructEnd,
            Token::String("End"),
            Token::UnitVariant { name: "Orientation", variant: $end.as_ref()},  
            Token::String("Std"),
            Token::UnitVariant { name: "Strand", variant: $strand.as_ref() },
            Token::String("Pos"),
            Token::TupleStruct{name: "Position", len: 1},
            Token::U32($pos),
            Token::TupleStructEnd,
            Token::String("C"),
            Token::U64($counts/2),
            Token::String("G"),
            Token::U64($counts/2),
            Token::String("C>T"),
            Token::U64(if $end == $crate::genome::Orientation::FivePrime  { (($counts/2) as f64 * $mis).floor() as u64 } else {0}),
            Token::String("G>A"),
            Token::U64(if $end == $crate::genome::Orientation::ThreePrime { (($counts/2) as f64 * $mis).floor() as u64 } else {0}),
            Token::StructEnd
            
        ]}
    }

    /// # Behavior
    /// Generate both a dummy [`MisincorporationRecord`] and an array of [`serde_test::Token`], and
    /// assert they both match using the [`serde_test::assert_de_tokens`] function.
    ///
    /// # Arguments:
    /// Same as those found in macros [`mis_record!`] and [`mis_tokens!`], which are :
    /// - `$chr`: A &str representation of a chromosome name 
    /// - `$strand`: A [`Strand`] enum variant
    /// - `$end`: An [`Orientation`] enum variant
    /// - `$pos`: A u32 representing the base pair position
    /// - `$counts`: The total number of nucleotide counts, found at `$pos`
    /// - `$mis`: the misincorporation frequency found at `$pos` (either `C>T` or `G>A`, depending on the value of `$end`)
    /// 
    /// # Example:
    /// ```
    /// use genome::{Strand, Orientation};
    /// assert_mis_de_tokens!("Y", Strand::Forward, Orientation::ThreePrime,  7, 124_816, 0.092);
    /// ```
    macro_rules! assert_mis_de_tokens {
        ($chr:expr, $strand:expr, $end:expr, $pos:expr, $counts:expr, $mis:expr) => {
            let record = mis_record!($chr, $strand, $end, $pos, $counts, $mis);
            let tokens = mis_tokens!($chr, $strand, $end, $pos, $counts, $mis);
            assert_de_tokens(&record, &tokens)
        }

    }

    #[test]
    fn deserialize() {
        assert_mis_de_tokens!("chr1", Strand::Forward, Orientation::FivePrime,  1, 10_000, 0.5);
        assert_mis_de_tokens!("chr1", Strand::Forward, Orientation::ThreePrime, 1, 10_000, 0.5);
        assert_mis_de_tokens!("chr1", Strand::Reverse, Orientation::FivePrime,  1, 10_000, 0.5);
        assert_mis_de_tokens!("chr1", Strand::Forward, Orientation::ThreePrime, 1, 10_000, 0.5);
        assert_mis_de_tokens!("chr1", Strand::Forward, Orientation::FivePrime,  1, 50_000, 0.0);
        assert_mis_de_tokens!("chr1", Strand::Forward, Orientation::FivePrime,  1, 50_000, 1.0);
        assert_mis_de_tokens!("MT",   Strand::Forward, Orientation::FivePrime, 10,   1000, 0.2);

    }

    #[test]
    fn g_to_a_freq() {
        let (counts, freq) = (2_000, 0.25);

        // ---- ThreePrime -> g_to_a_freq should be == 0.25 // c_to_a should be == 0 
        let record = mis_record!("X", Strand::Forward, Orientation::ThreePrime, 1, counts, freq);
        assert_eq!(record.g_to_a_freq(), freq as f32);

        // ---- FivePrime -> c_to_t_freq should be == 0.25 // g_to_a should be == 0 
        let record = mis_record!("X", Strand::Forward, Orientation::FivePrime, 1, counts, freq);
        assert_eq!(record.g_to_a_freq(), 0.0);

        // ---- Strand orientation should not affect this
        let record = mis_record!("X", Strand::Reverse, Orientation::ThreePrime, 1, counts, freq);
        assert_eq!(record.g_to_a_freq(), freq as f32);
        
        let record = mis_record!("X", Strand::Reverse, Orientation::FivePrime, 1, counts, freq);
        assert_eq!(record.g_to_a_freq(), 0.0);
    }

    #[test]
    fn c_to_t_freq() {
        let (counts, freq) = (2_000, 0.25);

        // ---- ThreePrime -> g_to_a_freq should be == 0.25 // c_to_a should be == 0 
        let record = mis_record!("X", Strand::Forward, Orientation::ThreePrime, 1, counts, freq);
        assert_eq!(record.c_to_t_freq(), 0.0);

        // ---- FivePrime -> c_to_t_freq should be == 0.25 // g_to_a should be == 0 
        let record = mis_record!("X", Strand::Forward, Orientation::FivePrime, 1, counts, freq);
        assert_eq!(record.c_to_t_freq(), freq as f32);

        // ---- Strand orientation should not affect this:
        let record = mis_record!("X", Strand::Reverse, Orientation::FivePrime, 1, counts, freq);
        assert_eq!(record.c_to_t_freq(), freq as f32);
        
        let record = mis_record!("X", Strand::Forward, Orientation::ThreePrime, 1, counts, freq);
        assert_eq!(record.c_to_t_freq(), 0.0);
    }

    #[test]
    fn target_freq() {
        let (counts, freq) = (2_000, 0.25);

        // ---- target_freq should return 0.25 no matter what, the only difference being :
        // - it returns G>A frequencies if the orientation is ThreePrime
        // - it returns C>T frequencies if the orientation is FivePrime
        let record = mis_record!("X", Strand::Forward, Orientation::ThreePrime, 1, counts, freq);
        assert_eq!(record.target_freq(), 0.25);
        assert_eq!(record.c_to_t, 0);
        assert_eq!(record.g_to_a, ((counts as f64/2.0) * freq).floor() as usize);

        let record = mis_record!("X", Strand::Forward, Orientation::FivePrime, 1, counts, freq);
        assert_eq!(record.target_freq(), 0.25);
        assert_eq!(record.g_to_a, 0);
        assert_eq!(record.c_to_t, ((counts as f64/2.0) * freq).floor() as usize);
    }

    #[test]
    fn display() {
        // ---- Ensure the function does not panic, or return an empty string, and that formatting is applied
        let record = mis_record!("X", Strand::Forward, Orientation::ThreePrime, 1, 654654, 0.12);
        assert!(!format!("{record}").is_empty());
        assert!(format!("{record:->1000}").starts_with("----------"));
    }
}
