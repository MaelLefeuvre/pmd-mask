use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;


use log::{warn, debug, trace};

pub mod entry;
pub use entry::MaskEntry;
pub use entry::MaskEntryError;

#[cfg(test)] pub use entry::dummy;

mod threshold;
pub use threshold::MaskThreshold;

mod error;
pub use error::MasksError;

use crate::misincorporation::Misincorporations;
use crate::genome::Orientation;


/// A [`HashMap`] collection of [`MaskThreshold`]s, mapped according to their respective [`MaskEntry`].
///
/// [`crate::mask::Masks`] contains all the entries contained within a [mapDamage-v2](https://github.com/ginolhac/mapDamage)'s
/// [`misincorporation.txt`](https://ginolhac.github.io/mapDamage/#a4) output file. These are structured through an internal
/// [`HashMap`] where:  
/// - keys are [`MaskEntry`] (themselves, containing the chromosome and Strand information of the entry)
/// - values are [`MaskThreshold`]s (themselves, containing the relative threshold positions for the 5p and 3p end of a read.)
/// 
#[derive(Debug)]
pub struct Masks {inner: HashMap<MaskEntry, MaskThreshold>}


impl TryFrom<&Misincorporations> for Masks {

    type Error = MasksError;

    /// Generate a structured [`Masks`] collection from a [`Misincorporations`] reference.
    /// 
    /// # Usage
    /// 
    /// ```
    /// use pmd_mask::misincorporation::Misincorporations;
    /// use pmd_mask::mask::Masks;
    /// 
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let reader = std::fs::File::open("tests/test-data/bam/dummy-MTonly/misincorporation.txt")?;
    ///     let misincorporations = Misincorporations::from_reader(reader, 0.01)?;
    ///     let masks = Masks::try_from(&misincorporations)?;
    ///     Ok(())
    /// }
    /// ```
    fn try_from(value: &Misincorporations) -> Result<Self, Self::Error> {
        // ---- Restructure Misincorporation records as a HashMap<MaskEntry, MaskThreshold>
        let mut masks = Self{ inner: HashMap::with_capacity(value.len()) };
        for position in value.iter() {
            let record = MaskEntry{chromosome: position.chromosome.clone(), strand: position.strand};
            masks.inner.entry(record)
                .or_insert(MaskThreshold::default())
                .set_threshold(position.end, position.position);
        }

        masks.validate()?;
        Ok(masks)
    }
}


impl Masks {

    /// Instantiate a [`Masks`] struct from a [mapDamage-v2](https://github.com/ginolhac/mapDamage)'s
    /// [`misincorporation.txt`](https://ginolhac.github.io/mapDamage/#a4) output file and a set, user-defined 
    /// threshold value.
    /// 
    /// # Usage
    /// ```
    /// use pmd_mask::mask::Masks;
    /// use std::error::Error;
    /// fn main() -> Result<(), Box<dyn Error>> {
    ///     let file  = "tests/test-data/bam/dummy-MTonly/misincorporation.txt";
    ///     let masks = Masks::from_path(&file, 0.01)?;
    ///     Ok(())
    /// }
    /// ```
    /// 
    /// # Errors
    /// 
    /// May return:  
    /// - a [`MasksError::OpenFile`] if the method failed to open the provided `misincorporations` file.
    /// - any [`MasksError`] spat out from the private [`Masks::from_reader()`](Masks::from_reader) function
    /// 
    pub fn from_path(misincorporations: impl AsRef<Path>, threshold: f32) -> Result<Self, MasksError> {
        let file = File::open(&misincorporations)
            .map_err(|e| MasksError::OpenFile{source: e})?;
        Self::from_reader(file, threshold)
    }

    /// Instantiate a [`Masks`] struct from a generic Reader and a set threshold. Used by [`Masks::from_path()`](Masks::from_path)
    /// 
    /// # Errors
    /// - May bubble out any errors arising from [`Misincorporations::from_reader()`](Misincorporations::from_reader)
    /// - May bubble out any errors arising from [`Masks::try_from::<&Misincorporations>()`]
    fn from_reader<R: std::io::Read>(misincorporations: R, threshold: f32) -> Result<Self, MasksError> {

        let mut threshold_positions = Misincorporations::from_reader(misincorporations, threshold)?;

        // ---- Validate misincorporation and issue warnings for any 'abnormal' frequency found.
        let mut abnormal_frequencies = threshold_positions
            .extrude_invalid_frequencies()
            .iter()
            .fold(String::new(), |abnormal_freqs, pos| abnormal_freqs + &format!("{pos}\n"));
        abnormal_frequencies.pop();

        if !abnormal_frequencies.is_empty() {
            warn!(
            "Found abnormal misincorporation frequencies (NaN, Infinite values, etc.). \
            Masking will apply along the full length of the sequence for those:\n{abnormal_frequencies}"
            );
        }

        debug!("Using the following positions as threshold for masking:\n{}",
            threshold_positions.iter().fold(String::new(), |acc, val| {
                acc + &format!("{val}\n")
            })
        );

        Masks::try_from(&threshold_positions)
    }

    /// Return the [`MaskThreshold`] of a provided [`MaskEntry`].
    /// 
    /// # Usage
    /// ```
    /// use std::error::Error;
    /// use pmd_mask::mask::{Masks, MaskEntry};
    /// use pmd_mask::genome::{ChrName, Strand};
    /// fn main() -> Result<(), Box<dyn Error>> {
    ///     let file  = "tests/test-data/bam/dummy-MTonly/misincorporation.txt";
    ///     let masks = Masks::from_path(&file, 0.01)?;
    ///     let entry = MaskEntry{ chromosome: ChrName::new("5"), strand: Strand::Reverse };
    /// 
    ///     let threshold = masks.get(&entry);
    ///     assert!(threshold.is_some());
    ///     Ok(())
    /// }
    /// ```
    pub fn get(&self, entry: &MaskEntry) -> Option<&MaskThreshold> {
        self.inner.get(entry)
    }

    /// Iteratively validate every contained [`MaskThreshold`]s found within the [`Masks`] struct, 
    /// i.e. each record must contain a [`HashMap`] w/ two keys: [`Orientation::FivePrime`] and [`Orientation::ThreePrime`]
    /// 
    /// # Errors
    /// May return a [`MasksError::ValidateThresholds`] if any threshold is invalid (i.e. [`MaskThreshold::validate()`] fails)
    fn validate(&self) -> Result<(), MasksError> {
        for (entry, threshold) in self.inner.iter() {
            trace!("Validating {entry} {threshold}");
            threshold.validate().map_err(|e|MasksError::ValidateThresholds(entry.clone(), e))?
        };
        Ok(())
    }

    /// Serialize the contents of this struct within a writer.
    /// Output is headed and Fields are '<Chr> <Std> <5p> <3p>'
    /// ```
    /// use std::error::Error;
    /// use pmd_mask::mask::{Masks, MaskEntry};
    /// use pmd_mask::genome::{ChrName, Strand};
    /// fn main() -> Result<(), Box<dyn Error>> {
    ///     let file  = "tests/test-data/bam/dummy-MTonly/misincorporation.txt";
    ///     let masks = Masks::from_path(&file, 0.01)?;
    ///
    ///     // ---- Write contents into cursor (or any output file)
    ///     let mut output = std::io::Cursor::new(Vec::new());
    ///     masks.write(&mut output)?; 
    ///     
    ///     // ---- Ensure utf8 validity
    ///     let raw_output = String::from_utf8(output.into_inner())?;
    ///     Ok(())
    /// }
    /// ```
    pub fn write(&self, writer: &mut impl Write) -> std::io::Result<()> {
        use threshold::ORIENTATIONS;

        // ---- Write header
        let header = format!("Chr\tStd\t{}\t{}\n", ORIENTATIONS[0], ORIENTATIONS[1]);
        writer.write_all(header.as_bytes())?;

        // ---- Sort entries according to Entry
        let mut sorted = self.inner.iter().collect::<Vec<_>>();
        sorted.sort_by(|a, b| a.0.cmp(b.0));

        // ---- Write entries in order.
        for (entry, threshold) in sorted.iter() {
            let mut output_string = format!("{}\t{}", entry.chromosome, entry.strand);
            for end in [Orientation::FivePrime, Orientation::ThreePrime] {
                let position = match threshold.get_threshold(&end).expect("Missing threshold").inner() {
                    usize::MAX => "\tNA".to_string(),
                    position   => format!("\t{position}"),
                };
                output_string.push_str(&position);
            }
            output_string.push('\n');
            writer.write_all(output_string.as_bytes())?;
        }
        writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::mask::threshold::MaskThresholdError;
    use crate::misincorporation::MisincorporationRecord;
    use crate::genome::{Strand::*, Orientation::*, ChrName, Position};


    use crate::mis_record;
    fn dummy_misincorporations() -> Misincorporations {
        Misincorporations::from_iter(vec![
            mis_record!("chr1", Forward, FivePrime,  10,  500, 0.01),
            mis_record!("chr1", Forward, ThreePrime, 12,  800, 0.01),
            mis_record!("chr1", Reverse, FivePrime,  8, 1000, 0.01),
            mis_record!("chr1", Reverse, ThreePrime, 6,  600, 0.01),

        ].into_iter())
    }

    #[test]
    pub fn try_from_misincorporations() {
        let misincorporations = dummy_misincorporations();
        let mask              = Masks::try_from(&misincorporations).expect("Invalid Masks");

        for i in 0..misincorporations.len() {
            let target_record = misincorporations.get(i).expect("Missing misincorporation record");
            let expected_mask = MaskEntry{chromosome: target_record.chromosome.clone(), strand: target_record.strand};
            let got           = &mask.inner[&expected_mask];

            assert_eq!(got.get_threshold(&target_record.end).expect("Missing Mask threshold"), &target_record.position);
        }
    }

    #[test]
    fn get_threshold() {
        let mut masks = Masks{inner: HashMap::new()};

        let entry = MaskEntry{chromosome: ChrName::new("chr1"), strand: Forward}; 
        masks.inner.insert(entry.clone(), MaskThreshold::default());

        assert_eq!(masks.get(&entry), Some(&MaskThreshold::default()));

        masks.inner.get_mut(&entry).expect("Missing entry").set_threshold(FivePrime, Position::new(10));
        assert_ne!(masks.get(&entry), Some(&MaskThreshold::default()));
    }

    #[test]
    fn validate() {
        let mut masks = Masks::try_from(&dummy_misincorporations()).expect("Invalid Misincorporations");
        match masks.validate() {
            Ok(()) => {},
            _ => panic!("Invalid mask")
        };

        let bad_entry       = MaskEntry{chromosome: ChrName::new("chr1"), strand: Reverse};
        let funky_threshold = masks.inner.get_mut(&bad_entry).expect("Entry not found");
        *funky_threshold    = MaskThreshold{inner: HashMap::new()};
        
        match masks.validate() {
            Ok(()) => panic!("Masks should be invalid at this point."),
            Err(e) => match e {
                MasksError::ValidateThresholds(entry, treshold_err) => {
                    assert_eq!(bad_entry, entry);
                    assert_eq!(treshold_err, MaskThresholdError::ValidateThresh { got: 0, want: 2 })
                },
                _ => panic!("Invalid error type.")
            }
        }
    }

}