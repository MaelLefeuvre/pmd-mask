use std::collections::HashMap;

use crate::misincorporation::Misincorporations;

use log::{trace};

mod entry;
pub use entry::MaskEntry;
pub use entry::MaskEntryError;

mod threshold;
pub use threshold::MaskThreshold;

mod error;
pub use error::MasksError;


/// A Hash collection of Masking thresholds, mapped according to their respective MaskEntries
/// (i.e. Chromosome name and strand information.)
#[derive(Debug)]
pub struct Masks {inner: HashMap<MaskEntry, MaskThreshold>}


impl TryFrom<&Misincorporations> for Masks {

    type Error = MasksError;

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
    pub fn get(&self, entry: &MaskEntry) -> Option<&MaskThreshold> {
        self.inner.get(entry)
    }

    // ---- Validate thresholds: Each record must contain a hashmap w/ two values (one for 3p, the other for 5p)
    pub fn validate(&self) -> Result<(), MasksError> {
        for (entry, threshold) in self.inner.iter() {
            trace!("Validating {entry:?} {threshold:?}");
            threshold.validate().map_err(|e|MasksError::ValidateThresholds(entry.clone(), e))?
        };
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
        assert_eq!(masks.validate(), Ok(()));

        let bad_entry       = MaskEntry{chromosome: ChrName::new("chr1"), strand: Reverse};
        let funky_threshold = masks.inner.get_mut(&bad_entry).expect("Entry not found");
        *funky_threshold    = MaskThreshold{inner: HashMap::new()};
        
        assert_eq!(masks.validate(), Err(MasksError::ValidateThresholds(bad_entry, MaskThresholdError::ValidateThresh{got:0, want:2})))
    }

}