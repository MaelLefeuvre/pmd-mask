use std::collections::HashMap;

use crate::misincorporation::Misincorporations;

use log::{trace};

mod entry;
pub use entry::MaskEntry;

mod threshold;
pub use threshold::MaskThreshold;

mod error;
pub use error::MasksError;


/// A Hash collection of Masking thresholds, mapped according to their respective MaskEntries
/// (i.e. Chromosome name and strand information.)
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