use std::{path::Path, ops::Deref};

use csv::ReaderBuilder;

mod error;
pub use error::MisincorporationsError;

mod record;
pub use record::MisincorporationRecord;

use crate::genome::{Strand, Orientation, ChrName};


pub struct Misincorporations{inner: Vec<MisincorporationRecord>}

impl Deref for Misincorporations {
    type Target = [MisincorporationRecord];

    fn deref(&self) -> &[MisincorporationRecord] {
        &self.inner
    }
}

impl Misincorporations {
    pub fn from_path(path: impl AsRef<Path>, threshold: f32) -> Result<Self, MisincorporationsError>{
        use MisincorporationsError::*;

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .comment(Some(b'#'))
            .from_path(&path)
            .map_err(|e| OpenMisincorporationFile(path.as_ref().display().to_string(), e))?; 

        let mut skip_chromosome: Option<(&ChrName, &Orientation, &Strand)> = None; 
        let mut threshold_positions = Vec::with_capacity(32 * 2 * 2);

        'nextline: for (line, result) in reader.deserialize::<MisincorporationRecord>().enumerate() {
            let record = result.map_err(|e|DeserializeRecord(line, e.into()))?;

            // Once we've found a position at which the threshold is met, we 
            // don't need to parse entries which belong from the same chromosome.
            if let Some(chromosome_to_skip) = skip_chromosome {
                if (&record.chromosome, &record.end, &record.strand) == chromosome_to_skip {
                    continue 'nextline
                }
            }

            // If we're below the requested treshold, keep that record!
            if record.target_freq() <= threshold {

                threshold_positions.push(record);
                let last_insert = threshold_positions.last().unwrap(); // We can unwrap here since we know we've just pushed a value
                skip_chromosome = Some((&last_insert.chromosome, &last_insert.end, &last_insert.strand)); 
            
            }
        }
        Ok(Self{ inner: threshold_positions }) 
    }

    pub fn extrude_invalid_frequencies(&mut self) -> Vec<MisincorporationRecord> {
        let mut invalid_positions = Vec::with_capacity(self.inner.len());

        self.inner.retain(|pos| {
            if pos.target_freq().is_normal() {
                true
            } else {
                invalid_positions.push(pos.clone());
                false
            }
        });

        invalid_positions
    
    }

}