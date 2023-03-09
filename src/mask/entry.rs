use std::fmt::{self, Display, Formatter};

use crate::genome::{ChrName, Strand};

#[derive(Debug, Hash, PartialEq, Eq)]
pub struct MaskEntry {
    pub chromosome: ChrName,
    pub strand    : Strand
}

impl Display for MaskEntry {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.chromosome, self.strand)
    }
}
