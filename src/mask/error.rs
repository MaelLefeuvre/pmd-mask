use thiserror::Error;

use crate::misincorporation::MisincorporationsError;

use super::{entry::MaskEntry, threshold::MaskThresholdError};

#[derive(Debug, Error)]
pub enum MasksError {
    #[error("Invalid MaskThreshold within MaskEntry {0}. Got {1}")]
    ValidateThresholds(MaskEntry, #[source] MaskThresholdError),

    #[error("Failed to open misincorporation file [{source}]")]
    OpenFile{#[source] source: std::io::Error},

    #[error("Failed to obtain Misincorporations from misincorporation file")]
    GenerateMisincorporations(#[from] MisincorporationsError)

}
