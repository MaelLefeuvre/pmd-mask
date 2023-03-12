use thiserror::Error;

use super::{entry::MaskEntry, threshold::MaskThresholdError};

#[derive(Debug, Error, PartialEq)]
pub enum MasksError {
    #[error("Invalid MaskThreshold within MaskEntry {0}. Got {1}")]
    ValidateThresholds(MaskEntry, #[source] MaskThresholdError)
}
