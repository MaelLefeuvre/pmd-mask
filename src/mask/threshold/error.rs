use thiserror::Error;

/// Error type enum for [`crate::mask::MaskThreshold`]
#[derive(Debug, Error, PartialEq)]
pub enum MaskThresholdError {
    #[error("Invalid number of masking threshold entries: expected {want}, found {got}")]
    ValidateThresh{got: usize, want: usize}
}
