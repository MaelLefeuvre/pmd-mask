use thiserror::Error;

#[derive(Debug, Error, PartialEq)]
pub enum MaskThresholdError {
    #[error("Invalid number of masking threshold entries: expected {want}, found {got}")]
    ValidateThresh{got: usize, want: usize}
}
