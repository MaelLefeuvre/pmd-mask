use thiserror::Error;

#[derive(Debug, Error)]
pub enum MaskEntryError {
    #[error("Failed to parse HtslibRecord into a valid Masking entry. Got {0}")]
    ParseFromHtslib(#[source] Box<dyn std::error::Error + Send + Sync> )
}