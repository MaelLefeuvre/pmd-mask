use thiserror::Error;

#[derive(Debug, Error)]
pub enum StrandError {
    #[error("Failed to parse string value '{0}' into a valid Strand representation")]
    ParseStrand(String)
}