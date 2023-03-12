use thiserror::Error;

#[derive(Debug, Error, PartialEq)]
pub enum StrandError {
    #[error("Failed to parse string value '{0}' into a valid Strand representation")]
    ParseStrand(String),

    #[error("Failed to construct a valid Strand from the provided Htslib Record")]
    ParseFromHtsLib(#[source] Box<Self>)
}