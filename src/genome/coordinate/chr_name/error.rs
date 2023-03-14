use thiserror::Error;

#[derive(Debug, Error, PartialEq)]
pub enum ChrNameError {
    #[error("Failed to convert Htslib record into a valid ChrName")]
    ParseFromHtsLib(#[from] std::str::Utf8Error)
}