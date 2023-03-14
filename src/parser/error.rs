
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CliError {
    #[error("Accepted values: 'Sam|Bam|Cram'")]
    InvalidBamOutputFmt,

    #[error("The provided value must either be 0, or a non negative integer. Got {0}")]
    InvalidThreadValue(#[source] std::num::ParseIntError)
}

