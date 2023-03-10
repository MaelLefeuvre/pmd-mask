use thiserror::Error;
#[derive(Debug, Error)]
pub enum MisincorporationsError {
    #[error("@line {0}: Failed to deserialize record in misincorporation file. Got {1}")]
    DeserializeRecord(usize, Box<dyn std::error::Error>),

    #[error("Failed to open {0}. Got {1}")]
    OpenMisincorporationFile(String, #[source] csv::Error)
}
