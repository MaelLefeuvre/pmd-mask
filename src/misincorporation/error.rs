use thiserror::Error;
#[derive(Debug, Error)]
pub enum MisincorporationsError {
    #[error("@line {0}: Failed to deserialize record in misincorporation file. Got {1}")]
    DeserializeRecord(usize, String),

    #[error("Failed to open {0}. Got {1}")]
    OpenFile(String, #[source] std::io::Error)
}
