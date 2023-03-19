use thiserror::Error;

#[derive(Debug, Error)]
pub enum RuntimeError {
    #[error(transparent)]
    HtsLibError(#[from] rust_htslib::errors::Error),

    #[error(transparent)]
    ParseMask(#[from] crate::mask::MaskEntryError),

    #[error(transparent)]
    ParseMisincorporation(#[from] crate::mask::MasksError),
    

    #[error("Both the output and input alignment files appears to be the same file! Exiting.")]
    InputIsOutput,

    #[error("Failed to open the requested threshold metrics file. [{0}]")]
    OpenMetrics(#[source] std::io::Error),
    
    #[error("Failed to write masking thresholds within the provided metrics file path. [{0}]")]
    WriteMasksMetrics(#[source] std::io::Error),

}

