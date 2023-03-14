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
    InputIsOutput

}

