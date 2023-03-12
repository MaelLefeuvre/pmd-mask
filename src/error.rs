use thiserror::Error;

#[derive(Debug, Error)]
pub enum RuntimeError {
    #[error(transparent)]
    HtsLibError(#[from] rust_htslib::errors::Error),

    #[error(transparent)]
    ParseMask(#[from] crate::mask::MaskEntryError)
}

