use thiserror::Error;

#[derive(Debug, Error)]
pub enum RuntimeError {
    HtsLibError(#[from] rust_htslib::errors::Error)
}