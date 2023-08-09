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

    #[error("Length of the retrieved reference sequence does not match the length of the read")]
    ReferenceOutOfIndexError,

    #[error(
    "Neither '--bam', nor the standard input is being sollicited at this point. \
    Please provide pmd-mask with an input to work from, either through piping, or with the --bam argument"
    )]
    NoStdin,

    #[error(
    "Failed to load fasta index into memory. \
    Ensure the provided fasta isn't corrupted, and is properly indexed with a companion '.fai' \
    file within the same directory."
    )]
    LoadFaidx,
}

