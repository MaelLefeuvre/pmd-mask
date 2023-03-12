use std::{fmt::{self, Display, Formatter}, str};

use serde::Deserialize;

use rust_htslib::bam::{HeaderView, Record};

/// Simple variant struct containing the raw string representation of a chromosome name.
/// ```
/// use genome::coordinate::ChrName;
/// 
/// let chr = ChrName::new("chrMT");
/// assert_eq!(chr.inner(), "chrMT");
/// ```
#[derive(Debug, Clone, Deserialize, PartialEq, Eq, Hash)]
pub struct ChrName(String);

use thiserror::Error;
#[derive(Debug, Error, PartialEq)]
pub enum ChrNameError {
    #[error("Failed to convert Htslib record into a valid ChrName")]
    ParseFromHtsLib(#[from] std::str::Utf8Error)
}

impl ChrName {
    pub fn new(name: &str) -> Self {
        ChrName(name.to_string())
    }

    pub fn from_htslib_record(header_view: &HeaderView, record: &Record) -> Result<Self, ChrNameError> {
        Ok(Self::new(str::from_utf8(header_view.tid2name(record.tid() as u32))?))
    }

    pub fn inner(&self) -> &str {
        &self.0
    }
}

impl Display for ChrName {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_test::Token;
    use serde_test::assert_de_tokens;

    #[test]
    fn display() {
        let chr = ChrName::new("chr1"); 
        assert_eq!("chr1", format!("{chr}"));
        assert_eq!("chr1-----", format!("{chr:-<9}"))
    }

    #[test]
    fn get_inner() {
        let chr = ChrName::new("chr1");
        assert_eq!(chr.inner(), "chr1")
    }

    #[test]
    fn deserialize() {
        assert_de_tokens(&ChrName::new("chr1"), &[
            Token::TupleStruct{name: "ChrName", len: 1},
            Token::String("chr1"),
            Token::TupleStructEnd
        ])
    }

    #[test]
    fn equality() {
        assert!(ChrName::new("chr1") == ChrName::new("chr1"));
        assert!(ChrName::new("chr1") != ChrName::new("chrY"));
    }
}
