use std::{fmt::{self, Display, Formatter}, str};

use serde::Deserialize;
use rust_htslib::bam::{HeaderView, Record};

mod error;
pub use error::ChrNameError;

/// Simple variant struct containing the raw string representation of a chromosome name.
/// 
/// # Example
/// ```
/// use pmd_mask::genome::ChrName;
/// 
/// let chr = ChrName::new("chrMT");
/// assert_eq!(chr.inner(), "chrMT");
/// ```
#[derive(Debug, Clone, Deserialize, PartialEq, Eq, Hash)]
pub struct ChrName(String);

impl ChrName {
    /// Instantiate a new ChrName from a raw string slice.
    /// 
    /// # Example
    /// ```
    /// use pmd_mask::genome::ChrName;
    /// let chr = ChrName::new("chrMT");
    /// ```
    pub fn new(name: &str) -> Self {
        ChrName(name.to_string())
    }

    /// Extract and instantiate a new ChrName from a [`rust_htslib::bam::Record`] and a [`rust_htslib::bam::HeaderView`]
    /// 
    /// # Example
    /// 
    /// ```no_run
    /// use rust_htslib::bam::{self, Read};
    /// use pmd_mask::genome::ChrName;
    /// 
    /// // Setup
    /// let mut reader = bam::Reader::from_stdin().expect("Failed to read BAM from stdin");
    /// let mut record = bam::Record::new();
    /// 
    /// // Read a single record...
    /// reader.read(&mut record);
    /// 
    /// // ... And easily extract the chromosome's name.
    /// let chr = ChrName::from_htslib_record(&reader.header(), &record).expect("Failed to retrieve this record's ChrName.");
    /// assert_eq!(chr.inner(), "chr1")
    /// 
    /// ```
    pub fn from_htslib_record(header_view: &HeaderView, record: &Record) -> Result<Self, ChrNameError> {
        Ok(Self::new(str::from_utf8(header_view.tid2name(record.tid() as u32))?))
    }

    /// Get the inner string representation back from this struct.
    /// # Example
    /// ```
    /// use pmd_mask::genome::ChrName;
    /// 
    /// let chr = ChrName::new("chrMT");
    /// assert_eq!(chr.inner(), "chrMT");
    /// ```
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
    use crate::{mask::dummy::dummy_bam, genome::Strand};

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

    #[test]
    fn from_htslib_record() {
        use rust_htslib::bam::{Record, Header, HeaderView};
    
        let (header, record): (Header, Record) = dummy_bam(Strand::Forward, 10_000);
        let chr = ChrName::from_htslib_record(&HeaderView::from_header(&header), &record).unwrap();
        
        assert_eq!(chr.inner(), "chr1");

    }
}
