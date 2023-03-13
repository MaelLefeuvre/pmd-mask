use std::{fmt::{self, Display, Formatter}, str::FromStr};

use serde::Deserialize;
use rust_htslib::bam::Record;

mod error;
pub use error::StrandError;

/// Strand orientation representation for paired-end sequencing reads.
#[derive(Debug, Clone, Copy, Deserialize, PartialEq, Eq, Hash)] 
pub enum Strand {
    #[serde(rename = "+")] Forward,
    #[serde(rename = "-")] Reverse
}


impl Strand {

    pub fn from_htslib_record(record: &mut Record) -> Result<Self, StrandError> {
        record.strand()
            .strand_symbol()
            .parse::<Strand>()
            .map_err(|e| StrandError::ParseFromHtsLib(Box::new(e)))
    }
}


impl AsRef<str> for Strand {
    fn as_ref(&self) -> &str {
        match self {
            Self::Forward => "+",
            Self::Reverse => "-",
        }
    }
}


impl Display for Strand {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl FromStr for Strand {
    type Err = StrandError;

    /// Attempt to convert a string sequence into a Strand representation.
    /// Valid values are either "+" or "-". Anything else will result in an error
    /// 
    /// # Errors
    /// 
    /// will return a [`StrandError::ParseStrand`] error upon encountering any character that is neither '+' nor '-'
    /// ```
    /// # use std::str::FromStr;
    /// # use pmd_mask::genome::{Strand, StrandError};
    /// 
    /// let forward_strand = "+".parse::<Strand>();
    /// assert_eq!(forward_strand, Ok(Strand::Forward));
    /// 
    /// let reverse_strand = Strand::from_str("-");
    /// assert_eq!(reverse_strand, Ok(Strand::Reverse));
    /// 
    /// let strange_strand = "x".parse::<Strand>();
    /// assert!(strange_strand.is_err())
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+"   => Ok(Self::Forward),
            "-"   => Ok(Self::Reverse),
            other => Err(Self::Err::ParseStrand(other.to_string()))

        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use serde_test::{Token, assert_de_tokens};

    #[test]
    fn display() {
        let forward = Strand::Forward;
        let reverse = Strand::Reverse;

        assert_eq!("+ ATCG\n- tagc", format!("{forward} ATCG\n{reverse} tagc"));
        assert_eq!("+    ATCG\n-    tagc", format!("{forward: <5}ATCG\n{reverse: <5}tagc"));
    }

    #[test]
    fn deserialize_forward() {
        assert_de_tokens(&Strand::Forward, &[
            Token::UnitVariant { name: "Strand", variant: "+" }
        ])
    }

    #[test]
    fn deserialise_reverse() {
        assert_de_tokens(&Strand::Reverse, &[
            Token::UnitVariant { name: "Strand", variant: "-" }
        ])
    }

    #[test]
    fn as_ref() {
        assert_eq!(Strand::Forward.as_ref(), "+");
        assert_eq!(Strand::Reverse.as_ref(), "-");
    }

    #[test]
    fn from_str() {
        assert_eq!(Strand::from_str("+"), Ok(Strand::Forward));
        assert_eq!(Strand::from_str("-"), Ok(Strand::Reverse));
        assert_eq!(Strand::from_str("*"), Err(StrandError::ParseStrand("*".to_string())));
    }

}