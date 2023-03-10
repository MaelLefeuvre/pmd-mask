use std::{fmt::{self, Display, Formatter}, str::FromStr};
use serde::Deserialize;

mod error;
pub use error::StrandError;

/// Strand orientation representation for paired-end sequencing reads.
#[derive(Debug, Clone, Copy, Deserialize, PartialEq, Eq, Hash)] 
pub enum Strand {
    #[serde(rename = "+")] Forward,
    #[serde(rename = "-")] Reverse
}

impl Strand {
    /// Return a string representation of this struct: `+` for reverse strands, `-` for forward strands.
    pub fn symbol(&self) -> &str {
        match self {
            Self::Forward => "+",
            Self::Reverse => "-",
        }
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.symbol().fmt(f)
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
    /// let forward_strand = "+".parse::<Strand>();
    /// assert_eq!(forward_strand, Ok(Strand::Forward));
    /// let reverse_strand = Strand::from_str("-");
    /// assert_eq!(reverse_strand, Ok(Strand::Reverse));
    /// let strange_strand = "x".parse::<Strand>();
    /// assert_eq!(strange_strand, StrandError::ParseStrand)
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+"   => Ok(Self::Forward),
            "-"   => Ok(Self::Reverse),
            other => Err(Self::Err::ParseStrand(other.to_string()))

        }
    }
}
