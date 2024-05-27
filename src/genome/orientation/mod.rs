use std::fmt::{self, Display, Formatter};
use serde::Deserialize;

/// Encodes the orientation of a relative bam record position. Two possible variants:  
/// - [`Orientation::ThreePrime`]|`'3p'`: 3'OH end of a fragment
/// - [`Orientation::FivePrime`]|`'5p'`: 5'P end of a fragment
/// 
/// Largely associated with [`crate::misincorporation::MisincorporationRecord`], 
/// to encode the relative position of a nucleotide within a read.
#[derive(Debug, Clone, Copy, Deserialize, PartialEq, Eq, Hash)]
pub enum Orientation {
    #[serde(rename = "3p")]
    ThreePrime,
    #[serde(rename = "5p")]
    FivePrime,
}

impl AsRef<str> for Orientation {
    /// Obtain the [`str`] representation of an [`Orientation`]
    /// ```
    /// use pmd_mask::genome::Orientation;
    /// 
    /// assert_eq!(Orientation::ThreePrime.as_ref(), "3p");
    /// assert_eq!(Orientation::FivePrime.as_ref(), "5p");
    /// 
    /// ```
    fn as_ref(&self) -> &str {
        match self {
            Self::ThreePrime => "3p",
            Self::FivePrime  => "5p",
        }
    }
}

impl Display for Orientation {
    /// Obtain a formatted  [`String`] representation of an [`Orientation`].
    /// 
    /// ```
    /// # use pmd_mask::genome::Orientation;
    /// assert_eq!(&format!("{: ^6}", Orientation::ThreePrime), "  3p  ");
    /// assert_eq!(&format!("{: <4}", Orientation::FivePrime), "5p  ");
    /// ```
    /// 
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use serde_test::{Token, assert_de_tokens};

    #[test]
    fn display() {
        let start = Orientation::FivePrime;
        let end   = Orientation::ThreePrime;

        assert_eq!("5p_ATCG_3p", format!("{start}_ATCG_{end}"));
        assert_eq!("5p____ATCG____3p", format!("{start:_<6}ATCG{end:_>6}"))
    }

    #[test]
    fn deserialize_5p() {
        assert_de_tokens(&Orientation::FivePrime, &[
            Token::UnitVariant { name: "Orientation", variant: "5p" }
        ])
    }


    #[test]
    fn deserialize_3p() {
        assert_de_tokens(&Orientation::ThreePrime, &[
            Token::UnitVariant { name: "Orientation", variant: "3p" }
        ])
    }
}
