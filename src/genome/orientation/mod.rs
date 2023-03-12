use std::{fmt::{self, Display, Formatter}};
use serde::Deserialize;

#[derive(Debug, Clone, Copy, Deserialize, PartialEq, Eq, Hash)]
pub enum Orientation {
    #[serde(rename = "3p")]
    ThreePrime,
    #[serde(rename = "5p")]
    FivePrime,
}

impl AsRef<str> for Orientation {
    fn as_ref(&self) -> &str {
        match self {
            Self::ThreePrime => "3p",
            Self::FivePrime  => "5p",
        }
    }
}

impl Display for Orientation {
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
