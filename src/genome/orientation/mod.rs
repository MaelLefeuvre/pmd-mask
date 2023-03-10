use std::{fmt::{self, Display, Formatter}};
use serde::Deserialize;

#[derive(Debug, Clone, Copy, Deserialize, PartialEq, Eq, Hash)]
pub enum Orientation {
    #[serde(rename = "3p")]
    ThreePrime,
    #[serde(rename = "5p")]
    FivePrime,
}

impl Display for Orientation {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let repr = match self {
            Self::ThreePrime => "3p",
            Self::FivePrime  => "5p",
        };
        repr.fmt(f)
    }
}
