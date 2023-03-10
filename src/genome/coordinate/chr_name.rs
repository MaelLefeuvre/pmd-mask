use std::{fmt::{self, Display, Formatter}};

use serde::Deserialize;

/// Contains the raw string representation of a chromosome name.
#[derive(Debug, Clone, Deserialize, PartialEq, Eq, Hash)]
pub struct ChrName(pub String);

impl Display for ChrName {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}
