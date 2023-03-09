use std::{fmt::{self, Display, Formatter}};
use serde::Deserialize;

/// Position within a chromosome or read in base pair
#[derive(Debug, Clone, Copy, Deserialize, Hash, PartialEq, Eq)]
pub struct Position(pub usize);

impl Display for Position {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}
