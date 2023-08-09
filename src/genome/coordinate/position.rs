use std::{fmt::{self, Display, Formatter}};
use serde::Deserialize;

/// Absolute or relative position within a chromosome or read (in base pairs).
/// 
/// This is just a struct containing a [`usize`].
#[derive(Debug, Clone, Copy, Deserialize, Hash, PartialEq, Eq)]
pub struct Position(usize);

impl Display for Position {
    /// Return a formatted [`String`] representation of a [`Position`].
    /// 
    /// # Example 
    /// ```
    /// use pmd_mask::genome::Position;
    /// let position: Position = 42.into();
    /// assert_eq!(&format!("{position: ^6}"), "  42  ");
    /// ```
    /// 
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl Position {
    /// Create a new [`Position`] from a usize
    /// 
    /// # Example
    /// ```
    /// use pmd_mask::genome::Position;
    /// 
    /// let position = Position::new(42usize);
    /// assert_eq!(position, 42.into());
    /// ```
    /// 
    pub fn new(value: usize) -> Self {
        Self(value)
    }

    //// Get the inner `usize` contained within a [`Position`] back
    /// 
    /// # Example 
    /// ```
    /// use pmd_mask::genome::Position;
    /// 
    /// assert_eq!(Position::from(42).inner(), 42);
    /// ```
    /// 
    pub fn inner(&self) -> usize {
        self.0
    }
}

impl From<usize> for Position {
    /// Convert a [`usize`] from and into a [`Position`]
    /// ```
    /// use pmd_mask::genome::Position;
    /// 
    /// let pos1: Position = 55.into();
    /// let pos2: Position = Position::from(55);
    /// assert_eq!(pos1, pos2);
    /// ```
    fn from(value: usize) -> Self {
        Position(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_test::Token;
    use serde_test::assert_de_tokens;

    #[test]
    fn display() {
        let chr = Position::new(10_000); 
        assert_eq!("10000", format!("{chr}"));
        assert_eq!("10000----", format!("{chr:-<9}"))
    }

    #[test]
    fn get_inner() {
        let chr = Position::new(15_654);
        assert_eq!(chr.inner(), 15_654)
    }

    #[test]
    fn deserialize() {
        assert_de_tokens(&Position::new(120_000_000), &[
            Token::TupleStruct{name: "Position", len: 1},
            Token::U32(120_000_000),
            Token::TupleStructEnd
        ])
    }

    #[test]
    fn equality() {
        assert!(Position::new(10) == Position::new(10));
        assert!(Position::new(10) != Position::new(11));
    }
}

