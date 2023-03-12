use std::{fmt::{self, Display, Formatter}};
use serde::Deserialize;

/// Position within a chromosome or read in base pair
#[derive(Debug, Clone, Copy, Deserialize, Hash, PartialEq, Eq)]
pub struct Position(usize);

impl Display for Position {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl Position {
    pub fn new(value: usize) -> Self {
        Self(value)
    }

    pub fn inner(&self) -> usize {
        self.0
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

