use std::{fmt::{self, Display, Formatter}, collections::HashMap};

use crate::genome::{Orientation, Position};

#[derive(Debug)]
pub struct MaskThreshold { inner: HashMap<Orientation, Position> }

impl Default for MaskThreshold {
    fn default() -> Self {
        use Orientation::*;
        let max_pos = Position(usize::MAX);
        Self { inner: HashMap::from_iter([(FivePrime, max_pos), (ThreePrime, max_pos)].into_iter())}
    }

}

impl MaskThreshold {
    
    pub fn set_threshold(&mut self, orientation: Orientation, position: Position) {
        self.inner.insert(orientation, position);
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.inner.len() != 2 {
            return Err("Invalid Threshold".to_string())
        }
        Ok(())
    }

    pub fn get_threshold(&self, orientation: &Orientation) -> &Position {
        self.inner.get(orientation).unwrap()
    }
}

impl Display for MaskThreshold {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut out = self.inner.iter().fold(String::new(), |acc, (orient, threshold)| {
            acc + &format!("({orient}: {threshold}bp) ")
        });
        out.pop();
        out.fmt(f)
    }
}
