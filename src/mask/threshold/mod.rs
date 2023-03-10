use std::{fmt::{self, Display, Formatter}, collections::HashMap};

use crate::genome::{Orientation, Position};

mod error;
pub use error::MaskThresholdError;

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

    pub fn validate(&self) -> Result<(), MaskThresholdError> {
        let expected_len = 2;
        if self.inner.len() != expected_len {
            return Err(MaskThresholdError::ValidateThresh{got: self.inner.len(), want: expected_len})
        }
        Ok(())
    }

    pub fn get_threshold(&self, orientation: &Orientation) -> Option<&Position> {
        self.inner.get(orientation)
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
