use std::{fmt::{self, Display, Formatter}, collections::HashMap};

use crate::genome::{Orientation, Position};

mod error;
pub use error::MaskThresholdError;

// Default orientations which MUST be contained within a MaskThreshold.
pub(crate) static ORIENTATIONS: [Orientation; 2] = [Orientation::FivePrime, Orientation::ThreePrime];

#[derive(Debug, PartialEq)]
pub struct MaskThreshold { pub(crate) inner: HashMap<Orientation, Position> }

impl Default for MaskThreshold {
    fn default() -> Self {
        let max_pos = Position::new(usize::MAX);
        let inner   = ORIENTATIONS.iter().map(|orientation| (*orientation, max_pos) );
        Self { inner: HashMap::from_iter(inner) }
    }

}


impl MaskThreshold {

    pub fn set_threshold(&mut self, orientation: Orientation, position: Position) {
        self.inner.insert(orientation, position);
    }

    pub fn validate(&self) -> Result<(), MaskThresholdError> {
        let expected_len = ORIENTATIONS.len();
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
        // ORIENTATIONS is used to ensure the key order remains the same (since our underlying data is a HashMap)
        //@SAFETY: We can expect the underlying hashmap has all the keys specified in order, since 
        //         it is also used during construction 
        let mut out = ORIENTATIONS.iter().fold(String::new(), |acc, orientation| {
            let threshold = self.inner.get(orientation).expect("Missing or invalid orientation value");
            acc + &format!("({orientation}: {threshold}bp) ")
        });
        out.pop();
        out.fmt(f)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn display() {
        let threshold = MaskThreshold::default();
        let want = format!("({}: {}bp) ({}: {}bp)", Orientation::FivePrime, usize::MAX, Orientation::ThreePrime, usize::MAX);

        assert_eq!(want, format!("{threshold}"));
        assert_eq!(format!("{want:-^80}"), format!("{threshold:-^80}"));
    }

    #[test]
    fn get_set_threshold() {
        let mut threshold = MaskThreshold::default();
        let target_orientation = ORIENTATIONS[0];
        let default_position   = Position::new(usize::MAX);

        // Ensure the target orientation has a default position
        assert_eq!(threshold.get_threshold(&target_orientation), Some(&default_position));

        // Update target_orientation.
        let new_position = Position::new(5);
        threshold.set_threshold(target_orientation, new_position);
        assert_eq!(threshold.get_threshold(&target_orientation), Some(&new_position));

        // Ensure other positions remain intact.
        for other_orientation in &ORIENTATIONS[1..] {
            assert_eq!(threshold.get_threshold(other_orientation), Some(&default_position));
        }
    }

    #[test]
    fn validate() {
        let default_threshold = MaskThreshold::default();
        assert_eq!(default_threshold.validate(), Ok(()));

        let mut funky_threshold = MaskThreshold::default();
        funky_threshold.inner.remove(&ORIENTATIONS[0]);
        assert_eq!(funky_threshold.validate(), Err(MaskThresholdError::ValidateThresh{got: 1, want: ORIENTATIONS.len()}));
    }
}
