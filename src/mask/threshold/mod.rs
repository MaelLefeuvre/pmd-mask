use std::{fmt::{self, Display, Formatter}, collections::HashMap};

use crate::genome::{Orientation, Position};

mod error;
pub use error::MaskThresholdError;

/// Default orientations which MUST be contained within a MaskThreshold.
pub(crate) static ORIENTATIONS: [Orientation; 2] = [Orientation::FivePrime, Orientation::ThreePrime];

/// Basic Building  block of a [`crate::mask::Masks`] struct. [`MaskThreshold`] contains 
/// the orientation and relative position of a [mapDamage-v2](https://github.com/ginolhac/mapDamage)'s
/// [`misincorporation.txt`](https://ginolhac.github.io/mapDamage/#a4) output file record.
/// 
/// Relative positions are encoded as [`crate::genome::Position`], while strand orientations are encoded
/// as [`crate::genome::Orientation`].
/// 
/// Internally, [`MaskThreshold`] is a [`HashMap`](std::collections::HashMap) containing two keys: 
/// [`Orientation::ThreePrime`] and [`Orientation::FivePrime`]. [`MaskThreshold`] will **always** 
/// start initialized with these two keys, and **cannot** contain any additional keys.
/// 
/// Relative [`Position`]s are on the other hand set as values.
/// 
/// [`MaskThreshold`]s are internally used by [`crate::mask::Masks`] structs as values within an internal
/// [`HashMap`](`std::collections::HashMap`)
/// 
#[derive(Debug, PartialEq)]
pub struct MaskThreshold { pub(crate) inner: HashMap<Orientation, Position> }

impl Default for MaskThreshold {

    /// Create a new, default [`MaskThreshold`], containing two [`HashMap`] keys: 
    /// [`Orientation::ThreePrime`] and [`Orientation::FivePrime`].
    /// Each [`Position`] value within the internal [`HashMap`] is set to [`usize::MAX`](`usize`)
    /// 
    /// # Usage
    /// ```
    /// use pmd_mask::mask::MaskThreshold;
    /// use pmd_mask::genome::{Orientation::*, Position};
    /// 
    /// let threshold = MaskThreshold::default();
    ///
    /// // Note threshold contains two Orientation, with positions set to usize::MAX
    /// let default_pos = Position::new(usize::MAX);
    /// assert_eq!(threshold.get_threshold(&ThreePrime), Some(&default_pos));
    /// assert_eq!(threshold.get_threshold(&FivePrime), Some(&default_pos));
    /// ```
    fn default() -> Self {
        let max_pos = Position::new(usize::MAX);
        let inner   = ORIENTATIONS.iter().map(|orientation| (*orientation, max_pos) );
        Self { inner: HashMap::from_iter(inner) }
    }

}


impl MaskThreshold {

    /// Set the threshold of a given `orientation` to a given `position`
    /// 
    /// # Usage
    /// ```
    /// use pmd_mask::mask::MaskThreshold;
    /// use pmd_mask::genome::{Orientation, Position};
    /// 
    /// let mut threshold = MaskThreshold::default();
    ///
    /// let new_position = Position::new(42);
    /// threshold.set_threshold(Orientation::ThreePrime, new_position);
    /// 
    /// # assert_eq!(threshold.get_threshold(&Orientation::ThreePrime), Some(&new_position))
    /// ```
    pub fn set_threshold(&mut self, orientation: Orientation, position: Position) {
        self.inner.insert(orientation, position);
    }

    /// Validate the integrity of a [`MaskThreshold`] i.e.: 
    /// - contains two keys: [`Orientation::FivePrime`] and [`Orientation::ThreePrime`]
    /// 
    /// # Errors
    /// Returns a [`MaskThresholdError::ValidateThresh`] if the number of keys contained within 
    /// the internal [`HashMap`] does not match the length of [`ORIENTATIONS`].
    /// 
    /// # @TODO: 
    /// - obvious code smells. set this validation as either a compile time check, or change the internal structure
    ///   to make the addition of keys an impossibility.
    pub fn validate(&self) -> Result<(), MaskThresholdError> {
        let expected_len = ORIENTATIONS.len();
        if self.inner.len() != expected_len {
            return Err(MaskThresholdError::ValidateThresh{got: self.inner.len(), want: expected_len})
        }
        Ok(())
    }

    /// Retrieve the threshold's [`Position`] for a given strand [`Orientation`].
    /// 
    /// # Usage 
    /// 
    /// ```
    /// use pmd_mask::mask::MaskThreshold;
    /// use pmd_mask::genome::Orientation;
    /// 
    /// let threshold = MaskThreshold::default();
    /// 
    /// let default_thresh = threshold.get_threshold(&Orientation::FivePrime);
    /// 
    /// # use pmd_mask::genome::Position;
    /// # assert_eq!(default_thresh, Some(&Position::new(usize::MAX)))
    /// ```
    pub fn get_threshold(&self, orientation: &Orientation) -> Option<&Position> {
        self.inner.get(orientation)
    }
}

impl Display for MaskThreshold {
    /// Return a formatted [`String`] representation of a [`MaskThreshold`] struct.
    /// 
    /// # Usage
    /// ```
    /// use pmd_mask::mask::MaskThreshold;
    /// let thresh = MaskThreshold::default();
    ///
    /// let thresh_str = format!("{thresh:_^100}");
    /// 
    /// # assert!(thresh_str.starts_with("_"));
    /// # assert!(thresh_str.ends_with("_"));
    /// ```
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
