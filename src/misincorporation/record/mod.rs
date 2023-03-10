use std::{fmt::{self, Display, Formatter}};
use crate::genome::{ChrName, Orientation, Strand, Position};

use serde::Deserialize;


#[derive(Debug, Clone, Deserialize)]
pub struct MisincorporationRecord {
    #[serde(rename(deserialize = "Chr"))] pub chromosome: ChrName,
    #[serde(rename(deserialize = "End"))] pub end       : Orientation,
    #[serde(rename(deserialize = "Std"))] pub strand    : Strand,
    #[serde(rename(deserialize = "Pos"))] pub position  : Position, 
    #[serde(rename(deserialize = "C"))]   pub c_counts  : usize,
    #[serde(rename(deserialize = "G"))]   pub g_counts  : usize,
    #[serde(rename(deserialize = "C>T"))] pub c_to_t    : usize, 
    #[serde(rename(deserialize = "G>A"))] pub g_to_a    : usize,
}

impl MisincorporationRecord {
    /// Return the relative C>T frequency 
    /// This is computed as the number of observed C>T, divided by the number of observed C.
    fn c_to_t_freq(&self) -> f32 {
        self.c_to_t as f32 / self.c_counts as f32
    }

    /// Return the relative G>A frequency.
    /// This is computed as the number of observed G>A, divided by the number of observed G.
    fn g_to_a_freq(&self) -> f32 {
        self.g_to_a as f32 / self.g_counts as f32
    }

    /// Return the frequency we're really interested in:
    /// If this entry is a 5p -> return C>T relative frequency (see [c_to_freq()])
    /// If this entry is a 3p -> return G>A relative frequency (see [g_to_a_freq()])
    pub fn target_freq(&self) -> f32 {
        match self.end {
            Orientation::FivePrime  => self.c_to_t_freq(),
            Orientation::ThreePrime => self.g_to_a_freq(),
        }
    }
}

impl Display for MisincorporationRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Chr: {:<15} End: {: <3} Strand: {: <2} Pos: {: <4} C: {: <9} G: {: <9} C>T: {: <9} ({: <9.6}) G>A: {: <9} ({: <9.6})", 
            self.chromosome,
            self.end,
            self.strand,
            self.position,
            self.c_counts,
            self.g_counts,
            self.c_to_t,
            self.c_to_t_freq(),
            self.g_to_a,
            self.g_to_a_freq(),
        )
    }
}
