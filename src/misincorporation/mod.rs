use std::{fs::File, path::Path, ops::Deref, io::Read};
use crate::genome::{Strand, Orientation, ChrName};
use csv::ReaderBuilder;

mod error;
pub use error::MisincorporationsError;

mod record;
pub use record::MisincorporationRecord;


/// A collection of *partially* deserialized CSV record from a [mapDamage-v2](https://github.com/ginolhac/mapDamage)'s
/// [`misincorporation.txt`](https://ginolhac.github.io/mapDamage/#a4) output file. 
/// 
/// Each row within the `misincorporation.txt` file is encoded as a [`MisincorporationRecord`]
#[derive(Debug)]
pub struct Misincorporations{inner: Vec<MisincorporationRecord>}

impl Deref for Misincorporations {
    type Target = [MisincorporationRecord];

    fn deref(&self) -> &[MisincorporationRecord] {
        &self.inner
    }
}


impl FromIterator<MisincorporationRecord> for  Misincorporations {
    fn from_iter<T: IntoIterator<Item = MisincorporationRecord>>(records: T) -> Self {
        let inner = Vec::from_iter(records);
        Self{inner}
    }
}

impl Misincorporations {
    /// Generate a [`Misincorporations`] struct from a [mapDamage-v2](https://github.com/ginolhac/mapDamage)'s
    /// [`misincorporation.txt`](https://ginolhac.github.io/mapDamage/#a4) output file.
    /// 
    /// # Usage
    /// ```
    /// use pmd_mask::misincorporation::Misincorporations;
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let file = "tests/test-data/bam/dummy-MTonly/misincorporation.txt";
    ///     let misincorporations = Misincorporations::from_path(&file, 0.01)?;
    ///     Ok(())
    /// }
    /// ```
    pub fn from_path(path: impl AsRef<Path>, threshold: f32) -> Result<Self, MisincorporationsError>{
        let file = File::open(&path)
            .map_err(|e| MisincorporationsError::OpenFile(path.as_ref().display().to_string(), e))?;
        Self::from_reader(file, threshold)
    }

    /// Private [`Misincorporations`] struct constructor from a generic Reader.
    /// See [`Misincorporations::from_path`](Misincorporations::from_path) for the public implementation
    pub(crate) fn from_reader<R: Read>(path: R, threshold: f32) -> Result<Self, MisincorporationsError>{
        use MisincorporationsError::*;

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .comment(Some(b'#'))
            .from_reader(path)
            ; 

        let mut skip_chromosome: Option<(&ChrName, &Orientation, &Strand)> = None; 
        let mut threshold_positions = Vec::with_capacity(32 * 2 * 2);

        'nextline: for (line, result) in reader.deserialize::<MisincorporationRecord>().enumerate() {
            let record = result.map_err(|e|DeserializeRecord(line, e.to_string()))?;
            // Once we've found a position at which the threshold is met, we 
            // don't need to parse entries which belong from the same chromosome.
            if let Some(chromosome_to_skip) = skip_chromosome {
                if (&record.chromosome, &record.end, &record.strand) == chromosome_to_skip { continue 'nextline }
            }

            // If we're below the requested treshold, keep that record!
            if record.target_freq() <= threshold {
                threshold_positions.push(record);
                let last_insert = threshold_positions.last().unwrap(); // We can unwrap here since we know we've just pushed a value
                skip_chromosome = Some((&last_insert.chromosome, &last_insert.end, &last_insert.strand)); 
            
            }
        }
        Ok(Self{ inner: threshold_positions }) 
    }

    /// Extrude invalid frequencies from the inner collection of [`MisincorporationRecord`] and return them
    /// into an owned [`Vec`].
    /// 
    /// Invalid [`MisincorporationRecord`]s are those whose [`target_freq()`](`MisincorporationRecord::target_freq`) is either:
    /// - a `NaN` value (see [`f32::is_nan()`](f32))
    /// - a `Inf` value (see [`f32::is_infinite()`](f32))
    /// - a negative float value (see [`f32::is_sign_negative()`](f32))
    /// 
    /// # Usage
    /// ```
    /// use pmd_mask::misincorporation::Misincorporations;
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let file = "tests/test-data/bam/dummy-MTonly/misincorporation.txt";
    ///     let mut misincorporations = Misincorporations::from_path(&file, 0.01)?;
    /// 
    ///     let invalid_freqs = misincorporations.extrude_invalid_frequencies();
    /// 
    ///     assert!(invalid_freqs.is_empty());
    ///     Ok(())
    /// }
    /// ```
    pub fn extrude_invalid_frequencies(&mut self) -> Vec<MisincorporationRecord> {
        let mut invalid_positions = Vec::with_capacity(self.inner.len());

        self.inner.retain(|pos| {
            if pos.target_freq().is_nan() || pos.target_freq().is_infinite() || pos.target_freq().is_sign_negative()  {
                invalid_positions.push(pos.clone());
                false
            } else {
                true
            }
        });

        invalid_positions
    }
}

#[cfg(test)]
macro_rules! mock_file {
    (header) => {
        "Chr	End	Std	Pos	A	C	G	T	Total	G>A	C>T	A>G	T>C	A>C	A>T	C>G	C>A	T>G	T>A	G>C	G>T	A>-	T>-	C>-	G>-	->A	->T	->C	->G	S\n"
    };

    (header, $chr:expr, $end:expr, $std:expr, $pos:expr, $n:expr, $freq:expr) => {
        format!("{}{}", mock_file!(header), mock_file!($chr, $end, $std, $pos, $n, $freq))
    };

    ($chr:expr, $end:expr, $std:expr, $pos:expr, $n:expr, $freq:expr) => {
        format!("{0}	{1}	{2}	{3}	{4}	{4}	{4}	{4}	{5}	{6}	{7}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}	{8}\n", 
            $chr, $end, $std, $pos, 
            $n/4, $n,
            if $end == "3p" { (($n as f64 /4.0) * $freq).floor() as usize } else {0},
            if $end == "5p" { (($n as f64 /4.0) * $freq).floor() as usize } else {0},
            0
        )
    };
}


#[cfg(test)]
mod test {
    use std::{ops::RangeInclusive, io::Cursor};

    use crate::genome::Position;

    use super::*;

    /// # Summary 
    /// Mock a misincorporation.txt file, given:
    /// - `chr_range`  : a range of chromosomes to simulate
    /// - `n`          : the number of read base pairs to consider
    /// - `count`      : the overall number of nucleotides counter @ each position
    /// - `start_freq` : the starting misincorporation frequency, @ position 1
    /// - `decay`      : the rate of misincorporation decay, in percent / position shift.
    /// - `threshold`  : the masking threshold
    /// # Behavior
    /// - for `5p` ends, `G>A` is counted as `(count / 4) * start_freq`, the rest (including `C>T`) is set to `0`  
    /// - for `3p` ends, `C>T` is counted as `(count / 4) * start_freq`, the rest (including `G>A`) is set to `0`  
    fn mock_misincorporation(chr_range: RangeInclusive<usize>, n: usize, count: usize, start_freq: f64, decay: f64, threshold: f32) -> Result<Misincorporations, MisincorporationsError> {
        let mut out = String::new();
        out += mock_file!(header);

        for chr in chr_range.map(|c| format!("chr{c}") ) {
            for end in ["3p", "5p"] {
                for std in ["+", "-"] {
                    let mut freq = start_freq;
                    for pos in 1..=n {
                        out += &mock_file!(chr, end, std, pos, count, freq);
                        freq *= decay;
                    }
                }
            }
        }

        Misincorporations::from_reader(Cursor::new(out), threshold)
    }

    #[test]
    fn test_threshold_long_decay() -> Result<(), MisincorporationsError> {
        let (start_mis, mis_decay, threshold) = (0.3, 0.88, 0.01);
        let misincorporations = mock_misincorporation(1..=22, 70, 100_000, start_mis, mis_decay, threshold)?;
        
        
        // ---- expected_pos = ln(threshold / start_mis) / ln(mis_decay) + 1 
        // Note the '+1', since we're messing between indices and bp positions. Careful with off by one errors..
        let expected = ( f64::ln(threshold as f64 / start_mis) / f64::ln(mis_decay) ).ceil() as usize + 1;
        println!("Expected: {expected}");

        assert!(!misincorporations.is_empty());
        for record in misincorporations.iter() {
            println!("got: {} | want: {}", record.position, expected);
            assert_eq!(record.position, Position::new(expected))
        }
        Ok(())
    }

    #[test]
    fn test_threshold_no_decay() -> Result<(), MisincorporationsError> {
        let (start_mis, mis_decay, threshold) = (0.3, 1.00, 0.01);
        let misincorporations = mock_misincorporation(1..=22, 70, 100_000, start_mis, mis_decay, threshold)?;

        // No decay, meaning we never attain the threshold
        // => Thus, we expect an empty list of Misincorporation Records.
        //    Which implies masking will apply on the whole sequence.
        assert!(misincorporations.is_empty());
        Ok(())
    }

    #[test]
    fn test_threshold_high_threshold() -> Result<(), MisincorporationsError> {
        let (start_mis, mis_decay, threshold) = (0.3, 0.88, 0.5);
        let misincorporations = mock_misincorporation(1..=22, 70, 100_000, start_mis, mis_decay, threshold)?;

        // Threshold is higher than the starting misincorporation, meaning we attain the threshold from the start
        // => Thus, we expect a full list of misincorporations
        //    Which implies masking will never apply

        let expected = ( f64::ln(threshold as f64 / start_mis) / f64::ln(mis_decay) ).ceil() as usize + 1;
        assert_eq!(expected, 1);
        println!("Expected: {expected}");

        assert!(!misincorporations.is_empty());
        for record in misincorporations.iter() {
            println!("got: {} | want: {}", record.position, expected);
            assert_eq!(record.position, Position::new(expected))
        }
        Ok(())
    }

    macro_rules! tamper_misincorporations {
        ($mis:expr, $field:ident, $pos:expr) => {{
            let funky_record = $mis.inner.get_mut($pos).expect("Empty misincorporations");
            funky_record.$field = 0;
            funky_record.clone()
        }
        }
    }

    #[test]
    fn extrude_invalid_frequencies() -> Result<(), MisincorporationsError> {
        let (start_mis, mis_decay, threshold) = (0.3, 0.88, 0.5);
        let mut misincorporations = mock_misincorporation(1..=1, 10, 100_000, start_mis, mis_decay, threshold)?;

        // Subsequent test rely on the assumption that the first misincorporation entry is ThreePrime
        assert_eq!(misincorporations.inner[0].end, Orientation::ThreePrime);
        assert_eq!(misincorporations.inner.len(), 4);

        // All frequencies should be valid at first.
        assert_eq!(misincorporations.extrude_invalid_frequencies().len(), 0);

        // ---- Let's tamper with the data...
        // Setting g_to_a to zero should not be considered invalid, as this makes sense if the data is simply not damaged.
        let _ = tamper_misincorporations!(misincorporations, g_to_a, 0);
        assert_eq!(misincorporations.extrude_invalid_frequencies().len(), 0);

        // Setting c_to_t and/or c_counts to zero should not be considered invalid: we're only interested in the validity
        // of the **target_freq(). i.e: g_to_a_freq() if orientation is ThreePrime. 
        let _            = tamper_misincorporations!(misincorporations, c_counts, 0);
        let funky_record = tamper_misincorporations!(misincorporations, c_to_t, 0);
        assert!(funky_record.c_to_t_freq().is_nan());                         // frequency is now invalid... 
        assert_eq!(misincorporations.extrude_invalid_frequencies().len(), 0); // ...but we don't care, since we're ThreePrime.

        // ---- Let's tamper with the data...
        // Setting g_counts to zero should return NaN in this case.
        //This is invalid. and extrude_invalid_frequencies should remove and return it.
        let funky_record = tamper_misincorporations!(misincorporations, g_counts, 0);
        println!("Funky record: {funky_record}");
        assert!(funky_record.g_to_a_freq().is_nan());
        assert_eq!(misincorporations.extrude_invalid_frequencies()[0], funky_record);
        assert_eq!(misincorporations.inner.len(), 3);

        Ok(())
    }

}