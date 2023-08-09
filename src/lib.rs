use std::{str, ops::Range};


pub mod misincorporation;
pub mod genome;
pub mod mask;
pub mod error;

use error::RuntimeError;
use genome::Orientation;
pub use mask::{Masks, MaskEntry, MaskThreshold};

use anyhow::{Result, Context};
use rust_htslib::{faidx, bam};
use log::{debug, trace, warn};


use rust_htslib::bam::ext::BamRecordExtensions;


/// Generic internal function intended to apply selective masking on either end of a raw `&mut [u8]` read.
/// Mainly used within [`mask_5p`] and [`mask_3p`].
/// 
/// # Parameters
/// - `range`: half-bounded [`Range`] of indices where masking should be applied within `seq`.
/// - `reference`: raw byte representation of reference sequence corresponding to the sequence to mask. 
///    Note that this sequence must contain and take into account any indel found within `seq`
/// - `seq`: raw byte representation of the sequence to mask. Note that this sequence must preserve any indel found within it.
/// - `quals`: raw byte representation of the phred-scores of `seq`.
/// - `positions`: matching positions found between `seq` and `ref`. These can be retrieved using 
///    [`rust_htslib::bam::Reader::aligned_pairs()`](`rust_htslib::bam::Reader`)

/// # Errors
/// - May emit a [`RuntimeError::ReferenceOutOfIndexError`] if the function ever fails to retrieve a reference nucleotide.
#[inline]
fn mask_sequence(range: Range<usize>, reference: &[u8], seq: &mut [u8], quals: &mut [u8], target_nucleotide: u8, positions: &[[usize; 2]]) -> Result<(), RuntimeError> {
    'mask: for [readpos, refpos] in positions.iter().copied().skip(range.start).take_while(|[readpos, _]| *readpos < range.end ) {
        if readpos >= seq.len() { break 'mask }
        let reference_nucleotide = reference.get(refpos).ok_or_else(|| RuntimeError::ReferenceOutOfIndexError)?;
        if *reference_nucleotide == target_nucleotide {
            seq[readpos] = b'N';
            quals[readpos] = 0;
        }
    }
    Ok(())
}

/// Apply selective masking from the [`Orientation::FivePrime`] end of a read.
/// 
/// # Parameters:  
/// - `threshold`: reference to a mask [`MaskThreshold`]. This is where the relative [`Position`](`crate::genome::Position`)
///   of the [`ThreePrime`](`Orientation::ThreePrime`) end is retrieved.
/// - `reference`: raw byte representation of reference sequence corresponding to the sequence to mask. 
///    Note that this sequence must contain and take into account any indel found within `seq`
/// - `seq`: raw byte representation of the read sequence to mask. Note that this sequence must preserve any indel found within it.
/// - `quals`: raw byte representation of the phred-scores of `seq`.
/// - `positions`: matching positions found between `seq` and `ref`. These can be retrieved using 
///    [`rust_htslib::bam::Reader::aligned_pairs()`](`rust_htslib::bam::Reader`)
/// 
/// # Errors: 
///  - May emit a [`RuntimeError::ReferenceOutOfIndexError`] (see [`mask_sequence`])
///
/// # @TODO:
/// Fix this horrible code stench: [`mask_3p`] and [`mask_5p`] quite identical, and may have to much responsibility
#[inline]
fn mask_5p(thresholds: &MaskThreshold, reference: &[u8], seq: &mut [u8], quals: &mut [u8], positions: &[[usize; 2]]) -> Result<(), RuntimeError> {
    // Unwrap cause we have previously validated the struct. [Code smell]
    let mask_5p_threshold = thresholds.get_threshold(&Orientation::FivePrime).unwrap().inner();
    let mask_5p_range     = 0..mask_5p_threshold -1;
    mask_sequence(mask_5p_range, reference, seq, quals, b'C', positions)
}

/// Apply selective masking from the [`Orientation::ThreePrime`] end of a read.
/// 
/// # Parameters:  
/// - `threshold`: reference to a mask [`MaskThreshold`]. This is where the relative [`Position`](`crate::genome::Position`)
///   of the [`ThreePrime`](`Orientation::ThreePrime`) end is retrieved.
/// - `reference`: raw byte representation of reference sequence corresponding to the sequence to mask. 
///    Note that this sequence must contain and take into account any indel found within `seq`
/// - `seq`: raw byte representation of the read sequence to mask. Note that this sequence must preserve any indel found within it.
/// - `quals`: raw byte representation of the phred-scores of `seq`.
/// - `positions`: matching positions found between `seq` and `ref`. These can be retrieved using 
///    [`rust_htslib::bam::Reader::aligned_pairs()`](`rust_htslib::bam::Reader`)
/// 
/// # Errors: 
///  - May emit a [`RuntimeError::ReferenceOutOfIndexError`] (see [`mask_sequence`])
///
/// # @TODO:
/// Fix this horrible code stench: [`mask_3p`] and [`mask_5p`] quite identical, and may have to much responsibility
#[inline]
fn mask_3p(thresholds: &MaskThreshold, reference: &[u8], seq: &mut [u8], quals: &mut [u8], positions: &[[usize; 2]]) -> Result<(), RuntimeError> {
    // Unwrap cause we have previously validated the struct. [Code smell]
    let mask_3p_threshold = thresholds.get_threshold(&Orientation::ThreePrime).unwrap().inner();
    let mask_3p_range     = seq.len().saturating_sub(mask_3p_threshold-1)..seq.len();
    mask_sequence(mask_3p_range, reference, seq, quals, b'G', positions)
}


/// Apply selective masking on any struct implementing [`rust_htslib::bam::Read`], using a reference genome and a
/// structured set of masking thresholds ([`Masks`]). Masked records and then written to the provided `writer`.
/// # Usage
/// ```
/// # use std::error::Error;
/// use rust_htslib::{bam::{self, Read}, faidx};
/// use pmd_mask::mask::Masks;
/// fn main() -> Result<(), Box<dyn Error>> {
/// 
///     // ---- Get an input bam, a reference genome, and a misincorpooration file.
///     let reference  = faidx::Reader::from_path("tests/test-data/reference/hs37d5-MTonly/hs37d5-MTonly.fa.gz")?;
///     let mut reader = bam::Reader::from_path("tests/test-data/bam/dummy-MTonly/dummy-MTonly-1000.bam")?;
///     let masks = Masks::from_path("tests/test-data/bam/dummy-MTonly/misincorporation.txt", 0.01)?;
///
///     // ----- Prepare an output
///     let header = bam::Header::from_template(reader.header());
///     let mut output = bam::Writer::from_stdout(&header, bam::Format::Sam)?;
/// 
///     // ---- Apply pmd-mask
///     pmd_mask::apply_pmd_mask(&mut reader, &reference, &masks, &mut output)?;
/// 
///     Ok(())
/// }
/// ```
#[inline]
pub fn apply_pmd_mask<B>(bam: &mut B, reference: &faidx::Reader, masks: &Masks, writer: &mut bam::Writer) -> Result<()>
where   B: bam::Read,
{
    // ---- Get header template
    let header            = bam::Header::from_template(bam.header());    
    let header_view       = bam::HeaderView::from_header(&header);

    let mut bam_record    = bam::Record::new(); // Input record buffer
    let mut out_record    = bam::Record::new(); // Output record buffer
    let default_threshold = MaskThreshold::default();
    // ---- Loop along input records
    while let Some(result) = bam.read(&mut bam_record) {
        result.unwrap();

        // ---- Get chromosome and strand info of this record.
        // EXPENSIVE: converting tid to string.
        // @ TODO: map Misincorporation chromosome names to tid. once, before looping.
        let current_record = MaskEntry::from_htslib_record(&header_view, &mut bam_record).map_err(RuntimeError::ParseMask)?;
        

        // ---- Get relevant misincorporation frequency:
        let relevant_thresholds = match masks.get(&current_record) {
            Some(threshold) => threshold,
            None => {
                debug!("{current_record} Not found in threshold dictionary. Setting default threshold {default_threshold}");
                &default_threshold
            }
        };


        // ---- Get the reference's position 
        // Memory leak here...
        let refseq = reference.fetch_seq(current_record.chromosome.inner(), bam_record.reference_start() as usize, bam_record.reference_end() as usize -1 )?;
        
        let aligned_pos = bam_record.aligned_pairs().map(|[readpos, refpos]| [readpos as usize, (refpos - bam_record.pos()) as usize]).collect::<Vec<[usize; 2]>>();
    
        trace!("-----------------------");
        trace!("---- Inspecting record: {current_record} {}", bam_record.pos());
        trace!("Relevant thresholds: {relevant_thresholds}");
        trace!("CIGAR              : {}", bam_record.cigar());
        let mut new_seq   = bam_record.seq().as_bytes(); // EXPENSIVE: Allocation
        let mut new_quals = bam_record.qual().to_vec();  // EXPENSIVE: Allocation


        // EXPENSIVE: converting sequence to utf8 for trace logging.
        // @ TODO: use this?: static DECODE_BASE: &[u8] = b"=ACMGRSVTWYHKDBN";
        let sequence = bam_record.seq().as_bytes();
        // @ SAFETY: Sequence is only parsed for TRACE logging, displaying potentially invalid UTF8-character is the whole point..
        trace!("Reference: {}", unsafe { str::from_utf8_unchecked(refseq) });
        trace!("Sequence : {}", unsafe { str::from_utf8_unchecked(&sequence) });

        const UNMAPPED_CONTEXT: &str = "(This is merely a warning because this record was already set as UnMapped (0x4)";
        let err_msg = |e: &RuntimeError , end: Orientation | {
            let position = bam_record.pos();
            format!("While attempting to mask the {end} end of record [{current_record} {position}]: {e}")
        };
        // ---- Mask 5p' positions
        if let Err(e) = mask_5p(relevant_thresholds, refseq, &mut new_seq, &mut new_quals, &aligned_pos) {
            let context = err_msg(&e, Orientation::FivePrime); 
            match bam_record.is_unmapped() {
                true => warn!("@ {context} {UNMAPPED_CONTEXT}"),
                false => return Err(e).with_context(|| format!("{current_record}"))
            }
        };

        // ---- Mask 3p' positions
        if let Err(e) = mask_3p(relevant_thresholds, refseq, &mut new_seq, &mut new_quals, &aligned_pos) {
            let context = err_msg(&e, Orientation::ThreePrime);
            match bam_record.is_unmapped() { 
                true  => warn!("{context} {UNMAPPED_CONTEXT}"), 
                false => return Err(e).context(context)
            }
        };

        // SAFETY: samtools performs UTF8 sanity checks on the raw sequence. So we're ok.
        trace!("Masked   : {}", unsafe{ std::str::from_utf8_unchecked(&new_seq) });

        // ---- Flush tampered record to the output.
        bam_record.clone_into(&mut out_record);
        out_record.set(bam_record.qname(), Some(&bam_record.cigar().take()), &new_seq, &new_quals);
        writer.write(&out_record)?;

        // Manually remove rust-htslib fetch_seq leak.
        unsafe {libc::free(refseq.as_ptr() as *mut std::ffi::c_void)}
    }
    Ok(())
}


/// Test examples to implement for INDELs
///
/// Relevant thresholds: (5p: 5bp) (3p: 6bp)
/// CIGAR              : 48M1I14M
/// Reference: CCAGCCTGGGCGACAGAGCAAGACTCTGTCTCAAAAAAAAAAAAAAAA TGGTGGGGTCGGGG
/// Sequence : CTAGCCTGGGCGACAGAGCAAGACTCTGTCTCAAAAAAAAAAAAAAAAGGGGTGGGGCCGGGG
/// Masked   : NNAGCCTGGGCGACAGAGCAAGACTCTGTCTCAAAAAAAAAAAAAAAAGGGGTGGGGCCNNNN
/// 
/// Relevant thresholds: (5p: 5bp) (3p: 6bp)
/// CIGAR              : 63M1D37M
/// Reference: TGTAGTGAGCTGAGATCGTGCCATTGCACTCCAGCCTGGGCAACAGGAGTGAAACTCTATCTCAAAAAAAAAAAAAAATTAAACAAAAACAAACCTGCCTC
/// Sequence : TGTAGTGAGCTGAGATCGTGCCATTGCACTCCAGCCTGGGCAACAGGAGTGAAACTCTATCTC AAAAAAAAAAAAAATTAAACAAAAACAAACCTGCCTC
/// Masked   : TGTAGTGAGCTGAGATCGTGCCATTGCACTCCAGCCTGGGCAACAGGAGTGAAACTCTATCTC AAAAAAAAAAAAAATTAAACAAAAACAAACCTNCCTC

#[cfg(test)]
mod test {
    use super::*;
    use crate::genome::Position;
    use anyhow::{anyhow, Result};

    //fn get_reference() -> faidx::Reader {
    //    const DUMMY_REFERENCE_PATH: &str = "tests/test-data/reference/hs37d5-MTonly/hs37d5-MTonly.fa.gz";
    //    faidx::Reader::from_path(DUMMY_REFERENCE_PATH).expect("Failed to open reference file with rust-htslib")
    //}

    fn dummy_threshold(position: usize ) -> MaskThreshold {
        let mut threshold = MaskThreshold::default();
        [Orientation::FivePrime, Orientation::ThreePrime].iter().for_each(|end| 
            threshold.set_threshold(*end, Position::new(position))
        );
        threshold
    }

    macro_rules! print_align {
        ($reference:expr, $seq:expr, $quals:expr) => {{
            println!("{}", std::str::from_utf8($reference).unwrap());
            println!("{}", std::str::from_utf8(&$seq).unwrap());
            println!("{}", std::str::from_utf8(&$quals.iter().map(|b| b + 33).collect::<Vec<u8>>()).unwrap());
        }};
    }
    


    #[allow(dead_code)]
    enum Cigar{
        Match(usize),
        Deletion(usize),
        Insertion(usize)
    }


    fn cigar2paired_indices(cigar: &[Cigar]) -> Vec<[usize; 2]> {
        use Cigar::*;
        let mut paired_indices: Vec<[usize; 2]> = Vec::new();

        macro_rules! push_pair {
            ($readincrement:expr, $refincrement:expr) => {
                [paired_indices.last().unwrap()[0] + $readincrement, paired_indices.last().unwrap()[1] + $refincrement]
            };

            ($times:expr, $readincrement:expr, $refincrement:expr) => {
                for _ in 0..$times {
                    paired_indices.push(if paired_indices.is_empty() { [1-$readincrement, 1-$refincrement] } else { push_pair!($refincrement, $readincrement)});
                }
            };
        }

        for cigarid in cigar {
            match cigarid {
                Match(len)     => push_pair!(*len, 1, 1) ,
                Insertion(len) =>  push_pair!(*len, 0, 1),
                Deletion(len)  =>  push_pair!(*len, 1, 0)
            }
        }
        paired_indices
    }



    /// 1. Takes an input reference, nucleotide sequence and PHRED qualities in their string representation,
    /// 2. Convert these to bytes,
    /// 2. Applies masking on these vector using [`mask_5p`] and [`mask_3p`]
    fn mask_and_validate(threshold_len: usize, reference: &str, seq: &str, quals: &str, cigar: &[Cigar]) -> Result<(String, String)> {
        let reference = reference.as_bytes();
        let mut seq   = seq.as_bytes().to_vec();
        let mut quals = quals.as_bytes().iter().map(|b| b - 33).collect::<Vec<u8>>();

        let threshold = dummy_threshold(threshold_len);
        
        // prepare missing errors
        let no_masking_err   = |i| Err(anyhow!("nucleotide {i} should be masked, but isn't"));
        let no_qual_err      = |i| Err(anyhow!("Phred quality {i} should be set to zero, but isn't"));
        let invalid_mask_err = |i| Err(anyhow!("nucleotide {i} has been masked, but should not have been"));

        // Mask 5p
        let pair_indices = cigar2paired_indices(cigar);
        println!("---- 5p masking with threshold set at {threshold_len}");
        mask_5p(&threshold, reference, &mut seq, &mut quals, &pair_indices)?;
        print_align!(reference, seq, quals);

        // ---- Validate 5p masking.
        // Since MaskThresholds are 1-based, and the position it provides is the first
        // position BELOW the threshold, we expect an off by one result. i.e. if threshold is 5, the first
        // four characters are masking candidates (if they are C). BUT, the fifth character should never equal b'N'
        for i in 0..threshold_len - 1 {
            if let Some(nuc) = reference.get(i) {
                match *nuc as char {
                    'C' => {
                        if seq[i] != b'N' {return no_masking_err(i)}
                        if quals[i] != 0  {return no_qual_err(i)}
                    },
                    _   => if seq[i] == b'N' { return invalid_mask_err(i) }, 
                };
            }
        }
        if let Some(nuc) = seq.get(threshold_len) {
            if *nuc == b'N' { return invalid_mask_err(threshold_len)}
        }

        // Mask 3p
        println!("---- 3p masking with threshold set at {threshold_len}");
        mask_3p(&threshold, reference, &mut seq, &mut quals, &pair_indices)?;
        print_align!(reference, seq, quals);

        // ---- Validate 3p masking
        // Since MaskThreshold is 1 based, AND the position it gives is the first position BELOW the threhold,
        // we expect an off by one result. i.e. if threshold is 5, the last 4 characters are candidates for maasking
        // (if they are 'G') BUT, the fifth to last character should never equal b'N'
        let start = seq.len().saturating_sub(threshold_len-1);
        for i in start..seq.len(){
            if let Some(nuc) = reference.get(i) {
                match *nuc as char {
                    'G' => {
                        if seq[i] != b'N' {return no_masking_err(i)}
                        if quals[i] != 0  {return no_qual_err(i)}
                    },
                    _   => if seq[i] == b'N' && *nuc != b'C' { return invalid_mask_err(i) }, 
                };
            }
        }
        if start > 0 {
            if let Some(nuc) = seq.get(start-1) {
                if *nuc == b'N' && reference[start-1] != b'C' { return invalid_mask_err(start-1) } 
            }
        }

        // Convert back to string and return
        let out_seq = std::str::from_utf8(&seq).unwrap().to_string();
        let out_quals = std::str::from_utf8(&quals.iter().map(|b| b + 33).collect::<Vec<u8>>()).unwrap().to_string();

        Ok((out_seq, out_quals))
    }

    #[test]
    fn mask_basic_sequence_threshold_le_seq_len() {
        let reference = "GCTCCTATTAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACGGGGGGGGG";
        let seq       = "GCTCCTATTAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATG";
        let quals     = "AAAAAEEEEEEEEEEEEEEEEEEEJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ";
        let cigar     = vec![Cigar::Match(seq.len())];
        for i in 2..reference.len() {
            match mask_and_validate(i, reference, seq, quals, &cigar) {
                Err(e) => panic!("{e}"),
                Ok((new_seq, new_quals)) => {
                    assert!(new_seq.contains(['C', 'G']));
                    assert!(new_quals.contains('!'));
                }
            };
        }
    }

    #[test]
    fn mask_basic_sequence_eq_seq_len() {
        let reference = "GTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACGGGGGGGGC";
        let seq       = "GTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATC";
        let quals     = "AAAAAEEEEEEEEEEEEEEEEEEEJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ";
        let cigar     = vec![Cigar::Match(seq.len())];
        match mask_and_validate(seq.len(), reference, seq, quals, &cigar) {
            Err(e) => panic!("{e}"),
            Ok((new_seq, new_quals)) => {
                // Threshold is exatly at the sequence length.
                // -> first and last nucleotites must remain untouched, even though they are candidates 
                //    in this instance.
                assert_eq!(new_seq.chars().next().unwrap(), 'G');
                assert_eq!(new_seq.chars().last().unwrap(), 'C');
                
                assert_ne!(new_quals.chars().next().unwrap(), '!');
                assert_ne!(new_quals.chars().last().unwrap(), '!');
                // However, anything between the first and last character must not contain any C / G
                assert!(!new_seq[1..new_seq.len()-1].contains(['C', 'G']));
            }
        }
    }

    #[test]
    fn mask_basic_sequence_threshold_gt_seq_len() {
        let reference = "GTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACGGGGGGGGG";
        let seq       = "GTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATC";
        let quals     = "AAAAAEEEEEEEEEEEEEEEEEEEJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ";
        let cigar     = vec![Cigar::Match(seq.len())];
        for i in reference.len()+1..reference.len()*2 {
            match mask_and_validate(i, reference, seq, quals, &cigar) {
                Err(e) => panic!("{e}"),
                Ok((sequence, _))  => {
                    assert!(!sequence.contains(['C', 'G'])) // Since our threshold goes beyond the length, everything should be masked.
                }
            };
        }

    }

    #[test]
    fn mask_basic_sequence_g_only() {
        let reference = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
        let seq       = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let quals     = "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ";
        let cigar     = vec![Cigar::Match(seq.len())];
        match mask_and_validate(1, reference, seq, quals, &cigar) {
            Err(e) => panic!("{e}"),
            Ok((seq, quals)) => {
                assert!(!seq.contains('N'));
                assert!(!quals.contains('!'));
            }
        }

        match mask_and_validate(seq.len()+1, reference, seq, quals, &cigar) {
            Err(e) => panic!("{e}"),
            Ok((seq, quals)) => {
                assert!(seq.chars().all(|c| c == 'N'));
                assert!(quals.chars().all(|c| c == '!'));
            }
        }
    }

    #[test]
    fn mask_basic_sequence_c_only() {
        let reference = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        let seq       = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let quals     = "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ";
        let cigar     = vec![Cigar::Match(seq.len())]; 
        match mask_and_validate(1, reference, seq, quals, &cigar) {
            Err(e) => panic!("{e}"),
            Ok((seq, quals)) => {
                assert!(!seq.contains('N'));
                assert!(!quals.contains('!'));
            }
        }

        match mask_and_validate(seq.len()+1, reference, seq, quals, &cigar) {
            Err(e) => panic!("{e}"),
            Ok((seq, quals)) => {
                assert!(seq.chars().all(|c| c == 'N'));
                assert!(quals.chars().all(|c| c == '!'));
            }
        }
    }

    #[test]
    fn mask_spurious_unmapped_sequence() {
        let reference = "N";
        let seq       = "GGATCACAGGTCTATCACCCTATTAACCACTCACGG";
        let quals     = "AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ";
        let cigar     = vec![Cigar::Match(seq.len())];
        println!("{:?}", mask_and_validate(10, reference, seq, quals, &cigar));
    }
}


