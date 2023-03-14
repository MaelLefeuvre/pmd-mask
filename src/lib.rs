use std::{str, ops::Range};


pub mod misincorporation;
pub mod genome;
pub mod mask;
pub mod error;

use error::RuntimeError;
use genome::Orientation;
use mask::{Masks, MaskEntry, MaskThreshold};

use anyhow::Result;
use rust_htslib::{faidx, bam};
use log::{debug, trace};

#[inline]
fn mask_sequence(range: Range<usize>, reference: &[u8], seq: &mut [u8], quals: &mut [u8], target_nucleotide: u8) {
    'mask: for i in range {
        if i >= seq.len() { break 'mask}
        if reference[i] == target_nucleotide {
            seq[i] = b'N';
            quals[i] = 0;
        }
    }
}

#[inline]
fn mask_5p(thresholds: &MaskThreshold, reference: &[u8], seq: &mut [u8], quals: &mut [u8]) {
    // Unwrap cause we have previously validated the struct. [Code smell]
    let mask_5p_threshold = thresholds.get_threshold(&Orientation::FivePrime).unwrap().inner();
    let mask_5p_range     = 0..mask_5p_threshold -1;
    mask_sequence(mask_5p_range, reference, seq, quals, b'C');
}

#[inline]
fn mask_3p(thresholds: &MaskThreshold, reference: &[u8], seq: &mut [u8], quals: &mut [u8]) {
    // Unwrap cause we have previously validated the struct. [Code smell]
    let mask_3p_threshold = thresholds.get_threshold(&Orientation::ThreePrime).unwrap().inner();
    let mask_3p_range     = seq.len().saturating_sub(mask_3p_threshold-1)..seq.len();
    mask_sequence(mask_3p_range, reference, seq, quals, b'G');
}

#[inline]
pub fn apply_pmd_mask<B>(bam: &mut B, reference: &faidx::Reader, masks: &Masks, writer: &mut bam::Writer) -> Result<(), RuntimeError>
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
        let position = bam_record.pos();
        let end = bam_record.seq_len();
        let reference = reference.fetch_seq(current_record.chromosome.inner(), position as usize, position as usize + end - 1)?;

        trace!("-----------------------");
        trace!("---- Inspecting record: {current_record} {position}");
        trace!("Relevant thresholds: {relevant_thresholds}");

        let mut new_seq   = bam_record.seq().as_bytes(); // EXPENSIVE: Allocation
        let mut new_quals = bam_record.qual().to_vec();  // EXPENSIVE: Allocation


        // EXPENSIVE: converting sequence to utf8 for trace logging.
        // @ TODO: use this?: static DECODE_BASE: &[u8] = b"=ACMGRSVTWYHKDBN";
        let sequence = bam_record.seq().as_bytes();
        // @ SAFETY: Sequence is only parsed for TRACE logging, displaying potentially invalid UTF8-character is the whole point..
        trace!("Reference: {}", unsafe { str::from_utf8_unchecked(reference) });
        trace!("Sequence : {}", unsafe { str::from_utf8_unchecked(&sequence) });

        // ---- Mask 5p' positions
        mask_5p(relevant_thresholds, reference, &mut new_seq, &mut new_quals);

        // ---- Mask 3p' positions
        mask_3p(relevant_thresholds, reference, &mut new_seq, &mut new_quals);

        // SAFETY: samtools performs UTF8 sanity checks on the raw sequence. So we're ok.
        trace!("Masked   : {}", unsafe{ std::str::from_utf8_unchecked(&new_seq) });

        // ---- Flush tampered record to the output.
        bam_record.clone_into(&mut out_record);
        out_record.set(bam_record.qname(), Some(&bam_record.cigar().take()), &new_seq, &new_quals);
        writer.write(&out_record)?;
    }
    Ok(())
}


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
    
    /// 1. Takes an input reference, nucleotide sequence and PHRED qualities in their string representation,
    /// 2. Convert these to bytes,
    /// 2. Applies masking on these vector using [`mask_5p`] and [`mask_3p`]
    fn mask_and_validate(threshold_len: usize, reference: &str, seq: &str, quals: &str) -> Result<(String, String)> {
        let reference = reference.as_bytes();
        let mut seq   = seq.as_bytes().to_vec();
        let mut quals = quals.as_bytes().iter().map(|b| b - 33).collect::<Vec<u8>>();

        let threshold = dummy_threshold(threshold_len);
        
        // prepare missing errors
        let no_masking_err   = |i| Err(anyhow!("nucleotide {i} should be masked, but isn't"));
        let no_qual_err      = |i| Err(anyhow!("Phred quality {i} should be set to zero, but isn't"));
        let invalid_mask_err = |i| Err(anyhow!("nucleotide {i} has been masked, but should not have been"));

        // Mask 5p
        println!("---- 5p masking with threshold set at {threshold_len}");
        mask_5p(&threshold, reference, &mut seq, &mut quals);
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
        mask_3p(&threshold, reference, &mut seq, &mut quals);
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
        for i in 2..reference.len() {
            match mask_and_validate(i, reference, seq, quals) {
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
        match mask_and_validate(seq.len(), reference, seq, quals) {
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
        for i in reference.len()+1..reference.len()*2 {
            match mask_and_validate(i, reference, seq, quals) {
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
        match mask_and_validate(1, reference, seq, quals) {
            Err(e) => panic!("{e}"),
            Ok((seq, quals)) => {
                assert!(!seq.contains('N'));
                assert!(!quals.contains('!'));
            }
        }

        match mask_and_validate(seq.len()+1, reference, seq, quals) {
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
        match mask_and_validate(1, reference, seq, quals) {
            Err(e) => panic!("{e}"),
            Ok((seq, quals)) => {
                assert!(!seq.contains('N'));
                assert!(!quals.contains('!'));
            }
        }

        match mask_and_validate(seq.len()+1, reference, seq, quals) {
            Err(e) => panic!("{e}"),
            Ok((seq, quals)) => {
                assert!(seq.chars().all(|c| c == 'N'));
                assert!(quals.chars().all(|c| c == '!'));
            }
        }
    }
}


