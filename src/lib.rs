use std::str;


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

        // EXPENSIVE: converting sequence to utf8
        // @ TODO: use this?: static DECODE_BASE: &[u8] = b"=ACMGRSVTWYHKDBN";
        let sequence = bam_record.seq().as_bytes(); // That's an expensive operation.
        let sequence = str::from_utf8(&sequence)?;

        // ---- Get the reference's position
        let position = bam_record.pos();
        let end = bam_record.seq_len();
        let reference = reference.fetch_seq(current_record.chromosome.inner(), position as usize, position as usize + end - 1)?;

        trace!("-----------------------");
        trace!("---- Inspecting record: {current_record} {position}");
        trace!("Relevant thresholds: {relevant_thresholds}");

        let mut new_seq   = bam_record.seq().as_bytes(); // EXPENSIVE: Allocation
        let mut new_quals = bam_record.qual().to_vec();  // EXPENSIVE: Allocation

        trace!("Reference: {}", std::str::from_utf8(reference)?);
        trace!("Sequence : {sequence}");

        // ---- Mask 5p' positions
        // Unwrap cause we have previously validated the struct. [Code smell]
        let mask_5p_threshold = relevant_thresholds.get_threshold(&Orientation::FivePrime).unwrap().inner();
        'mask5p: for i in 0..mask_5p_threshold - 1 {

            // ---- Bail if our index has gone past the read length.
            if i >= bam_record.seq_len() { break 'mask5p}

            if reference[i] == b'C' {
                new_seq[i] = b'N';
                new_quals[i] = 0;
            }
        }

        // ---- Mask 3p' positions
        // Unwrap cause we have previously validated the struct. [Code smell]
        let mask_3p_threshold = relevant_thresholds.get_threshold(&Orientation::ThreePrime).unwrap().inner();
        'mask3p: for i in (bam_record.seq_len().saturating_sub(mask_3p_threshold - 1)..bam_record.seq_len()).rev() {
            // ---- Bail if our index has gone past the read length.
            if i == 0 { break 'mask3p }
            
            if reference[i] == b'G' {
                new_seq[i] = b'N';
                new_quals[i] = 0;
            }

        }

        // SAFETY: samtools performs UTF8 sanity checks on the raw sequence. So we're ok.
        trace!("Masked   : {}", unsafe{ std::str::from_utf8_unchecked(&new_seq) });

        // ---- Flush tampered record to the output.
        bam_record.clone_into(&mut out_record);
        out_record.set(bam_record.qname(), Some(&bam_record.cigar().take()), &new_seq, &new_quals);
        writer.write(&out_record)?;
    }
    Ok(())
}
