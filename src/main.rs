use std::str;

use rust_htslib::{faidx, bam, bam::Read, tpool::ThreadPool};

use num_cpus::{self};

use clap::Parser;

use anyhow::Result;

mod logger;
use logger::Logger;

mod parser;
use parser::Cli;

mod misincorporation;
use misincorporation::Misincorporations;

mod genome;
use genome::{ChrName, Orientation, Strand};

mod mask;
use mask::{Masks, MaskEntry, MaskThreshold};


use log::{info, warn, debug, trace};



// What's gonna happen:
// 1. Open a misincorporation file and map each relative position to a given probability.
// 2. Within this misincorporation file, find the relative position at which the probability of
//    either C->T (forward) or G->A (reverse) misincorporation is lower than the user-provided threshold
// 3. Open and read a bam file. 
//    For each line of the bam file, look through the sequence read.
//    - for forward strands, if at any point there is a C in the **reference** sequence AND we've not attained the treshold, mask it.
//    - for reverse strands, if at any point there is a G in the **reference** sequence AND we've not attained the treshold, mask it.

// ---- Clap:
// - input (maybe stdin)
// - misincorporation-file
// - reference
// - treshold
// - output (maybe stdout)

// ---- Sanity checks:
// Ensure Provided reference matches the reference in the bam header.
// [DONE] Warn user if misincorporation frequencies are Inf, Nan, or were computed from empty C.

fn main() -> Result<()> {
    // ---- Parse arguments
    let args = Cli::parse();

    // ---- initialize logger
    Logger::init(if args.quiet {0} else {args.verbose + 1} );
    debug!("Provided Command line Arguments:{args:#?}");

    // ---- define threads:
    let threads = if args.threads == 0 {num_cpus::get() as u32 } else {args.threads};
    info!("Setting threads to {threads}");
    let pool = ThreadPool::new(threads)?;


    // ---- Read Misincorporation file as a tsv file.
    let mut threshold_positions = Misincorporations::from_path(args.misincorporation, args.threshold).unwrap();

    // ---- Validate misincorporation and issue warnings for any 'abnormal' frequency found.
    let mut abnormal_frequencies = threshold_positions.extrude_invalid_frequencies().iter().fold(String::new(), |abnormal_freqs, pos| {
        abnormal_freqs + &format!("{pos}\n")
    });
    abnormal_frequencies.pop();

    if !abnormal_frequencies.is_empty() {
        warn!("Found abnormal misincorporation frequencies (NaN, Infinite values, etc.). Masking will apply along the full length of the sequence for those:\n{abnormal_frequencies}");
    }

    debug!("Using the following positions as threshold for masking:\n{}",
        threshold_positions.iter().fold(String::new(), |acc, val| {
            acc + &format!("{val}\n")
        })
    );

    // ---- Convert to a more usable format. This is a crash test after all...
    let thresholds = Masks::try_from(&threshold_positions)?;

    // ---- Open Reference File
    let reference = faidx::Reader::from_path(&args.reference)?;


    // ---- Open bam file
    let mut bam = match args.bam {
        Some(path) => bam::Reader::from_path(path),
        None       => bam::Reader::from_stdin(),
    }?;

    bam.set_reference(&args.reference)?;
    bam.set_thread_pool(&pool)?;

    // ---- Get header template
    let header = bam::Header::from_template(bam.header());    
    let header_view = bam::HeaderView::from_header(&header);
    let mut record = bam::Record::new();


    // ---- Define output format
    let out_fmt = if let Some(fmt) = args.output_fmt {
        match fmt.to_ascii_uppercase().as_str() {
            "B" | "BAM" => bam::Format::Bam,
            "S" | "SAM" => bam::Format::Sam,
            "C" | "CRAM" => bam::Format::Cram,
            _ => panic!("Invalid output format"),
        }
    } else {
        bam::Format::Sam // It'd be nice to set this to the same format as the input...
    }; 
    
    // ---- Prepare Bam Writer
    let mut writer = match args.output {
        Some(path) => bam::Writer::from_path(path, &header, out_fmt),
        None       => bam::Writer::from_stdout(&header, out_fmt),
    }?;
    writer.set_compression_level(bam::CompressionLevel::Level(args.compress_level))?;
    writer.set_thread_pool(&pool)?;
    writer.set_reference(args.reference)?;


    let mut out_record    = bam::Record::new();
    let default_threshold = MaskThreshold::default();
    // ---- Loop along input records
    while let Some(result) = bam.read(&mut record) {
        result.unwrap();

        // ---- Get chromosome and strand info of this record.
        // EXPENSIVE: converting tid to string.
        // @ TODO: map Misincorporation chromosome names to tid. once, before looping.
        let chromosome = ChrName(str::from_utf8(header_view.tid2name(record.tid() as u32))?.to_string());
        let strand = record.strand().strand_symbol().parse::<Strand>()?;
        let current_record = MaskEntry{chromosome, strand};
        

        // ---- Get relevant misincorporation frequency:
        let relevant_thresholds = match thresholds.get(&current_record) {
            Some(threshold) => threshold,
            None => {
                debug!("{current_record} Not found in threshold dictionary. Setting default threshold {default_threshold}");
                &default_threshold
            }
        };

        // EXPENSIVE: converting sequence to utf8
        // @ TODO: use this?: static DECODE_BASE: &[u8] = b"=ACMGRSVTWYHKDBN";
        let sequence = record.seq().as_bytes(); // That's an expensive operation.
        let sequence = str::from_utf8(&sequence)?;

        // ---- Get the reference's position
        let position = record.pos();
        let end = record.seq_len();
        let reference = reference.fetch_seq(&current_record.chromosome.0, position as usize, position as usize + end - 1)?;

        trace!("-----------------------");
        trace!("---- Inspecting record: {current_record} {position}");
        trace!("Relevant thresholds: {relevant_thresholds}");

        let mut new_seq   = record.seq().as_bytes(); // EXPENSIVE: Allocation
        let mut new_quals = record.qual().to_vec();  // EXPENSIVE: Allocation

        trace!("Reference: {}", std::str::from_utf8(reference)?);
        trace!("Sequence : {sequence}");

        // ---- Mask 5p' positions
        // Unwrap cause we have previously validated the struct. [Code smell]
        let mask_5p_threshold = relevant_thresholds.get_threshold(&Orientation::FivePrime).unwrap().0;
        'mask5p: for i in 0..mask_5p_threshold {

            // ---- Bail if our index has gone past the read length.
            if i >= record.seq_len() { break 'mask5p}

            if reference[i] == b'C' {
                new_seq[i] = b'N';
                new_quals[i] = 0;
            }
        }

        // ---- Mask 3p' positions
        // Unwrap cause we have previously validated the struct. [Code smell]
        let mask_3p_threshold = relevant_thresholds.get_threshold(&Orientation::ThreePrime).unwrap().0;
        'mask3p: for i in (record.seq_len().saturating_sub(mask_3p_threshold)..record.seq_len()).rev() {
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
        record.clone_into(&mut out_record);
        out_record.set(record.qname(), Some(&record.cigar().take()), &new_seq, &new_quals);
        writer.write(&out_record)?;
    }

    Ok(())
}
