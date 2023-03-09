use std::collections::HashMap;

use rust_htslib::{bam, bam::Read};
use rust_htslib::faidx;


use csv::{ReaderBuilder};
use num_cpus::{self};

use clap::Parser;

mod logger;

mod parser;
use parser::Cli;

mod misincorporation;
use misincorporation::MisincorporationRecord;

mod genome;
use genome::{ChrName, Orientation, Strand};

mod mask;
use mask::{MaskEntry, MaskThreshold};


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

fn main() {
    // ---- Parse arguments
    let args = Cli::parse();

    // ---- initialize logger
    logger::Logger::init(if args.quiet {0} else {args.verbose + 1} );
    debug!("Provided Command line Arguments:{args:#?}");

    // ---- define threads:
    let threads = if args.threads == 0 {num_cpus::get() as u32 } else {args.threads};
    info!("Setting threads to {threads}");
    let pool = rust_htslib::tpool::ThreadPool::new(threads).unwrap();


    // ---- Read Misincorporation file as a tsv file.
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .comment(Some(b'#'))
        .from_path(args.misincorporation)
        .unwrap(); 

    let mut skip_chromosome: Option<(&ChrName, &Orientation, &Strand)> = None; 
    let mut threshold_positions = Vec::with_capacity(32 * 2 * 2);

    'nextline: for result in reader.deserialize::<MisincorporationRecord>() {
        let record = result.unwrap();

        // Once we've found a position at which the threshold is met, we 
        // don't need to parse entries which belong from the same chromosome.
        if let Some(chromosome_to_skip) = skip_chromosome {
            if (&record.chromosome, &record.end, &record.strand) == chromosome_to_skip {
                continue 'nextline
            }
        }

        // If we're below the requested treshold, keep that record!
        if record.target_freq() <= args.threshold {
            threshold_positions.push(record);
            let last_insert = threshold_positions.last().unwrap(); // We can unwrap here since we know we've just pushed a value
            skip_chromosome = Some((&last_insert.chromosome, &last_insert.end, &last_insert.strand)); 
            
        }
    }

    // ---- Validate misincorporation and issue warnings for any 'abnormal' frequency found.
    let mut abnormal_frequencies = String::new();
    let threshold_positions: Vec<MisincorporationRecord> = threshold_positions.into_iter().filter_map(|pos| {
        if pos.target_freq().is_normal() {
            Some(pos)
        } else {
            abnormal_frequencies += &format!("{pos}\n");
            None
        }
    }).collect();
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
    let mut thresholds = HashMap::new();
    for position in threshold_positions {
        let record = MaskEntry{chromosome: position.chromosome, strand: position.strand};
        thresholds.entry(record)
            .or_insert(MaskThreshold::default())
            .set_threshold(position.end, position.position);
    }

    // ---- Validate thresholds: Each record must contain a hashmap w/ two values (one for 3p, the other for 5p)
    for (record, threshold) in thresholds.iter() {
        trace!("Validating {record:?} {threshold:?}");
        threshold.validate().unwrap()
    };

    // ---- Open Reference File
    let reference = faidx::Reader::from_path(&args.reference).unwrap();


    // ---- Open bam file
    let mut bam = match args.bam {
        Some(path) => rust_htslib::bam::Reader::from_path(path),
        None       => rust_htslib::bam::Reader::from_stdin(),
    }.unwrap();

    bam.set_reference(&args.reference).unwrap();
    bam.set_thread_pool(&pool).unwrap();

    // ---- Get header template
    let header = rust_htslib::bam::Header::from_template(bam.header());    
    let header_view = rust_htslib::bam::HeaderView::from_header(&header);
    let mut record = rust_htslib::bam::Record::new();


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
        Some(path) => rust_htslib::bam::Writer::from_path(path, &header, out_fmt).unwrap(),
        None       => rust_htslib::bam::Writer::from_stdout(&header, out_fmt).unwrap(),
    };
    writer.set_compression_level(rust_htslib::bam::CompressionLevel::Level(args.compress_level)).unwrap();
    writer.set_thread_pool(&pool).unwrap();
    writer.set_reference(args.reference).unwrap();


    let mut out_record    = rust_htslib::bam::Record::new();
    let default_threshold = MaskThreshold::default();
    // ---- Loop along input records
    'record: while let Some(result) = bam.read(&mut record) {
        result.unwrap();

        // ---- Get chromosome and strand info of this record.
        // EXPENSIVE: converting tid to string.
        // @ TODO: map Misincorporation chromosome names to tid. once, before looping.
        let chromosome = ChrName(std::str::from_utf8(header_view.tid2name(record.tid() as u32)).unwrap().to_string());
        let strand = record.strand().strand_symbol().parse::<Strand>().unwrap();
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
        let sequence = std::str::from_utf8(&sequence).unwrap();

        // ---- Get the reference's position
        let position = record.pos();
        let end = record.seq_len();
        let reference = reference.fetch_seq(&current_record.chromosome.0, position as usize, position as usize + end - 1).unwrap();

        trace!("-----------------------");
        trace!("---- Inspecting record: {current_record} {position}");
        trace!("Relevant thresholds: {relevant_thresholds}");

        let mut new_seq   = record.seq().as_bytes(); // EXPENSIVE: Allocation
        let mut new_quals = record.qual().to_vec();  // EXPENSIVE: Allocation

        trace!("Reference: {}", std::str::from_utf8(reference).unwrap());
        trace!("Sequence : {sequence}");

        // ---- Mask 5p' positions
        'mask5p: for i in 0..relevant_thresholds.get_threshold(&Orientation::FivePrime).0 {
            if i >= record.seq_len() {
                break 'mask5p
            }

            //let decoded_nuc = unsafe { record.seq().decoded_base_unchecked(i) };
            if reference[i] == b'C' {
                new_seq[i] = b'N';
                new_quals[i] = 0;
            }
        }

        // ---- Mask 3p' positions
        'mask3p: for i in (record.seq_len()-relevant_thresholds.get_threshold(&Orientation::ThreePrime).0..record.seq_len()).rev() {
            if i == 0 {
                break 'mask3p
            }
            
            //let decoded_nuc = unsafe { record.seq().decoded_base_unchecked(i) };
            if reference[i] == b'G' {
                new_seq[i] = b'N';
                new_quals[i] = 0;
            }

        }

        trace!("Masked   : {}", std::str::from_utf8(&new_seq).unwrap());


        record.clone_into(&mut out_record);
        out_record.set(record.qname(), Some(&record.cigar().take()), &new_seq, &new_quals);
        //println!("{record:#?}");
        writer.write(&out_record).unwrap();
    }

}
