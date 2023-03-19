//! # `pmd-mask` Main binary.
//! 
//! # What `pmd-mask` needs:
//! 1. A sam/bam/cram file
//! 2. A `misincorporation.txt` file, obtained from the same input bam file, using MapDamage
//! 3. A reference genome
//! # What `pmd-mask` does:
//! 
//! 1. Opens a MapDamage misincorporation file and maps each relative position to a given probability of misincorporation.
//! 
//! 2. Within this misincorporation file, finds the relative position at which the probability of either C->T (forward) or 
//!    G->A (reverse) misincorporation is lower than the user-provided threshold (the default is 1%)
//! 
//!    - This relative position is defined in base pairs, starting from the relevant read end, (i.e: we start counting from
//!      the `5p` end for `C>T` misincorporations | `3p` end, for `G>A` misincorporations. 
//! 
//! 3. Open and read through the input bam file. 
//!    For each line of the bam file, look through the sequence read.
//!    - starting at the 5p end, if at any point there is a C in the **reference** sequence AND we've not attained the treshold, mask it.
//!    - starting at the 3p end, if at any point there is a G in the **reference** sequence AND we've not attained the treshold, mask it.
//! 
//!    By masking, we mean *hard masking*. i.e. downscale the Phred base-quality to 0, and set the nucleotide value to 'N'.
//! 
//! 4. Stream the masked sequences to either an output file, or the standard output. The output format and compression level 
//!    can be specified by the used. 

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use pmd_mask::apply_pmd_mask;
use pmd_mask::mask::Masks;
use pmd_mask::error::RuntimeError;

mod logger;
use logger::Logger;

mod parser;
use parser::Cli;

use clap::Parser;
use anyhow::Result;
use rust_htslib::{faidx, bam, bam::Read, tpool::ThreadPool};
use rust_htslib::errors::Error as HtslibError;

use log::{error, info, debug};




/// Open a bam, from either a file, or from standard input.
fn open_bam_reader(maybe_file: &Option<impl AsRef<Path>>) -> Result<bam::Reader, HtslibError> {
    match maybe_file {
        Some(ref path) => { info!("Opening {}", &path.as_ref().display()); bam::Reader::from_path(path) },
        None           => { info!("Reading bam from standard input"); bam::Reader::from_stdin()},
    }
}


/// Open a bam writer, either to a file, or to standard output.
fn open_bam_writer<F>(maybe_file: &Option<F>, header: &bam::Header, output_fmt: bam::Format) -> Result<bam::Writer, HtslibError>
where   F: AsRef<Path> 
{
    match maybe_file {
        Some(ref path) => {info!("Writing output to {}", &path.as_ref().display()); bam::Writer::from_path(path, header, output_fmt)},
        None           => {info!("Writing output to standard output"); bam::Writer::from_stdout(header, output_fmt)},
    }
}


fn run(args: &Cli) -> Result<(), RuntimeError> {

    // ---- Ensure Input bam and output bam are not the same.
    if let Some(ref input_file) = args.bam{
        if let Some(ref output_file) = args.output {
            if *input_file == *output_file {
                return Err(RuntimeError::InputIsOutput)
            }
        }
    }



    // ---- define threadpool if required:
    let thread_pool = match args.threads {
        1    => None ,
        more => {debug!("Firing up threadpool..."); Some(ThreadPool::new(more)?) }
    };

    // ---- Read Misincorporation file as a tsv file and obtain a list of Masking threshold
    //      for each chromosome, strand, and orientation.
    info!("Computing masking positions from {}, using {} as threshold", &args.misincorporation.display(), args.threshold);
    let thresholds = Masks::from_path(&args.misincorporation, args.threshold)?;

    if let Some(ref file) = args.metrics_file {
        info!("Writing masking thresholds to {}", file.display());
        let mut metrics_writer = BufWriter::new(File::create(file).map_err(RuntimeError::OpenMetrics)?);
        thresholds.write(&mut metrics_writer).map_err(RuntimeError::WriteMasksMetrics)?;
    }

    // ---- Open Reference File
    info!("Opening reference file {}", &args.reference.display());
    let reference = faidx::Reader::from_path(&args.reference)?;

    // ---- Open bam file
    let mut bam = open_bam_reader(&args.bam)?;
    
    // ---- IndexedRead: Does not seem to provide with any underlying async/multithread support either..
    //let mut bam = bam::IndexedReader::from_path(args.bam.unwrap()).unwrap();
    //bam.fetch(bam::FetchDefinition::All);
    

    // ---- Define an output format if the user never specified it.
    // @TODO: It'd be nice to set this to the same format as the input... rust_htslib might have a way to access the header's magic number
    let output_format = args.output_fmt.unwrap_or(bam::Format::Sam);

    // ---- Prepare Bam Writer
    let output_header = bam::Header::from_template(bam.header());
    let mut writer = open_bam_writer(&args.output, &output_header, output_format)?;

    // ---- Set output compression level for BAM/CRAM output.
    writer.set_compression_level(bam::CompressionLevel::Level(args.compress_level))?;

    // ---- Set reference for CRAM files. 
    bam.set_reference(&args.reference)?;
    writer.set_reference(&args.reference)?;

    // ---- Set thread pool if the user requested multi-threading
    if let Some(ref pool) = thread_pool { 
        debug!("Allocating threadpool to Reader and Writer");
        bam.set_thread_pool(pool)?;
        writer.set_thread_pool(pool)?;
    };

    info!("Applying PMD-masking...");
    apply_pmd_mask(&mut bam, &reference, &thresholds, &mut writer)?;
    info!("Done");
    Ok(())
}

fn main() {
    // ---- Parse arguments
    let args = Cli::parse();

    // ---- Initialize logger
    Logger::init(if args.quiet {0} else {args.verbose + 1} );

    debug!("Provided Command line Arguments:{args:#?}");

    // ---- Run main process
    if let Err(e) = run(&args) {
        error!("{e}");
        std::process::exit(1)
    }
}
