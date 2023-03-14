use std::path::PathBuf;

mod error;
use error::CliError;

use clap::{Parser, ArgAction, ColorChoice};
use num_cpus::{self};
use rust_htslib::bam;
use log::info;


/// Convert the user provided output format string to a htslib-friendly enum
/// 
/// # Behaviour:
///
/// - Accepts initials or fully-specified, case-insensitive file formats for BAM, SAM, CRAM
/// 
/// i.e: b, B, Bam, BAM, bam, BaM, bAm, etc. will work., 
/// d, V, vam, BA,  bame, BLAM, abam, etc. will return an error. 
/// 
/// 
/// # Errors
/// returns a [`CliError::InvalidBamOutputFmt`] in case the function fails to pattern match the provided string with a [`bam::Format`] enum variant.
fn parse_output_fmt(s: &str) -> Result<bam::Format, CliError> {
    match s.to_ascii_uppercase().as_str() {
        "S" | "SAM"  => Ok(bam::Format::Sam),
        "B" | "BAM"  => Ok(bam::Format::Bam),
        "C" | "CRAM" => Ok(bam::Format::Cram),
        _ => Err(CliError::InvalidBamOutputFmt),
    }
}


// parses the user-provided string into a u32, specifying the number of allocated
// 
fn parse_threads(s: &str) -> Result<u32, CliError> {
    Ok(match s.parse::<u32>().map_err(CliError::InvalidThreadValue)? {
        0 => {
            let available_cores = num_cpus::get() as u32;
            info!("Setting threadpool to all available cores ({available_cores})");
            available_cores 
        }
        other => {
            match other {
                1    => info!("Requesting a single thread for computation"),
                more => info!("Setting threadpool to {more} additional worker threads."),
            }
            other
        }
    })
}


/// pmd-mask: Perform hard selective masking of ancient DNA deamination patterns, using the output misincorporation frequency estimates of MapDamage (see: https://github.com/ginolhac/mapDamage.git).
#[derive(Parser, Debug)]
#[command(name="pmd-mask", author, version, about, long_about = None, color=ColorChoice::Always)]
#[clap(propagate_version = true)]
pub struct Cli {
    /// Misincorporation frequency threshold.
    /// 
    /// This is the misincorporation frequency threshold at which masking is applied on reads.
    /// pmd-mask will search through the provided misincorporation file, and obtain the base-pair position at which this threshold is mask.  
    /// 
    /// - At the 5p end, pmd-mask will mask all encountered reference Cytosines starting at the start of the read until the threshold is met.  
    /// - At the 3p end, pmd-mask will mask all encountered reference Guanines starting at the end of the read, until the threshold is met.   
    #[arg(short, long, default_value("0.01"))]
    pub threshold: f32,

    /// Input alignment file (SAM|BAM|CRAM)
    /// 
    /// Input bam file, on which pmd-masking should be performed. When unspecified, the pmd-mask will look for standard input.
    #[arg(short, long, required(false))]
    pub bam: Option<PathBuf>,

    /// Output file (SAM|BAM|CRAM. See --output-fmt).
    /// 
    /// Output bam file If not specified, print to stdout.
    #[arg(short, long, required(false))]
    pub output: Option<PathBuf>,

    /// Output format. (SAM|BAM|CRAM). 
    /// 
    /// If unspecified: Will output to SAM if stdout is targeted. Or use the original input format.
    #[arg(short='O', long, required(false), value_parser(parse_output_fmt))]
    pub output_fmt: Option<bam::Format>,

    /// Output compression level. 
    /// 
    /// Set the compression for BAM|CRAM output files. 9 if unspecified
    ///
    /// Note that specifying this value is meaningless for SAM files.
    #[arg(long, default_value("9"))]
    pub compress_level: u32,

    /// Reference genome. (fasta|fa)[.gz]
    /// 
    /// Path to a reference genome, in fasta file format. The provided file must be indexed, and the corresponding .fai file should be located at the same directory, and carry the same name.
    /// 
    /// The provided reference genome must, of course, be the same as the one used to align the input sequences. At this stage, pmd-mask does **NOT** ensure the provided reference is consistent with the input's header information, and thus using a different reference may result in undefined behavior.
    #[arg(short='f', long)]
    pub reference: PathBuf,

    /// Input misincorporation file
    /// 
    /// Path leading to an input MapDamage-v2 misincorporation file, obtained from the input alignment file. This file is usually located in the output folder of MapDamage and is simply named 'misincorporations.txt' 
    /// This file provides with strand-specific PMD frequency estimates, which are then used by pmd-mask to compute the pb thresholds at which masking should be performed.
    /// 
    /// Note that this file MUST have been obtained using the same input bam file as the one used with this program. Applying pmd-mask using a misincorporation file from a different sample may result with imprecise thresholds estimates, and thus either create (over|under)correction. Note that pmd-mask does not, and most probably cannot check that the two files are consistent.
    #[arg(short, long)]
    pub misincorporation: PathBuf,

    /// Set the verbosity level (-v|-vv|-vvv)
    /// 
    /// Set the verbosity level of this program. Multiple levels available, depending on the number of calls to this argument.  
    /// 
    /// Levels: -v (Info) | -vv (Debug) | -vvv (Trace)  
    /// 
    /// Warn level is active no matter what, but can be disabled by using the --quiet argument
    #[arg(short='v', long, action = ArgAction::Count)]
    pub verbose: u8,

    /// Disable warnings.
    /// 
    /// By defaults, warnings are emmited and redirected to stderr no matter the verbosity level.
    /// Use this argument to disable Warnings. Only errors will be displayed. Note that using this argument will also have the effect of removing verbosity.
    #[clap(short='q', long)]
    pub quiet: bool,

    /// Set additional worker threads
    ///
    /// Set the number of additional (de)compression threads. Setting this value to zero will preempt all the available cores.
    /// 
    /// Note that this requires instantiating a Threadpool, which also lives in its own thread. Thus setting this to 8 will actually instantiate 10 threads: (1 runtime thread, 8 worker threads, 1 threadpool manager)
    /// 
    /// # General guidelines:  
    /// 
    ///     - Set this value to 2 if your output is in SAM format (1 thread for the reader, 1 for the writer). Setting this value to anything greater than two is pointless in this case as htslib does not require compression.
    /// 
    ///     - If your output is in BAM/CRAM format, set this value to your liking. Note that anything greater than 8 threads is rarely beneficial for compressing bam files., especially for low coverage samples.
    /// 
    /// 
    #[clap(short='@', long, default_value("1"), value_parser(parse_threads))]
    pub threads: u32
}




#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn threads_parser() {

        // 0 means all cores
        match parse_threads("0") {
            Ok(nprocs) => assert!(nprocs == num_cpus::get() as u32),
            Err(e)     => panic!("{e}")
        }

        // 1..u32::MAX means the provided number
        for n in 1..256 {
            match parse_threads(n.to_string().as_str()) {
                Ok(nprocs) => assert!(nprocs == n),
                Err(e)     => panic!("{e}")
            }
        }

        // Negative values and such are considered errors.
        for n in -256..0 {
            assert!(parse_threads(n.to_string().as_str()).is_err())
        }
    }

    #[test]
    fn output_format_parser() {

        for sam in ["SAM", "Sam", "sam", "S", "s"] {
            assert!(matches!(parse_output_fmt(sam), Ok(bam::Format::Sam)));
        }
        for bam in ["BAM", "Bam", "bam", "B", "b"] {
            assert!(matches!(parse_output_fmt(bam), Ok(bam::Format::Bam)));
        }
        
        for cram in ["CRAM", "Cram", "cram", "C", "c"] {
            assert!(matches!(parse_output_fmt(cram), Ok(bam::Format::Cram)));
        }

        for invalid  in ["CRAME", "BLAM", "Same", "vam", "v", "1", "2", "3", "ba", "sa", "Cra"] {
            assert!(parse_output_fmt(invalid).is_err())
        }
    }


}

