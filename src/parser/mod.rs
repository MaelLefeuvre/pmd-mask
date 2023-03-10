use std::{path::PathBuf};

use clap::{Parser, ArgAction, ColorChoice};


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
    #[arg(short='O', long, required(false))]
    pub output_fmt: Option<String>,

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
    #[clap(short='@', long, default_value("1"))]
    pub threads: u32
}

