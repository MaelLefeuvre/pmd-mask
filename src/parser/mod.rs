use std::{path::PathBuf};

use clap::{Parser, ArgAction, ColorChoice};



#[derive(Parser, Debug)]
#[command(name="pmd-mask", author, version, about, long_about = None, color=ColorChoice::Always)]
#[clap(propagate_version = true)]
pub struct Cli {
    /// Optional name to operate on
    #[arg(short, long, default_value("0.01"))]
    pub threshold: f32,

    /// Input bam file, on which pmd-masking should be performed. When unspecified, the program will look for standard input.
    #[arg(short, long, required(false))]
    pub bam: Option<PathBuf>,

    /// Output bam file. If not specified, print to stdout.
    #[arg(short, long, required(false))]
    pub output: Option<PathBuf>,

    /// Output format. (Sam, Bam, CRAM). If unspecified:
    /// - Sam if stdout
    /// - original 
    #[arg(short='O', long, required(false))]
    pub output_fmt: Option<String>,

    /// Set the compression level. 9 if unspecified
    #[arg(long, default_value("9"))]
    pub compress_level: u32,

    /// Path to a fasta indexed reference genome.
    #[arg(short='f', long)]
    pub reference: PathBuf,

    /// Path to a misincorporation file
    #[arg(short, long)]
    pub misincorporation: PathBuf,

    /// Set the verbosity level (-v -vv -vvv -vvv)
    /// 
    /// Set the verbosity level of this program. Multiple levels available, depending on the number of calls to this argument.
    /// -v (Info) | -vv (Debug) | -vvv (Trace)
    /// Warn level is active no matter what, but can be disabled by using the --quiet argument
    #[arg(short='v', long, action = ArgAction::Count)]
    pub verbose: u8,

    /// Disable warnings.
    /// 
    /// By defaults, warnings are emmited and redirected to stderr no matter the verbosity level.
    /// Use this argument to disable warnings. Only errors will be displayed
    /// Note that using this argument will also have the effect of removing verbosity.
    #[clap(short='q', long)]
    pub quiet: bool,

    /// Set threads threads 
    ///
    /// Set the number of additional decompression threads. Setting this value to zero will preempt all the available cores.
    #[clap(short='@', long, default_value("0"))]
    pub threads: u32
}

