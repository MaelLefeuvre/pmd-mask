use std::collections::HashMap;
use std::fmt::{Display, Formatter, self};
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::path::PathBuf;
use std::str::FromStr;

use rust_htslib::{bam, bam::Read};
use rust_htslib::faidx;

use clap::{Parser, Args};
use serde::Deserialize;
use csv::{Error, ReaderBuilder};
use num_cpus;

mod logger;

use log::{info, warn, debug, trace};

#[derive(Parser, Debug)]
#[command(name="pmd-mask", author, version, about, long_about = None)]
#[clap(propagate_version = true)]
struct Cli {
    /// Optional name to operate on
    #[arg(short, long, default_value("0.01"))]
    threshold: f32,

    /// Input bam file, on which pmd-masking should be performed. When unspecified, the program will look for standard input.
    #[arg(short, long, required(false))]
    bam: Option<PathBuf>,

    /// Output bam file. If not specified, print to stdout.
    #[arg(short, long, required(false))]
    output: Option<PathBuf>,

    /// Output format. (Sam, Bam, CRAM). If unspecified:
    /// - Sam if stdout
    /// - original 
    #[arg(short='O', long, required(false))]
    output_fmt: Option<String>,

    /// Set the compression level. 9 if unspecified
    #[arg(long, default_value("9"))]
    compress_level: u32,

    /// Path to a fasta indexed reference genome.
    #[arg(short='f', long)]
    reference: PathBuf,

    /// Path to a misincorporation file
    #[arg(short, long)]
    misincorporation: PathBuf,

    /// Set the verbosity level (-v -vv -vvv -vvv)
    /// 
    /// Set the verbosity level of this program. Multiple levels available, depending on the number of calls to this argument.
    /// -v (Info) | -vv (Debug) | -vvv (Trace)
    /// Warn level is active no matter what, but can be disabled by using the --quiet argument
    #[arg(short='v', long, action = clap::ArgAction::Count)]
    verbose: u8,

    /// Disable warnings.
    /// 
    /// By defaults, warnings are emmited and redirected to stderr no matter the verbosity level.
    /// Use this argument to disable warnings. Only errors will be displayed
    /// Note that using this argument will also have the effect of removing verbosity.
    #[clap(short='q', long)]
    quiet: bool,

    /// Set threads threads 
    ///
    /// Set the number of additional decompression threads. Setting this value to zero will preempt all the available cores.
    #[clap(short='@', long, default_value("1"))]
    threads: u32
}


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
// Warn user if misincorporation frequencies are Inf, Nan, or were computed from empty C.


#[derive(Debug, Deserialize, PartialEq, Eq, Hash)] 
enum Strand {
    #[serde(rename = "+")]
    Forward,
    #[serde(rename = "-")]
    Reverse
}

impl Strand {
    pub fn symbol(&self) -> &str {
        match self {
            Self::Forward => "+",
            Self::Reverse => "-",
        }
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.symbol().fmt(f)
    }
}

impl FromStr for Strand {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Self::Forward),
            "-" => Ok(Self::Reverse),
            other => Err(format!("Incalid Strand. Got {other}"))

        }
    }
}

#[derive(Debug, Deserialize, PartialEq, Eq, Hash)]
enum Orientation {
    #[serde(rename = "3p")]
    ThreePrime,
    #[serde(rename = "5p")]
    FivePrime,
}

impl Display for Orientation {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let repr = match self {
            Self::ThreePrime => "3p",
            Self::FivePrime  => "5p",
        };
        repr.fmt(f)
    }
}

#[derive(Debug, Deserialize, PartialEq, Eq, Clone, Hash)]
struct Chromosome(String);

impl Display for Chromosome {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Debug, Deserialize, Hash, PartialEq, Eq)]
struct Position(usize);

impl Display for Position {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[derive(Debug, Deserialize)]
struct MisincorporationEntry {
    #[serde(rename(deserialize = "Chr"))] chromosome: Chromosome,
    #[serde(rename(deserialize = "End"))] end: Orientation,
    #[serde(rename(deserialize = "Std"))] strand: Strand,
    #[serde(rename(deserialize = "Pos"))] position: Position, 
    #[serde(rename(deserialize = "C"))]   c_counts: usize,
    #[serde(rename(deserialize = "G"))]   g_counts: usize,
    #[serde(rename(deserialize = "C>T"))] c_to_t: usize, 
    #[serde(rename(deserialize = "G>A"))] g_to_a: usize,
}

impl MisincorporationEntry {
    /// Return the relative C>T frequency 
    /// This is computed as the number of observed C>T, divided by the number of observed C.
    fn c_to_t_freq(&self) -> f32 {
        self.c_to_t as f32 / self.c_counts as f32
    }

    /// Return the relative G>A frequency.
    /// This is computed as the number of observed G>A, divided by the number of observed G.
    fn g_to_a_freq(&self) -> f32 {
        self.g_to_a as f32 / self.g_counts as f32
    }

    /// Return the frequency we're really interested in:
    /// If this entry is a 5p -> return C>T relative frequency (see [c_to_freq()])
    /// If this entry is a 3p -> return G>A relative frequency (see [g_to_a_freq()])
    fn target_freq(&self) -> f32 {
        match self.end {
            Orientation::FivePrime  => self.c_to_t_freq(),
            Orientation::ThreePrime => self.g_to_a_freq(),
        }
    }
}

impl Display for MisincorporationEntry {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Chr: {:<15} End: {: <3} Strand: {: <2} Pos: {: <4} C: {: <9} G: {: <9} C>T: {: <9} ({: <9.6}) G>A: {: <9} ({: <9.6})", 
            self.chromosome,
            self.end,
            self.strand,
            self.position,
            self.c_counts,
            self.g_counts,
            self.c_to_t,
            self.c_to_t_freq(),
            self.g_to_a,
            self.g_to_a_freq(),
        )
    }
}

#[derive(Debug, Hash, PartialEq, Eq)]
struct MisincorporationRecord {
    chromosome: Chromosome,
    strand: Strand
}

impl Display for MisincorporationRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.chromosome, self.strand)
    }
}



#[derive(Debug)]
struct MisincorporationThreshold { inner: HashMap<Orientation, Position>}

impl Default for MisincorporationThreshold {
    fn default() -> Self {
        Self { inner: HashMap::with_capacity(2) }
    }

}

impl MisincorporationThreshold {
    
    fn set_threshold(&mut self, orientation: Orientation, position: Position) {
        self.inner.insert(orientation, position);
    }

    fn validate(&self) -> Result<(), String> {
        if self.inner.len() != 2 {
            return Err("Invalid Threshold".to_string())
        }
        Ok(())
    }

    fn get_threshold(&self, orientation: &Orientation) -> &Position {
        self.inner.get(orientation).unwrap()
    }

}


fn main() {
    // ---- parse arguments
    let args = Cli::parse();

    // ---- initialize logger
    logger::Logger::init(if args.quiet {0} else {args.verbose + 1} );
    debug!("Provided Command line Arguments:{args:#?}");

    // ---- define threads:
    let threads = if args.threads == 0 {num_cpus::get() as u32 } else {args.threads};
    info!("Setting threads to {threads}");
    let pool = rust_htslib::tpool::ThreadPool::new(threads).unwrap();

    
    // ---- Read Misincorporation file
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .comment(Some(b'#'))
        .from_path(args.misincorporation)
        .unwrap(); 

    let mut skip_chromosome: Option<(&Chromosome, &Orientation, &Strand)> = None; 
    let mut threshold_positions = Vec::with_capacity(32);
    'nextline: for result in reader.deserialize::<MisincorporationEntry>() {
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

    // ---- Validate misincorporation and issue warnings
    let mut abnormal_frequencies = String::new();
    let threshold_positions: Vec<MisincorporationEntry> = threshold_positions.into_iter().filter_map(|pos| {
        if pos.target_freq().is_normal() {
            Some(pos)
        } else {
            abnormal_frequencies += &format!("{pos}\n");
            None
        }
    }).collect();
    abnormal_frequencies.pop();


    warn!("Found abnormal misincorporation frequencies (NaN, Infinite values, etc.). Masking will be disabled for those:\n{abnormal_frequencies}");
    debug!("Using the following positions as threshold for masking:\n{}",
        threshold_positions.iter().fold(String::new(), |acc, val| {
            acc + &format!("{val}\n")
        })
    );

    // ---- Convert to a more usable format. This is a crash test after all...
    let mut thresholds = HashMap::new(); //HashMap<MisincorporationRecord, MisincorporationThreshold> = HashMap::new();
    for position in threshold_positions {
        let record = MisincorporationRecord{chromosome: position.chromosome, strand: position.strand};
        thresholds.entry(record).or_insert(MisincorporationThreshold::default()).set_threshold(position.end, position.position);
    }

    // ---- Validate thresholds: Each one must contain a hashmap w/ two values (one for 3p, the other for 5p)
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


    let mut out_record = rust_htslib::bam::Record::new();
    // ---- Loop along input records
    'record: while let Some(result) = bam.read(&mut record) {
        result.unwrap();

        // ---- Get chromosome and strand info of this record.
        // EXPENSIVE: converting tid to string.
        // @ TODO: map Misincorporation chromosome names to tid. once, before looping.
        let chromosome = Chromosome(std::str::from_utf8(header_view.tid2name(record.tid() as u32)).unwrap().to_string());
        let strand = Strand::from_str(record.strand().strand_symbol()).unwrap();
        let current_record = MisincorporationRecord{chromosome, strand};
        

        // ---- Get relevant misincorporation frequency:
        let Some(relevant_thresholds) = thresholds.get(&current_record) else {
            trace!("{current_record} Not found in threshold dictionary. Skipping");
            continue 'record
        };

        // EXPENSIVE: converting sequence to utf8
        // @ TODO: use this?: static DECODE_BASE: &[u8] = b"=ACMGRSVTWYHKDBN";
        let sequence = record.seq().as_bytes(); // That's an expensive operation.
        let sequence = std::str::from_utf8(&sequence).unwrap();

        // ---- Get the reference's position
        let position = record.pos();
        let end = record.seq_len();
        let reference = reference.fetch_seq(&current_record.chromosome.0, position as usize, position as usize + end - 1).unwrap();
        trace!("Relevant thresholds: {relevant_thresholds:#?}");
        trace!("{} {position} {} {sequence:?}", current_record.chromosome, current_record.strand);
        
        let mut new_seq = record.seq().as_bytes(); // EXPENSIVE: Allocation
        let mut new_quals = record.qual().to_vec(); // EXPENSIVE: Allocation

        trace!("Reference: {}", std::str::from_utf8(reference).unwrap());
        trace!("Sequence : {sequence}");

        // ---- Mask 5p' positions
        'mask5p: for i in 0..relevant_thresholds.get_threshold(&Orientation::FivePrime).0 {
            if i >= record.seq_len() {
                break 'mask5p
            }

            let decoded_nuc = unsafe { record.seq().decoded_base_unchecked(i) };
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
            
            let decoded_nuc = unsafe { record.seq().decoded_base_unchecked(i) };
            if reference[i] == b'G' {
                new_seq[i] = b'N';
                new_quals[i] = 0;
            }

        }

        trace!("Masked   : {}", std::str::from_utf8(&new_seq).unwrap());


        let mut out_record = record.clone();
        out_record.set(record.qname(), Some(&record.cigar().take()), &new_seq, &new_quals);
        //println!("{record:#?}");
        writer.write(&out_record).unwrap();
    }

}
