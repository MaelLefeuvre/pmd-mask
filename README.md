# PMD-mask

Perform hard selective masking of ancient DNA deamination patterns, using the output misincorporation frequency estimates of [MapDamage](https://github.com/ginolhac/mapDamage.git).

![Build-Ubuntu](https://github.com/MaelLefeuvre/pmd-mask/workflows/Ubuntu/badge.svg) | ![Build-Ubuntu](https://github.com/MaelLefeuvre/pmd-mask/workflows/MacOS/badge.svg)


# Installation

## Dependencies

### Cargo
This project is written in [Rust](https://www.rust-lang.org/), and thus requires [cargo](https://crates.io/) for source compilation.

To install cargo:
```Bash
curl --proto '=https' --tlsv1.2 https://sh.rustup.rs -sSf | sh
```

### Compilation

1. Clone this repository
```Bash
git clone git@github.com:MaelLefeuvre/pmd-mask.git
```

2. Run the test suite from the repository's root
```Bash
cd pmd-mask && cargo test
```

2. Compile and install
```Bash
RUSTFLAGS="-Ctarget-cpu=native" cargo install --path .
```

3. `pmd-mask` should be located within `~/.cargo/bin/` and included in your PATH
```Bash
pmd-mask --help
```

# Usage
## Data requirements:

The following inputs are required to use PMD-mask:
1. An input bam file (SAM|BAM|CRAM formats are accepted). pmd-mask can either read from a file (using `-b`|`--bam`) or from the standard input, through shell piping.
2. A [MapDamage-v2](https://github.com/ginolhac/mapDamage.git)  `misinscorporation.txt` file. This file provides strand-specific PMD frequency estimates, which are used to compute the threshold at which masking should be performed. Use `-m`|`--misincorporation` to specify this input. Of course, this file must have been obtained from your input bam file to provide with a sound estimate.
3. A reference genome. This genome must of course be the same as the one used to align the aforementionned bam file. Use `-f`|`--reference` to specify the path to your reference

```Bash
pmd-mask --reference data/GRCh37/Homo_sapiens.GRCh37.dna.primary_assemby.fa --misincorporation test-sample-MD-folder/misincorporation.txt --bam ./test-sample.srt.rmdup.bam 
```

## Optional parameters:

- The PMD-frequency threshold used to apply masking can be specified with the `-t`|`--threshold` parameter (Default: `0.01`)
- The name of the output can be specified using `-o`|`--output`. When unspecified, pmd-mask will flush results to the standard output.
- The output format can be specified using `-O`|`--output-fmt`. (SAM|BAM|CRAM format accepted). When using `BAM` or `CRAM`, the compression level can be specified using `--compress-level`.
- Use `-@`|`--threads` to allocate additional cores to the program. This can speed-up the (de)compression rate of your input and output files.
- Add `-v`|`--verbose` flags to increase the verbosity. Multiple levels: `-v`: Will output general information (INFO) `-vv`: will output general information (DEBUG) `-vvv`: will output detailled debug information (TRACE). Not that warnings are still emitted, no matter the verbosity level. This behavior can be disabled using the `-q`|`--quiet` flag, which will inhibit all logging.

## A more detailled example:

1. Filter autosomes using samtools (notice the `-h` flags on this command, which is required to communicate the header information to `pmd-mask`) 
2. Apply PMD-masking with a threshold of 2%
3. Pileup this sample, while targeting the Reich 1240K Compendium positions.
4. Run [grups-rs](https://github.com/MaelLefeuvre/grups)  to compute an approximate estimate of the individual's average heterozygocity.
```
samtools view -h ./MT23/MT23.srt.rmdup.rescaled.bam {1..22} | pmd-mask -f ./hs37d5.fa -m ./MT23/misincorporation.txt -Ob --threshold 0.02 --quiet | samtools mpileup -RB -q25 -Q25 -f ./hs37d5.fa.gz -l ./v52.2_1240K_public.bed - | grups pwd-from-stin --samples 0 --self-comparison --sample-name MT23
```

