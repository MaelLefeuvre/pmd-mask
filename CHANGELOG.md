# version 0.1.3 (2023-03-14)
## Features
**A More trustworthy codebase :** ~91.86% line coverage at this point.
- Basic suite of Unit tests (PRs #4 #6 #9)
- Basic suite of integration tests (PRs #6 and #9)
- Basic Continuous Integration workflow for Ubuntu and MacOS (PRs #5 #7)
- Some very basic doctests within certain user libraries. 
- Better separation between libraries and the main logic.

## Known bugs
- [Uncaught rust_htslib exception when the specified --reference points to an invalid reference file.](https://github.com/MaelLefeuvre/pmd-mask/issues/8)

# version 0.1.2 (2023-03-10)
## Features
 - `--help` is now a bit more verbose.
## Bugfixes
 - Fixed threadpool construction and cpu_count (#2)

# version 0.1.1 (2023-03-10)
## Features
- Restructured code architecture
- Implemented basic error handling using `thiserror` and `anyhow`. Might attempt to go from `anyhow` to `miette` in the future, just for the 'gist of it. 
- Masking is now more aggressive, and will be applied on the whole read length in case of an invalid or missing MaskingThreshold.
- Command line interface `--help` message now supports linewrap.

## Bugfixes
- Fixed threshold validation. (#1)
- Fixed overflowing error when the `3p G>A` threshold length was greater than the read length

## Known Bugs:
- (#2) Threadpool construction and cpu count is broken + thread preemption appears under-utilized.
## Documentation
- **very** light documentation. Need to take some time to add doctests, now that we have a semblance of structure within the source directory.

# version 0.1.0 (2023-03-06)
## Features
- An (Absolutely barbaric) prototype.

## Known bugs
- (#1) Threshold validation falsely triggers if the user defined threshold is to low. This is because we currently have no way of distinguishing 'Invalid thresholds' from those which were not found (and thus implies we should apply masking on the whole sequence.)
