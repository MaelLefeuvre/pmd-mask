# version 0.1.1 (2023-03-10)
## Features
- Restructured code architecture
- Implemented basic error handling using `thiserror` and `anyhow`. Might attempt to go from `anyhow` to `miette` in the future, just for the 'gist of it. 
- Masking is now more aggressive, and will be applied on the whole read length in case of an invalid or missing MaskingThreshold.

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