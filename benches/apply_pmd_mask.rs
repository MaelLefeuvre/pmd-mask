use criterion::{black_box, criterion_group, criterion_main, Criterion};

use pmd_mask::{Masks, apply_pmd_mask};


use rust_htslib::bam::{Read, Header, Format};
use rust_htslib::bam::Reader as BamReader;
use rust_htslib::bam::Writer as BamWriter;
use rust_htslib::faidx::Reader as FaReader;

use assert_fs::NamedTempFile;

macro_rules! test_dir {
    (bam $path:literal) => { concat!("tests/test-data/bam/dummy-MTonly/", $path) };
    (fa $path:literal)  => { concat!("tests/test-data/reference/", $path) };
}

fn bench_apply_pmd_mask(c: &mut Criterion) {


    let mut bam     = black_box(BamReader::from_path(test_dir!(bam "dummy-MTonly-1000.bam")).expect("Can't open bam Reader"));
    
    let mut fixture = black_box(NamedTempFile::new("output.bam").expect("Failed to create fixture for output bam"));    
    let mut writer  = black_box(BamWriter::from_path(&mut fixture, &Header::from_template(bam.header()), Format::Sam).expect("Can't open Sam Writer"));

    let reference   = black_box(FaReader::from_path(test_dir!(fa "/hs37d5-MTonly/hs37d5-MTonly.fa.gz")).expect("Can't open reference file"));


    let mut bench_masks = Vec::new();
    for threshold in [0.0, 0.01, 0.05, 0.5, 1.0] {
        let masks       = black_box(Masks::from_path(test_dir!(bam "misincorporation.txt"), black_box(threshold)).expect("Failed to open misincorporation file"));
        bench_masks.push((threshold, masks));
    }

    for (threshold, masks) in bench_masks {
        let bench_name = format!("apply_pmd_mask-{threshold}");
        c.bench_function(&bench_name, |bench| bench.iter(|| {        
            apply_pmd_mask(&mut bam, &reference, &masks, &mut writer).unwrap();
        }));
    }


}

criterion_group!(pmd_mask_lib, bench_apply_pmd_mask);
criterion_main!(pmd_mask_lib);