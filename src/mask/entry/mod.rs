use std::fmt::{self, Display, Formatter};

//use genome::{ChrName, Strand};

use crate::genome::{ChrName, Strand};

use rust_htslib::bam::{HeaderView, Record};

mod error;
pub use error::MaskEntryError;

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct MaskEntry {
    pub chromosome: ChrName,
    pub strand    : Strand
}

impl Display for MaskEntry {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let out = format!("{} {}", self.chromosome, self.strand);
        out.fmt(f)
    }
}

impl MaskEntry {
    pub fn from_htslib_record(header_view: &HeaderView, record: &mut Record) -> Result<Self, MaskEntryError> {
        use MaskEntryError::ParseFromHtslib;
        let chromosome = ChrName::from_htslib_record(header_view, record).map_err(|e| ParseFromHtslib(Box::new(e)))?;
        let strand     = Strand::from_htslib_record(record).map_err(|e| ParseFromHtslib(Box::new(e)))?;
        Ok(MaskEntry{chromosome, strand})
    }
}

 
#[cfg(test)]
pub mod dummy {
    use super::*;
    use rust_htslib::bam::{Header, header::HeaderRecord};

    /// Generates a dummy bam record and header, using the provided Strand and position.
    pub fn dummy_bam(strand: Strand, pos: usize) -> (Header, Record){
        let mut header = Header::new();
        let mut chr_record = HeaderRecord::new("SQ".as_bytes());
        chr_record.push_tag("SN".as_bytes(), &"chr1".to_string());
        chr_record.push_tag("LN".as_bytes(), &249250621);
        header.push_record(&chr_record);

        let mut record = Record::new();
        record.set("*".as_bytes(), None, &[b'A', b'T', b'C', b'G'], &[37, 37, 37, 37]);
        record.set_tid(0);
        record.set_pos(pos as i64);
        if strand == Strand::Reverse {
            record.set_reverse();
        };

        println!("Record strand symbol: {}", record.strand().strand_symbol());

        (header.clone(), record.clone())
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use super::dummy::dummy_bam;
    #[test]
    fn from_htslib_record() {
        let (header, mut record) = dummy_bam(Strand::Forward, 10_000);
        let header_view          = HeaderView::from_header(&header);
        let mask                 = MaskEntry::from_htslib_record(&header_view, &mut record);

        let want = MaskEntry {chromosome: ChrName::new("chr1"), strand: Strand::Forward};
        assert_eq!(mask.expect("Invalid Mask"), want)

    }

    #[test]
    fn display() {
    
    let (header, mut record) = dummy_bam(Strand::Forward, 10_000);
    let header_view          = HeaderView::from_header(&header);

    let mask = MaskEntry::from_htslib_record(&header_view, &mut record).expect("Invalid Header");

    let want = "chr1 +";
    assert_eq!(want, format!("{mask}"));
    assert_eq!(format!("{want:-^15}"), format!("{mask:-^15}"));
    }
}