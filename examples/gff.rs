use seq_io::fasta::Record;

extern crate mincer;
extern crate seq_io;

fn main() -> std::io::Result<()> {
    let crispr_scanner = mincer::Scanner::default();

    let path = std::env::args().skip(1).next().ok_or(std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        "Missing required path",
    ))?;
    let reader = std::fs::File::open(&path)
        .map(seq_io::fasta::Reader::new)
        .unwrap();
    for result in reader.into_records() {
        let record = result.unwrap();

        let id = record.id().unwrap();
        let seq = std::str::from_utf8(record.seq()).unwrap();

        for (i, crispr) in crispr_scanner.scan(seq).enumerate() {
            if i == 0 {
                println!("##gff-version 3");
            }
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t.\t.\tID=CRISPR{};rpt_type=direct;rpt_family=CRISPR;rpt_unit_seq={}",
                id,
                "mincer:0.1.0",
                "repeat_region",
                crispr.start() + 1,
                crispr.end(),
                crispr.len(),
                i + 1,
                crispr.repeat(1),
            );
        }
    }

    Ok(())
}
