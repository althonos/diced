#![allow(non_snake_case)]

fn test(builder: &mincer::ScannerBuilder, gff_path: &str) {
    let mut reader = std::fs::File::open("tests/data/Aquifex_aeolicus_VF5.fna")
        .map(std::io::BufReader::new)
        .map(noodles_fasta::Reader::new)
        .unwrap();
    let record = reader.records().next().unwrap().unwrap();
    let seq = std::str::from_utf8(record.sequence().as_ref()).unwrap();

    let regions = builder.scan(seq).collect::<Vec<_>>();

    let gff = std::fs::File::open(gff_path)
        .map(std::io::BufReader::new)
        .map(noodles_gff::Reader::new)
        .unwrap()
        .records()
        .map(Result::unwrap)
        .collect::<Vec<_>>();

    let repeat_regions = gff
        .iter()
        .filter(|record| record.ty() == "repeat_region")
        .collect::<Vec<_>>();

    assert_eq!(regions.len(), repeat_regions.len());
    for (actual_region, expected_region) in regions.iter().zip(repeat_regions) {
        assert_eq!(actual_region.start() + 1, expected_region.start().get());
        assert_eq!(actual_region.end(), expected_region.end().get());
        assert_eq!(
            actual_region.len(),
            expected_region.score().unwrap_or(0.0) as usize
        );
        let unit_seq = expected_region.attributes().get("rpt_unit_seq").unwrap();
        assert_eq!(actual_region.repeat(1), unit_seq.as_string().unwrap());
    }
}

#[test]
fn default() {
    test(
        &mincer::ScannerBuilder::new(),
        "tests/data/Aquifex_aeolicus_VF5.gff",
    )
}

#[test]
fn minSL12() {
    test(
        mincer::ScannerBuilder::new().min_spacer_length(12),
        "tests/data/Aquifex_aeolicus_VF5.minSL12.gff",
    )
}

#[test]
fn maxSL30() {
    test(
        mincer::ScannerBuilder::new().max_spacer_length(30),
        "tests/data/Aquifex_aeolicus_VF5.maxSL30.gff",
    )
}

#[test]
fn minRL30_maxRL40() {
    test(
        mincer::ScannerBuilder::new()
            .min_repeat_length(30)
            .max_repeat_length(40),
        "tests/data/Aquifex_aeolicus_VF5.minRL30.maxRL40.gff",
    )
}

#[test]
fn maxSL34_minSL20_minNR4_minRL20() {
    test(
        mincer::ScannerBuilder::new()
            .max_spacer_length(34)
            .min_spacer_length(20)
            .min_repeat_count(4)
            .min_repeat_length(20),
        "tests/data/Aquifex_aeolicus_VF5.maxSL34.minSL20.minNR4.minRL20.gff",
    )
}

#[test]
fn maxSL34_minSL20_minNR4_minRL20_maxRL40() {
    test(
        mincer::ScannerBuilder::new()
            .max_spacer_length(34)
            .min_spacer_length(20)
            .min_repeat_count(4)
            .min_repeat_length(20)
            .max_repeat_length(40),
        "tests/data/Aquifex_aeolicus_VF5.maxSL34.minSL20.minNR4.minRL20.maxRL40.gff",
    )
}
