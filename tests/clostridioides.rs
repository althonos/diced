#![allow(non_snake_case)]

fn test(builder: &mincer::ScannerBuilder, gff_path: &str) {
    let mut reader = std::fs::File::open("tests/data/NZ_CP019870.1.fna")
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
        // println!();
        // for repeat in 0..actual_region.len() {
        //     println!("{}", actual_region.repeat(1).as_str());
        // }

        assert_eq!(
            actual_region.len(),
            expected_region.score().unwrap_or(0.0) as usize
        );

        assert_eq!(
            actual_region.start() + 1,
            expected_region.start().get(),
            "start coordinates differ"
        );
        assert_eq!(
            actual_region.end(),
            expected_region.end().get(),
            "end coordinates differ"
        );

        let unit_seq = expected_region.attributes().get("rpt_unit_seq").unwrap();
        assert_eq!(
            actual_region.repeat(1).as_str(),
            unit_seq.as_string().unwrap()
        );
    }
}

#[test]
fn default() {
    test(
        &mincer::ScannerBuilder::new(),
        "tests/data/NZ_CP019870.1.gff",
    )
}

#[test]
fn maxRL36() {
    test(
        &mincer::ScannerBuilder::new().max_repeat_length(36),
        "tests/data/NZ_CP019870.1.maxRL36.gff",
    )
}

#[test]
fn minRL30() {
    test(
        &mincer::ScannerBuilder::new().min_repeat_length(30),
        "tests/data/NZ_CP019870.1.minRL30.gff",
    )
}

#[test]
fn minRL30_minSL10_maxSL60() {
    test(
        &mincer::ScannerBuilder::new()
            .min_repeat_length(30)
            .min_spacer_length(10)
            .max_spacer_length(60),
        "tests/data/NZ_CP019870.1.minRL30.minSL10.maxSL60.gff",
    )
}

#[test]
fn minRL36_minSL10_maxSL60() {
    test(
        &mincer::ScannerBuilder::new()
            .min_repeat_length(36)
            .min_spacer_length(10)
            .max_spacer_length(60),
        "tests/data/NZ_CP019870.1.minRL36.minSL10.maxSL60.gff",
    )
}

#[test]
fn minSL10_maxSL60() {
    test(
        &mincer::ScannerBuilder::new()
            .min_spacer_length(10)
            .max_spacer_length(60),
        "tests/data/NZ_CP019870.1.minSL10.maxSL60.gff",
    )
}
