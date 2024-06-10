# 🔪🧅 Diced [![Star me](https://img.shields.io/github/stars/althonos/mincer?style=social&label=Star&maxAge=3600)](https://github.com/althonos/diced/stargazers)

*A Rust re-implementation of the [MinCED](https://github.com/ctSkennerton/minced) algorithm to Detect Instances of [CRISPRs](https://en.wikipedia.org/wiki/CRISPR) in Environmental Data.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/diced/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/diced/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/diced?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/diced/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/diced.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/diced)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/diced?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/diced)
[![AUR](https://img.shields.io/aur/version/python-diced?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-diced)
[![Wheel](https://img.shields.io/pypi/wheel/diced.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/diced/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/diced.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/diced/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/diced.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/diced/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/diced/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/diced.svg?style=flat-square&maxAge=600)](https://github.com/althonos/diced/issues)
[![Docs](https://img.shields.io/readthedocs/diced/latest?style=flat-square&maxAge=600)](https://diced.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/diced/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/diced?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/diced)

## 🗺️ Overview

MinCED is a method developed by [Connor T. Skennerton](https://github.com/ctSkennerton) 
to identify [Clustered Regularly Interspaced Short Palindromic Repeats (CRISPRs)](https://en.wikipedia.org/wiki/CRISPR) 
in isolate and metagenomic-assembled genomes. It was derived from the CRISPR 
Recognition Tool [\[1\]](#ref1). It uses a fast scanning algorithm to identify
candidate repeats, combined with an extension step to find maximally spanning
regions of the genome that feature a CRISPR repeat.

Diced is a Rust reimplementation of the MinCED method, using the original
Java code as a reference. It produces exactly the same results as MinCED,
corrects some bugs, and is much faster. The Diced implementation is available 
as a Rust library for convenience.

### 📋 Features

- **library interface**: The Rust implementation is written as library to facilitate 
  reusability in other projects. It is used to implement a Python library using
  PyO3 to generate a native extension.
- **single dependency**: Pyrodigal is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  Prodigal binary being present on the end-user machine.
- **zero-copy**: The `Scanner` which iterates over candidate CRISPRs is zero-copy if 
  provided with a simple `&str` reference, but it also supports data behind smart 
  pointers such as `Rc<str>` or `Arc<str>`.
- **fast string matching**: The Java implementation uses a handwritten implementation 
  of the [Boyer-Moore algorithm](https://en.wikipedia.org/wiki/Boyer%E2%80%93Moore_string-search_algorithm)[\[2\]](#ref2), while the Rust 
  implementation uses the `str::find` method of the standard library, which 
  implements the [Two-way algorithm](https://en.wikipedia.org/wiki/Two-way_string-matching_algorithm)[\[3\]](#ref3). In addition, the [`memchr`](https://crates.io/crates/memchr) crate can be used as a fast SIMD-capable 
  implementation of the `memmem` function.

## 💡 Example

Diced supports any sequence in string format.

```python
import Bio.SeqIO
import diced

record = Bio.SeqIO.read("diced/tests/data/Aquifex_aeolicus_VF5.fna", "fasta")
sequence = str(record.seq)

for crispr in diced.scan(sequence):
    print(
        crispr.start,
        crispr.end,
        len(crispr.repeats),
        crispr.repeats[0],
    )
```

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/diced/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

<!-- ### 🏗️ Contributing

Contributions are more than welcome! See [`CONTRIBUTING.md`](https://github.com/althonos/diced/blob/master/CONTRIBUTING.md) for more details. -->

## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/diced/blob/master/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.

## ⚖️ License

This library is provided under the open-source
[GPLv3 license](https://choosealicense.com/licenses/gpl-3.0/), or later. 
The code for this implementation was derived from the 
[MinCED source code](https://github.com/ctSkennerton/minced), which is 
available under the GPLv3 as well.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [original MinCED authors](https://github.com/ctSkennerton). It was developed 
by [Martin Larralde](https://github.com/althonos/) during his PhD project at 
the [Leiden University Medical Center](https://www.lumc.nl/en/) in the 
[Zeller team](https://github.com/zellerlab).*

## 📚 References

- <a id="ref1">\[1\]</a> Bland, C., Ramsey, T. L., Sabree, F., Lowe, M., Brown, K., Kyrpides, N. C., & Hugenholtz, P. (2007). 'CRISPR recognition tool (CRT): a tool for automatic detection of clustered regularly interspaced palindromic repeats'. BMC bioinformatics, 8, 209. [PMID:17577412](https://pubmed.ncbi.nlm.nih.gov/17577412/) [doi:10.1186/1471-2105-8-209](https://doi.org/10.1186/1471-2105-8-209).
- <a id="ref2">\[2\]</a> Boyer, R. S. and & Moore, J. S. (1977). 'A fast string searching algorithm'. Commun. ACM 20, 10 762–772. [doi:10.1145/359842.359859](https://doi.org/10.1145/359842.359859)
- <a id="ref3">\[3\]</a> Crochemore, M. & Perrin, D. (1991). 'Two-way string-matching'. J. ACM 38, 3, 650–674. [doi:10.1145/116825.116845](https://doi.org/10.1145/116825.116845)

  


  


