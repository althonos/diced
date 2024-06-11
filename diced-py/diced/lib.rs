#![doc = include_str!("../README.md")]

extern crate diced;
extern crate pyo3;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use pyo3::types::PySlice;
use pyo3::types::PyString;

/// A sequence region.
#[pyclass(module = "diced.lib", frozen, subclass)]
pub struct Region {
    region: diced::Region<PyBackedStr>,
}

#[pymethods]
impl Region {
    #[new]
    pub fn __new__<'py>(
        py: Python<'py>,
        sequence: PyBackedStr,
        start: usize,
        end: usize,
    ) -> PyResult<PyClassInitializer<Self>> {
        if start > end || start > sequence.len() || end > sequence.len() {
            let s = PySlice::new_bound(py, start as isize, end as isize, 1);
            return Err(PyIndexError::new_err(s.to_object(py)));
        }
        Ok(Region {
            region: diced::Region::new(sequence, start, end),
        }
        .into())
    }

    /// `int`: The start coordinate of the region (zero-based).
    #[getter]
    pub fn start(&self) -> usize {
        self.region.start()
    }

    /// `int`: The end coordinate of the region (zero-based, exclusive).
    #[getter]
    pub fn end(&self) -> usize {
        self.region.end()
    }

    /// Get the sequence region as a string.
    pub fn __str__<'py>(&self, py: Python<'py>) -> Bound<'py, PyString> {
        PyString::new_bound(py, self.region.as_str())
    }
}

/// A CRISPR repeat.
#[pyclass(module="diced.lib", extends=Region)]
pub struct Repeat {}

#[pymethods]
impl Repeat {
    #[new]
    pub fn __new__<'py>(
        py: Python<'py>,
        sequence: PyBackedStr,
        start: usize,
        end: usize,
    ) -> PyResult<PyClassInitializer<Self>> {
        Region::__new__(py, sequence, start, end).map(|r| r.add_subclass(Repeat {}))
    }
}

/// A list of repeats inside a CRISPR region.
#[pyclass(module = "diced.lib", sequence)]
pub struct Repeats {
    crispr: Py<Crispr>,
}

#[pymethods]
impl Repeats {
    pub fn __len__<'py>(&self, py: Python<'py>) -> usize {
        self.crispr.borrow(py).crispr.len()
    }

    pub fn __getitem__<'py>(&self, py: Python<'py>, index: usize) -> PyResult<Py<Repeat>> {
        self.crispr
            .bind(py)
            .borrow()
            .crispr
            .repeats()
            .nth(index)
            .ok_or(PyIndexError::new_err(index))
            .and_then(|region| {
                Py::new(
                    py,
                    PyClassInitializer::from(Region { region }).add_subclass(Repeat {}),
                )
            })
    }
}

/// A CRISPR spacer.
#[pyclass(module="diced.lib", extends=Region)]
pub struct Spacer {}

#[pymethods]
impl Spacer {
    #[new]
    pub fn __new__<'py>(
        py: Python<'py>,
        sequence: PyBackedStr,
        start: usize,
        end: usize,
    ) -> PyResult<PyClassInitializer<Self>> {
        Region::__new__(py, sequence, start, end).map(|r| r.add_subclass(Spacer {}))
    }
}

/// A list of spacers inside a CRISPR region.
#[pyclass(module = "diced.lib", sequence)]
pub struct Spacers {
    crispr: Py<Crispr>,
}

#[pymethods]
impl Spacers {
    pub fn __len__<'py>(&self, py: Python<'py>) -> usize {
        self.crispr.borrow(py).crispr.len().saturating_sub(1)
    }

    pub fn __getitem__<'py>(&self, py: Python<'py>, index: usize) -> PyResult<Py<Spacer>> {
        self.crispr
            .bind(py)
            .borrow()
            .crispr
            .spacers()
            .nth(index)
            .ok_or(PyIndexError::new_err(index))
            .and_then(|region| {
                Py::new(
                    py,
                    PyClassInitializer::from(Region { region }).add_subclass(Spacer {}),
                )
            })
    }
}

/// A CRISPR region in a nucleotide sequence.
#[pyclass(module = "diced.lib")]
pub struct Crispr {
    crispr: diced::Crispr<PyBackedStr>,
}

#[pymethods]
impl Crispr {
    /// `int`: The start coordinate of the CRISPR region (zero-based).
    #[getter]
    pub fn start(&self) -> usize {
        self.crispr.start()
    }

    /// `int`: The end coordinate of the CRISPR region (zero-based, exclusive).
    #[getter]
    pub fn end(&self) -> usize {
        self.crispr.end()
    }

    /// `~diced.Repeats`: The list of repeats inside the CRISPR region.
    #[getter]
    pub fn repeats(slf: Py<Self>) -> Repeats {
        Repeats { crispr: slf }
    }

    /// `~diced.Spacers`: The list of spacers inside the CRISPR region.
    #[getter]
    pub fn spacers(slf: Py<Self>) -> Spacers {
        Spacers { crispr: slf }
    }

    pub fn __len__(&self) -> usize {
        self.crispr.len()
    }

    pub fn __str__<'py>(&self, py: Python<'py>) -> Bound<'py, PyString> {
        PyString::new_bound(py, self.crispr.to_region().as_str())
    }
}

/// A scanner for iterating on the CRISPR regions of a genome.
#[pyclass(module = "diced.lib")]
pub struct Scanner {
    scanner: diced::Scanner<PyBackedStr>,
}

#[pymethods]
impl Scanner {
    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    /// Return the next CRISPR region, if any.
    ///
    /// Returns:
    ///     `~diced.Crispr`: The next CRISPR region in the sequence.
    ///
    /// Raises:
    ///     `StopIteration`: When the end of the sequence has been reached
    ///         without finding new CRISPR regions.
    ///
    fn __next__<'py>(&mut self, py: Python<'py>) -> PyResult<Option<Crispr>> {
        match py.allow_threads(move || self.scanner.next()) {
            Some(crispr) => Ok(Some(Crispr { crispr })),
            None => Ok(None),
        }
    }

    /// `str`: The genomic sequence being scanned.
    #[getter]
    fn sequence<'py>(&self, py: Python<'py>) -> Py<PyAny> {
        self.scanner.sequence().clone().to_object(py)
    }
}

/// Scan a genome sequence for CRISPRs repeats.
///
/// Arguments:
///     sequence (`str`): A string containing the genomic sequence to build
///         a scanner for.
///
/// Returns:
///     `~diced.Scanner`: A scanner yielding CRISPRs in the given contig.
///
#[pyfunction]
pub fn scan(sequence: PyBackedStr) -> PyResult<Scanner> {
    let builder = diced::ScannerBuilder::new();
    let scanner = builder.scan(sequence);
    Ok(Scanner { scanner })
}

/// PyO3 bindings to ``diced``, a library for CRISPRs detection.
///
/// Diced is re-implementation of MinCED, a method developed by
/// `Connor T. Skennerton <https://github.com/ctSkennerton>`_ to identify
/// CRISPRs in isolate and metagenomic-assembled genomes. It was derived
/// from the CRISPR recognition tool developed by Charles Bland *et al.*.
///
/// Example:
///     Load a genome from a FASTA file using Biopython::
///
///         >>> import Bio.SeqIO
///         >>> record = Bio.SeqIO.read("Aquifex_aeolicus_VF5.fna", "fasta")
///
///     Detect CRISPR regions with Diced using the default parameters::
///
///         >>> import diced
///         >>> for crispr in diced.scan(str(record.seq[:300000])):
///         ...     print(
///         ...         crispr.start,
///         ...         crispr.end,
///         ...         len(crispr.repeats),
///         ...         crispr.repeats[0],
///         ...     )
///         156459 156767 5 GTTCCTAATGTACCGTGTGGAGTTGAAACC
///         244560 244791 4 GTTTCAACTCCACACGGTACATTAGGAAC
///         279263 279555 5 GTTTTAACTCCACACGGTACATTAGAAAC
///
#[pymodule]
#[pyo3(name = "lib")]
pub fn init(_py: Python, m: Bound<PyModule>) -> PyResult<()> {
    m.add("__package__", "diced")?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", env!("CARGO_PKG_AUTHORS").replace(':', "\n"))?;

    m.add_class::<Crispr>()?;
    m.add_class::<Region>()?;
    m.add_class::<Scanner>()?;
    m.add_class::<Repeat>()?;
    m.add_class::<Repeats>()?;
    m.add_class::<Spacer>()?;
    m.add_class::<Spacers>()?;

    m.add_function(wrap_pyfunction!(scan, &m)?)?;

    Ok(())
}
