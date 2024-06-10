extern crate diced;
extern crate pyo3;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use pyo3::types::PyString;

/// A sequence region.
#[pyclass]
pub struct Region {
    region: diced::Region<PyBackedStr>,
}

#[pymethods]
impl Region {
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

    pub fn __str__<'py>(&self, py: Python<'py>) -> Bound<'py, PyString> {
        PyString::new_bound(py, self.region.as_str())
    }
}

#[pyclass]
pub struct Repeats {
    crispr: Py<Crispr>,
}

#[pymethods]
impl Repeats {
    pub fn __len__<'py>(&self, py: Python<'py>) -> usize {
        self.crispr.borrow(py).crispr.len()
    }

    pub fn __getitem__<'py>(&self, py: Python<'py>, index: usize) -> PyResult<Region> {
        self.crispr
            .bind(py)
            .borrow()
            .crispr
            .repeats()
            .nth(index)
            .ok_or(PyIndexError::new_err(index))
            .map(|region| Region { region })
    }
}

#[pyclass]
pub struct Spacers {
    crispr: Py<Crispr>,
}

#[pymethods]
impl Spacers {
    pub fn __len__<'py>(&self, py: Python<'py>) -> usize {
        self.crispr.borrow(py).crispr.len().saturating_sub(1)
    }

    pub fn __getitem__<'py>(&self, py: Python<'py>, index: usize) -> PyResult<Region> {
        self.crispr
            .bind(py)
            .borrow()
            .crispr
            .spacers()
            .nth(index)
            .ok_or(PyIndexError::new_err(index))
            .map(|region| Region { region })
    }
}

/// A CRISPR region in a nucleotide sequence.
#[pyclass]
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
///
/// Attributes:
///     sequence (`str` or `bytes`): The sequence of the genome being
///         scanned by the `Scanner`.
///
#[pyclass]
pub struct Scanner {
    scanner: diced::Scanner<PyBackedStr>,
}

#[pymethods]
impl Scanner {
    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __next__<'py>(&mut self, py: Python<'py>) -> PyResult<Option<Crispr>> {
        match py.allow_threads(move || self.scanner.next()) {
            Some(crispr) => Ok(Some(Crispr { crispr })),
            None => Ok(None),
        }
    }

    #[getter]
    fn sequence<'py>(&self, py: Python<'py>) -> Py<PyAny> {
        self.scanner.sequence().clone().to_object(py)
    }
}

/// Scan a genome sequence for CRISPRs repeats.
///
/// Arguments:
///     sequence (`str` or `bytes`): A string containing the genomic
///         sequence to build a scanner for.
///
/// Returns:
///     `~diced.Scanner`: A scanner returning CRISPRs in the given genome.
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
    m.add_class::<Repeats>()?;
    m.add_class::<Spacers>()?;

    m.add_function(wrap_pyfunction!(scan, &m)?)?;

    Ok(())
}