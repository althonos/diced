extern crate diced;
extern crate pyo3;

use pyo3::{prelude::*, pybacked::PyBackedStr};

#[pyclass]
pub struct Crispr {
    crispr: diced::Crispr<PyBackedStr>,
}

#[pymethods]
impl Crispr {
    #[getter]
    pub fn start(&self) -> usize {
        self.crispr.start()
    }

    #[getter]
    pub fn end(&self) -> usize {
        self.crispr.end()
    }
}

#[pyclass]
pub struct Scanner {
    scanner: diced::Scanner<PyBackedStr>,
}

#[pymethods]
impl Scanner {
    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __next__(&mut self) -> PyResult<Option<Crispr>> {
        match self.scanner.next() {
            Some(crispr) => Ok(Some(Crispr { crispr })),
            None => Ok(None),
        }
    }
}

#[pyfunction]
pub fn scan(s: PyBackedStr) -> PyResult<Scanner> {
    let builder = diced::ScannerBuilder::new();
    let scanner = builder.scan(s);
    Ok(Scanner { scanner })
}

/// PyO3 bindings to ``diced``, a library for CRISPRs detection.
#[pymodule]
#[pyo3(name = "lib")]
pub fn init(_py: Python, m: Bound<PyModule>) -> PyResult<()> {
    m.add("__package__", "diced")?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", env!("CARGO_PKG_AUTHORS").replace(':', "\n"))?;

    m.add_class::<Crispr>()?;
    m.add_class::<Scanner>()?;

    m.add_function(wrap_pyfunction!(scan, &m)?)?;

    Ok(())
}
