extern crate diced;
extern crate pyo3;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use pyo3::types::PyString;

#[pyclass]
pub struct Region {
    region: diced::Region<PyBackedStr>,
}

#[pymethods]
impl Region {
    #[getter]
    pub fn start(&self) -> usize {
        self.region.start()
    }

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

    #[getter]
    pub fn repeats(slf: Py<Self>) -> Repeats {
        Repeats { crispr: slf }
    }

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

    pub fn repeat(&self, index: isize) -> PyResult<Region> {
        let mut _index = index;
        if _index < 0 {
            _index += self.crispr.len() as isize;
        }
        if _index < 0 || _index >= self.crispr.len() as isize {
            return Err(PyIndexError::new_err(index));
        }
        Ok(Region {
            region: self.crispr.repeat(_index as usize),
        })
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

    fn __next__<'py>(&mut self, py: Python<'py>) -> PyResult<Option<Crispr>> {
        match py.allow_threads(move || self.scanner.next()) {
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
    m.add_class::<Region>()?;
    m.add_class::<Scanner>()?;
    m.add_class::<Repeats>()?;
    m.add_class::<Spacers>()?;

    m.add_function(wrap_pyfunction!(scan, &m)?)?;

    Ok(())
}
