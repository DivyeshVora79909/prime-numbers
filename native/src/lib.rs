use pyo3::prelude::*;

mod arithmetic;
mod multiply;
mod storage;
mod wire;

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    multiply::register(m)?;
    arithmetic::register(m)?;
    storage::register(m)?;
    wire::register(m)
}
