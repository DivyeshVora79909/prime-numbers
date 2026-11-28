//! Exact arithmetic on ordinary-coordinate signed base-two rows.
//! Grouping nearby powers before convolution saves work while keeping every
//! intermediate bounded. No coefficient ever holds the represented integer.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule};

const MODULUS: u32 = 998_244_353; // 119 × 2^23 + 1
const GENERATOR: u32 = 3;
const MAX_SPAN: usize = 1 << 23;
const DIRECT_WIDTH: usize = 24;

fn mod_mul(a: u32, b: u32) -> u32 {
    (u64::from(a) * u64::from(b) % u64::from(MODULUS)) as u32
}

fn mod_pow(mut base: u32, mut exponent: u32) -> u32 {
    let mut result = 1;
    while exponent > 0 {
        if exponent & 1 != 0 {
            result = mod_mul(result, base);
        }
        base = mod_mul(base, base);
        exponent >>= 1;
    }
    result
}

fn transform(values: &mut [u32], roots: &[u32]) {
    let n = values.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            values.swap(i, j);
        }
    }
    let mut width = 2;
    while width <= n {
        let stride = n / width;
        for block in values.chunks_exact_mut(width) {
            let (left, right) = block.split_at_mut(width / 2);
            for (i, (lo, hi)) in left.iter_mut().zip(right).enumerate() {
                let a = *lo;
                let b = mod_mul(*hi, roots[i * stride]);
                let sum = a + b;
                *lo = if sum >= MODULUS { sum - MODULUS } else { sum };
                *hi = if a >= b { a - b } else { a + MODULUS - b };
            }
        }
        width *= 2;
    }
}

fn pack(row: &[u8], width: usize) -> Vec<i32> {
    row.chunks(width)
        .map(|block| {
            block
                .iter()
                .enumerate()
                .fold(0, |value, (i, &d)| value + (i32::from(d as i8) << i))
        })
        .collect()
}

fn direct_product(left: &[i32], right: &[i32]) -> Vec<i64> {
    let mut coefficients = vec![0; left.len() + right.len() - 1];
    for (i, &a) in left.iter().enumerate() {
        if a != 0 {
            for (slot, &b) in coefficients[i..].iter_mut().zip(right) {
                *slot += i64::from(a) * i64::from(b);
            }
        }
    }
    coefficients
}

fn transform_product(left: &[i32], right: &[i32]) -> Vec<i64> {
    let length = left.len() + right.len() - 1;
    let n = length.next_power_of_two();
    let omega = mod_pow(GENERATOR, (MODULUS - 1) / n as u32);
    let mut roots = Vec::with_capacity(n / 2);
    let mut root = 1;
    for _ in 0..n / 2 {
        roots.push(root);
        root = mod_mul(root, omega);
    }
    let forward = |row: &[i32]| {
        let mut values = vec![0; n];
        for (slot, &d) in values.iter_mut().zip(row) {
            *slot = d.rem_euclid(MODULUS as i32) as u32;
        }
        transform(&mut values, &roots);
        values
    };
    let mut values = forward(left);
    for (a, b) in values.iter_mut().zip(forward(right)) {
        *a = mod_mul(*a, b);
    }
    // Applying the same transform twice gives n times the reversed vector.
    transform(&mut values, &roots);
    values[1..].reverse();
    let scale = mod_pow(n as u32, MODULUS - 2);
    values[..length]
        .iter()
        .map(|&v| {
            let x = mod_mul(v, scale);
            if x > MODULUS / 2 {
                i64::from(x) - i64::from(MODULUS)
            } else {
                i64::from(x)
            }
        })
        .collect()
}

fn canonical_row(coefficients: &[i64], width: usize, sign: i64, max_len: usize) -> Vec<u8> {
    let mask = (1_i64 << width) - 1;
    let mut output = vec![0; max_len];
    let mut carry = 0;
    let mut triple_carry = 0;
    let mut block = 0;
    let mut used = 0;
    while block < coefficients.len() || carry != 0 || triple_carry != 0 {
        let value = sign * coefficients.get(block).copied().unwrap_or(0) + carry;
        let digit = value & mask;
        carry = value >> width;
        // For a positive magnitude x, its NAF digit at i is bit(i+1,3x)
        // minus bit(i+1,x). Read both strings a bounded block at a time.
        let triple = 3 * digit + triple_carry;
        triple_carry = triple >> width;
        let mut occupied = (digit ^ triple) & mask;
        while occupied != 0 {
            let bit = occupied.trailing_zeros() as usize;
            used = block * width + bit;
            let signed_digit = if triple & (1_i64 << bit) != 0 {
                sign
            } else {
                -sign
            };
            output[used - 1] = signed_digit as u8;
            occupied &= occupied - 1;
        }
        block += 1;
    }
    output.truncate(used);
    output
}

fn trimmed(row: &[u8]) -> Result<&[u8], &'static str> {
    if row.iter().any(|&d| !matches!(d, 0 | 1 | 255)) {
        return Err("signed digits must be encoded as 0, 1 or 255");
    }
    let len = row.iter().rposition(|&d| d != 0).map_or(0, |i| i + 1);
    Ok(&row[..len])
}

fn combine(left: &[u8], right: &[u8], subtract: bool) -> Result<Vec<u8>, &'static str> {
    let a = trimmed(left)?;
    let b = trimmed(right)?;
    let span = a.len().max(b.len());
    if span > MAX_SPAN {
        return Err("signed addition exceeds the 2^23-coordinate limit");
    }
    let sign = if subtract { -1 } else { 1 };
    let raw =
        |i| a.get(i).copied().unwrap_or(0) as i8 + sign * b.get(i).copied().unwrap_or(0) as i8;
    let mut row = Vec::with_capacity(span + 2);
    let mut carry = 0;
    let mut coefficient = raw(0);
    while row.len() < span || carry != 0 {
        let following = raw(row.len() + 1);
        let value = carry + coefficient;
        // Choose the odd digit so the next carried coefficient is even.
        let digit = if value & 1 != 0 {
            2 - (value + 2 * following).rem_euclid(4)
        } else {
            0
        };
        row.push(digit as u8);
        carry = (value - digit) / 2;
        coefficient = following;
    }
    row.truncate(row.iter().rposition(|&d| d != 0).map_or(0, |i| i + 1));
    Ok(row)
}

#[pyfunction]
fn add_rows<'py>(
    py: Python<'py>,
    left: &[u8],
    right: &[u8],
    subtract: bool,
) -> PyResult<Bound<'py, PyBytes>> {
    let result = py
        .allow_threads(|| combine(left, right, subtract))
        .map_err(PyValueError::new_err)?;
    Ok(PyBytes::new(py, &result))
}

pub(crate) fn multiply(left: &[u8], right: &[u8]) -> Result<Vec<u8>, &'static str> {
    let mut a = trimmed(left)?;
    let mut b = trimmed(right)?;
    if a.len() > b.len() {
        std::mem::swap(&mut a, &mut b);
    }
    if a.is_empty() {
        return Ok(Vec::new());
    }
    let span = a
        .len()
        .checked_add(b.len() - 1)
        .ok_or("signed row is too wide")?;
    if span > MAX_SPAN {
        return Err("signed convolution exceeds the 2^23-coordinate limit");
    }
    let sign = if a.last() == b.last() { 1 } else { -1 };
    let short = a.len().div_ceil(DIRECT_WIDTH);
    let long = b.len().div_ceil(DIRECT_WIDTH);
    let (width, coefficients) = if short <= 32 || short <= 16_384 / long {
        // The direct branch has at most 128 summands per coefficient, each < 2^48.
        (
            DIRECT_WIDTH,
            direct_product(&pack(a, DIRECT_WIDTH), &pack(b, DIRECT_WIDTH)),
        )
    } else {
        // Choose the widest signed blocks for which centered modular residues
        // recover every convolution coefficient exactly, even for adjacent digits.
        let width = (4..=14)
            .rev()
            .find(|&w| {
                let mask = (1_u64 << w) - 1;
                a.len().div_ceil(w) as u64 * mask * mask <= u64::from(MODULUS / 2)
            })
            .expect("four-bit blocks fit the supported span");
        (width, transform_product(&pack(a, width), &pack(b, width)))
    };
    Ok(canonical_row(&coefficients, width, sign, span + 2))
}

#[pyfunction]
fn multiply_rows<'py>(py: Python<'py>, left: &[u8], right: &[u8]) -> PyResult<Bound<'py, PyBytes>> {
    let result = py
        .allow_threads(|| multiply(left, right))
        .map_err(PyValueError::new_err)?;
    Ok(PyBytes::new(py, &result))
}

pub(crate) fn normalize(row: &[u8]) -> Result<Vec<u8>, &'static str> {
    multiply(row, &[1])
}

#[pyfunction]
fn normalize_row<'py>(py: Python<'py>, row: &[u8]) -> PyResult<Bound<'py, PyBytes>> {
    let result = py
        .allow_threads(|| normalize(row))
        .map_err(PyValueError::new_err)?;
    Ok(PyBytes::new(py, &result))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("MAX_MULTIPLY_SPAN", MAX_SPAN)?;
    module.add_function(wrap_pyfunction!(add_rows, module)?)?;
    module.add_function(wrap_pyfunction!(normalize_row, module)?)?;
    module.add_function(wrap_pyfunction!(multiply_rows, module)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_edges_and_cancellation() {
        assert_eq!(multiply(&[255], &[255]).unwrap(), vec![1]);
        assert_eq!(multiply(&[1, 1], &[1, 255, 0]).unwrap(), vec![1, 0, 255]);
        assert!(multiply(&[], &[0]).unwrap().is_empty());
        assert!(multiply(&[0, 2], &[1]).is_err());
    }

    #[test]
    fn signed_addition_handles_cancellation_and_noncanonical_rows() {
        assert_eq!(combine(&[1, 1, 1], &[1], false).unwrap(), [0, 0, 0, 1]);
        assert_eq!(
            combine(&[255, 255, 255], &[1], true).unwrap(),
            [0, 0, 0, 255]
        );
        assert_eq!(combine(&[], &[1, 1], true).unwrap(), [1, 0, 255]);
        assert!(combine(&[1, 255], &[1, 255, 0], true).unwrap().is_empty());
        assert!(combine(&[2], &[], false).is_err());
    }

    #[test]
    fn wide_transform_matches_direct_coefficients() {
        let row = |n: usize, period: usize| -> Vec<u8> {
            (0..n)
                .map(|i| if i % period < 2 { 255 } else { 1 })
                .collect()
        };
        let left = row(5_015, 5);
        let right = row(4_729, 4);
        let coefficients = direct_product(&pack(&left, 16), &pack(&right, 16));
        let sign = if left.last() == right.last() { 1 } else { -1 };
        let expected = canonical_row(&coefficients, 16, sign, left.len() + right.len() + 1);
        assert_eq!(multiply(&left, &right).unwrap(), expected);
    }
}
