//! Bounded-remainder arithmetic on signed base-two coefficients, low position first.
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{ToPrimitive, Zero};
use pyo3::exceptions::{PyNotImplementedError, PyOverflowError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

const WORD: usize = 64;

fn orientation(row: &[u8]) -> PyResult<i8> {
    let mut sign = 1;
    for &d in row {
        match d {
            0 => {}
            1 | 255 => sign = d as i8,
            _ => return Err(PyValueError::new_err("invalid signed binary coefficient")),
        }
    }
    Ok(sign)
}

/// Read unsigned magnitude one word at a time, without materializing the magnitude.
/// A lower signed prefix borrows exactly when its highest occupied position is negative.
fn words(row: &[u8], sign: i8) -> impl Iterator<Item = u64> + '_ {
    let mut lower = row.len().checked_sub(1);
    row.chunks(WORD).enumerate().rev().map(move |(i, chunk)| {
        if lower.is_some_and(|j| j >= i * WORD) {
            lower = row[..i * WORD].iter().rposition(|&d| d != 0);
        }
        let borrow = lower.is_some_and(|j| row[j] as i8 != sign);
        let value = chunk.iter().rev().fold(0u64, |v, &d| {
            (v << 1).wrapping_add_signed(i64::from(d as i8 * sign))
        });
        value.wrapping_sub(u64::from(borrow))
    })
}

/// Supply at most one quotient word per input word. With a remainder-only query,
/// the compiler removes quotient extraction and the callback.
fn divide_words<const QUOTIENT: bool>(
    row: &[u8],
    sign: i8,
    modulus: &BigUint,
    mut emit: impl FnMut(usize, u64),
) -> BigUint {
    let offsets = (0..row.len()).step_by(WORD).rev();
    let input = offsets.zip(words(row, sign));
    if let Some(m) = modulus.to_u64() {
        let m = u128::from(m);
        let mut remainder = 0;
        for (offset, word) in input {
            let current = (remainder << WORD) + u128::from(word);
            if QUOTIENT {
                emit(offset, (current / m) as u64);
            }
            remainder = current % m;
        }
        return BigUint::from(remainder);
    }

    let mut remainder = BigUint::zero();
    for (offset, word) in input {
        remainder <<= WORD;
        remainder += word;
        if QUOTIENT {
            // Previous remainder < modulus, hence the quotient fits in one word.
            let (q, r) = remainder.div_rem(modulus);
            emit(offset, q.to_u64().expect("single quotient word"));
            remainder = r;
        } else {
            remainder %= modulus;
        }
    }
    remainder
}

fn nonzero(data: &[u8]) -> PyResult<BigUint> {
    let value = BigUint::from_bytes_le(data);
    if value.is_zero() {
        Err(PyValueError::new_err(
            "expected a positive divisor or modulus",
        ))
    } else {
        Ok(value)
    }
}

#[pyfunction]
fn dense_residue<'py>(
    py: Python<'py>,
    row: &[u8],
    modulus_le: &[u8],
) -> PyResult<Bound<'py, PyBytes>> {
    let modulus = nonzero(modulus_le)?;
    let bytes = py.allow_threads(|| -> PyResult<Vec<u8>> {
        let sign = orientation(row)?;
        let mut r = divide_words::<false>(row, sign, &modulus, |_, _| {});
        if sign < 0 && !r.is_zero() {
            r = modulus - r;
        }
        Ok(r.to_bytes_le())
    })?;
    Ok(PyBytes::new(py, &bytes))
}

fn naf_digit(bit: u8, next: u8, carry: &mut u8) -> i8 {
    let sum = bit + *carry;
    if sum == 1 {
        *carry = next;
        return 1 - 2 * next as i8;
    }
    *carry = sum / 2;
    0
}

#[pyfunction]
#[pyo3(signature = (row, divisor_le, max_terms=None))]
fn dense_divide<'py>(
    py: Python<'py>,
    row: &[u8],
    divisor_le: &[u8],
    max_terms: Option<usize>,
) -> PyResult<Option<Bound<'py, PyBytes>>> {
    let divisor = nonzero(divisor_le)?;
    let result = py.allow_threads(|| -> PyResult<Option<Vec<u8>>> {
        let sign = orientation(row)?;
        let mut digits = vec![0; row.len() + 2];
        let r = divide_words::<true>(row, sign, &divisor, |offset, mut word| {
            while word != 0 {
                digits[offset + word.trailing_zeros() as usize] = 1;
                word &= word - 1;
            }
        });
        let exact = r.is_zero();
        if !exact && max_terms.is_none() {
            return Ok(None);
        }
        let mut carry = 0;
        let mut ceiling_carry = 1;
        let mut weight = 0;
        let mut ceiling_weight = 0;
        for i in 0..digits.len() - 1 {
            let bit = digits[i];
            let next = digits[i + 1];
            let d = naf_digit(bit, next, &mut carry);
            if !exact {
                let upper = naf_digit(bit, next, &mut ceiling_carry);
                ceiling_weight += usize::from(upper != 0);
            }
            weight += usize::from(d != 0);
            digits[i] = (d * sign) as u8;
        }
        if let Some(limit) = max_terms {
            // For an inexact division the symbolic greedy prefix ends within one
            // unit of |n/d|. Its nonadjacent digits therefore encode floor or ceil.
            if limit < weight && (exact || limit < ceiling_weight) {
                return Err(PyOverflowError::new_err(
                    "quotient exceeds its signed-digit budget",
                ));
            }
            if !exact && limit < weight.max(ceiling_weight) {
                // The budget distinguishes the neighboring candidates. Delegate
                // their precise tie decision to the original symbolic algorithm.
                return Err(PyNotImplementedError::new_err(
                    "boundary of inexact quotient budget",
                ));
            }
        }
        if !exact {
            return Ok(None);
        }
        while digits.last() == Some(&0) {
            digits.pop();
        }
        Ok(Some(digits))
    })?;
    Ok(result.map(|digits| PyBytes::new(py, &digits)))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(dense_residue, module)?)?;
    module.add_function(wrap_pyfunction!(dense_divide, module)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounded_words_match_signed_definition_across_boundaries() {
        for span in [1, 63, 64, 65, 128, 129, 256, 257] {
            let mut row = vec![0; span];
            row[0] = 255;
            row[span - 1] = 1;
            let m = (BigUint::from(1u8) << 130) + 13u8;
            let n = if span == 1 {
                BigUint::from(1u8)
            } else {
                (BigUint::from(1u8) << (span - 1)) - 1u8
            };
            let mut q = BigUint::zero();
            let r = divide_words::<true>(&row, 1, &m, |offset, word| {
                q += BigUint::from(word) << offset;
            });
            assert_eq!(&q * &m + &r, n);
            assert_eq!(r, divide_words::<false>(&row, 1, &m, |_, _| {}));
        }
    }
}
