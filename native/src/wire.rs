//! H2N1 fast path for roots with bounded ordinary exponent coordinates.
//! Only exponent values are materialized. The root always remains a signed row.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

use crate::multiply::normalize;
use crate::storage::{naf_weight, weighted_digits};

const MAX_SPAN: usize = 1 << 23;
const UNDEFINED: u32 = u32::MAX;

fn field_size(x: usize) -> usize {
    ((usize::BITS - x.leading_zeros()) as usize)
        .max(1)
        .div_ceil(7)
}

fn push_field(out: &mut Vec<u8>, mut x: usize) {
    while x >= 128 {
        out.push((x as u8) | 128);
        x >>= 7;
    }
    out.push(x as u8);
}

/// Number exponent definitions in the same depth-first postorder as Circuit.encode.
fn intern(position: usize, indices: &mut [u32], count: &mut u32, definitions: &mut Vec<u8>) -> u32 {
    if indices[position] != UNDEFINED {
        return indices[position];
    }
    let mut refs = [0; usize::BITS as usize];
    let mut arity = 0;
    let mut value = position;
    let mut exponent = 0;
    while value != 0 {
        let gap = value.trailing_zeros() as usize;
        exponent += gap;
        value >>= gap;
        let negative = value & 3 == 3;
        let child = intern(exponent, indices, count, definitions);
        refs[arity] = 2 * child as usize + usize::from(negative);
        arity += 1;
        value = if negative { value + 1 } else { value - 1 };
    }
    let id = *count;
    *count += 1;
    indices[position] = id;
    push_field(definitions, arity);
    for &edge in &refs[..arity] {
        push_field(definitions, edge);
    }
    id
}

fn encode(row: &[u8]) -> Vec<u8> {
    let mut indices = vec![UNDEFINED; row.len()];
    let mut count = 0;
    let mut definitions = Vec::new();
    let mut arity = 0;
    for (e, &digit) in row.iter().enumerate() {
        if digit != 0 {
            intern(e, &mut indices, &mut count, &mut definitions);
            arity += 1;
        }
    }
    let mut data =
        Vec::with_capacity(definitions.len() + arity * field_size(2 * count as usize) + 16);
    data.extend_from_slice(b"H2N1");
    push_field(&mut data, count as usize + 1);
    data.extend_from_slice(&definitions);
    push_field(&mut data, arity);
    for (e, &digit) in row.iter().enumerate() {
        if digit != 0 {
            push_field(
                &mut data,
                2 * indices[e] as usize + usize::from(digit == 255),
            );
        }
    }
    data
}

fn canonical(row: &[u8]) -> bool {
    row.len() <= MAX_SPAN
        && row.last() != Some(&0)
        && row.iter().all(|d| matches!(d, 0 | 1 | 255))
        && row.windows(2).all(|pair| pair[0] == 0 || pair[1] == 0)
}

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl Reader<'_> {
    fn field(&mut self) -> Result<usize, &'static str> {
        let mut value = 0usize;
        let bound = self.data.len().saturating_mul(2);
        for shift in (0..usize::BITS).step_by(7) {
            let byte = *self
                .data
                .get(self.pos)
                .ok_or("truncated signed circuit integer")?;
            self.pos += 1;
            let payload = usize::from(byte & 127);
            if payload > bound >> shift {
                return Err("signed circuit integer exceeds format bounds");
            }
            value |= payload << shift;
            if byte < 128 {
                if shift != 0 && byte == 0 {
                    return Err("overlong signed circuit integer");
                }
                return Ok(value);
            }
        }
        Err("signed circuit integer exceeds format bounds")
    }
}

// None asks the general decoder to handle a tower or a sparse range.
fn decode(data: &[u8]) -> Result<Option<Vec<u8>>, &'static str> {
    if !data.starts_with(b"H2N1") {
        return Err("invalid signed circuit encoding");
    }
    let mut r = Reader { data, pos: 4 };
    let count = r.field()?;
    if count == 0 || count > data.len() - r.pos {
        return Err("invalid signed circuit node count");
    }
    if count - 1 > MAX_SPAN {
        return Ok(None); // The ordinary path needs distinct coordinates below MAX_SPAN.
    }
    let mut values = Vec::<u32>::with_capacity(count - 1);
    let mut defined = vec![false; count.min(MAX_SPAN)];
    let mut referenced = vec![false; count];
    let mut root = Vec::new();
    for row in 0..count {
        let arity = r.field()?;
        if arity > (data.len() - r.pos).saturating_sub(count - row - 1) {
            return Err("invalid signed circuit node arity");
        }
        let last = row + 1 == count;
        let mut previous = None;
        let mut value = 0i64;
        for _ in 0..arity {
            let edge = r.field()?;
            let child = edge / 2;
            let e = *values
                .get(child)
                .ok_or("signed circuit reference must point to an earlier node")?
                as usize;
            if previous.is_some_and(|p| p >= e) {
                return Err("signed circuit exponents must be strictly increasing");
            }
            previous = Some(e);
            referenced[child] = true;
            let digit = if edge & 1 == 0 { 1i64 } else { -1 };
            if last {
                if e >= MAX_SPAN || e >= arity.saturating_mul(16).max(256) {
                    return Ok(None);
                }
                root.resize(e + 1, 0);
                root[e] = digit as u8;
            } else {
                if e > 30 {
                    return Ok(None);
                }
                value += digit << e;
            }
        }
        if !last {
            if value < 0 {
                return Err("negative exponent in signed circuit");
            }
            if value as usize >= MAX_SPAN {
                return Ok(None);
            }
            let index = value as usize;
            if index >= defined.len() {
                defined.resize(index + 1, false);
            }
            if defined[index] {
                return Err("duplicate value in signed circuit");
            }
            defined[index] = true;
            values.push(value as u32);
        }
    }
    if r.pos != data.len() {
        return Err("trailing data after signed circuit");
    }
    if referenced[..count - 1].contains(&false) {
        return Err("unreachable node in signed circuit");
    }
    let root = normalize(&root)?;
    // All preceding rows are nonnegative coordinates below MAX_SPAN.
    if root.len() <= 25 {
        let value: i64 = root
            .iter()
            .enumerate()
            .map(|(e, &d)| i64::from(d as i8) << e)
            .sum();
        if value >= 0 && defined.get(value as usize) == Some(&true) {
            return Err("duplicate value in signed circuit");
        }
    }
    Ok(Some(root))
}

/// Weighted root candidate and whole-graph lower bound for wide signed roots.
fn pack(row: &[u8]) -> Result<(Vec<u8>, usize), &'static str> {
    let minimum_terms = row.iter().filter(|&&d| d != 0).count();
    if !canonical(row) || minimum_terms < 8192 {
        return Err("expected a wide canonical signed row");
    }
    let sign = *row.last().unwrap() as i8;
    let mut bits = Vec::with_capacity(row.len() + 1);
    let mut carry = 0i8;
    for &digit in row {
        carry += sign * digit as i8;
        bits.push((carry & 1) as u8);
        carry >>= 1;
    }
    bits.push(0);
    let costs: Vec<_> = (0..bits.len()).map(|e| 4 + naf_weight(e)).collect();
    let digits = weighted_digits(&bits, &costs).ok_or("unrepresentable root")?;
    let mut alternative = vec![0u8; bits.len()];
    let mut cost = 0usize;
    for (e, d) in digits {
        alternative[e] = d as u8;
        cost += usize::from(costs[e]);
    }
    let k = (usize::BITS - bits.len().leading_zeros()) as usize;
    // For an odd local prefix P of width w, 2B - (2^w - 1) = P.
    // Writing +/-1 at each position according to B references low nodes early.
    if let Some(start) = alternative[..=k].iter().position(|&d| d != 0) {
        let low = &mut alternative[start..=k];
        let prefix: i64 = low
            .iter()
            .enumerate()
            .map(|(i, &d)| i64::from(d as i8) << i)
            .sum();
        let positive = (prefix + (1i64 << low.len()) - 1) / 2;
        for (i, digit) in low.iter_mut().enumerate() {
            *digit = if positive >> i & 1 != 0 { 1 } else { 255 };
        }
    }
    if sign < 0 {
        for d in &mut alternative {
            *d = d.wrapping_neg();
        }
    }
    let bound = |weight: usize| {
        4 + field_size(minimum_terms) + field_size(minimum_terms + 1) + weight - 64 - 8192
            + minimum_terms.saturating_sub(1 << 20)
    };
    let lower = bound(cost);
    let mut best = encode(&alternative);
    // Fixing the root to its canonical digits gives a stronger bound for
    // that one candidate. Serialize it only if it could win (including ties).
    let canonical_weight = row
        .iter()
        .zip(&costs)
        .filter(|(digit, _)| **digit != 0)
        .map(|(_, weight)| usize::from(*weight))
        .sum();
    if best.len() >= bound(canonical_weight) {
        let canonical = encode(row);
        if canonical.len() <= best.len() {
            best = canonical;
        }
    }
    Ok((best, lower))
}

#[pyfunction]
fn encode_dense<'py>(py: Python<'py>, row: &[u8]) -> PyResult<Bound<'py, PyBytes>> {
    if !canonical(row) {
        return Err(PyValueError::new_err("expected a canonical signed row"));
    }
    let bytes = py.allow_threads(|| encode(row));
    Ok(PyBytes::new(py, &bytes))
}

#[pyfunction]
fn decode_dense<'py>(py: Python<'py>, data: &[u8]) -> PyResult<Option<Bound<'py, PyBytes>>> {
    let row = py
        .allow_threads(|| decode(data))
        .map_err(PyValueError::new_err)?;
    Ok(row.map(|x| PyBytes::new(py, &x)))
}

#[pyfunction]
fn pack_dense<'py>(py: Python<'py>, row: &[u8]) -> PyResult<(Bound<'py, PyBytes>, usize)> {
    let (bytes, bound) = py
        .allow_threads(|| pack(row))
        .map_err(PyValueError::new_err)?;
    Ok((PyBytes::new(py, &bytes), bound))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(encode_dense, m)?)?;
    m.add_function(wrap_pyfunction!(decode_dense, m)?)?;
    m.add_function(wrap_pyfunction!(pack_dense, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonicalizes_alternative_root_and_validates_all_references() {
        // 3 = 2^0 + 2^1, using the preceding definitions of 0 and 1.
        assert_eq!(
            decode(b"H2N1\x03\x00\x01\x00\x02\x00\x02").unwrap(),
            Some(vec![255, 0, 1])
        );
        for data in [
            &b"H2N1\x81\x00\x00"[..],        // overlong integer
            b"H2N1\x01\x01\x00",             // forward edge
            b"H2N1\x02\x00\x02\x00\x00",     // repeated exponent
            b"H2N1\x03\x00\x01\x00\x01\x00", // duplicate value
            b"H2N1\x03\x00\x01\x00\x01\x01", // unreachable definition of 1
            b"H2N1\x01\x00\x00",             // trailing bytes
        ] {
            assert!(decode(data).is_err());
        }
    }

    #[test]
    fn preserves_both_signs_in_dense_wire_round_trips() {
        for n in [0, 1, 64, 8192] {
            let mut row = vec![0; n];
            for (i, d) in row.iter_mut().rev().step_by(3).enumerate() {
                *d = if i % 2 == 0 { 1 } else { 255 };
            }
            for _ in 0..2 {
                assert_eq!(
                    decode(&encode(&row)).unwrap().as_deref(),
                    Some(row.as_slice())
                );
                for d in &mut row {
                    *d = d.wrapping_neg();
                }
            }
        }
    }
}
