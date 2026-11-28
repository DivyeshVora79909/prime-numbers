//! Carry dynamic programming for dense, ordinary exponent coordinates.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

const INF: u64 = u64::MAX;

pub(crate) fn naf_weight(position: usize) -> u8 {
    // (3x XOR x) has one set bit per nonadjacent signed digit of x.
    let x = position as u128;
    ((3 * x) ^ x).count_ones() as u8
}

/// Minimum suffix costs when the incoming binary carry is zero or one.
/// The optional trace records the same positive-on-tie decision as Python.
fn root_cost(bits: &[u8], cost: impl Fn(usize) -> u64, mut negative: Option<&mut [u8]>) -> u64 {
    let (mut clear, mut carry) = (0, INF);
    for e in (0..bits.len()).rev() {
        if let Some(trace) = negative.as_deref_mut() {
            trace[e] = (carry < clear) as u8;
        }
        let odd = clear.min(carry).saturating_add(cost(e));
        if bits[e] == 0 {
            carry = odd;
        } else {
            clear = odd;
        }
    }
    clear
}

pub(crate) fn weighted_digits(bits: &[u8], costs: &[u8]) -> Option<Vec<(usize, i8)>> {
    if bits.len() != costs.len() {
        return None;
    }
    let mut negative = vec![0; bits.len()];
    if root_cost(bits, |e| u64::from(costs[e]), Some(&mut negative)) == INF {
        return None;
    }
    let mut carry = 0;
    let mut digits = Vec::new();
    for (e, (&bit, &negative)) in bits.iter().zip(&negative).enumerate() {
        let coefficient = bit + carry;
        if coefficient == 1 {
            digits.push((e, if negative == 1 { -1 } else { 1 }));
            carry = negative;
        } else {
            carry = coefficient / 2;
        }
    }
    Some(digits)
}

#[pyfunction]
fn naf_weights(py: Python<'_>, limit: usize) -> PyResult<Bound<'_, PyBytes>> {
    let count = limit
        .checked_add(1)
        .filter(|&n| n <= isize::MAX as usize)
        .ok_or_else(|| PyValueError::new_err("weight table too large"))?;
    PyBytes::new_with(py, count, |data| {
        for (position, weight) in data.iter_mut().enumerate() {
            *weight = naf_weight(position);
        }
        Ok(())
    })
}

#[pyfunction]
fn weighted_root_digits<'py>(
    py: Python<'py>,
    bits: &[u8],
    costs: &[u8],
) -> PyResult<Bound<'py, PyBytes>> {
    if bits.len() != costs.len() || bits.iter().any(|&bit| bit > 1) {
        return Err(PyValueError::new_err(
            "expected one weight per binary column",
        ));
    }
    let digits = py
        .allow_threads(|| weighted_digits(bits, costs))
        .ok_or_else(|| PyValueError::new_err("weighted row exceeds native cost range"))?;
    // Build Python's (position, sign) tuples with struct.iter_unpack in one
    // pass, avoiding a separate Python C API conversion for every object.
    PyBytes::new_with(py, digits.len() * 9, |buffer| {
        for (field, (position, sign)) in buffer.chunks_exact_mut(9).zip(digits) {
            field[..8].copy_from_slice(&(position as u64).to_le_bytes());
            field[8] = sign as u8;
        }
        Ok(())
    })
}

// Higher positions use the ordinary minimum weight. At restricted positions
// a missing bit prohibits both signs of the corresponding digit.
fn restricted_weights(ordinary: &[u8], positions: usize, possible: usize) -> Vec<u64> {
    let limit = ordinary.len() - 1;
    let mut above: Vec<_> = ordinary.iter().map(|&w| u64::from(w)).collect();
    for position in (0..positions).rev() {
        let scale = 1usize << position;
        let available = possible & scale != 0;
        let largest = limit.div_ceil(scale);
        above = (0..=largest)
            .map(|x| {
                let half = x / 2;
                if x % 2 == 0 {
                    above[half]
                } else if available {
                    above[half].min(above[half + 1]).saturating_add(1)
                } else {
                    INF
                }
            })
            .collect();
    }
    above
}

/// Score one branch of the existing proof search. Undecided low nodes serve
/// descendant rows for free in both scores. Only the completion must actually
/// pay for all their definitions. The Python caller controls the state budget.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn storage_subset_costs(
    py: Python<'_>,
    bits: &[u8],
    ordinary: &[u8],
    low_count: usize,
    restricted_count: usize,
    present: usize,
    possible: usize,
    row_header: u64,
    ceiling: u64,
) -> PyResult<(u64, u64)> {
    if low_count >= usize::BITS as usize
        || restricted_count > low_count
        || low_count > ordinary.len()
        || bits.len() != ordinary.len()
        || bits.len() < 2
        || bits.iter().any(|&b| b > 1)
        || present & !possible != 0
        || possible >> low_count != 0
    {
        return Err(PyValueError::new_err("invalid ordinary storage subset"));
    }
    Ok(py.allow_threads(|| {
        let weights = restricted_weights(ordinary, restricted_count, possible);
        let mut costs: Vec<_> = weights
            .iter()
            .map(|w| w.saturating_add(row_header).saturating_add(1))
            .collect();
        let mut base = 0u64;
        for e in 0..low_count {
            let bit = 1 << e;
            if present & bit != 0 {
                base = base.saturating_add(row_header).saturating_add(weights[e]);
                costs[e] = 1;
            } else if possible & bit == 0 {
                costs[e] = INF;
            }
        }
        if base >= ceiling {
            return (ceiling, ceiling);
        }
        let optimistic = base.saturating_add(root_cost(bits, |e| costs[e], None));
        if optimistic >= ceiling || present == possible {
            return (optimistic.min(ceiling), optimistic.min(ceiling));
        }
        for e in 0..low_count {
            if (possible & !present) & (1 << e) != 0 {
                base = base.saturating_add(row_header).saturating_add(weights[e]);
                costs[e] = 1;
            }
        }
        let completion = base.saturating_add(root_cost(bits, |e| costs[e], None));
        (optimistic, completion.min(ceiling))
    }))
}

/// The subgraph needed to score a change to the first 64 wire IDs. The
/// adjacency index is reused for every candidate within one layout pass.
#[pyclass(frozen)]
struct PrefixLayout {
    rows: Vec<Vec<usize>>,
    counts: Vec<usize>,
}

#[pymethods]
impl PrefixLayout {
    #[new]
    fn new(rows: Vec<Vec<usize>>, counts: Vec<usize>) -> PyResult<Self> {
        if rows.is_empty()
            || rows.len() != counts.len()
            || rows.iter().flatten().any(|&i| i >= rows.len())
        {
            return Err(PyValueError::new_err("invalid prefix layout"));
        }
        Ok(Self { rows, counts })
    }

    fn credit(&self, changed: usize, replacement: Vec<usize>) -> PyResult<usize> {
        if changed >= self.rows.len() || replacement.iter().any(|&i| i >= self.rows.len()) {
            return Err(PyValueError::new_err("invalid replacement row"));
        }
        // DFS postorder gives the IDs. A change can bring fresh nodes into
        // the prefix, so walk the proposed row just as the encoder would.
        let mut state = vec![0u8; self.rows.len()];
        let mut stack = vec![(0, 0usize)];
        state[0] = 1;
        let mut remaining = 64;
        let mut credit: usize = 0;
        while let Some((node, index)) = stack.last_mut() {
            let node = *node;
            let row = if node == changed {
                &replacement
            } else {
                &self.rows[node]
            };
            if let Some(&child) = row.get(*index) {
                *index += 1;
                match state[child] {
                    0 => {
                        state[child] = 1;
                        stack.push((child, 0));
                    }
                    1 => return Err(PyValueError::new_err("cyclic prefix layout")),
                    _ => {}
                }
            } else {
                state[node] = 2;
                credit = credit
                    .checked_add(self.counts[node])
                    .ok_or_else(|| PyValueError::new_err("prefix reference count overflow"))?;
                remaining -= 1;
                stack.pop();
                if remaining == 0 {
                    break;
                }
            }
        }
        let in_prefix = |i: &&usize| state[**i] == 2;
        let added = replacement.iter().filter(in_prefix).count();
        let removed = self.rows[changed].iter().filter(in_prefix).count();
        credit
            .checked_add(added)
            .and_then(|n| n.checked_sub(removed))
            .ok_or_else(|| PyValueError::new_err("inconsistent prefix reference counts"))
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(naf_weights, m)?)?;
    m.add_function(wrap_pyfunction!(weighted_root_digits, m)?)?;
    m.add_function(wrap_pyfunction!(storage_subset_costs, m)?)?;
    m.add_class::<PrefixLayout>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimum_weight_matches_recurrence() {
        let mut expected = vec![0, 1];
        for x in 2..8194 {
            let half = x / 2;
            let w = if x % 2 == 0 {
                expected[half]
            } else {
                1 + expected[half].min(expected[half + 1])
            };
            expected.push(w);
        }
        assert!(expected
            .iter()
            .enumerate()
            .all(|(e, w)| naf_weight(e) == *w));
    }

    #[test]
    fn weighted_traceback_prefers_positive_on_ties() {
        assert_eq!(
            weighted_digits(&[1, 1, 0], &[1, 1, 1]),
            Some(vec![(0, 1), (1, 1)])
        );
        assert_eq!(
            weighted_digits(&[1, 1, 0], &[1, 8, 1]),
            Some(vec![(0, -1), (2, 1)])
        );
    }

    #[test]
    fn forbidden_low_digit_makes_odd_exponents_unavailable() {
        let ordinary: Vec<_> = (0..=100).map(naf_weight).collect();
        let restricted = restricted_weights(&ordinary, 1, 0);
        for (value, &cost) in restricted.iter().enumerate() {
            let expected = if value % 2 == 0 {
                u64::from(ordinary[value / 2])
            } else {
                INF
            };
            assert_eq!(cost, expected);
        }
    }
}
