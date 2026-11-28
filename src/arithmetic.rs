//! Block arithmetic on canonical positive, little-endian magnitudes.

use std::{fs::File, io, os::unix::fs::FileExt};

use bytemuck::{cast_slice, cast_slice_mut};
use rug::{Assign, Integer, integer::Order};

const BLOCK: usize = 64 * 1024;

pub(crate) fn label_len(input: &File) -> io::Result<u64> {
    let metadata = input.metadata()?;
    if !metadata.is_file() || metadata.len() == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "expected a nonempty label file",
        ));
    }
    let mut high = [0];
    input.read_exact_at(&mut high, metadata.len() - 1)?;
    if high[0] == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "label has a zero most significant byte",
        ));
    }
    Ok(metadata.len())
}

// Use aligned limbs in GMP even when an input or output ends partway through a limb.
fn read_block(
    file: &File,
    offset: u64,
    n: usize,
    words: &mut [u64],
    value: &mut Integer,
) -> io::Result<()> {
    let words = &mut words[..n.div_ceil(8)];
    words[words.len() - 1] = 0;
    file.read_exact_at(&mut cast_slice_mut(words)[..n], offset)?;
    value.assign_digits(words, Order::LsfLe);
    Ok(())
}

fn write_block(
    file: &File,
    offset: u64,
    n: usize,
    words: &mut [u64],
    value: &Integer,
) -> io::Result<()> {
    let words = &mut words[..n.div_ceil(8)];
    value.write_digits(words, Order::LsfLe);
    file.write_all_at(&cast_slice(words)[..n], offset)
}

/// Long division walks from the most significant block down. Each new
/// dividend is previous_remainder * block_radix + current_block.
fn scan(input: &File, divisor: &Integer, output: Option<&File>) -> io::Result<Integer> {
    if divisor <= &0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "expected a positive divisor",
        ));
    }
    let mut end = label_len(input)?;
    if output.is_none() && divisor == &1 {
        return Ok(Integer::new());
    }
    let mut words = vec![0; (end.min(BLOCK as u64) as usize).div_ceil(8)];
    let (mut block, mut rem, mut quotient) = (Integer::new(), Integer::new(), Integer::new());
    let mut output_len = 0;
    while end > 0 {
        let n = (1 + (end - 1) % BLOCK as u64) as usize;
        end -= n as u64;
        read_block(input, end, n, &mut words, &mut block)?;
        rem <<= 8 * n;
        block += &rem;
        if let Some(out) = output {
            (&mut quotient, &mut rem).assign(block.div_rem_ref(divisor));
            let count = if output_len == 0 {
                let significant = quotient.significant_digits::<u8>();
                if significant == 0 {
                    continue;
                }
                output_len = end + significant as u64; // bounded by this input block's end
                significant
            } else {
                n
            };
            write_block(out, end, count, &mut words, &quotient)?;
        } else {
            rem.assign(&block % divisor);
        }
    }
    if let Some(out) = output {
        if rem != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "division is not exact",
            ));
        }
        out.set_len(output_len)?;
    }
    Ok(rem)
}

pub(crate) fn remainder(input: &File, modulus: &Integer) -> io::Result<Integer> {
    scan(input, modulus, None)
}

/// The caller owns publication of the scratch output, including on inexact division.
pub(crate) fn divide(input: &File, divisor: &Integer, output: &File) -> io::Result<()> {
    scan(input, divisor, Some(output)).map(|_| ())
}

/// Low-to-high multiplication permits the unpublished input and output to be
/// the same file. Capture the source length once; the final carry can extend it.
pub(crate) fn multiply(input: &File, factor: &Integer, output: &File) -> io::Result<()> {
    if factor <= &0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "expected a positive factor",
        ));
    }
    let len = label_len(input)?;
    let mut words = vec![0; (len.min(BLOCK as u64) as usize).div_ceil(8)];
    let (mut block, mut carry) = (Integer::new(), Integer::new());
    let mut offset = 0;
    while offset < len {
        let n = (len - offset).min(BLOCK as u64) as usize;
        read_block(input, offset, n, &mut words, &mut block)?;
        block *= factor;
        block += &carry;
        carry.assign(&block >> (8 * n));
        block.keep_bits_mut((8 * n) as u32);
        write_block(output, offset, n, &mut words, &block)?;
        offset += n as u64;
    }
    let tail_len = carry.significant_digits::<u8>();
    let final_len = len.checked_add(tail_len as u64).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "product exceeds file offset range",
        )
    })?;
    if tail_len != 0 {
        let tail = carry.to_digits::<u64>(Order::LsfLe);
        output.write_all_at(&cast_slice(&tail)[..tail_len], len)?;
    }
    output.set_len(final_len)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempfile;

    fn file_with(bytes: &[u8]) -> File {
        let file = tempfile().unwrap();
        file.write_all_at(bytes, 0).unwrap();
        file
    }

    fn contents(file: &File) -> Vec<u8> {
        let mut bytes = vec![0; file.metadata().unwrap().len() as usize];
        file.read_exact_at(&mut bytes, 0).unwrap();
        bytes
    }

    #[test]
    fn boundaries_and_in_place_product_match_full_integer_arithmetic() {
        let factor = (Integer::from(1) << 132u32) + 17;
        for n in [1, 7, 8, 9, BLOCK - 1, BLOCK, BLOCK + 1, 2 * BLOCK + 13] {
            let mut bytes = vec![255; n];
            bytes[n / 2] = 3;
            let src = file_with(&bytes);
            let value = Integer::from_digits(&bytes, Order::LsfLe);
            assert_eq!(
                remainder(&src, &factor).unwrap(),
                Integer::from(&value % &factor)
            );
            multiply(&src, &factor, &src).unwrap();
            let product = Integer::from(&value * &factor);
            assert_eq!(contents(&src), product.to_digits::<u8>(Order::LsfLe));
            let dest = file_with(&vec![0xaa; n + 90]);
            divide(&src, &factor, &dest).unwrap();
            assert_eq!(contents(&dest), bytes);
        }
    }

    #[test]
    fn batch_factor_can_exceed_a_block() {
        let value = file_with(&[3]);
        let large = (Integer::from(1) << (8 * BLOCK + 37)) + 5;
        let scratch = tempfile().unwrap();
        multiply(&value, &large, &scratch).unwrap();
        assert_eq!(remainder(&scratch, &large).unwrap(), 0);
        divide(&scratch, &large, &scratch).unwrap();
        assert_eq!(contents(&scratch), [3]);
    }

    #[test]
    fn malformed_inputs_and_inexact_division_return_errors() {
        let out = tempfile().unwrap();
        let two = Integer::from(2);
        for bytes in [&[][..], &[0], &[1, 0]] {
            let src = file_with(bytes);
            assert!(remainder(&src, &two).is_err());
            assert!(multiply(&src, &two, &out).is_err());
            assert!(divide(&src, &two, &out).is_err());
        }
        let three = file_with(&[3]);
        assert!(divide(&three, &two, &out).is_err());
        assert!(remainder(&three, &Integer::new()).is_err());
    }
}
