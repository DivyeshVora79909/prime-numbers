use rug::{Integer, integer::Order, ops::Pow};
use std::{ffi::c_void, fmt, io, str::FromStr};

unsafe extern "C" {
    fn prime_label_is_prime(input: *const c_void) -> std::ffi::c_int;
}

/// An entity ID whose primality has been proved by FLINT.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct PrimeId(Integer);

impl PrimeId {
    fn checked(value: Integer) -> Option<Self> {
        // SAFETY: rug owns this initialized GMP integer for the call. The C
        // shim copies it to a local FLINT integer and retains no pointers.
        (unsafe { prime_label_is_prime(value.as_raw().cast()) } == 1).then_some(Self(value))
    }

    /// Generate a uniformly sampled 40-digit prime using the operating system RNG.
    pub fn generate() -> io::Result<Self> {
        let low = Integer::from(10).pow(39u32);
        let high = Integer::from(&low * 10u32);
        loop {
            let mut bytes = [0u8; 17];
            getrandom::fill(&mut bytes).map_err(|error| io::Error::other(error.to_string()))?;
            // Each odd 133-bit candidate has the same number of preimages.
            bytes[16] &= 0x1f;
            bytes[0] |= 1;
            let candidate = Integer::from_digits(&bytes, Order::Lsf);
            if candidate >= low
                && candidate < high
                && let Some(prime) = Self::checked(candidate)
            {
                return Ok(prime);
            }
        }
    }

    pub fn value(&self) -> &Integer {
        &self.0
    }

    pub(crate) fn into_integer(self) -> Integer {
        self.0
    }
}

impl FromStr for PrimeId {
    type Err = io::Error;

    fn from_str(text: &str) -> io::Result<Self> {
        if text.is_empty() || !text.bytes().all(|digit| digit.is_ascii_digit()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "expected a decimal prime ID",
            ));
        }
        let value = Integer::from_str_radix(text, 10)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
        Self::checked(value)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "entity ID must be prime"))
    }
}

impl fmt::Display for PrimeId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn proves_decimals_and_rejects_strong_pseudoprime() {
        let member = "3589358312608496142804550037657486194871";
        assert_eq!(member.parse::<PrimeId>().unwrap().to_string(), member);
        assert_eq!("002".parse::<PrimeId>().unwrap().to_string(), "2");
        for invalid in [
            "",
            "0",
            "1",
            "4",
            "-3",
            "+7",
            "3.0",
            " 7",
            "7\0",
            "318665857834031151167461",
        ] {
            assert!(invalid.parse::<PrimeId>().is_err(), "{invalid:?}");
        }
    }

    #[test]
    fn issued_id_is_a_proved_40_digit_prime() {
        let prime = PrimeId::generate().unwrap();
        let text = prime.to_string();
        assert_eq!(text.len(), 40);
        assert_eq!(text.parse::<PrimeId>().unwrap(), prime);
    }
}
