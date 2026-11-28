#include <gmp.h>
#include <flint/fmpz.h>

int prime_label_is_prime(mpz_srcptr input)
{
    fmpz_t value;
    fmpz_init(value);
    fmpz_set_mpz(value, input);
    int result = fmpz_is_prime(value);
    fmpz_clear(value);
    return result;
}
