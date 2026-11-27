import unittest
from random import Random

from symbolic_primes import Circuit, decode_storage, multiply_fast, multiply_storage, optimize_storage


class StorageTests(unittest.TestCase):
    def test_stored_factors_and_direct_integer_have_the_same_proved_optimum(self):
        first = 943433720027638044163815075978910058741
        second = 705383844795881156131636713278929157563
        source = Circuit()
        left = optimize_storage(source, source.integer(first))
        right = optimize_storage(source, source.integer(second))
        self.assertTrue(left.optimal and right.optimal)
        self.assertEqual((len(left.data), len(right.data)), (233, 213))

        calculated = multiply_storage(Circuit(), left.data, right.data)
        direct = Circuit()
        integer = direct.integer(first * second)
        self.assertEqual(calculated, optimize_storage(direct, integer))
        self.assertIs(decode_storage(direct, calculated.data), integer)
        self.assertTrue(calculated.optimal)
        self.assertEqual(len(calculated.data), 496)

    def test_multiplying_with_a_budget_returns_proof_status(self):
        circuit = Circuit()
        stored = circuit.encode(circuit.integer(-11))
        result = multiply_storage(circuit, stored, circuit.integer(48), max_states=0)
        self.assertFalse(result.optimal)
        self.assertIs(decode_storage(circuit, result.data), circuit.integer(-528))

    def test_joint_expansion_can_save_an_exponent_node(self):
        circuit = Circuit()
        value = circuit.integer(40)

        result = optimize_storage(circuit, value)
        self.assertTrue(result.optimal)
        self.assertEqual(len(circuit.encode(value)), 19)
        self.assertEqual(len(result.data), 18)

        restored = Circuit()
        value = decode_storage(restored, result.data)
        self.assertIs(value, restored.integer(40))
        self.assertIs(
            multiply_fast(restored, value, restored.integer(3)), restored.integer(120)
        )

    def test_state_limit_retains_an_exact_encoding(self):
        circuit = Circuit()
        value = circuit.integer(-40)

        result = optimize_storage(circuit, value, max_states=0)
        self.assertFalse(result.optimal)
        self.assertEqual(result.states, 0)
        self.assertEqual(result.data, circuit.encode(value))
        self.assertIs(decode_storage(circuit, result.data), value)

    def test_quality_bound_and_determinism_on_wide_input(self):
        number = (1 << 1023) | Random(130259).getrandbits(1023) | 1
        circuit = Circuit()
        value = circuit.integer(number)
        result = optimize_storage(circuit, value)
        self.assertIs(decode_storage(circuit, result.data), value)
        self.assertLessEqual(len(result.data), len(circuit.encode(value)))
        self.assertLessEqual(100 * result.max_extra_bytes, result.lower_bound_bytes)
        self.assertLessEqual(result.states, 2048)
        fresh = Circuit()
        self.assertEqual(result, optimize_storage(fresh, fresh.integer(number)))

    def test_gap_target_and_budget_against_known_optimum(self):
        circuit = Circuit()
        value = circuit.integer(3)
        exact = optimize_storage(circuit, value, max_states=None)
        self.assertTrue(exact.optimal)
        self.assertEqual(exact.max_extra_bytes, 0)
        for budget in (0, 1, 5, 32):
            with self.subTest(budget=budget):
                result = optimize_storage(circuit, value, max_states=budget)
                self.assertLessEqual(result.states, budget)
                self.assertLessEqual(result.lower_bound_bytes, len(exact.data))
                self.assertLessEqual(len(result.data) - len(exact.data), result.max_extra_bytes)
        loose = optimize_storage(circuit, value, target_gap_percent=100)
        self.assertFalse(loose.optimal)
        self.assertEqual(loose.max_extra_bytes, len(loose.data) - len(exact.data))

    def test_sparse_tower_with_repeated_and_irregular_terms(self):
        circuit = Circuit()
        height = circuit.power(circuit.integer(65_536))
        block = circuit.integer((1 << 12) - (1 << 5) - 1)
        shifts = circuit.integer(sum(1 << (40 * index) for index in range(6)))
        repeated = multiply_fast(circuit, block, shifts)
        high = multiply_fast(circuit, circuit.power(height), repeated)
        high = circuit.add(high, circuit.power(circuit.add(height, circuit.integer(251))))
        value = circuit.add(high, circuit.integer((1 << 13) - 1))
        result = optimize_storage(circuit, value)
        self.assertTrue(result.optimal)
        self.assertEqual(len(result.data), 148)
        self.assertIs(decode_storage(circuit, result.data), value)

    def test_tower_storage_can_be_proved_without_search(self):
        circuit = Circuit()
        value = circuit.zero
        for _ in range(1_250):
            value = circuit.power(value)

        result = optimize_storage(circuit, value, max_states=0)
        self.assertTrue(result.optimal)
        self.assertEqual(result.states, 0)
        self.assertEqual(result.data, circuit.encode(value))
        self.assertIs(decode_storage(circuit, result.data), value)

    def test_invalid_signed_graph_reference(self):
        circuit = Circuit()
        with self.assertRaises(ValueError):
            decode_storage(circuit, b"H2N1\x01\x01\x00")
