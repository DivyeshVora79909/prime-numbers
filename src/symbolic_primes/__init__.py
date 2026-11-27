"""Exact arithmetic on recursively represented signed base-2 integers."""

from .adaptive import Adaptive, Node
from .quotient_kernel import QuotientWorkLimit, divide_fast
from .similarity_kernel import multiply_fast
from .storage import StorageResult, decode_storage, multiply_storage, optimize_storage
from .symbolic import Circuit, ExactQuotient, Expression, Geometric, Signed, exact_divide, residue

__all__ = (
    "Adaptive",
    "Circuit",
    "ExactQuotient",
    "Expression",
    "Geometric",
    "Node",
    "QuotientWorkLimit",
    "Signed",
    "StorageResult",
    "decode_storage",
    "divide_fast",
    "exact_divide",
    "multiply_fast",
    "multiply_storage",
    "optimize_storage",
    "residue",
)
