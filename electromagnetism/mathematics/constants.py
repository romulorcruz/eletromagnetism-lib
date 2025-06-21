"""Cosntants Module.

This module contains all constants used in this eletromagnetism package.
"""
from numpy import pi

# Defining constants for clarity's sake
X, Y, Z, BX, BY, BZ, B = 0, 1, 2, 3, 4, 5, 6

# This is MU0/4pi in T.m/A
MU0_PRIME = 1e-7

MU0 = MU0_PRIME * 4 * pi
