#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
#
# Copyright (c) 2026 Matthew Badger
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Verification of the computations in

    M. Badger, "An exact formula for optimal axis-parallel packing of
    eight squares".

All arithmetic is exact (Python integers).  The script uses only the
standard library.  Run

    python3 verify_eight_squares.py                 # all checks
    python3 verify_eight_squares.py --tex paper.tex # also compare the data
                                                    # below with the LaTeX source

and see README.md for a description of each check.  The exit status is 0
exactly when every check passes.

Conventions.  Side vectors are nonincreasing tuples s = (s_0,...,s_7) of
nonnegative integers.  A side sum is written as a string of indices, so
"047" means s_0 + s_4 + s_7.  Coordinates of lower left corners are linear
forms in s_0,...,s_7 and the container side L, written as in Table 3.1.

Section, table, figure, and equation numbers refer to the paper.
"""

import argparse
import functools
import itertools
import random
import re
import sys

N8 = 8
GREEK = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta",
         "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron",
         "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi"]
PIECES = {
    "alpha": ['047', '056', '127', '135'],
    "beta": ['037', '056', '136', '2567'],
    "gamma": ['037', '056', '126'],
    "delta": ['045', '127', '136'],
    "epsilon": ['056', '134', '2467'],
    "zeta": ['067', '124', '3456'],
    "eta": ['036', '127'],
    "theta": ['027', '056', '136'],
    "iota": ['047', '056', '135', '2567'],
    "kappa": ['056', '124', '3467'],
    "lambda": ['045', '136', '2567'],
    "mu": ['056', '135', '1467'],
    "nu": ['067', '137', '1456'],
    "xi": ['135', '1367', '2456'],
    "omicron": ['047', '125', '3456'],
    "pi": ['056', '1367'],
    "rho": ['056', '136', '1457'],
    "sigma": ['1467', '2356'],
    "tau": ['1456', '2367'],
    "upsilon": ['1567', '2346'],
    "phi": ['125', '1267', '3456'],
    "chi": ['0567', '136'],
    "psi": ['14567'],
}
COORDS = {
    "alpha": [('0', '0'), ('L-s1', 's5'), ('0', 'L-s2'), ('L-s3', 'L-s3'), ('s2', 'L-s4'), ('L-s5', '0'), ('s0', '0'), ('s2', 's0')],
    "beta": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 'L-s4'), ('L-s5', 's1'), ('L-s5-s6', 's1'), ('L-s5-s6-s7', 's0')],
    "gamma": [('0', '0'), ('L-s1', '0'), ('s3', 'L-s2'), ('0', 'L-s3'), ('L-s4', 'L-s4'), ('L-s5', 's1'), ('L-s5-s6', 's1'), ('0', 's0')],
    "delta": [('0', '0'), ('L-s1', '0'), ('s4', 'L-s2'), ('L-s3', 'L-s3'), ('0', 'L-s4'), ('0', 's0'), ('L-s6', 's1'), ('L-s6-s7', 's1')],
    "epsilon": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('L-s3', 's1'), ('L-s4', 'L-s4'), ('s2', 's0'), ('s2+s7', 'L-s6'), ('s2', 'L-s7')],
    "zeta": [('0', '0'), ('L-s1', '0'), ('L-s2', 's1'), ('s6', 'L-s3'), ('L-s4', 'L-s4'), ('s3+s6', 'L-s5'), ('0', 'L-s6'), ('0', 's0')],
    "eta": [('0', '0'), ('L-s1', '0'), ('s3', 'L-s2'), ('0', 'L-s3'), ('L-s4', 'L-s4'), ('L-s5', 's1'), ('0', 's0'), ('L-s5-s7', 's1')],
    "theta": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 'L-s4'), ('L-s5', 's1'), ('L-s5-s6', 's1'), ('0', 's0')],
    "iota": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('L-s3', 'L-s3'), ('s2', 'L-s4'), ('L-s5', 's1'), ('L-s5-s6', 's1'), ('s2', 's0')],
    "kappa": [('0', '0'), ('L-s1', '0'), ('L-s2', 's1'), ('0', 'L-s3'), ('L-s4', 'L-s4'), ('s3', 's0'), ('s3', 'L-s6'), ('s3+s6', 'L-s7')],
    "lambda": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('L-s3', 'L-s3'), ('s2', 'L-s4'), ('s2', 's0'), ('L-s6', 's1'), ('L-s6-s7', 's1')],
    "mu": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 'L-s4'), ('L-s5-s6', 's1'), ('L-s6', 's1'), ('L-s7', 's1+s6')],
    "nu": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 'L-s4'), ('L-s5', 's1+s6'), ('L-s6', 's1'), ('L-s6-s7', 's1')],
    "xi": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('L-s3', 's1'), ('s2', 'L-s4'), ('s2+s4', 'L-s5'), ('L-s6', 'L-s6'), ('L-s7', 's1+s3')],
    "omicron": [('0', '0'), ('L-s1', '0'), ('L-s2', 's1'), ('0', 'L-s3'), ('s3', 'L-s4'), ('L-s5', 'L-s5'), ('s3+s4', 'L-s6'), ('s3', 'L-s4-s7')],
    "pi": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 'L-s4'), ('L-s5', 's1'), ('L-s5-s6', 's1'), ('L-s5-s7', 's1+s6')],
    "rho": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 'L-s4'), ('L-s5', 's1'), ('L-s5-s6', 's1'), ('L-s7', 's1+s5')],
    "sigma": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 's1'), ('s2+s3', 'L-s5'), ('L-s6', 'L-s6'), ('L-s7', 's1+s4')],
    "tau": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 's1'), ('L-s5', 's1+s4'), ('L-s6', 'L-s6'), ('s2+s3', 'L-s7')],
    "upsilon": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('s2+s3', 'L-s4'), ('L-s5', 's1'), ('L-s6', 'L-s6'), ('L-s7', 's1+s5')],
    "phi": [('0', '0'), ('L-s1', '0'), ('L-s2', 's1'), ('0', 'L-s3'), ('s3', 'L-s4'), ('s3+s4', 'L-s5'), ('L-s6', 'L-s6'), ('L-s7', 's1+s2')],
    "chi": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 'L-s4'), ('L-s5', 's1'), ('L-s5-s6', 's1'), ('L-s5-s6-s7', 's1')],
    "psi": [('0', '0'), ('L-s1', '0'), ('0', 'L-s2'), ('s2', 'L-s3'), ('L-s4', 's1'), ('L-s5', 's1+s4'), ('L-s6', 's1+s4+s5'), ('L-s7', 'L-s7')],
}
CERT = {
    "alpha": ((21,18,15,13,12,10,9,8), 41),
    "beta": ((21,19,16,13,11,10,8,7), 41),
    "gamma": ((20,18,15,14,12,11,9,8), 42),
    "delta": ((85,83,69,67,53,52,39,38), 190),
    "epsilon": ((61,38,33,32,31,20,19,18), 101),
    "zeta": ((690,413,412,279,278,277,269,144), 1103),
    "eta": ((164,162,131,128,98,97,66,65), 358),
    "theta": ((71,69,54,52,50,37,35,32), 157),
    "iota": ((50,47,38,30,29,21,20,19), 98),
    "kappa": ((952,580,579,478,477,342,341,340), 1636),
    "lambda": ((26,23,20,17,13,12,10,9), 51),
    "mu": ((38,27,23,22,20,17,10,9), 66),
    "nu": ((195,194,162,161,66,65,64,6), 389),
    "xi": ((71,32,31,27,26,24,23,22), 104),
    "omicron": ((225,166,128,100,98,97,96,68), 391),
    "pi": ((39,27,24,23,19,18,9,8), 67),
    "rho": ((38,29,26,25,16,15,11,8), 68),
    "sigma": ((73,34,31,28,27,25,24,23), 108),
    "tau": ((527,212,211,181,177,176,175,173), 740),
    "upsilon": ((77,38,31,30,29,27,26,25), 116),
    "phi": ((1034,466,433,369,368,367,366,236), 1501),
    "chi": ((98,87,66,65,54,35,27,26), 186),
    "psi": ((265,264,231,230,68,67,66,65), 530),
}
ILLUSTRATIONS = {
    "alpha": ((21,18,15,13,12,10,9,8), 41),
    "beta": ((20,18,15,13,11,10,8,7), 40),
    "gamma": ((20,18,15,14,12,11,9,8), 42),
    "delta": ((23,21,18,16,13,12,10,9), 48),
    "epsilon": ((24,17,15,13,12,9,8,7), 42),
    "zeta": ((28,18,17,13,12,11,10,9), 47),
    "eta": ((24,22,19,16,14,13,10,9), 50),
    "theta": ((24,22,19,17,15,14,12,9), 52),
    "iota": ((25,22,18,15,14,11,10,9), 48),
    "kappa": ((24,16,15,14,13,10,9,8), 44),
    "lambda": ((26,23,20,17,13,12,10,9), 51),
    "mu": ((74,51,47,46,32,29,22,21), 126),
    "nu": ((51,50,46,45,22,21,20,19), 114),
    "xi": ((25,13,12,11,10,9,8,7), 39),
    "omicron": ((37,27,22,18,17,15,14,11), 65),
    "pi": ((47,33,28,27,25,22,11,10), 81),
    "rho": ((98,69,66,65,36,35,32,28), 168),
    "sigma": ((30,17,14,13,12,11,10,9), 48),
    "tau": ((32,17,16,15,12,11,10,9), 50),
    "upsilon": ((35,21,16,15,14,13,12,11), 57),
    "phi": ((42,21,20,18,17,16,12,11), 64),
    "chi": ((35,28,24,23,16,14,8,7), 64),
    "psi": ((41,40,35,34,12,11,10,9), 82),
}

# ---------------------------------------------------------------------------
# Further data from the paper
# ---------------------------------------------------------------------------

# Equation (1.1): L_6^* is the maximum of these side sums.
L6_TERMS = ["01", "234", "145"]

# Equation (2.5) and Table 4.1: recurring side sums.
NAMES = {
    "W": "3456", "H": "125", "C0": "056", "C3": "136", "R": "1456",
    "D": "137", "E": "047", "P": "14567",
    "G6": "036", "G7": "037", "K": "027", "T": "067", "U": "045",
    "H4": "134", "H5": "135", "J4": "124", "J6": "126",
    "J7": "127", "J67": "1267", "M": "1467", "Ppi": "1367",
    "Prho": "1457", "Pchi": "0567", "V": "2567", "X": "3467", "Y": "2467",
}

# Propositions 4.1-4.11: the nine lower bounds b_i, as minima of side sums.
# "KK" stands for the balanced-partition bound \mathcal K of Proposition 4.1.
BOUNDS = {
    1: ["067", "KK"],
    2: ["137", "14567"],
    3: ["056", "047", "124", "1267"],
    4: ["127", "2567"],
    5: ["037", "045", "135", "1367", "1457", "0567"],
    6: ["027", "036", "045", "135", "126", "2567"],
    7: ["047", "3467"],
    8: ["047", "124", "2467"],
    9: ["047", "134", "1467", "0567"],
}

# Table 5.1: the same bounds written with the names of Table 4.1.
BOUNDS_BY_NAME = {
    1: ["T", "KK"], 2: ["D", "P"], 3: ["C0", "E", "J4", "J67"],
    4: ["J7", "V"], 5: ["G7", "U", "H5", "Ppi", "Prho", "Pchi"],
    6: ["K", "G6", "U", "H5", "J6", "V"], 7: ["E", "X"],
    8: ["E", "J4", "Y"], 9: ["E", "H4", "M", "Pchi"],
}

# Lemma 5.1: the five balanced partitions (part containing 1, other part).
FIVE_PARTITIONS = {
    "xi": ("137", "245"), "sigma": ("147", "235"), "tau": ("145", "237"),
    "upsilon": ("157", "234"), "phi": ("127", "345"),
}
# Lemma 5.1, proof: each remaining partition (named by its part containing 1)
# is dominated by the indicated one of the five.
PARTITION_DOMINATION = {"123": "phi", "124": "phi", "125": "phi",
                        "134": "xi", "135": "xi"}

# The three seven-square strategies of the companion paper (its Table 5.1),
# with their costs from equation (1.2).
SEVEN_COORDS = {
    "top": [("0", "0"), ("L-s1", "0"), ("L-s2", "s1"), ("0", "L-s3"),
            ("s3", "L-s4"), ("L-s5", "L-s5"), ("s3+s4", "L-s6")],
    "center": [("0", "0"), ("L-s1", "0"), ("0", "L-s2"), ("s2", "L-s3"),
               ("L-s4", "L-s4"), ("L-s5", "s1"), ("L-s5-s6", "s1")],
    "right": [("0", "0"), ("L-s1", "0"), ("0", "L-s2"), ("s2", "L-s3"),
              ("L-s4", "L-s4"), ("L-s5", "s1"), ("L-s6", "s1+s5")],
}
SEVEN_COSTS = {"top": L6_TERMS + ["3456", "125"],
               "center": L6_TERMS + ["056", "136"],
               "right": L6_TERMS + ["1456"]}

# Section 3, after (3.1): restrictions with the cost of a seven-square strategy.
RESTRICTION_CLASS = {
    "beta": "center", "theta": "center", "pi": "center", "rho": "center",
    "chi": "center", "omicron": "top", "phi": "top",
    "nu": "right", "tau": "right", "psi": "right",
}

# Figures with explicit packings: side vector, container side, and corners.
FIGURE_1_1 = {
    "a": ("top", (22, 15, 14, 11, 10, 9, 8), 38,
          {0: (0, 0), 1: (23, 0), 2: (24, 15), 3: (0, 27), 4: (11, 28), 5: (29, 29), 6: (21, 30)}),
    "b": ("center", (19, 16, 13, 12, 10, 9, 8), 36,
          {0: (0, 0), 1: (20, 0), 2: (0, 23), 3: (13, 24), 4: (26, 26), 5: (27, 16), 6: (19, 16)}),
    "c": ("right", (26, 25, 21, 20, 10, 9, 8), 52,
          {0: (0, 0), 1: (27, 0), 2: (0, 31), 3: (21, 32), 4: (42, 42), 5: (43, 25), 6: (44, 34)}),
    "d": ("beta", (20, 18, 15, 13, 11, 10, 8, 7), 40,
          {0: (0, 0), 1: (22, 0), 2: (0, 25), 3: (15, 27), 4: (29, 29), 5: (30, 18), 6: (22, 18), 7: (15, 20)}),
    "e": ("alpha", (21, 18, 15, 13, 12, 10, 9, 8), 41,
          {0: (0, 0), 1: (23, 10), 2: (0, 26), 3: (28, 28), 4: (15, 29), 5: (31, 0), 6: (21, 0), 7: (15, 21)}),
}

FIGURE_2_1 = {  # Lemma 2.11, m = 6
    "a": ((12, 8, 7, 6, 5, 4, 3), 25,
          {0: (0, 0), 1: (17, 0), 2: (0, 12), 3: (11, 12), 4: (17, 8), 5: (7, 12), 6: (18, 13)}),
    "b": ((12, 8, 7, 6, 5, 4, 3), 20,
          {0: (0, 0), 1: (12, 0), 2: (0, 12), 3: (11, 13), 4: (13, 8), 5: (7, 12), 6: (17, 13)}),
}

FIGURE_4_1 = {  # Proposition 4.2
    "a": ((53, 29, 28, 20, 19, 17, 15, 11), 82,
          {0: (0, 0), 1: (53, 0), 2: (54, 29), 3: (0, 62), 4: (20, 63), 5: (39, 65), 6: (67, 67), 7: (56, 71)}),
    "b": ((44, 25, 23, 17, 15, 14, 12, 11), 69,
          {0: (0, 0), 1: (44, 0), 2: (46, 25), 3: (0, 52), 4: (17, 54), 5: (43, 55), 6: (57, 57), 7: (32, 58)}),
}

FIGURE_4_2 = {  # Lemma 4.4, Propositions 4.5 and 4.7
    "a": ((39, 32, 28, 25, 18, 15, 9, 8), 71,
          {0: (0, 0), 1: (39, 0), 2: (0, 43), 3: (28, 46), 4: (53, 53), 5: (56, 32), 6: (39, 32), 7: (48, 32)}),
    "b": ((33, 29, 26, 20, 16, 14, 10, 5), 62,
          {0: (0, 0), 1: (33, 0), 2: (0, 36), 3: (26, 42), 4: (46, 46), 5: (48, 29), 6: (33, 29), 7: (43, 29)}),
}

FIGURE_4_3 = {  # Lemma 4.6
    "": ((28, 26, 24, 22, 8, 7, 6, 4), 54,
         {0: (22, 0), 1: (24, 28), 2: (0, 30), 3: (0, 8), 4: (0, 0), 5: (8, 0), 6: (15, 0), 7: (50, 0)}),
}

FIGURE_4_4 = {  # Lemma 4.8
    "a": ((37, 20, 19, 18, 17, 6, 5, 4), 57,
          {0: (0, 0), 1: (37, 0), 2: (38, 20), 3: (39, 39), 4: (0, 40), 5: (17, 51), 6: (23, 52), 7: (28, 53)}),
    "b": ((37, 20, 19, 18, 17, 6, 5, 4), 57,
          {0: (0, 0), 1: (37, 0), 2: (38, 20), 3: (0, 39), 4: (18, 40), 5: (35, 39), 6: (52, 39), 7: (53, 44)}),
    "c": ((46, 30, 29, 22, 21, 12, 11, 10), 76,
          {0: (0, 0), 1: (46, 0), 2: (47, 30), 3: (0, 54), 4: (22, 55), 5: (43, 59), 6: (55, 59), 7: (66, 59)}),
}

# ---------------------------------------------------------------------------
# Side sums, linear forms, and the coefficient test
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def idx(term):
    """'047' -> (0, 4, 7)."""
    return tuple(int(c) for c in term)


def vec(term, n=N8):
    v = [0] * n
    for i in idx(term):
        v[i] += 1
    return v


def ssum(s, term):
    return sum(s[i] for i in idx(term))


def prefix_nonneg(c):
    """Coefficient test (2.4): sum_i c_i s_i >= 0 for every nonincreasing
    nonnegative s if and only if every prefix sum of c is nonnegative."""
    acc = 0
    for x in c:
        acc += x
        if acc < 0:
            return False
    return True


def dominates(big, small, n=N8):
    """big >= small on every nonincreasing nonnegative tuple (side sums)."""
    return prefix_nonneg([a - b for a, b in zip(vec(big, n), vec(small, n))])


def parse_form(expr, n):
    """'L-s5-s6' -> coefficient list of length n+1; the last entry is L."""
    v = [0] * (n + 1)
    e = expr.replace(" ", "")
    if e == "0":
        return v
    pos = 0
    for m in re.finditer(r"([+-]?)(L|s\d)", e):
        if m.start() != pos:
            raise ValueError("cannot parse " + expr)
        pos = m.end()
        sign = -1 if m.group(1) == "-" else 1
        if m.group(2) == "L":
            v[n] += sign
        else:
            v[int(m.group(2)[1:])] += sign
    if pos != len(e):
        raise ValueError("cannot parse " + expr)
    return v


def eval_form(f, s, L):
    n = len(f) - 1
    return sum(f[i] * s[i] for i in range(n)) + f[n] * L


def certify(f, cost_terms, n):
    """Certify a*L + c.s >= 0 for every nonincreasing nonnegative s and every
    L >= max of cost_terms(s).  Returns a description of the certificate, or
    None.  As in Section 2, the certificate adds a times one cost term and
    applies the coefficient test."""
    a, c = f[n], f[:n]
    if a < 0:
        return None
    if a == 0:
        return "prefix" if prefix_nonneg(c) else None
    for t in cost_terms:
        tv = vec(t, n)
        if prefix_nonneg([a * tv[i] + c[i] for i in range(n)]):
            return t
    return None


def reduce_terms(terms, n):
    """Remove side sums dominated by another one; the maximum is unchanged."""
    terms = sorted(set(terms), key=lambda t: (len(t), t))
    keep = []
    for a in terms:
        if not any(b != a and dominates(b, a, n) for b in terms):
            keep.append(a)
    return sorted(keep)


# ---------------------------------------------------------------------------
# The formulas
# ---------------------------------------------------------------------------

def L6(s):
    return max(ssum(s, t) for t in L6_TERMS)


def L7(s):
    """Equation (1.2), equivalently (2.6)."""
    n = {k: ssum(s, v) for k, v in NAMES.items() if k in ("W", "H", "C0", "C3", "R")}
    return max(L6(s), min(max(n["W"], n["H"]), max(n["C0"], n["C3"]), n["R"]))


def piece(s, w):
    return max(ssum(s, t) for t in PIECES[w])


def ell(s, w):
    """Strategy cost ell_omega = max(L_6^*, omega) (Theorem 1.1)."""
    return max(L6(s), piece(s, w))


def F8(s):
    """Theorem 1.1, first form: max(L_7^*, min_omega omega)."""
    return max(L7(s), min(piece(s, w) for w in GREEK))


PARTITIONS = []
for _I in itertools.combinations((1, 2, 3, 4, 5, 7), 3):
    if 1 in _I:
        _J = tuple(i for i in (1, 2, 3, 4, 5, 7) if i not in _I)
        PARTITIONS.append(("".join(map(str, _I)), "".join(map(str, _J))))


def KK(s):
    """Balanced-partition bound \\mathcal K of Proposition 4.1."""
    return min(max(s[6] + ssum(s, I), s[6] + ssum(s, J)) for I, J in PARTITIONS)


def nonincreasing(s):
    return all(s[i] >= s[i + 1] for i in range(len(s) - 1)) and s[-1] >= 0


# ---------------------------------------------------------------------------
# Packings and strategies
# ---------------------------------------------------------------------------

def check_packing(s, L, X):
    """Exact containment and disjointness of interiors.  Returns an error
    string or None."""
    for i, (x, y) in X.items():
        if not (0 <= x <= L - s[i] and 0 <= y <= L - s[i]):
            return "containment fails for Q_%d" % i
    for i, j in itertools.combinations(sorted(X), 2):
        (xi, yi), (xj, yj) = X[i], X[j]
        if not (xi + s[i] <= xj or xj + s[j] <= xi or
                yi + s[i] <= yj or yj + s[j] <= yi):
            return "Q_%d and Q_%d overlap" % (i, j)
    return None


def left(s, X, i, j):
    return X[i][0] + s[i] <= X[j][0]


def below(s, X, i, j):
    return X[i][1] + s[i] <= X[j][1]


def xov(s, X, i, j):
    return not left(s, X, i, j) and not left(s, X, j, i)


def yov(s, X, i, j):
    return not below(s, X, i, j) and not below(s, X, j, i)


def hchain(s, X, I):
    return all(not xov(s, X, i, j) for i, j in itertools.combinations(I, 2))


def vchain(s, X, I):
    return all(not yov(s, X, i, j) for i, j in itertools.combinations(I, 2))


def coords_at(table_row, s, L):
    n = len(table_row)
    return {i: (eval_form(parse_form(a, n), s, L), eval_form(parse_form(b, n), s, L))
            for i, (a, b) in enumerate(table_row)}


def transitive_closure(R):
    R = set(R)
    changed = True
    while changed:
        changed = False
        for (a, b) in list(R):
            for (c, d) in list(R):
                if b == c and (a, d) not in R:
                    R.add((a, d))
                    changed = True
    return R


def uniform_strategy(table_row, cost_terms):
    """Lemma 2.6 for the coordinates in table_row and L >= max(cost_terms).
    Certifies containment and, for each ordered pair, the separations that
    hold uniformly.  Returns (H, V, problems)."""
    n = len(table_row)
    forms = [(parse_form(a, n), parse_form(b, n)) for a, b in table_row]
    problems = []
    for i, (fx, fy) in enumerate(forms):
        for k, f in enumerate((fx, fy)):
            lo = f
            hi = [-x for x in f]
            hi[n] += 1
            hi[i] -= 1
            if certify(lo, cost_terms, n) is None or certify(hi, cost_terms, n) is None:
                problems.append("containment of Q_%d (%s) not certified" % (i, "xy"[k]))
    H, V = set(), set()
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for k, R in ((0, H), (1, V)):
                f = [forms[j][k][t] - forms[i][k][t] for t in range(n + 1)]
                f[i] -= 1
                if certify(f, cost_terms, n) is not None:
                    R.add((i, j))
    for i, j in itertools.combinations(range(n), 2):
        if not ({(i, j), (j, i)} & (H | V)):
            problems.append("no certified separation for Q_%d, Q_%d" % (i, j))
    return transitive_closure(H), transitive_closure(V), problems


def chains(H, V, idxs):
    """All nonempty chains of the strategy (H, V) within idxs."""
    out = []
    for r in range(1, len(idxs) + 1):
        for I in itertools.combinations(idxs, r):
            for R in (H, V):
                if all((a, b) in R or (b, a) in R for a, b in itertools.combinations(I, 2)):
                    out.append("".join(map(str, I)))
                    break
    return out


def is_chain(H, V, term):
    I = idx(term)
    return any(all((a, b) in R or (b, a) in R for a, b in itertools.combinations(I, 2))
               for R in (H, V))


def sorted_tuples(n, top):
    """All nonincreasing n-tuples with entries in {0,...,top}."""
    for c in itertools.combinations_with_replacement(range(top, -1, -1), n):
        yield c

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

class Report:
    def __init__(self, verbose=False):
        self.failures = 0
        self.verbose = verbose
        self.section = None

    def begin(self, title):
        self.section = title
        print("\n== " + title)

    def ok(self, cond, msg):
        if not cond:
            self.failures += 1
            print("  FAIL  " + msg)
        elif self.verbose:
            print("  ok    " + msg)
        return cond

    def note(self, msg):
        print("  " + msg)


# ---------------------------------------------------------------------------
# Check A: consistency of the formulas and of the summary tables
# ---------------------------------------------------------------------------

def check_formulas(rep):
    rep.begin("A. Formulas and summary tables (Table 1.1, Table 4.1, Table 5.1)")
    rep.ok(sorted(PIECES) == sorted(GREEK), "Table 1.1 lists the twenty-three maxima")
    seven_names = {"W": "3456", "H": "125", "C0": "056", "C3": "136", "R": "1456"}
    rep.ok(all(NAMES[k] == v for k, v in seven_names.items()),
           "names W, H, C_0, C_3, R agree with (2.5) and (1.2)")
    for i in range(1, 10):
        exp = sorted("KK" if t == "KK" else NAMES[t] for t in BOUNDS_BY_NAME[i])
        rep.ok(exp == sorted(BOUNDS[i]), "Table 5.1 row b_%d agrees with Proposition 4.x" % i)
    rep.ok(len(PARTITIONS) == 10, "ten balanced 3+3 partitions of {1,2,3,4,5,7}")
    rep.note("checked %d recurring names and nine lower bounds" % len(NAMES))


# ---------------------------------------------------------------------------
# Check B: the constructions (Proposition 3.1, Table 3.1)
# ---------------------------------------------------------------------------

def strategy_report(rep, label, row, cost_terms, piece_terms, n):
    H, V, problems = uniform_strategy(row, cost_terms)
    for p in problems:
        rep.ok(False, "%s: %s" % (label, p))
    # every side sum in the maximum is a chain; the six-square sums are
    # chains or dominated by a side sum of the maximum
    for t in piece_terms:
        rep.ok(is_chain(H, V, t), "%s: side sum %s is a chain" % (label, t))
    for t in L6_TERMS:
        ok = is_chain(H, V, t) or any(dominates(u, t, n) for u in piece_terms)
        rep.ok(ok, "%s: six-square sum %s is a chain or dominated" % (label, t))
    # every chain sum is dominated by a single cost term (upper bound on cost)
    for c in chains(H, V, list(range(n))):
        ok = any(dominates(t, c, n) for t in cost_terms)
        rep.ok(ok, "%s: chain %s dominated by a cost term" % (label, c))
    return H, V, not problems


def check_constructions(rep):
    rep.begin("B. Constructions: Table 3.1 and Proposition 3.1")
    strategies = {}
    for w in GREEK:
        H, V, _ = strategy_report(rep, w, COORDS[w], L6_TERMS + PIECES[w], PIECES[w], N8)
        strategies[w] = (H, V)
    rep.note("23 coordinate rows: containment and a separation for every pair "
             "certified by the coefficient test; cost = max(L_6^*, omega)")
    for k in SEVEN_COORDS:
        extra = [t for t in SEVEN_COSTS[k] if t not in L6_TERMS]
        strategy_report(rep, "seven-square " + k, SEVEN_COORDS[k], SEVEN_COSTS[k], extra, 7)
    rep.note("seven-square strategies top, center, right: costs as in (1.2)")
    return strategies


# ---------------------------------------------------------------------------
# Check C: the illustrations in Figures 3.1 and 3.2
# ---------------------------------------------------------------------------

def tight_terms(s, w, L):
    return sorted(t for t in PIECES[w] if ssum(s, t) == L)


def check_illustrations(rep):
    rep.begin("C. Figures 3.1 and 3.2")
    for w in GREEK:
        s, L = ILLUSTRATIONS[w]
        rep.ok(nonincreasing(s), "%s: side vector is nonincreasing" % w)
        rep.ok(F8(s) == L, "%s: L_8^*(s) = %d" % (w, L))
        rep.ok(L7(s) < L, "%s: L_8^*(s) > L_7^*(s)" % w)
        others = [v for v in GREEK if v != w and ell(s, v) <= L]
        rep.ok(ell(s, w) == L and not others,
               "%s: the only one of the 23 strategy costs equal to L_8^*" % w)
        X = coords_at(COORDS[w], s, L)
        err = check_packing(s, L, X)
        rep.ok(err is None, "%s: coordinates give a packing (%s)" % (w, err or "ok"))
        for t in tight_terms(s, w, L):
            ok = hchain(s, X, idx(t)) or vchain(s, X, idx(t))
            rep.ok(ok, "%s: tight side sum %s is a chain of the drawn packing" % (w, t))
    rep.note("23 illustrated side vectors: value, uniqueness, packing, tight chains")


# ---------------------------------------------------------------------------
# Check D: restrictions to the seven largest squares (Section 1, (3.1))
# ---------------------------------------------------------------------------

def restricted_terms(w):
    """Equation (3.1): side sums of max(L_6^*, omega(s_0,...,s_6,0))."""
    terms = list(L6_TERMS) + [t.replace("7", "") for t in PIECES[w]]
    return reduce_terms([t for t in terms if t], 7)


def restricted_cost(s, w):
    return max(ssum(s, t) for t in restricted_terms(w))


def check_restrictions(rep, strategies):
    rep.begin("D. Restrictions to Q_0,...,Q_6: equation (3.1) and Section 1")
    canon = {k: reduce_terms(v, 7) for k, v in SEVEN_COSTS.items()}
    gens = [[1] * k + [0] * (7 - k) for k in range(1, 8)]
    small = []
    for c in itertools.product(range(4), repeat=7):
        small.append(tuple(sum(c[j] * gens[j][i] for j in range(7)) for i in range(7)))
    count_exceed = 0
    for w in GREEK:
        # (3.1) agrees with the restriction of the uniform strategy
        H, V = strategies[w]
        H7 = {p for p in H if 7 not in p}
        V7 = {p for p in V if 7 not in p}
        from_coords = reduce_terms(chains(H7, V7, list(range(7))), 7)
        rep.ok(from_coords == restricted_terms(w),
               "%s: (3.1) agrees with the restricted strategy of Table 3.1" % w)
        match = [k for k in canon if canon[k] == restricted_terms(w)]
        expected = RESTRICTION_CLASS.get(w)
        rep.ok(match == ([expected] if expected else []),
               "%s: restriction has the cost of %s" % (w, expected or "none of the three"))
        if not expected:
            for k in canon:
                wit = next((v for v in small
                            if max(ssum(v, t) for t in canon[k]) != restricted_cost(v, w)), None)
                rep.ok(wit is not None, "%s: restricted cost differs from %s (witness %s)" % (w, k, wit))
        s, L = ILLUSTRATIONS[w]
        r, l7 = restricted_cost(s[:7], w), L7(s)
        if not expected or w == "tau":
            rep.ok(r == l7 + 1 == L, "%s: restricted cost = L_7^* + 1 = L_8^* at the illustration" % w)
            count_exceed += 1
        else:
            rep.ok(r == l7, "%s: restricted cost = L_7^* at the illustration" % w)
    rep.ok(count_exceed == 14, "fourteen restrictions exceed L_7^* at their illustrations")
    # the alpha example in Section 1 and Figure 1.1(e)
    rep.ok(restricted_terms("alpha") == reduce_terms(L6_TERMS + ["056", "135"], 7),
           "alpha: restricted cost = max(L_6^*, s0+s5+s6, s1+s3+s5)")
    s = ILLUSTRATIONS["alpha"][0]
    rep.ok((restricted_cost(s[:7], "alpha"), L7(s), F8(s)) == (41, 40, 41),
           "alpha example: 41 = L_8^*(s) and L_7^*(s) = 40")
    H, V = strategies["alpha"]
    rep.ok({(5, 1), (1, 3)} <= V and ssum(s, "135") == 41,
           "alpha example: 5, 1, 3 is a vertical chain with side sum 41")
    rep.ok(COORDS["beta"][:7] == SEVEN_COORDS["center"],
           "beta: deleting Q_7 leaves the arrangement center")
    rep.note("ten restrictions with a seven-square cost, thirteen without; tau as stated")


# ---------------------------------------------------------------------------
# Check E: the minimality certificate (Table 5.2, Theorem 1.1)
# ---------------------------------------------------------------------------

def check_minimality(rep):
    rep.begin("E. Minimality certificate: Table 5.2 and (5.2)")
    for w in GREEK:
        p, L = CERT[w]
        rep.ok(nonincreasing(p), "%s: p_omega is nonincreasing" % w)
        rep.ok(F8(p) == L, "%s: L_8^*(p_omega) = %d" % (w, L))
        others = [v for v in GREEK if v != w and ell(p, v) <= L]
        rep.ok(ell(p, w) == L and not others, "%s: only strategy with side L_8^*(p_omega)" % w)
    diffs = []
    for v, w in itertools.combinations(GREEK, 2):
        p, q = CERT[v][0], CERT[w][0]
        d = F8(tuple(a + b for a, b in zip(p, q))) - F8(p) - F8(q)
        diffs.append(d)
        rep.ok(d > 0, "pair %s, %s is strictly superadditive (difference %d)" % (v, w, d))
    rep.ok(len(diffs) == 253, "253 pairs")
    rep.ok(min(diffs) == 1 and max(diffs) == 135, "differences range from 1 to 135")
    pa, pb = CERT["alpha"][0], CERT["beta"][0]
    rep.ok(F8(tuple(a + b for a, b in zip(pa, pb))) == 83, "L_8^*(p_alpha + p_beta) = 83")
    rep.note("253 differences, minimum %d, maximum %d" % (min(diffs), max(diffs)))


# ---------------------------------------------------------------------------
# Check F: arithmetic of Section 5 (exact parts and a finite sanity check)
# ---------------------------------------------------------------------------

def check_arithmetic(rep, grid, samples, seed):
    rep.begin("F. Section 5: Lemma 5.1 exactly; Theorem 5.3 and Theorem 1.1 on test vectors")
    # Lemma 5.1, exact: each other partition is dominated by one of the five
    rep.ok(sorted(I for I, J in PARTITIONS) ==
           sorted(list(PARTITION_DOMINATION) + [FIVE_PARTITIONS[w][0] for w in FIVE_PARTITIONS]),
           "the ten partitions are the five of Lemma 5.1 and the five dominated ones")
    for I, w in PARTITION_DOMINATION.items():
        I2, J2 = FIVE_PARTITIONS[w]
        ok = dominates(I, I2) and dominates(I, J2)
        rep.ok(ok, "partition with part %s is dominated by the %s partition" % (I, w))
    for w, (I, J) in FIVE_PARTITIONS.items():
        parts = {"".join(sorted(I + "6")), "".join(sorted(J + "6"))}
        rep.ok(parts <= set(PIECES[w]),
               "maximum %s contains the two part-sums of its partition" % w)
    # Finite checks.  These are consistency checks, not proofs: the proofs
    # are in the paper.
    rng = random.Random(seed)
    vectors = list(sorted_tuples(N8, grid))
    for _ in range(samples):
        top = rng.choice([10, 100, 1000, 10 ** 6])
        vectors.append(tuple(sorted((rng.randint(0, top) for _ in range(N8)), reverse=True)))
    bad_forms = bad_cov = bad_lb = bad_trunc = 0
    for s in vectors:
        l6, l7, k = L6(s), L7(s), KK(s)
        pw = [piece(s, w) for w in GREEK]
        f = max(l7, min(pw))
        if f != min(max(l6, x) for x in pw):
            bad_forms += 1
        bs = [min(k if t == "KK" else ssum(s, t) for t in BOUNDS[i]) for i in range(1, 10)]
        lb = max([l7] + bs)
        if lb < f:
            bad_cov += 1
        if lb > f:
            bad_lb += 1
        T = ssum(s, "067")
        if l7 < T and min(f, T) != max(l7, min(k, T)):
            bad_trunc += 1
    rep.ok(bad_forms == 0, "Theorem 1.1: the two forms agree (%d failures)" % bad_forms)
    rep.ok(bad_cov == 0, "Theorem 5.3: max(L_7^*, b_1,...,b_9) >= F_8 (%d failures)" % bad_cov)
    rep.ok(bad_lb == 0, "each b_i <= F_8, consistent with Theorem 1.1 (%d failures)" % bad_lb)
    rep.ok(bad_trunc == 0, "Lemma 5.2, identity (5.1) (%d failures)" % bad_trunc)
    rep.note("%d test vectors (all with entries at most %d, plus %d random)"
             % (len(vectors), grid, samples))


# ---------------------------------------------------------------------------
# Check G: the other figures with explicit packings
# ---------------------------------------------------------------------------

def arms(s, X, idxs):
    R = {i for i in idxs if X[i][0] >= s[0]}
    U = {i for i in idxs if X[i][1] >= s[0]}
    return R, U


def check_other_figures(rep):
    rep.begin("G. Figures 1.1, 2.1, 4.1, 4.2, 4.3, 4.4")
    S = ssum
    # Figure 1.1
    for key, (name, s, L, X) in FIGURE_1_1.items():
        err = check_packing(s, L, X)
        rep.ok(err is None, "Fig 1.1(%s): packing (%s)" % (key, err or "ok"))
        if len(s) == 7:
            rep.ok(L7(s) == L, "Fig 1.1(%s): L_7^*(s) = %d" % (key, L))
            rep.ok(X == coords_at(SEVEN_COORDS[name], s, L), "Fig 1.1(%s): follows %s" % (key, name))
        else:
            rep.ok(F8(s) == L and (s, L) == ILLUSTRATIONS[name], "Fig 1.1(%s): L_8^*(s) = %d" % (key, L))
            rep.ok(X == coords_at(COORDS[name], s, L), "Fig 1.1(%s): follows sigma_%s" % (key, name))
    # Figure 2.1: arms share Q_6, covering chains as in the caption
    for key, (s, L, X) in FIGURE_2_1.items():
        err = check_packing(s, L, X)
        rep.ok(err is None, "Fig 2.1(%s): packing (%s)" % (key, err or "ok"))
        R, U = arms(s, X, range(1, 7))
        rep.ok(X[0] == (0, 0) and R | U == set(range(1, 7)) and R & U == {6},
               "Fig 2.1(%s): Q_0 cornered, arms cover 1..6 and share 6" % key)
        rep.ok(all(xov(s, X, i, j) for i, j in itertools.combinations(R, 2)) and vchain(s, X, R)
               and all(yov(s, X, i, j) for i, j in itertools.combinations(U, 2)) and hchain(s, X, U),
               "Fig 2.1(%s): arm structure of Lemma 2.10" % key)
        if key == "a":
            C1, C2 = R, (U - {6}) | {4}
            rep.ok(left(s, X, 3, 4), "Fig 2.1(a): r = 4 right of u = 3")
        else:
            C1, C2 = (R - {6}) | {3}, U
            rep.ok(below(s, X, 4, 3), "Fig 2.1(b): u = 3 above r = 4")
        rep.ok(vchain(s, X, C1) and hchain(s, X, C2) and C1 | C2 == set(range(1, 7))
               and len(C1 & C2) == 1 and 6 not in C1 & C2,
               "Fig 2.1(%s): covering chains with common index other than 6" % key)
    # Figure 4.1 (Proposition 4.2)
    stated = {"a": (85, 83, 83), "b": (70, 70, 71)}
    for key, (s, L, X) in FIGURE_4_1.items():
        rep.ok(check_packing(s, L, X) is None, "Fig 4.1(%s): packing" % key)
        vals = (S(s, "056"), S(s, "047"), S(s, "1267"))
        rep.ok(vals == stated[key] and min(vals) > L >= S(s, "124"),
               "Fig 4.1(%s): C_0, E, J_67 = %s exceed L = %d, and J_4 <= L" % (key, vals, L))
        R, U = arms(s, X, range(1, 8))
        t, u = 2, 5
        ok = {1, 2} <= R - U and below(s, X, 1, 2) and {3, 4} <= U - R
        ok = ok and left(s, X, 3, t) and left(s, X, 4, t) and 6 in R and below(s, X, t, 6)
        ok = ok and u in U - R and below(s, X, t, u) and X[u][0] < s[0] <= X[t][0] < X[u][0] + s[u]
        ok = ok and yov(s, X, 5, 6)
        if key == "a":
            ok = ok and X[7][0] >= s[0] and xov(s, X, 7, 1) and xov(s, X, 7, 2) \
                and below(s, X, 2, 7) and hchain(s, X, [3, 4, 5, 6, 7])
        else:
            ok = ok and xov(s, X, 7, 0) and below(s, X, 0, 7) and not xov(s, X, 7, t) \
                and hchain(s, X, [3, 4, 7, t])
        rep.ok(ok, "Fig 4.1(%s): configuration described in the caption" % key)
    # Figure 4.2 (Lemma 4.4)
    lemma_families = ["123", "015", "135", "017", "2345", "045", "1457"]
    stated = {"a": {"G7": 72, "U": 72, "H5": 72, "Ppi": 74, "Prho": 73, "Pchi": 71},
              "b": {"K": 64, "G6": 63, "U": 63, "H5": 63, "J6": 65, "V": 55}}
    below_L = {"a": "Pchi", "b": "V"}
    for key, (s, L, X) in FIGURE_4_2.items():
        rep.ok(check_packing(s, L, X) is None, "Fig 4.2(%s): packing" % key)
        rep.ok(min(S(s, t) for t in lemma_families) > L,
               "Fig 4.2(%s): hypotheses of Lemma 4.4 hold" % key)
        vals = {k: S(s, NAMES[k]) for k in stated[key]}
        rep.ok(vals == stated[key] and all(v > L for k, v in vals.items() if k != below_L[key])
               and vals[below_L[key]] <= L,
               "Fig 4.2(%s): values in the caption" % key)
        c, p, b, t = 3, 2, 5, 4
        ok = X[c][0] < s[0] <= X[c][1] and left(s, X, p, c) and xov(s, X, 4, 5) and below(s, X, b, t)
        ok = ok and all(left(s, X, c, m) and X[m][1] >= s[1] for m in (4, 5)) and X[b][1] < s[0] <= X[b][0]
        for e in (6, 7):
            ok = ok and X[e][1] < s[0] <= X[e][0] and xov(s, X, e, c) and below(s, X, e, c) \
                and yov(s, X, e, b) and not yov(s, X, e, t) and left(s, X, e, b)
        ok = ok and yov(s, X, 6, 7) and yov(s, X, b, c) and all(X[j][1] >= s[1] for j in range(2, 8))
        bar = [0, 6, 7, 5] if key == "a" else [2, 6, 7, 5]
        ok = ok and hchain(s, X, bar) and sum(s[i] for i in bar) <= L
        rep.ok(ok, "Fig 4.2(%s): right-hand stack and the final chain" % key)
    # Figure 4.3 (Lemma 4.6)
    s, L, X = FIGURE_4_3[""]
    rep.ok(check_packing(s, L, X) is None, "Fig 4.3: packing")
    vals = {"K": S(s, "027"), "G6": S(s, "036"), "H5": S(s, "135"), "J6": S(s, "126"),
            "123": S(s, "123"), "U": S(s, "045"), "V": S(s, "2567")}
    rep.ok(vals == {"K": 56, "G6": 56, "H5": 55, "J6": 56, "123": 72, "U": 43, "V": 41}
           and min(vals[k] for k in ("K", "G6", "H5", "J6", "123")) > L >= max(vals["U"], vals["V"]),
           "Fig 4.3: values in the caption")
    ok = X[0][1] == 0 and left(s, X, 3, 0) and yov(s, X, 3, 0) and left(s, X, 0, 7) and yov(s, X, 0, 7)
    ok = ok and all(xov(s, X, i, 0) and below(s, X, 0, i) for i in (1, 2)) and left(s, X, 2, 1)
    ok = ok and xov(s, X, 3, 2) and below(s, X, 3, 2)
    for m in (4, 5):
        ok = ok and yov(s, X, m, 0) and left(s, X, m, 0) and xov(s, X, m, 3) \
            and below(s, X, m, 3) and below(s, X, m, 2)
    ok = ok and yov(s, X, 4, 5) and hchain(s, X, [4, 5, 0])
    rep.ok(ok, "Fig 4.3: double blockage described in the caption")
    # Figure 4.4 (Lemma 4.8)
    core = (1, 2, 3, 4)
    for key, (s, L, X) in FIGURE_4_4.items():
        rep.ok(check_packing(s, L, X) is None, "Fig 4.4(%s): packing" % key)
        C = min(S(s, "".join(map(str, sorted((0, b, j))))) for b in core for j in range(1, 8) if j != b)
        Q = min(sum(s[i] for i in I) for I in itertools.combinations(range(1, 8), 4)
                if len(set(I) & set(core)) >= 3)
        rep.ok(C > L and Q > L, "Fig 4.4(%s): exclusions (C) and (Q) hold" % key)
        R, U = arms(s, X, range(1, 8))
        if key == "a":
            rep.ok((S(s, "047"), S(s, "2347")) == (58, 58), "Fig 4.4(a): values in the caption")
            k, t = 4, 3
            ok = {1, 2, 3} <= R and vchain(s, X, [1, 2, 3]) and {5, 6, 7} <= U - R and k in U
            ok = ok and all(not xov(s, X, a, b) for a in (k, t) for b in (5, 6, 7)) and not xov(s, X, k, t)
            rep.ok(ok and vchain(s, X, [1, 2, 3]), "Fig 4.4(a): Alternative I with t = 3, k = 4")
        else:
            b_, t, p, q = 1, 2, 3, 4
            ok = {b_, t} <= R - U and {p, q} <= U - R and left(s, X, q, t) and hchain(s, X, [p, q, t])
            ok = ok and all(xov(s, X, e, t) and not yov(s, X, e, b_) and not yov(s, X, e, t)
                            for e in (5, 6, 7))
            if key == "b":
                rep.ok((S(s, "047"), S(s, "2347")) == (58, 58), "Fig 4.4(b): values in the caption")
            else:
                rep.ok((S(s, "047"), S(s, "2347"), S(s, "1267")) == (77, 82, 80) and S(s, "1267") > L,
                       "Fig 4.4(c): values in the caption")
                ok = ok and hchain(s, X, [5, 6, 7]) and all(yov(s, X, a, b) for a, b in
                                                             itertools.combinations((5, 6, 7), 2))
                ok = ok and all(below(s, X, t, e) for e in (5, 6, 7))
                ok = ok and all(left(s, X, a, e) for a in (p, q) for e in (6, 7))
                ok = ok and hchain(s, X, [3, 4, 5, 6, 7])
            rep.ok(ok, "Fig 4.4(%s): Alternative II as in the caption" % key)
    rep.note("all drawn packings are exact packings with the stated properties")

# ---------------------------------------------------------------------------
# Check H (optional): the data above agree with the LaTeX source
# ---------------------------------------------------------------------------

TEX_GREEK = {"epsilon": "varepsilon", "phi": "varphi"}


def tex_name(w):
    return TEX_GREEK.get(w, w)


def tex_sums(text):
    """'s_0+s_4+s_7, s_1+s_2' -> ['047', '12']"""
    return ["".join(re.findall(r"s_(\d)", z)) for z in re.findall(r"s_\d(?:\+s_\d)*", text)]


def tex_float(t, label):
    k = t.index(r"\label{%s}" % label)
    env = "figure" if label.startswith("fig") else "table"
    i = t.rindex(r"\begin{%s}" % env, 0, k)
    j = t.index(r"\end{%s}" % env, k)
    return t[i:j]


def tex_panels(fig):
    """Squares, container side, bars, and caption of each TikZ panel."""
    out = []
    for p in re.split(r"\\begin\{minipage\}", fig)[1:]:
        sq = re.findall(r"\\filldraw\[fill=black!\d+,draw=black,line width=[.\d]+pt\] "
                        r"\(([-\d.]+),([-\d.]+)\) rectangle \(([-\d.]+),([-\d.]+)\);\s*"
                        r"\\node\[[^\]]*\] at \([-\d.,]+\) \{\$(\d)\$\}", p)
        box = re.findall(r"\\draw\[line width=[.\d]+pt\] \(0,0\) rectangle \(([\d.]+),([\d.]+)\)", p)
        bars = re.findall(r"\\filldraw\[fill=(hc|vc)!\d+,draw=[^\]]*\] \(([-\d.]+),([-\d.]+)\) "
                          r"rectangle \(([-\d.]+),([-\d.]+)\);\s*\\node\[[^\]]*\] at \([-\d.,]+\) "
                          r"\{\$(\d)\$\}", p)
        X = {int(i): (float(x0), float(y0)) for x0, y0, x1, y1, i in sq}
        sides = {int(i): float(x1) - float(x0) for x0, y0, x1, y1, i in sq}
        groups = {}
        for kind, x0, y0, x1, y1, i in bars:
            key = (kind, y0 if kind == "hc" else x0)
            groups.setdefault(key, []).append(int(i))
        out.append(dict(X=X, sides=sides, L=float(box[0][0]) if box else None,
                        bars=groups, text=p))
    return out


def as_int_map(X):
    return {i: (int(x), int(y)) for i, (x, y) in X.items() if x == int(x) and y == int(y)}


def check_tex(rep, path):
    rep.begin("H. Comparison with the LaTeX source " + path)
    t = open(path, encoding="utf-8").read()
    # Table 1.1
    tab = t[t.index(r"\label{tab:pieces}"):]
    tab = tab[:tab.index(r"\end{longtable}")]
    for w in GREEK:
        m = re.search(r"\$\\%s\$&\$(.*?)\$\\\\" % tex_name(w), tab)
        rep.ok(m is not None and tex_sums(m.group(1)) == PIECES[w], "Table 1.1 row %s" % w)
    # Table 3.1
    ct = t[t.index(r"\label{tab:coords8}"):]
    ct = ct[:ct.index(r"\end{longtable}")]
    for w in GREEK:
        m = re.search(r"\$\\%s\$\s*&(.*?)\\\\\*?\s*\n\s*&(.*?)\\\\" % tex_name(w), ct, re.S)
        cells = re.findall(r"\$\((.*?)\)\$", m.group(1) + m.group(2)) if m else []
        row = [tuple(x.replace("s_", "s").strip() for x in c.split(r",\,")) for c in cells]
        rep.ok(row == [tuple(p) for p in COORDS[w]], "Table 3.1 row %s" % w)
    # Table 5.2
    cb = tex_float(t, "tab:cert")
    found = {w: (tuple(map(int, v.split(","))), int(L))
             for g, v, L in re.findall(r"\$\\(\w+)\$ & \$\(([0-9,]+)\)\$ & \$(\d+)\$", cb)
             for w in GREEK if tex_name(w) == g}
    rep.ok(found == {w: (tuple(CERT[w][0]), CERT[w][1]) for w in GREEK}, "Table 5.2")
    # Table 4.1
    nt = t[t.index(r"\label{tab:names}"):]
    nt = nt[:nt.index(r"\end{tabular}")]
    texnames = {}
    for name, rhs in re.findall(r"\$\\?([A-Za-z]+(?:_\{?\d+\}?|_\\(?:pi|rho|chi))?)=([s_\d+]+)\$", nt):
        key = name.replace("_", "").replace("{", "").replace("}", "").replace("\\", "")
        texnames[key] = tex_sums(rhs)[0]
    rep.ok(texnames == NAMES, "Table 4.1 (%d names)" % len(texnames))
    # Table 5.1
    bt = t[t.index(r"\label{tab:bounds}"):]
    bt = bt[:bt.index(r"\end{tabular}")]
    for i in range(1, 10):
        m = re.search(r"\$b_%d=\\min\((.*?)\)\$" % i, bt)
        names = [x.strip().replace("_", "").replace("{", "").replace("}", "").replace("\\", "")
                 for x in m.group(1).split(",")] if m else []
        names = ["KK" if x == "KK" else x for x in names]
        rep.ok(names == BOUNDS_BY_NAME[i], "Table 5.1 row b_%d" % i)
    # Propositions 4.1-4.11 and the definition of b_2
    for i in range(1, 10):
        m = re.search(r"b_%d\(s\)(?::=|=)\\min(.*?)(?:\.|\\end|\$)" % i, t, re.S)
        sums = tex_sums(m.group(1)) if m else []
        if i == 1:
            sums = sums + ["KK"]
        rep.ok(sorted(sums) == sorted(BOUNDS[i]), "definition of b_%d" % i)
    # Figures 3.1 and 3.2: labels, drawn packings, and bars
    for label in ("fig:a", "fig:b"):
        for p in tex_panels(tex_float(t, label)):
            m = re.search(r"\$\\(\w+)\$,\\quad \$L_8\^\*=(\d+)\$\\\\ \$s=\(([0-9,]+)\)", p["text"])
            w = next(v for v in GREEK if tex_name(v) == m.group(1))
            s, L = tuple(map(int, m.group(3).split(","))), int(m.group(2))
            rep.ok((s, L) == tuple(ILLUSTRATIONS[w]), "%s: panel %s label" % (label, w))
            rep.ok(as_int_map(p["X"]) == coords_at(COORDS[w], s, L) and p["L"] == L
                   and all(p["sides"][i] == s[i] for i in range(N8)),
                   "%s: panel %s drawn as Table 3.1 at L = L_8^*" % (label, w))
            drawn = sorted("".join(map(str, sorted(v))) for v in p["bars"].values())
            rep.ok(drawn == tight_terms(s, w, L),
                   "%s: panel %s bars are the side sums of %s equal to L_8^*" % (label, w, w))
    # Figures with explicit packings
    figs = [("fig:insert", [FIGURE_1_1[k][1:] for k in "abcde"]),
            ("fig:cover8", list(FIGURE_2_1.values())),
            ("fig:bridge", list(FIGURE_4_1.values())),
            ("fig:stack", list(FIGURE_4_2.values())),
            ("fig:pincer", list(FIGURE_4_3.values())),
            ("fig:frame8", list(FIGURE_4_4.values()))]
    for label, data in figs:
        panels = tex_panels(tex_float(t, label))
        ok = len(panels) == len(data)
        for p, (s, L, X) in zip(panels, data):
            ok = ok and as_int_map(p["X"]) == X and p["L"] == L
            ok = ok and all(p["sides"][i] == s[i] for i in X)
        rep.ok(ok, "%s: drawn packings agree with the data in this script" % label)
    rep.note("tables, bound definitions, and drawn packings agree with the source")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Verify the computations in M. Badger, "An exact formula for '
                    'optimal axis-parallel packing of eight squares".')
    ap.add_argument("--tex", help="LaTeX source of the paper; compare its tables and figures with the data here")
    ap.add_argument("--grid", type=int, default=12,
                    help="check F uses all nonincreasing 8-tuples with entries at most GRID (default 12)")
    ap.add_argument("--samples", type=int, default=20000,
                    help="number of additional random test vectors in check F (default 20000)")
    ap.add_argument("--seed", type=int, default=2026, help="random seed for check F")
    ap.add_argument("-v", "--verbose", action="store_true", help="print every individual check")
    args = ap.parse_args(argv)

    rep = Report(verbose=args.verbose)
    check_formulas(rep)
    strategies = check_constructions(rep)
    check_illustrations(rep)
    check_restrictions(rep, strategies)
    check_minimality(rep)
    check_arithmetic(rep, args.grid, args.samples, args.seed)
    check_other_figures(rep)
    if args.tex:
        check_tex(rep, args.tex)
    print()
    if rep.failures:
        print("FAILED: %d check(s) failed." % rep.failures)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
