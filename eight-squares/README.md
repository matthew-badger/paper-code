# Eight squares

This repository contains code that checks the explicit computations in

> Matthew Badger, *An exact formula for optimal axis-parallel packing of eight squares*.

The paper determines L\*\_8, the least side of a square containing eight axis-parallel squares with prescribed side lengths, as a minimum of twenty-three strategy costs. It also proves that no minimum of fewer convex functions gives L\*\_8. The lower bounds are proved by hand in the paper. The code checks the constructions, the tables, the figures, and the integer certificates on which the proofs rely.

## Running the checks

The script needs Python 3.8 or later and uses only the standard library. All arithmetic is exact integer arithmetic.

```
python3 verify_eight_squares.py
python3 verify_eight_squares.py --tex 8-square.tex   # also compare with the LaTeX source
python3 verify_eight_squares.py -v                    # list every individual check
```

A full run takes about 15 seconds. The script prints a summary for each group of checks and ends with `All checks passed.` The exit status is 0 if every check passes and 1 otherwise.

The options `--grid N` (default 12), `--samples K` (default 20000), and `--seed S` control the test vectors in check F.

## What is checked

Numbers refer to sections, tables, figures, and equations of the paper.

**A. Formulas and summary tables.** Table 1.1 lists twenty-three maxima. The names in Table 4.1 agree with (2.5), and each row of Table 5.1 agrees with the definition of the corresponding bound b\_i in Section 4.

**B. Constructions (Table 3.1, Proposition 3.1).** For each of the twenty-three rows of Table 3.1, the script parses the coordinates as linear forms in s\_0, …, s\_7 and L. It then certifies, for every nonincreasing nonnegative s and every L ≥ max(L\*\_6(s), ω(s)):

- containment of every square;
- a separation, horizontal or vertical, for every pair of squares.

Each certificate is the coefficient (prefix-sum) test (2.4), applied after adding a multiple of one side sum of the cost, as described in Section 2. The certified relations define the strategy σ\_ω. The script also checks that:

- every side sum in ω is a chain of σ\_ω;
- each of the three sums in L\*\_6 is a chain or is dominated by a side sum in ω;
- every chain of σ\_ω is dominated by a single side sum of max(L\*\_6, ω).

Together these give cost(σ\_ω) = max(L\*\_6, ω). The same checks are run for the three seven-square strategies top, center, and right of the companion paper.

**C. Figures 3.1 and 3.2.** At each illustrated side vector, the script checks four things:

- L\*\_8(s) is the stated value and exceeds L\*\_7(s);
- σ\_ω is the only one of the twenty-three strategies with cost L\*\_8(s);
- the coordinates of Table 3.1 at L = L\*\_8(s) form an exact packing;
- every side sum of ω equal to L\*\_8(s) is a chain of that packing.

**D. Restrictions to the seven largest squares (Section 1, equation (3.1)).** The script checks:

- Equation (3.1) agrees with the restriction to Q\_0, …, Q\_6 of the strategy obtained from Table 3.1.
- Exactly ten restrictions have the cost of top, center, or right, as listed after (3.1).
- Each of the other thirteen differs from all three. The script finds an explicit integer side vector for each comparison.
- At the illustrated side vectors, the fourteen restrictions named in Section 1 cost L\*\_7 + 1 = L\*\_8. These are the thirteen and σ\_τ.
- The α example of Section 1 and Figure 1.1(e) holds as stated.

**E. Minimality certificate (Table 5.2).** At each p\_ω, the script checks the stated value of L\*\_8(p\_ω) and that σ\_ω is the only strategy attaining it. It then checks all 253 pairwise inequalities L\*\_8(p\_ω + p\_ω′) > L\*\_8(p\_ω) + L\*\_8(p\_ω′), whose differences range from 1 to 135.

**F. Arithmetic of Section 5.** Lemma 5.1 is checked exactly: each of the five omitted balanced partitions is dominated, by the coefficient test, by the partition named in its proof. In addition, on every nonincreasing 8-tuple with entries at most 12 and on 20000 random nonincreasing tuples, the script checks:

- the two forms of the formula in Theorem 1.1 agree;
- Theorem 5.3 holds;
- each b\_i is at most the formula;
- identity (5.1) of Lemma 5.2 holds.

These finite checks are consistency tests, not proofs. The proofs are in the paper.

**G. Figures 1.1, 2.1, 4.1, 4.2, 4.3, and 4.4.** Every drawn configuration is an exact packing with the properties stated in its caption. These properties include:

- the arm structure and covering chains in Figure 2.1;
- the inequalities between L and the named side sums that the captions quote;
- the chains recorded by the bars.

**H. Comparison with the LaTeX source (optional, `--tex`).** The data in the script are compared with the paper's source:

- Tables 1.1, 3.1, 4.1, 5.1, and 5.2;
- the definitions of b\_1, …, b\_9;
- the label, coordinates, and chain bars of every panel in Figures 3.1 and 3.2;
- the coordinates of every square drawn in Figures 1.1, 2.1, and 4.1–4.4.

This guards against transcription errors between the paper and the code.

## What is not checked

The nine lower bounds of Section 4 are proved by geometric arguments in the paper, and this code does not verify those arguments. Check F tests the arithmetic that combines them (Theorem 5.3) on many side vectors. The proof of Theorem 5.3 itself is the argument in Section 5.1.

## License

MIT License; see [`LICENSE`](../LICENSE) at the top of the repository. The license text is also included at the top of `verify_eight_squares.py`.
