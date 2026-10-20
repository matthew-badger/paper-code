# Eight squares

This folder contains code that checks the explicit computations in

> Matthew Badger, *An exact formula for optimal axis-parallel packing of eight squares*.

The paper determines L\*\_8, the least side of a square containing eight axis-parallel squares with prescribed side lengths, as a minimum of the costs of twenty-three packing strategies. It also proves that no minimum of fewer convex functions gives L\*\_8. The lower bounds are proved by hand in the paper. The code checks the constructions, the tables, the figures, and the integer certificates on which the proofs rely.

## Version

The current version of `verify_eight_squares.py` is **1.1.0** (October 20, 2026). Run `python3 verify_eight_squares.py --version` to print it. Every run also prints the version on its first line.

- **1.1.0** (October 20, 2026). Matches the current manuscript: eight lower bounds b\_1, …, b\_8 together with the column bound of the seven-square paper (Lemma 2.4), Figure 3.1 set over two pages, and the current numbering of sections, tables, figures, and equations. The comparison with the LaTeX source now checks the general statement of Lemma 2.4, every side-sum definition in the text, and both dimensions of every drawn square. Adds the version number.
- **1.0.0** (October 16, 2026). First release.

## Files

- `verify_eight_squares.py`: the verification script.
- `8-square.tex`: the LaTeX source of the paper (version of October 20, 2026), for use with `--tex`.

## Running the checks

The script needs Python 3.8 or later and uses only the standard library. All arithmetic is exact integer arithmetic.

```
python3 verify_eight_squares.py
python3 verify_eight_squares.py --tex 8-square.tex   # also compare with the LaTeX source
python3 verify_eight_squares.py -v                    # list every individual check
python3 verify_eight_squares.py --version
```

A full run takes about 10 to 15 seconds. The script prints a summary for each group of checks and ends with `All checks passed.` The exit status is 0 if every check passes and 1 otherwise.

The options `--grid N` (default 12), `--samples K` (default 20000), and `--seed S` control the test vectors in check F.

## What is checked

Numbers refer to sections, tables, figures, and equations of the paper.

**A. Formulas and summary tables.** Table 1.1 lists twenty-three maxima. The names in Table 4.1 agree with (1.2) and (5.1). Each row of Table 5.1 agrees with the definition of the corresponding bound b\_i in Section 4. The column bound of Lemma 2.4 at n = 8 is min(D, P) in the names of Table 4.1.

**B. Constructions (Table 3.1, Proposition 3.1).** For each of the twenty-three rows of Table 3.1, the script parses the coordinates as linear forms in s\_0, …, s\_7 and L. It then certifies, for every nonincreasing nonnegative s and every L ≥ max(L\*\_6(s), ω(s)):

- containment of every square;
- a separation, horizontal or vertical, for every pair of squares.

Each certificate is the coefficient test (3.1), applied after adding a multiple of one side sum of the cost, as in the proof of Proposition 3.1. The certified relations define the strategy σ\_ω. The script also checks that:

- every side sum in ω is a chain of σ\_ω;
- each of the three sums in L\*\_6 is a chain or is dominated by a side sum in ω;
- every chain of σ\_ω is dominated by a single side sum of max(L\*\_6, ω).

Together these give cost(σ\_ω) = max(L\*\_6, ω). The same checks are run for the three seven-square strategies top, center, and right of the companion paper, with the costs in (1.2).

**C. Figure 3.1.** At each of the twenty-three illustrated side vectors, the script checks four things:

- L\*\_8(s) is the stated value and exceeds L\*\_7(s);
- σ\_ω is the only one of the twenty-three strategies with cost L\*\_8(s);
- the coordinates of Table 3.1 at L = L\*\_8(s) form an exact packing;
- every side sum of ω equal to L\*\_8(s) is a chain of that packing.

**D. Restrictions to the seven largest squares (Sections 1 and 3, equation (3.2)).** The script checks:

- Equation (3.2) agrees with the restriction to Q\_0, …, Q\_6 of the strategy obtained from Table 3.1.
- Exactly ten restrictions have the cost of top, center, or right, as listed after (3.2).
- Each of the other thirteen differs from all three. The script finds an explicit integer side vector for each comparison.
- At the illustrated side vectors, the fourteen restrictions named in Section 3 cost L\*\_7 + 1 = L\*\_8. These are the thirteen and σ\_τ.
- The α example of Section 1 and Figure 1.1(e) holds as stated, and deleting Q\_7 from σ\_β leaves the arrangement center.

**E. Minimality certificate (Table 5.2).** At each p\_ω, the script checks the stated value of L\*\_8(p\_ω) and that σ\_ω is the only strategy attaining it. It then checks all 253 pairwise inequalities L\*\_8(p\_ω + p\_ω′) > L\*\_8(p\_ω) + L\*\_8(p\_ω′) of (5.5), whose differences range from 1 to 135.

**F. Arithmetic of Section 5.** Lemma 5.1 is checked exactly: each of the five omitted balanced partitions is dominated, by the coefficient test, by the partition named in its proof. In addition, on every nonincreasing 8-tuple with entries at most 12 and on 20000 random nonincreasing tuples, the script checks:

- the two forms of the formula in Theorem 1.1 agree;
- Theorem 5.3 holds: max(L\*\_7, min(D, P), b\_1, …, b\_8) ≥ F\_8;
- each b\_i is at most the formula;
- identity (5.2) in the proof of Lemma 5.2 holds.

These finite checks are consistency tests, not proofs. The proofs are in the paper.

**G. Figures 1.1, 2.1, 4.1, 4.2, 4.3, and 4.4.** Every drawn configuration is an exact packing with the properties stated in its caption. These properties include:

- the arm structure and covering chains in Figure 2.1;
- the inequalities between L and the named side sums that the captions quote;
- the chains recorded by the bars.

**H. Comparison with the LaTeX source (optional, `--tex`).** The data in the script are compared with the paper's source:

- Tables 1.1, 3.1, 4.1, 5.1, and 5.2;
- the general statement of Lemma 2.4 and its instance for n = 8 in Section 5.1;
- the definitions of b\_1, …, b\_8;
- every equation in the text that defines a recurring side sum;
- the label, coordinates, and chain bars of every panel in Figure 3.1, and the order of the panels;
- the coordinates and both dimensions of every square drawn in Figures 1.1, 2.1, and 4.1–4.4.

This guards against transcription errors between the paper and the code.

## What is not checked

The eight lower bounds of Section 4 are proved by geometric arguments in the paper, and this code does not verify those arguments. The column bound (Lemma 2.4) is cited from the seven-square paper. Check F tests the arithmetic that combines these bounds (Theorem 5.3) on many side vectors. The proof of Theorem 5.3 itself is the argument in Section 5.1.

## License

MIT License; see [`LICENSE`](../LICENSE) at the top of the repository. The license text is also included at the top of `verify_eight_squares.py`.
