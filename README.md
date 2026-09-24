# FCC + L1₂ Multi-Principal-Element Superalloy Designer

A rule-based (ML-derived) design and screening toolkit for **FCC + L1₂ dual-phase
multi-principal-element superalloys (MPESAs)** — high-entropy analogues of Ni-base
superalloys combining an FCC solid-solution matrix with coherent L1₂ (γ′-type)
nanoprecipitates.

The design rules implemented here convert machine-learning model knowledge into
transparent, human-readable criteria (Tao et al., 2025), so no black-box model is
required: the "model" is fully shareable as rules + code.

## 1. Design rules (ML-derived)

| # | Criterion | Threshold |
|---|-----------|-----------|
| 1 | VEC | > 8 |
| 2 | ΔH_mix | −16.0 < ΔH_mix < −9.7 **kJ/mol** |
| 3 | Average melting point T̄m | 1671 < T̄m < 1822 K |
| 4 | Co content (if present) | 31 – 72 at.% |
| 5 | Fe content (if present) | ≤ 10 at.% |
| 6 | Al content (if present) | ≤ 9 at.% |

> **Unit note:** the source paper prints rule 2 as J·mol⁻¹·K⁻¹; ΔH_mix is in
> **kJ/mol** (consistent with Table 2 of the paper and with this implementation).
> Favorable (not enforced) ranges from the paper: δr ≈ 3–5.4 %, Δχ < 0.12,
> ΔS_mix ≈ 5.3–13.4 J/(mol·K).

## 2. Element database

15 elements: **9 refractory** (Nb Mo Ta W V Cr Ti Zr Hf) + **6 FCC-family**
(Ni Co Fe Cu Mn Al). Per element: atomic radius, VEC, e/a, Pauling χ, Tm,
bulk modulus B, thermal expansion α, and DFT total energies of BCC/FCC/HCP
(Wang et al. 2004 lattice-stability set).

Binary mixing enthalpies ΔH_ij (kJ/mol), two sources:
- Refractory↔refractory pairs: Takeuchi & Inoue (2005) semi-empirical table.
- Pairs involving Al and the 3d FCC family: Miedema/de Boer values as used in
  the HEA literature (e.g. Ni-Al = −22, Co-Al = −19, Fe-Al = −11, Cu-Al = −1).
  These differ slightly from the Takeuchi–Inoue table; both are standard sources.

## 3. The five programs

| File | Role | Method | Input | Output |
|---|---|---|---|---|
| `fcc_l12_ml_designer.py` | v1 scanner (basic ML rules) | exhaustive grid, 5% step | element set | valid compositions, parameter ranges, best candidate (lowest E_FCC) |
| `fcc_l12_ml_designer_v2.py` | v2 scanner (strict) | exhaustive grid, 2% step | element set | valid compositions + detrimental-phase filters + safe concentration ranges |
| `fcc_l12_composition_evaluator.py` | single-alloy verifier | rule evaluation | user-defined composition | pass/fail per rule + final verdict |
| `hea_fcc_calculator.py` | general SS + FCC evaluator | empirical criteria | 3–6 elements + fractions | MC1–MC7 solid-solution verdicts + FCC structural criteria |
| `hea_fcc_range_dft_scanner.py` | SS + FCC range scanner | exhaustive grid + DFT | element set | valid ranges incl. DFT BCC/FCC/HCP energies |

Workflow analogy to the source paper: the scanners implement the *design strategy*
(rules → generate candidates), while the evaluators play the role of *verification*
(check one composition against all rules).

## 4. Installation

Python 3.8+. No hard dependencies (`tabulate` is optional — all scripts fall back
to plain-text tables if it is not installed):

```bash
pip install tabulate
```

## 5. Usage

```bash
python fcc_l12_ml_designer.py            # v1 scan, 5% step
python fcc_l12_ml_designer_v2.py         # strict scan, 2% step + phase filters
python fcc_l12_composition_evaluator.py  # check ONE composition
python hea_fcc_calculator.py             # general SS + FCC check
python hea_fcc_range_dft_scanner.py      # SS+FCC ranges with DFT energies
```

Example: `fcc_l12_ml_designer_v2.py` with `Ni Co Cr Al Ti` scans ~2×10⁵ compositions
(few seconds) and reports the valid FCC+L1₂ window and the candidate with the lowest
DFT E_FCC.

## 6. Scientific caveats

1. The DFT energies are per-element reference data (Wang et al. 2004); the
   composition-weighted alloy average mixes elements with different references —
   treat E_mix values as *relative indicators*, not absolute thermodynamics.
2. The rules are correlations extracted from a Co/Ni/Fe/Cr/Al/Ti-rich MPESA
   dataset; validity outside that chemical space is unverified.
3. Currently enforced element filters: Co 31–72 %, Fe ≤ 10 %, Al ≤ 9 %.
   The paper additionally recommends Cr < 15 at.%, Al ≥ 1.5 at.% and Ti 1–9 at.%
   (planned update). Selecting element sets without Al/Co may pass all rules while
   L1₂ formation is actually unlikely.
4. These rules predict *phase presence*, not mechanical performance.
5. Grid scanners evaluate discrete steps (5% / 2%); narrow windows can be missed.

## 7. References

1. Q. Tao, X. Yang, L. Bao, Y. Zhou, T. Yang, Y. Zhao, R. Shi, Z. Yao, X. Liu,
   *Transforming machine learning model knowledge into material insights for
   multi-principal-element superalloy phase design*, npj Computational Materials
   11, 99 (2025). https://doi.org/10.1038/s41524-025-01578-6
2. Y. Wang et al., *Ab initio lattice stability in comparison with CALPHAD lattice
   stability*, CALPHAD 28, 79–90 (2004).
3. A. Takeuchi, A. Inoue, *Classification of bulk metallic glasses...*, Materials
   Transactions 46, 2817–2829 (2005).
4. A.R. Miedema / F.R. de Boer et al., cohesive-energy data as tabulated in the
   HEA literature.
5. Y. Zhang et al., Adv. Eng. Mater. 10, 534 (2008); X. Yang & Y. Zhang, Mater.
   Chem. Phys. 132, 233 (2012); S. Guo et al., Intermetallics 41, 96 (2013) and
   J. Appl. Phys. 109, 103505 (2011); M.G. Poletti & L. Battezzati, Acta Mater.
   75, 297 (2014); Z. Wang et al., Scr. Mater. 94, 28 (2014).

## License

MIT (see LICENSE)