# mcp-verifier

A collection of R scripts for mathematically verifying statistical results
and generating reproducible analysis code.

Currently focused on **network meta-analysis (NMA)** verification, with the
intention of growing into a broader library of verification functions for
evidence synthesis methods.

Part of the [MCP tools for evidence synthesis](https://biostatistics.med.auth.gr/mcp/)
project at Biostatistics, Aristotle University of Thessaloniki.

---

## Contents

| Script | Purpose |
|--------|---------|
| `verify_fixed_effect_mle.R` | Verify a fixed-effect NMA league table is the unique MLE |
| `verify_random_effect_reml.R` | Verify a random-effects NMA satisfies REML optimality conditions |
| `generate_r_code.R` | Generate a standalone reproducible R script for any NMA |

---

## Installation

```r
install.packages(c("netmeta", "igraph", "jsonlite"))
remotes::install_github("tpapak/multiarmvars")

# Pinned netmeta version
remotes::install_github("guido-s/netmeta@5ecfc1d7739c3df360a694d60af0563bc43d68ea")
```

---

## Usage

### Fixed-effect MLE verification

```r
source("verify_fixed_effect_mle.R")

library(netmeta)
data(Senn2013)

nma <- netmeta(
  TE = Senn2013$TE, seTE = Senn2013$seTE,
  treat1 = Senn2013$treat1, treat2 = Senn2013$treat2,
  studlab = Senn2013$studlab,
  sm = "MD", common = TRUE, random = FALSE,
  reference.group = "plac"
)

v <- verify_netmeta_mle(nma)
print_verification(v)
```

Output:

```
=== Fixed-Effect NMA MLE Verification ===

Observations: 28  |  Studies: 26 (1 multi-arm)  |  Parameters: 9
Reference: plac   |  Tolerance: 1.00e-06

[PASS] Score equations = 0
[PASS] Information matrix positive definite
[PASS] Unique solution
[PASS] Matches WLS solution
[PASS] Multi-arm variance consistency

Result: 5/5 tests passed
CONCLUSION: Output is the UNIQUE maximum likelihood estimate.
```

### Verify output from any NMA software

```r
source("verify_fixed_effect_mle.R")

league_table <- matrix(
  c( 0, -0.5, -0.8,
     0.5,  0, -0.3,
     0.8,  0.3,  0),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("A","B","C"), c("A","B","C"))
)

v <- verify_fixed_effect_mle(
  TE = c(0.5, 0.8, 0.3), seTE = c(0.1, 0.15, 0.12),
  treat1 = c("A","A","B"), treat2 = c("B","C","C"),
  studlab = c("S1","S2","S3"),
  league_table = league_table, reference = "C"
)
print_verification(v)
```

### Random-effect REML verification

```r
source("verify_random_effect_reml.R")

nma <- netmeta(
  TE = Senn2013$TE, seTE = Senn2013$seTE,
  treat1 = Senn2013$treat1, treat2 = Senn2013$treat2,
  studlab = Senn2013$studlab,
  sm = "MD", common = FALSE, random = TRUE,
  method.tau = "REML", reference.group = "plac"
)

v <- verify_netmeta_reml(nma)
print_reml_verification(v)
```

### Generate reproducible R code

```r
source("generate_r_code.R")

nma <- netmeta(...)

code <- generate_nma_code_from_netmeta(nma)
cat(code)
save_nma_code(code, "my_analysis.R")
```

---

## Running tests

```bash
Rscript tests/test_verify_fixed_effect_mle.R
Rscript tests/test_verify_fixed_effect_multiarm.R
Rscript tests/test_verify_random_effect_reml.R
Rscript tests/test_generate_r_code.R
```

---

## Theory

### Fixed-effect MLE

The MLE for a fixed-effect NMA must satisfy:

1. **Score = 0** &mdash; X'W(y &minus; X&theta;&#x0302;) = 0
2. **Positive definite Hessian** &mdash; X'WX &succ; 0
3. **Full rank** &mdash; unique solution
4. **WLS equivalence** &mdash; &theta;&#x0302; = (X'WX)&minus;&sup1; X'Wy
5. **Multi-arm consistency** &mdash; contrast variances decompose into non-negative arm variances

### Random-effect REML

Seven conditions are checked: fixed-effects score, REML score for &tau;&sup2;,
positive-definite information matrix, &tau;&sup2; &ge; 0, uniqueness, GLS
equivalence, and a perturbation test confirming the log-likelihood is at
a maximum.

### Multi-arm studies

Multi-arm studies have correlated comparisons. The
[multiarmvars](https://github.com/tpapak/multiarmvars) package is used to
recover arm-level variances from contrast variances and construct the correct
off-diagonal covariance elements.

**Note:** netmeta uses an approximate variance inflation method for multi-arm
studies. Verification of netmeta output will therefore show small discrepancies
from the exact MLE (typically &lt; 0.05); this is expected and documented.

---

## MCP server

The scripts in this repository are the backend for the
**netmeta-verify MCP server**, which exposes them as tools callable from
any AI assistant:

| Tool | Script |
|------|--------|
| `verify_fixed_effect_mle` | `verify_fixed_effect_mle.R` |
| `verify_random_effect_reml` | `verify_random_effect_reml.R` |
| `generate_nma_code` | `generate_r_code.R` |

Endpoint: `https://biostatistics.med.auth.gr/mcp/netmeta-verify`

Setup instructions: [biostatistics.med.auth.gr/mcp/](https://biostatistics.med.auth.gr/mcp/)

---

## License

GNU Lesser General Public License v3.0 &mdash; see [LICENSE](LICENSE).

LGPL allows these scripts to be used as a library by other software without
requiring that software to be open-source, while keeping modifications to
these scripts themselves open.

## Authors

Thodoris Papakonstantinou  
Biostatistics, Aristotle University of Thessaloniki

## References

- R&uuml;cker G, et al. (2024). netmeta: Network Meta-Analysis using Frequentist Methods. R package.
- Papakonstantinou T. multiarmvars: Arm variance decomposition for multi-arm studies.
- Lu G, Ades AE (2004). Combination of direct and indirect evidence in mixed treatment comparisons. *Statistics in Medicine*, 23(20), 3105&ndash;3124.
