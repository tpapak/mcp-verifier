# mcp-verifier

R scripts for mathematical verification of network meta-analysis (NMA) results
and reproducible R code generation.

Used as a submodule by [tpapak/mcp-netmeta](https://github.com/tpapak/mcp-netmeta)
and served via the **netmeta-verify MCP server** at:

```
https://biostatistics.med.auth.gr/mcp/netmeta-verify
```

---

## Scripts

| Script | Purpose |
|--------|---------|
| `verify_fixed_effect_mle.R` | Verify a fixed-effect NMA league table is the unique MLE |
| `verify_random_effect_reml.R` | Verify a random-effects NMA satisfies REML optimality conditions |
| `generate_r_code.R` | Generate a standalone reproducible R script for any NMA |

---

## Installation

```r
# Required R packages
install.packages(c("netmeta", "igraph", "jsonlite"))
remotes::install_github("tpapak/multiarmvars")

# Install netmeta at the pinned version used by the MCP server
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
  TE       = Senn2013$TE,
  seTE     = Senn2013$seTE,
  treat1   = Senn2013$treat1,
  treat2   = Senn2013$treat2,
  studlab  = Senn2013$studlab,
  sm       = "MD",
  common   = TRUE,
  random   = FALSE,
  reference.group = "plac"
)

v <- verify_netmeta_mle(nma)
print_verification(v)
```

Example output:

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

### Verify any software output (not just netmeta)

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
  TE           = c(0.5, 0.8, 0.3),
  seTE         = c(0.1, 0.15, 0.12),
  treat1       = c("A", "A", "B"),
  treat2       = c("B", "C", "C"),
  studlab      = c("S1", "S2", "S3"),
  league_table = league_table,
  reference    = "C"
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

data(Senn2013)
nma <- netmeta(...)   # run your analysis first

code <- generate_nma_code_from_netmeta(nma)
cat(code)             # print to console
save_nma_code(code, "my_analysis.R")   # save to file
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

### Fixed-effect MLE conditions

For a fixed-effect NMA with normal contrast-level data the MLE must satisfy:

1. **Score = 0** — X'W(y − Xθ̂) = 0
2. **Positive definite Hessian** — X'WX ≻ 0
3. **Full rank** — rank(X'WX) = number of parameters
4. **WLS equivalence** — θ̂ = (X'WX)⁻¹ X'Wy
5. **Multi-arm consistency** — contrast variances decompose into non-negative arm variances

### Random-effect REML conditions

Seven conditions are checked, including the REML score for τ² and a
perturbation test confirming the log-likelihood is at a maximum.

Full theoretical background: [`docs/verify.md`](https://github.com/tpapak/mcp-netmeta/blob/paper/docs/verify.md)

---

## Multi-arm studies

Multi-arm studies have correlated comparisons. The verifier uses the
[multiarmvars](https://github.com/tpapak/multiarmvars) package to recover
arm-level variances from contrast variances and construct the correct
off-diagonal covariance elements.

**Note on netmeta's approximation:** netmeta inflates variances to achieve
approximate independence between contrasts. This is valid but gives slightly
different estimates from the exact MLE. Verification of netmeta output with
multi-arm data will therefore fail the exact MLE test; the discrepancy is
typically small (< 0.05).

---

## MCP server

These scripts are the backend for the **netmeta-verify MCP server**, which
exposes them as tools callable from any AI assistant:

- `verify_fixed_effect_mle`
- `verify_random_effect_reml`
- `generate_nma_code`

See [biostatistics.med.auth.gr/mcp/](https://biostatistics.med.auth.gr/mcp/)
for setup instructions.

---

## License

GNU Lesser General Public License v3.0 — see [LICENSE](LICENSE).

The LGPL allows these scripts to be used as a library by other software
(including proprietary software) without requiring that software to be
open-source, provided that modifications to these scripts themselves are
shared under the same licence.

## Authors

Thodoris Papakonstantinou
Biostatistics, Aristotle University of Thessaloniki

## References

- Rücker G, et al. (2024). netmeta: Network Meta-Analysis using Frequentist Methods. R package.
- Papakonstantinou T, et al. multiarmvars: Arm variance decomposition for multi-arm studies.
- Lu G, Ades AE (2004). Combination of direct and indirect evidence in mixed treatment comparisons. *Statistics in Medicine*, 23(20), 3105–3124.
