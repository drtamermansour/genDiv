# MIGRATION.md — per-group popRefs rollout

This document is the contract between the upstream `genDiv` pipeline and the downstream GPA consumer for requested refactor. It covers:

1. The filename convention upstream now guarantees.
2. Schema invariants GPA reads by fixed column index.
3. Numerical equivalence vs genuine semantic change per file.
4. A checklist for the GPA-side `create_popFiles.sh` update.
5. How to use the shared `validate_popRefs.sh` on both sides.

Keep this document authoritative — if upstream or GPA diverges from it, fix whichever side drifted rather than the doc.

---

## 1. Filename convention

Every per-group reference file follows `<stem>.${rg}.<ext>` with `rg ∈ {wholePop, Trotter, Pacer}`. No file that GPA consumes uses a bare whole-pop name anymore; `wholePop` is spelled out in the filename just like `Trotter` and `Pacer`.

Upstream source paths (after Step 7) and the flat popFiles/ targets:

| Logical | Upstream path | popFiles target | Produced by |
|---|---|---|---|
| Per-group afreq | `${OUTPUT_DIR}/LD_pruned/pruned.${rg}.afreq` | `popFiles/pruned.${rg}.afreq` | `genDiversity_per_group.sh` §1 |
| Per-group tabix AF | `${OUTPUT_DIR}/divStats/freqs.${rg}.tab.gz` (+ `.tbi`) | `popFiles/freqs.${rg}.tab.gz` (+ `.tbi`) | `genDiversity_per_group.sh` §6 |
| Per-group het (F_SNP) | `${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.${rg}.het` | `popFiles/filtered.LD_prune.het_stats.${rg}.het` | `genDiversity_per_group.sh` §4 |
| Per-group F_ROH summary | `${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.${rg}.txt` | `popFiles/roh_summary_by_RG_L3_Froh.${rg}.txt` | `genDiversity_per_group.sh` §8 |
| Per-group Inbreeding_Comparison | `${OUTPUT_DIR}/rep_ROHRM/roh_1Mb.Threshold_3SD/Inbreeding_Comparison.${rg}.csv` | `popFiles/Inbreeding_Comparison.${rg}.csv` | `genDiversity_per_group.sh` §9 |
| Per-group Pairwise_Differences | `${OUTPUT_DIR}/rep_ROHRM/roh_1Mb.Threshold_3SD/Pairwise_Differences.${rg}.csv` | `popFiles/Pairwise_Differences.${rg}.csv` | `genDiversity_per_group.sh` §9 |
| Per-group consensus ROH BED | `${OUTPUT_DIR}/divStats/roh.L3.consensus_25pct.merged.${rg}.smoothed.bed` | `popFiles/roh.L3.consensus_25pct.merged.${rg}.smoothed.bed` | `genDiversity_per_group.sh` §7 |
| Per-group ROH_sh intersect | `${OUTPUT_DIR}/divStats/roh.L3.perSample_intersect_${rg}_consensus_25pct.summary.txt` | `popFiles/roh.L3.perSample_intersect_${rg}_consensus_25pct.summary.txt` | `genDiversity_per_group.sh` §7 |
| Sample → group mapping | `${OUTPUT_DIR}/preprocess/sample_groups.tsv` | `popFiles/sample_groups.tsv` | `genDiversity_shared.sh` |

### Scope note on `freqs.${rg}.tab.gz`

The `freqs.${rg}.tab.gz` files **are** produced. Upstream's own `bcftools roh` run uses `--estimate-AF -` on the group-subset VCF (no external AF file needed in-pipeline), but GPA's downstream `bcftools roh` runs on a 1–2 animal mate-pair VCF that cannot estimate AF from its own samples — it requires an external tabix-indexed AF table via `--AF-file`. That's what `freqs.${rg}.tab.gz` + `.tbi` is for: downstream consumption, one per group, so each GPA report tab can call ROH on the candidate animal under the same AF prior that the tab's reference distribution was built with. GPA should **stop** building `freqs.tab.gz` locally in its `create_popFiles.sh` and instead copy the three upstream-produced per-group files.

### 1.1 Secondary outputs (PCA artifacts)

The nine files above are the per-group reference contract. The PLINK2 `--pca` outputs are a second family that upstream produces but does **not** promise as a fixed-shape interface — PC count, eigenvalue count, and the SNP-row count of `eigenvec.allele` all drift with the input data. Treat the *names* as a commitment and the *shapes* as the consumer's concern (see §5).

**GPA consumes exactly one of them: `eigenvec.allele`.** It is the PLINK2 `--score` weight file for the gait classifier's PCA projection, used on both sides of that classifier. `create_popFiles.sh` reads it from `divStats/filtered.LD_prune.wholePop.pca.eigenvec.allele` to project the reference cohort into `popFiles/filtered.LD_prune.projected.sscore` (which `gait_model.py` trains on to produce `popFiles/gait_model.pkl`), then copies it to the flat, rg-less `popFiles/eigenvec.allele`. `GPA.sh` and the `Snakefile` read that flat copy as `pca_allele_wt` to project the *candidate* into the same coordinate frame. The older path — training directly off `.eigenvec.wGait` — is commented out at `create_popFiles.sh:94`, so the overlay files are no longer on GPA's read path either. Nothing downstream reads `eigenvec`, `eigenval`, the `wBook_Size` / `wCOI` overlays, or `pca_pairwise_euclidean.dist`, and nothing reads the Trotter / Pacer PCA family at all.

**Rename status: already applied downstream.** The GPA-side read was updated to the infixed name on GPA's own `per-group-references` branch (`create_popFiles.sh`, commit `bd14464`). The tables below are therefore a record of what moved, not an outstanding work item — see §4 step 4 for the verification that remains.

The refactor inserted the same `.${rg}.` infix into the PCA outputs that previously carried bare whole-pop names. Source: `genDiversity_per_group.sh` §2 builds `pca_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune${rg_tag}.pca"` with `rg_tag=".${rg}"`.

**Renamed (wholePop) — 9 files.** Only `eigenvec.allele` is consumed downstream; the other eight moved for naming-convention consistency.

| pre-refactor (main) | per-group-references (current) |
|---|---|
| `divStats/filtered.LD_prune.pca.eigenval` | `divStats/filtered.LD_prune.wholePop.pca.eigenval` |
| `divStats/filtered.LD_prune.pca.eigenvec` | `divStats/filtered.LD_prune.wholePop.pca.eigenvec` |
| `divStats/filtered.LD_prune.pca.eigenvec.allele` | `divStats/filtered.LD_prune.wholePop.pca.eigenvec.allele` |
| `divStats/filtered.LD_prune.pca.eigenvec.wBook_Size` | `divStats/filtered.LD_prune.wholePop.pca.eigenvec.wBook_Size` |
| `divStats/filtered.LD_prune.pca.eigenvec.wCOI` | `divStats/filtered.LD_prune.wholePop.pca.eigenvec.wCOI` |
| `divStats/filtered.LD_prune.pca.eigenvec.wGait` | `divStats/filtered.LD_prune.wholePop.pca.eigenvec.wGait` |
| `divStats/filtered.LD_prune.pca.eigenvec.wSex` | `divStats/filtered.LD_prune.wholePop.pca.eigenvec.wSex` |
| `divStats/filtered.LD_prune.pca.log` | `divStats/filtered.LD_prune.wholePop.pca.log` |
| `divStats/filtered.LD_prune.pca.pca_pairwise_euclidean.dist` | `divStats/filtered.LD_prune.wholePop.pca.pca_pairwise_euclidean.dist` |

**New per-group PCA artifacts — Trotter / Pacer.** The refactor also added per-gait PCA runs. Trotter and Pacer use plain `--pca` (no `'allele-wts'`), so they produce a strict subset of wholePop's PCA file family — no `eigenvec.allele`, no `eigenvec.wGait` (single-gait by construction), no `eigenvec.wSex` (the sex overlay is wholePop-only).

```
divStats/filtered.LD_prune.${rg}.pca.eigenval                 # rg ∈ {Trotter, Pacer}
divStats/filtered.LD_prune.${rg}.pca.eigenvec
divStats/filtered.LD_prune.${rg}.pca.eigenvec.wBook_Size
divStats/filtered.LD_prune.${rg}.pca.eigenvec.wCOI
divStats/filtered.LD_prune.${rg}.pca.log
divStats/filtered.LD_prune.${rg}.pca.pca_pairwise_euclidean.dist
```

These reflect per-group population structure: Trotter PC1 is the strongest axis of variation *within Trotters*, not the Trotter-vs-Pacer axis that wholePop's PC1 captures.

**Dropped — 2 files no longer produced.** Neither is read by GPA (verified by `git grep` over the GPA repo), so the replacements below are recorded for completeness rather than as a required migration.

| file | replacement |
|---|---|
| `divStats/filtered.LD_prune.pca.correlation_plot_PCA_EUCLIDEAN_DIST_vs_KINSHIP_PLINK.png` | `divStats/${rg}.relatedness_correlation.{correlation_heatmap,pairplot}.png` (per-group, `genDiversity_per_group.sh` §14). |
| `divStats/filtered.LD_prune.pca.pca_pairwise_euclidean.dist.withKIN0` | The `.withKIN0` augmentation step is gone; §14's cross-method correlation reads `divStats/filtered.LD_prune.king_gait.${rg}.kin0.withIBS` and `…pca.pca_pairwise_euclidean.dist` separately and merges at consumption time. |

---

## 2. Schema invariants

GPA reads by fixed column index (0-based in comments below; `awk -F` uses 1-based). Upstream guarantees these schemas; any upstream change that reorders columns is a breaking change and requires a coordinated GPA update.

### `pruned.${rg}.afreq` (PLINK2 `--freq`)
Tab-separated. Header row starts with `#CHROM`.
```
#CHROM   ID   REF   ALT   PROVISIONAL_REF?   ALT_FREQS   OBS_CT
  0      1     2     3           4                5         6
```
GPA reads `ALT_FREQS` at 0-based column 5.

### `freqs.${rg}.tab.gz` + `.tbi` (bcftools-compatible tabix AF)
bgzip-compressed, tabix-indexed (`tabix -s1 -b2 -e2`). No header line. Four tab-separated columns:
```
CHROM   POS   REF,ALT   AF
  0      1       2      3
```
Built via `bcftools +fill-tags $group_vcf -- -t AF | bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n'`. Consumed by `bcftools roh --AF-file` only — not by PLINK or R scripts. Schema is identical to GPA's legacy self-built `freqs.tab.gz` (the whole-pop file it produced from `popVCF` via the same `+fill-tags`/`query` pipeline).

### `filtered.LD_prune.het_stats.${rg}.het` (PLINK2 `--het` + `--read-freq`)
Tab-separated. Header starts with `#FID` (PLINK2 output under `--het cols=fid,hom,het,nobs,f`, 8 columns).
```
#FID   IID   O(HOM)   E(HOM)   O(HET)   E(HET)   OBS_CT   F
 0     1      2         3         4        5        6      7
```
GPA reads `F` (F_SNP) at 0-based column 7.

### `roh_summary_by_RG_L3_Froh.${rg}.txt`
Tab-separated, header `IID  NSEG  KB  KBAVG  F_ROH`.
```
IID   NSEG   KB   KBAVG   F_ROH
 0     1      2     3        4
```
GPA reads `F_ROH` at 0-based column 4.

### `Inbreeding_Comparison.${rg}.csv`
Comma-separated. Header `IID,D_STD,D_ROH,Phenotype`.
```
IID, D_STD, D_ROH, Phenotype
 0     1      2       3
```
GPA reads `D_STD` (= D_SNP) at column 1 and `D_ROH` at column 2 (0-based).

### `Pairwise_Differences.${rg}.csv`
Comma-separated. Header `ID1,ID2,Pheno1,Pheno2,Kinship_Std,Kinship_ROH,Difference` (the pipeline also appends `centered_Kinship_diff` as column 7).
```
ID1, ID2, Pheno1, Pheno2, Kinship_Std, Kinship_ROH, Difference, centered_Kinship_diff
 0    1     2       3          4             5           6               7
```
GPA reads `Kinship_Std` (= G_SNP) at column 4 and `Kinship_ROH` (= G_ROH) at column 5 (0-based). The proposal's existing `awk '{print $5}'` / `'{print $6}'` (1-based) covers columns 4/5 here.

### `roh.L3.consensus_25pct.merged.${rg}.smoothed.bed`
Tab-separated BED, no header.
```
chr   start   end   mean_coverage   size_Mb
 0      1      2          3              4
```

### `roh.L3.perSample_intersect_${rg}_consensus_25pct.summary.txt`
Tab-separated. Header `IID  Total_ROH_in_Consensus_region(bp)  Percent_of_Consensus_ROH`.
```
IID   Total_ROH_in_Consensus_region(bp)   Percent_of_Consensus_ROH
 0                      1                                  2
```
GPA reads `Percent_of_Consensus_ROH` (= ROH_sh) at column 2.

### `sample_groups.tsv`
Tab-separated. Header `IID  group`. `group ∈ {Trotter, Pacer, wholePop}`. Samples without a gait label are recorded as `group=wholePop` explicitly. wholePop membership is implicit for every sample regardless of the `group` column value.

### PCA artifacts (§1.1 — secondary, not validator-enforced)

The schemas below are described so a consumer can read them by column index, but unlike the files above they are **not** guaranteed stable across runs: the PC count, the eigenvalue count, and the `eigenvec.allele` row count follow the input data. Header layout is stable; row/column counts are not.

Of these, only `eigenvec.allele` is on GPA's read path today (§1.1). The rest are documented because they are the files a future consumer is most likely to reach for, and because their shapes explain why the family is excluded from the validator (§5).

#### `filtered.LD_prune.${rg}.pca.eigenvec` (PLINK2 `--pca`)
Tab-separated. Header starts with `#FID`. The pipeline does not override PLINK2's default of 10 PCs, so 12 columns total.
```
#FID   IID   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8   PC9   PC10
 0     1      2     3     4     5     6     7     8     9    10    11
```
GPA reads PCs by the 0-based indices above (PC1 at 2 … PC10 at 11). Row count is `NGRP + 1`.

#### `filtered.LD_prune.${rg}.pca.eigenval`
No header, one eigenvalue per row, 10 rows by default. GPA reads the column-0 floats and divides by their sum to derive variance-explained percentages (`var_explained` in `pca_plots.R:29`). No row-count assertion — PLINK2's PC count is the only source of truth.

#### `filtered.LD_prune.wholePop.pca.eigenvec.allele` (wholePop only — `--pca 'allele-wts'`)
Tab-separated. One row per LD-pruned SNP plus header.
```
#CHROM   ID   REF   ALT   A1   PC1   ...   PC10
 0       1     2     3     4     5    ...   14
```
Trotter / Pacer do **not** produce this file — they use plain `--pca`.

#### `filtered.LD_prune.${rg}.pca.eigenvec.{wBook_Size,wCOI}`
Tab-separated. Identical to `eigenvec` plus one trailing column carrying the overlay attribute (`Book_Size` or `COI`), appended by inline `awk` in `genDiversity_per_group.sh` §2 / §5. Consumed only by R plot scripts, which re-derive the overlay column name from the header.

#### `filtered.LD_prune.${rg}.pca.pca_pairwise_euclidean.dist`
Tab-separated. Header `FID1 IID1 FID2 IID2 PCA_EUCLIDEAN_DIST DIST_KINSHIP`. Row count is `(NGRP × (NGRP − 1)) / 2 + 1`. Built inline in `genDiversity_per_group.sh` §11.

---

## 3. Numerical equivalence vs genuine change

This table exists to avoid surprise when comparing the new per-group files to the original pipeline's single whole-pop output.

### ⚠️ Important baseline caveat — the SNP set has drifted

The claims in this section are about **pipeline logic equivalence on a fixed input**, *not* byte-equality against whichever older `popFiles/` bundle GPA last imported. The upstream filtered/phased VCF and the LD-pruned SNP set can change between pipeline runs because (a) the source `SNPdata_iScan_Standardbred/` on the shared Google Drive may be refreshed, and (b) LD pruning has stochastic tie-breaking. An empirical check between GPA's `popFiles/` snapshot from `2026-03-28` and upstream run `results_20260421_200003` on the **same 560 samples** showed:

- **Filtered SNP set:** 57,829 → 58,411 (+590 new, 8 dropped, 57,821 common).
- **LD-pruned SNP set:** 45,576 → 46,096 (+2,080 new, 1,560 dropped, 44,016 common).
- **Sample set:** unchanged (560 ↔ 560, exact IID overlap).

Because the SNP set drifts, every wholePop file downstream of the filtered / pruned SNP sets drifts numerically as well. The magnitude is small for most metrics but non-zero. The table below states the claim for pipeline-logic equivalence (what a re-run on the *same* filtered VCF would produce) and, separately, the measured drift against the legacy popFiles bundle so GPA does not build a byte-equality regression test by mistake.

| File | wholePop (pipeline-logic equivalence, same input) | wholePop measured drift vs legacy popFiles (2026-03-28) | Trotter / Pacer |
|---|---|---|---|
| `pruned.${rg}.afreq` | Identical to the old `pl1_pruned.freq_stats.afreq` when run on the same pruned SNP set. | On the 44,016 common SNPs: `ALT_FREQS` exactly equal (max \|Δ\|=0). File as a whole differs by 1,560 legacy-only + 2,080 current-only rows. | **New content.** AF from group members only; a locus common in whole-pop but rare in Trotter gets a different number. |
| `freqs.${rg}.tab.gz` (+ `.tbi`) | Byte-identical to GPA's legacy self-built `freqs.tab.gz` when built on the same whole-pop phased VCF via the same `bcftools +fill-tags \| bcftools query` pipeline. | On the 57,806 common (CHROM, POS) pairs: mean \|ΔAF\|=3e-6, RMS=6e-5, max=3.6e-3 (floating-point noise). Plus 605 current-only + 23 legacy-only sites from the filtered-set drift. | **New content.** AF computed over group members only via `bcftools +fill-tags` on the group-subset VCF. Changes `bcftools roh` calibration for any downstream ROH call that uses `--AF-file freqs.${rg}.tab.gz` on a candidate animal. |
| `filtered.LD_prune.het_stats.${rg}.het` | Numerically equivalent to the old `filtered.LD_prune.het_stats.het` on the same pruned set. We pass `--read-freq pruned.wholePop.afreq`, which equals PLINK2's internal default for the whole-pop case. | On the 560 common IIDs: mean \|ΔF\|=0.0019, RMS=0.0053, **max=0.060**. Inherits the pruned-set drift via `--read-freq`. | **New content.** `F` reflects inbreeding relative to group AF, not whole-pop AF. |
| `roh_summary_by_RG_L3_Froh.${rg}.txt` | Identical to the old bare whole-pop F_ROH summary given the same phased VCF. `bcftools roh --estimate-AF -` on the wholePop subset == on the full VCF. | On the 560 common IIDs: mean \|ΔF_ROH\|=0.0005, RMS=0.0008, max=0.0047. Drift is tiny because bcftools roh is robust to the filtered-SNP-set churn. | **New content.** `bcftools roh` is called on the group-subset VCF, so ROH segments themselves differ — this is the "option (ii)" decision from the refactor brainstorm. |
| `Inbreeding_Comparison.${rg}.csv` / `Pairwise_Differences.${rg}.csv` | Identical to the primary-cutoff (1.0 Mb) whole-pop output from the old Section-6 ROHRM loop given the same inputs. | Current wholePop has 560 rows / 156,520 pairs (all samples). Legacy popFiles has 542 rows because pre-fix `scripts/analysis_comparison.py` silently dropped the 18 gait-less samples — fixed in the per-group-references PR; resolved via `ISSUE_wholePop_ROHRM_sample_drop.md`. On the 542 common IIDs: D_STD mean \|Δ\|=0.0017 max=0.095; D_ROH mean \|Δ\|=0.0003 max=0.001. D_STD drift tracks the pruned-set churn; D_ROH is essentially stable. The 18 newly-included samples carry `Phenotype=Unknown`. | **New content.** ROHRM run on the group-subset phased VCF; GRM on the group-subset BED with `--read-freq pruned.${rg}.afreq`. |
| `roh.L3.consensus_25pct.merged.${rg}.smoothed.bed` | Identical to the old per-subpop consensus BED for wholePop given the same ROH calls. | Region count differs (legacy 67 → current 68); boundaries drift because the underlying per-base ROH coverage changes with the filtered-set churn. | **New content** (option ii). Consensus is built from the group's own ROH calls, so Trotter's islands reflect Trotter-specific selection / drift, not whole-pop. |
| `roh.L3.perSample_intersect_${rg}_consensus_25pct.summary.txt` | Identical to today's wholePop ROH_sh table given the same consensus BED. | On 560 common IIDs: mean \|Δ ROH_sh\|=0.20 pct-points, max=1.50 pct-points. Inherits the consensus-BED drift. | **New content.** Derived from the group's own consensus BED above. |

**Practical implication for GPA reports:**
- **Tab-shape expectations:** tabs that used to show identical distributions (whole-pop everywhere) will now show different distributions on Trotter / Pacer tabs for F_SNP, F_ROH, D_SNP, D_ROH, G_SNP, G_ROH, and ROH_sh. The wholePop tab should *shape*-match the legacy report — distribution histograms and percentile lines will look the same to the eye — but individual sample values will shift by the drift magnitudes above.
- **Do not build a byte-equality regression test** against the March-28 popFiles bundle. The SNP set has moved on since. A useful regression test is a shape/percentile check (e.g. P25/P50/P75 of each distribution within ±0.01 of legacy) or a re-run of this very comparison on a freshly-refreshed upstream `OUTPUT_DIR`.
- **`ISSUE_wholePop_ROHRM_sample_drop.md` is now resolved** in `scripts/analysis_comparison.py` (the bug-line filter at L48 was replaced with an empty-pheno_map guard). The drift table above already reflects the fixed-vs-legacy comparison.

---

## 4. Checklist for the GPA-side `create_popFiles.sh` PR

**Branch name convention.** Open the GPA-side PR on a branch named **`per-group-references`** (no prefix). The upstream `genDiv` work lives on the same branch name there — keeping the two in lockstep makes the coordinated merge easier to reason about in PR URLs and chat threads. Neither repo's branch is merged yet as of this doc's writing; pair-merge them when both sides are green.

Mechanical steps. Run through in order:

1. **Pin the upstream pipeline version** in the GPA repo's README or a constants file so it's obvious which `genDiv` commit this popFiles format corresponds to.
2. **Update the copy block** in `create_popFiles.sh` to iterate over `rg ∈ {wholePop, Trotter, Pacer}` and copy the nine per-group files listed in §1's table. Template:
   ```bash
   for rg in wholePop Trotter Pacer; do
       cp "${OUTPUT_DIR}/LD_pruned/pruned.${rg}.afreq"                                             popFiles/.
       cp "${OUTPUT_DIR}/divStats/freqs.${rg}.tab.gz"                                              popFiles/.
       cp "${OUTPUT_DIR}/divStats/freqs.${rg}.tab.gz.tbi"                                          popFiles/.
       cp "${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.${rg}.het"                           popFiles/.
       cp "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.${rg}.txt"                             popFiles/.
       cp "${OUTPUT_DIR}/rep_ROHRM/roh_1Mb.Threshold_3SD/Inbreeding_Comparison.${rg}.csv"          popFiles/.
       cp "${OUTPUT_DIR}/rep_ROHRM/roh_1Mb.Threshold_3SD/Pairwise_Differences.${rg}.csv"           popFiles/.
       cp "${OUTPUT_DIR}/divStats/roh.L3.consensus_25pct.merged.${rg}.smoothed.bed"                popFiles/.
       cp "${OUTPUT_DIR}/divStats/roh.L3.perSample_intersect_${rg}_consensus_25pct.summary.txt"    popFiles/.
   done
   cp "${OUTPUT_DIR}/preprocess/sample_groups.tsv"                                                 popFiles/.
   ```
3. **Remove any reads of the old bare names** (`filtered.LD_prune.het_stats.het`, `roh_summary_by_RG_L3_Froh.txt`, `Inbreeding_Comparison.csv`, `Pairwise_Differences.csv`, `pruned.freq_stats.afreq`). `git grep` for each bare name to confirm no leftover references in GPA's R scripts or `generate_report.py`.
4. **PCA-artifact reads (§1.1) — done; verify only.** `create_popFiles.sh` already reads the infixed `divStats/filtered.LD_prune.wholePop.pca.eigenvec.allele` and flattens it to `popFiles/eigenvec.allele`. Confirm no bare-name read crept back in — `git grep -n 'LD_prune\.pca\.'` over the GPA repo should return nothing — and leave the rest of the PCA family alone; GPA reads none of it.
5. **Delete the in-repo `freqs.tab.gz` self-build** in `create_popFiles.sh` — specifically the `bcftools +fill-tags $popVCF | bcftools query … > popFiles/freqs.tab.gz` block and the follow-up `tabix`. That file is now shipped per-group by upstream (step 2 above).
6. **Teach `GPA.sh` to pick the right `--AF-file` per report tab.** At present it runs `bcftools roh -G30 --AF-file $allele_freqs $vcf_filtered` exactly once with a single whole-pop AF file, then the three tabs all read from the same outputs. For the per-group refactor to do real work, GPA needs to loop `bcftools roh` three times — once per `rg` — each with `--AF-file popFiles/freqs.${rg}.tab.gz`, producing per-group `roh_out.${rg}.txt` → L1/L2/L3 → `roh_summary_by_RG_L3_Froh.${rg}.txt` → F_ROH on the candidate, and pipe each into the matching tab's histogram. This is the only non-mechanical GPA change in this migration.
7. **Add a validator call** at the end of `create_popFiles.sh` (see §5 below).
8. **Update the report generator** (`generate_report.py` and the `plot_*.R` set in `../GPA/`) to read `${rg}.` suffixed filenames for all three tabs. Column indices unchanged; only filenames change.
9. **Regenerate test fixtures** if GPA has recorded-output tests pinned to the old bare-named files.
10. **Before merge:** run a full GPA report end-to-end on a populated `popFiles/` dir and confirm each tab (wholePop, Trotter, Pacer) renders with a distinct distribution for F_SNP, F_ROH, D_SNP, G_SNP, D_ROH, G_ROH, and ROH_sh. (The wholePop tab should look like today's single report; the other two should be new.) Expected direction for F_ROH on the candidate animal: on the Trotter tab its F_ROH should shift relative to the legacy whole-pop calibration, and the shift's sign should match the Fst landscape between Trotter and wholePop at the variants in that animal's ROH segments.

---

## 5. Using `validate_popRefs.sh` in the GPA repo

The upstream pipeline ships a validator at `scripts/validate_popRefs.sh` that can run in two modes.

### On the upstream side (automatic)

Run after `bash genDiversity.sh` to catch upstream drift:

```bash
bash scripts/validate_popRefs.sh --mode upstream --root "${OUTPUT_DIR}"
```

Exits 0 if every per-group file exists, is non-empty, has the right header, and meets its row-count lower bound. Exits 1 with a specific failure list otherwise.

### On the GPA side (copy into GPA repo)

**Installation:** copy `scripts/validate_popRefs.sh` from this repo into the GPA repo (for example as `GPA/scripts/validate_popRefs.sh`). Commit it. Future upstream edits to the script can be pulled manually when the popRefs schema evolves.

**Invocation:** append to the bottom of `GPA/create_popFiles.sh`:

```bash
bash scripts/validate_popRefs.sh --mode popFiles --root "$gpa_root"
```

where `$gpa_root` is the parent directory containing `popFiles/` (typically the GPA repo root or whatever `$results_dir` the GPA wrapper uses). The script appends `popFiles/` itself, mirroring how `--mode upstream` appends each per-file subdir under `$OUTPUT_DIR`. The popFiles tree is flat — every file has the same basename as upstream and lives directly under `$gpa_root/popFiles/`. Exits 1 on any missing / malformed file, so the GPA pipeline stops before a broken report is generated.

**Adapting it:** the list of files checked lives in the `FILE_SPECS` array at the top of the script. Each row is `key|upstream_subdir|filename_template|header_regex|min_rows_formula`. If GPA doesn't need a particular file, delete that row in the GPA copy. If GPA adds a new file GPA wants validated, add a new row. The `upstream_subdir` column is ignored in `--mode popFiles`, so only the filename template and schema fields matter for the GPA copy.

**Row-count formulas** use `NGRP` as the group size (looked up from `samples.${rg}.txt` if present, else `0`). Use `NGRP+1` for one-row-per-sample files, `(NGRP*(NGRP-1))/2+1` for pairwise files, and the literal string `any` to skip the row-count check.

**Header regexes** are extended-regex (fed to `grep -E`). Keep them anchored (`^...$`) wherever possible; a passing match on a header is evidence the schema hasn't silently drifted. If upstream genuinely changes a header, update the script in both repos in the same PR pair — the script is literally the schema contract.

**Scope — the PCA artifacts are deliberately out of the upstream validator.** Upstream's `FILE_SPECS` covers the nine per-group reference files of §1 and nothing else. The PCA family of §1.1 / §2 is *not* validated there, by design: PC count, eigenvalue count, and the `eigenvec.allele` SNP-row count all move with the input data, so a row-count assertion would fail on ordinary reruns.

Existence checks for the one PCA file that matters are the consumer's concern, and the GPA copy already does this. Note that the GPA copy has **diverged from upstream by adding a second array**, `POPFILES_INVARIANT_SPECS` — 11 rows covering the rg-invariant bundle files that exist only in the flat `popFiles/` tree: `EquCab3_map`, the four SNP-list files, `pop.vcf.gz` (+ `.tbi`), `effective_autosomal_genome_length.txt`, `gait_model.pkl`, and the PCA row `pca_allele_wt|eigenvec.allele||any`. That array has no upstream counterpart — upstream never produces the flattened names — so when pulling a newer `validate_popRefs.sh` from this repo into GPA, port the `FILE_SPECS` changes and keep the local `POPFILES_INVARIANT_SPECS` block rather than overwriting the file wholesale.

If a future change makes a PCA file's schema part of the contract, add it to upstream `FILE_SPECS` at that point.

---

## 6. When upstream has to change a file

If upstream needs to rename a file, add/remove a column, or change a file's content semantics:

1. Update `FILE_SPECS` in `scripts/validate_popRefs.sh` upstream.
2. Update `§1` (paths) and `§2` (schemas) of this document in the same commit.
3. Open a coordinated GPA-side PR that updates `create_popFiles.sh` and any downstream reads, plus copies the new `validate_popRefs.sh` over.
4. Merge both within the same window. If that's not possible, add transitional symlinks on the upstream side pointing old names at new names, remove them once the GPA PR lands.
