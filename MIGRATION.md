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
Tab-separated. Header starts with `#FID` (PLINK2 output).
```
#FID   IID   O(HOM)   E(HOM)   OBS_CT   F
 0     1      2         3         4      5
```
GPA reads `F` (F_SNP) at 0-based column 5.

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

---

## 3. Numerical equivalence vs genuine change

This table exists to avoid surprise when comparing the new per-group files to the original pipeline's single whole-pop output.

| File | wholePop | Trotter / Pacer |
|---|---|---|
| `pruned.${rg}.afreq` | Identical to the old `pl1_pruned.freq_stats.afreq` (PLINK2 `--freq` on the same samples). | **New content.** AF from group members only; a locus common in whole-pop but rare in Trotter gets a different number. |
| `freqs.${rg}.tab.gz` (+ `.tbi`) | Byte-identical to GPA's legacy self-built `freqs.tab.gz` when built on the same whole-pop phased VCF via the same `bcftools +fill-tags \| bcftools query` pipeline. | **New content.** AF computed over group members only via `bcftools +fill-tags` on the group-subset VCF. Changes `bcftools roh` calibration for any downstream ROH call that uses `--AF-file freqs.${rg}.tab.gz` on a candidate animal. |
| `filtered.LD_prune.het_stats.${rg}.het` | Numerically equivalent to the old `filtered.LD_prune.het_stats.het`. We pass `--read-freq pruned.wholePop.afreq`, which is identical to PLINK2's internal default AF for the whole-pop case. | **New content.** `F` reflects inbreeding relative to group AF, not whole-pop AF. |
| `roh_summary_by_RG_L3_Froh.${rg}.txt` | Identical to the old bare whole-pop F_ROH summary. `bcftools roh --estimate-AF -` on the wholePop subset == on the full VCF. | **New content.** `bcftools roh` is called on the group-subset VCF, so ROH segments themselves differ — this is the "option (ii)" decision from the refactor brainstorm. |
| `Inbreeding_Comparison.${rg}.csv` / `Pairwise_Differences.${rg}.csv` | Identical to the primary-cutoff (1.0 Mb) whole-pop output from the old Section-6 ROHRM loop. | **New content.** ROHRM run on the group-subset phased VCF; GRM on the group-subset BED with `--read-freq pruned.${rg}.afreq`. |
| `roh.L3.consensus_25pct.merged.${rg}.smoothed.bed` | Identical to the old per-subpop consensus BED for wholePop. | **New content** (option ii). Consensus is built from the group's own ROH calls, so Trotter's islands reflect Trotter-specific selection / drift, not whole-pop. |
| `roh.L3.perSample_intersect_${rg}_consensus_25pct.summary.txt` | Identical to today's wholePop ROH_sh table. | **New content.** Derived from the group's own consensus BED above. |

**Practical implication for GPA reports:** tabs that used to show identical distributions (whole-pop everywhere) will now show different distributions on Trotter / Pacer tabs for F_SNP, F_ROH, D_SNP, D_ROH, G_SNP, G_ROH, and ROH_sh. The wholePop tab should continue to match the current report's numbers.

---

## 4. Checklist for the GPA-side `create_popFiles.sh` PR

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
4. **Delete the in-repo `freqs.tab.gz` self-build** in `create_popFiles.sh` — specifically the `bcftools +fill-tags $popVCF | bcftools query … > popFiles/freqs.tab.gz` block and the follow-up `tabix`. That file is now shipped per-group by upstream (step 2 above).
5. **Teach `GPA.sh` to pick the right `--AF-file` per report tab.** At present it runs `bcftools roh -G30 --AF-file $allele_freqs $vcf_filtered` exactly once with a single whole-pop AF file, then the three tabs all read from the same outputs. For the per-group refactor to do real work, GPA needs to loop `bcftools roh` three times — once per `rg` — each with `--AF-file popFiles/freqs.${rg}.tab.gz`, producing per-group `roh_out.${rg}.txt` → L1/L2/L3 → `roh_summary_by_RG_L3_Froh.${rg}.txt` → F_ROH on the candidate, and pipe each into the matching tab's histogram. This is the only non-mechanical GPA change in this migration.
6. **Add a validator call** at the end of `create_popFiles.sh` (see §5 below).
7. **Update the report generator** (`generate_report.py` and the `plot_*.R` set in `../GPA/`) to read `${rg}.` suffixed filenames for all three tabs. Column indices unchanged; only filenames change.
8. **Regenerate test fixtures** if GPA has recorded-output tests pinned to the old bare-named files.
9. **Before merge:** run a full GPA report end-to-end on a populated `popFiles/` dir and confirm each tab (wholePop, Trotter, Pacer) renders with a distinct distribution for F_SNP, F_ROH, D_SNP, G_SNP, D_ROH, G_ROH, and ROH_sh. (The wholePop tab should look like today's single report; the other two should be new.) Expected direction for F_ROH on the candidate animal: on the Trotter tab its F_ROH should shift relative to the legacy whole-pop calibration, and the shift's sign should match the Fst landscape between Trotter and wholePop at the variants in that animal's ROH segments.

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
bash scripts/validate_popRefs.sh --mode popFiles --root "$popFiles_dir"
```

where `$popFiles_dir` is the target directory the copy step just populated (typically `popFiles/`). The script treats `popFiles/` as flat — every file is expected directly under `--root`, without subdirectories. Exits 1 on any missing / malformed file, so the GPA pipeline stops before a broken report is generated.

**Adapting it:** the list of files checked lives in the `FILE_SPECS` array at the top of the script. Each row is `key|upstream_subdir|filename_template|header_regex|min_rows_formula`. If GPA doesn't need a particular file, delete that row in the GPA copy. If GPA adds a new file GPA wants validated, add a new row. The `upstream_subdir` column is ignored in `--mode popFiles`, so only the filename template and schema fields matter for the GPA copy.

**Row-count formulas** use `NGRP` as the group size (looked up from `samples.${rg}.txt` if present, else `0`). Use `NGRP+1` for one-row-per-sample files, `(NGRP*(NGRP-1))/2+1` for pairwise files, and the literal string `any` to skip the row-count check.

**Header regexes** are extended-regex (fed to `grep -E`). Keep them anchored (`^...$`) wherever possible; a passing match on a header is evidence the schema hasn't silently drifted. If upstream genuinely changes a header, update the script in both repos in the same PR pair — the script is literally the schema contract.

---

## 6. When upstream has to change a file

If upstream needs to rename a file, add/remove a column, or change a file's content semantics:

1. Update `FILE_SPECS` in `scripts/validate_popRefs.sh` upstream.
2. Update `§1` (paths) and `§2` (schemas) of this document in the same commit.
3. Open a coordinated GPA-side PR that updates `create_popFiles.sh` and any downstream reads, plus copies the new `validate_popRefs.sh` over.
4. Merge both within the same window. If that's not possible, add transitional symlinks on the upstream side pointing old names at new names, remove them once the GPA PR lands.
