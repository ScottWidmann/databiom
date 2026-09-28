# Data Access Notes

## `breast_cancer.txt` (pending removal)

`breast_cancer.txt` was a real clinical mutation-annotation (MAF-format) file
used in the mutation/functional-enrichment portion of the bioinformatics
lecture. It contains real tumor/normal sample identifiers
(`Tumor_Sample_Barcode`, `Matched_Norm_Sample_Barcode`,
`Tumor_Sample_UUID`, etc.) consistent with an MSK-IMPACT-style clinical
sequencing cohort (e.g., data distributed via
[cBioPortal](https://www.cbioportal.org)). This is third-party clinical
research data, and redistributing it directly in this repository is not
appropriate without confirming the originating study's data-use terms.

**This file is being removed from the repository.** To reproduce the same
kind of analysis:

1. Go to [cbioportal.org](https://www.cbioportal.org) and search for a
   breast cancer study using MSK-IMPACT clinical sequencing (patient IDs in
   the `P-0000000-T00-IM0` format).
2. Download the mutation data (MAF) for that study directly from
   cBioPortal, following their terms of use and citing the study.
3. Cite cBioPortal itself in any derived materials:
   - Cerami et al., *Cancer Discovery*, 2012.
   - Gao et al., *Sci. Signal.*, 2013.

The lecture slide that previously showed a `wget` command pulling this file
directly from this repository (`Lectures/Introduction to Bioinformatics
DATABIOM.pdf`, "Using the wget command to download from internet") has had
that specific example redacted, since it would otherwise 404 once the file
is removed. It should be replaced with a public-data example (e.g., a GEO/
SRA accession, matching the rest of the course) the next time that deck is
updated.
