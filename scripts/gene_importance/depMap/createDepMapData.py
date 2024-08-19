import pandas as pd
import re

CERES_columns_of_interest = [
    "chrom",
    "chromStart",
    "chromEnd",
    "pValue",
    "annotation",
    "geneChr",
    "geneStart",
    "geneEnd",
    "geneLength",
    "geneStrand",
    "geneSymbol",
    "distanceToTSS",
    "DHS_length",
    "DHS_prop_repeat",
    "DHS_prop_GC",
    "ploidyZhou",
    "LossHetZhou",
    "SV_Zhou",
    "n_SNV_Zhou",
    "SNV_Zhou",
    "n_SNV_Zhou_per_bp",
    "probIntolerantLoF",
    "probIntolerantLoF_gt_0.9",
    "numTKOHits_Hart",
    "anyTKOHits_Hart",
    "HartEssential",
    "OGEE_n_Essential",
    "OGEE_n_NonEssential",
    "OGEE_n",
    "OGEE_prop_Essential",
    "OGEE_prop_NonEssential",
    "cancer_census_tier",
    "cancer_census_tissue_type",
    "cancer_census_role",
    "n_conserved_LindbladToh",
    "n_conserved_LindbladToh_per_bp",
    "chromHMM_cat_longest",
    "segway_cat_longest",
    "DNase_CPM_per_1kbp",
    "H3K27ac_CPM_per_1kbp",
    "wgCERES_score",
    "dhs_0_1_wg",
]

CERES_df = pd.read_csv("given_data/K562-DHS-v6.noChIP.csv")[CERES_columns_of_interest]

depMap_df = (
    pd.read_csv("given_data/CRISPRGeneDependency.csv")
    .set_index("ModelID")
    .transpose()
    .rename(columns={"ACH-000551": "dependency"})["dependency"]
    .reset_index()
    .rename(columns={"index": "geneSymbol"})
    .assign(
        geneSymbol=lambda df: df["geneSymbol"].apply(
            lambda x: re.sub(r"\s*\(\d+\)", "", x)
        )
    )
)

merged_df = (
    pd.merge(CERES_df, depMap_df, on="geneSymbol", how="inner")
    .assign(
        elem_count=lambda df: df.groupby("geneSymbol")["geneSymbol"].transform("count")
    )
    .drop(["geneSymbol"], axis=1)
    .reset_index(drop=True)
    .pipe(
        lambda df: pd.get_dummies(
            df, columns=df.select_dtypes(include=["object"]).columns
        )
    )
    .dropna()
)

merged_df.to_csv("created_data/gene_importance/DepMapData.csv")
