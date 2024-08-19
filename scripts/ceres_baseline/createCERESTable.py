import pandas as pd
import matplotlib.pyplot as plt

k562_columns_of_interest = [
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


k562_df = (
    pd.read_csv("given_data/K562-DHS-v6.noChIP.csv")[k562_columns_of_interest]
    .pipe(
        lambda df: pd.get_dummies(
            df,
            columns=df.select_dtypes(include=["object"]).columns,
            prefix_sep="_",
            drop_first=False,
        )
    )
    .dropna()
    .reset_index(drop=True)
)


k562_df.to_csv(
    "created_data/ceres_baseline/CERES.csv",
    index=False,
)
