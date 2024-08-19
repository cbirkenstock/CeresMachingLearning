import pandas as pd
import pyranges as pr

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)  # Display all columns

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

CERES_df = pd.read_csv("given_data/K562-DHS-v6.noChIP.csv")[
    CERES_columns_of_interest
].rename(columns={"chrom": "Chromosome", "chromStart": "Start", "chromEnd": "End"})


phyloP_df = (
    pd.read_csv("given_data/K562.selection.csv")
    .assign(
        **{
            "Chromosome": lambda df: df["chromosome"].str.split(":|-", expand=True)[0],
            "Start": lambda df: df["chromosome"].str.split(":|-", expand=True)[1],
            "End": lambda df: df["chromosome"].str.split(":|-", expand=True)[2],
        }
    )
    .drop(["chromosome", "gene"], axis=1)
)

CERES_gr = pr.PyRanges(CERES_df)
phyloP_gr = pr.PyRanges(phyloP_df)

merged_df = CERES_gr.join(phyloP_gr, strandedness=False).df.pipe(
    lambda df: pd.get_dummies(
        df,
        columns=df.select_dtypes(include=["object", "category"]).columns,
        prefix_sep="_",
        drop_first=False,
    )
    .dropna()
    .reset_index(drop=True)
)

merged_df.to_csv("created_data/phyloP/phyloP.csv", index=False)
