import pandas as pd
import pyranges as pr


def classify_tad_relation(row):
    dhs_start, dhs_end = row["Start"], row["End"]
    gene_start, gene_end = row["gene_start"], row["gene_end"]
    tad_start, tad_end = row["Start_b"], row["End_b"]

    dhsStartOverlap = dhs_start > tad_start
    dhsEndOverlap = dhs_end < tad_end
    dhsRelation = ""

    if dhsStartOverlap and dhsEndOverlap:
        dhsRelation = "full"
    else:
        dhsRelation = "partial"

    geneStartOverlap = gene_start > tad_start
    geneEndOverlap = gene_end < tad_end
    geneRelation = ""

    if geneStartOverlap and geneEndOverlap:
        geneRelation = "full"
    elif geneStartOverlap or geneEndOverlap:
        geneRelation = "partial"
    else:
        geneRelation = "off"

    return dhsRelation, geneRelation


tads_df = pd.read_csv(
    "C:/Users/Ictinike/Documents/WrayLab/raw_data/OCRs_TADdomains_int.bed",
    delimiter="\t",
    header=None,
)

tads_df.columns = ["_", "_", "_", "chrom", "Start", "End", "TAD_ID"]
tads_df = tads_df[["chrom", "Start", "End", "TAD_ID"]]
tads = pr.PyRanges(tads_df)

DHS_gene_df = pd.read_csv("given_data/K562-DHS-v6.noChIP.csv").rename(
    columns={"chrom": "Chromosome", "chromStart": "Start", "chromEnd": "End"}
)

DHSs = pr.PyRanges(DHS_gene_df)

dhs_tads_df = DHSs.join(tads).as_df()

dhs_tads_df["dhsRelation"], dhs_tads_df["geneRelation"] = dhs_tads_df.apply(
    classify_tad_relation, axis=1
)
