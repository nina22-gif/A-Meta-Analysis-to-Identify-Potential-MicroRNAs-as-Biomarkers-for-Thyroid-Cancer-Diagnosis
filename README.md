1. Introduction

This project aims to identify potential microRNAs (miRNAs) as diagnostic biomarkers for thyroid cancer through an integrated analysis of publicly available gene expression datasets.

The study focuses on identifying differentially expressed miRNAs (DEMs) in papillary thyroid carcinoma (PTC) and medullary thyroid carcinoma (MTC), as well as exploring potential biomarkers that may distinguish between these two thyroid cancer subtypes.

2. Objectives

Identify differentially expressed miRNAs in PTC and MTC compared with normal thyroid tissues.

Integrate results from multiple GEO datasets to identify consistently dysregulated miRNAs.

Explore potential miRNA–mRNA regulatory interactions.

Identify candidate miRNA biomarkers specific to each thyroid cancer subtype.

3. Data Sources

The analysis uses publicly available datasets from the Gene Expression Omnibus (GEO):

GSE40807

GSE73182

GSE97070

GSE113629

GSE104006

4. Methods

Differential expression analysis using the limma package in R.

Identification of differentially expressed miRNAs.

Meta-analysis and integration of results across datasets.

Functional enrichment analysis.

miRNA–mRNA interaction network analysis.

5. Repository Structure

Data/: Input datasets and processed data.

R_Scripts/: R scripts used for data analysis.

Results/: Figures, tables, and analysis results.

References/: References and supporting documents.

6. Software and Tools

R and RStudio

limma

GEO database

Cytoscape

STRING

multiMiR
