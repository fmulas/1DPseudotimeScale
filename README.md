# 1DPseudotimeScale
This is a collection of R and Python Scripts and data files used to analyze 

single-cell data related to the Manuscript "Pseudotemporal ordering of single 

cells reveals metabolic control of postnatal beta cell proliferation".


FOLDERS:

- "data": contains single cell data and bulk data described in the manuscript, as 

well as gene lists used in the analyses.

- "res": contains results produced by other tools and will store results from 

the scripts included in the folder.



SCRIPTS:

1_batch_eff_removal_svaseq.R: identifies negative control genes and remove batch effect using SVASeq. Selects most variant genes using MAD and saves the data in Orange format.


2_1D_Pseudotime_Scale.py: builds the 1D Pseudotime Scale starting from normalized data. Saves cell 

projections and plots their values. Also saves projections of external data on a scale previously 

defined and saves predicted cell ordering on bootstrap samples for method validation.

Requires several Python and Orange (https://orange.biolab.si/) modules. 

Data files have to be in the Orange format. See the script for details.

A widget-based implementation of PCA data projection, related to a previous manuscript 
(DOI: 10.1093/bioinformatics/btr422) can be downloaded at: https://bitbucket.org/biolab/orange-differentiation).



3_ordering_validation_POS.R: compute CI of POS score based on bootstrap distribution of cell coordinates

inferred with 1D Pseudotime Scale; computes POS Score on pseudotime coordinates obtained from 

different unsupervised cell ordering methods.



4_ordering_validation_roughness.py: Computes distances between pairs of contiguously placed cells 

and builds a null distribution of the distances, for several distance metrics; measures 

"Roughness" distance for coordinates obtained with 1DPseudotime and the supervised method 

Delorean. 



5_save_pseudotime_binned_data.R: Saves data with binned points, i.e. with the same number of samples 

for each time point, but ordered according to pseudotime.



6_differential_analysis_binned_and_average.py: Computes fold change and permutation-based p-value 

for consecutive points, both from time-ordered and from pseudotime-ordered data.



7_denovo_genesets_GSEA: clusters cells with a similar pattern along pseudotime and runs continuous 

GSEA to score gene sets correlated with pseudotime coordinates. Requires utils.py and Java 

library to run GSEA (here gsea2-2.2.2.jar was used).



8_corr_individual_genes.R: computes correlation of individual genes with pseudotime coordinates 

and scores each value through a permutation-based null distribution.



9_network_propagation.py: saves a network from STRING in gpickle format for analysis; propagate 

interest in a set of select genes (e.g. TFs) by increasing scores of their neighbors. Related

to manuscript with DOI: 10.3233/978-1-61499-240-0-182. Requires several libraries of the 

Orange software.



10_GSEA_gof_data.R: runs GSEA to score differential expression of selected gene sets in Srf 

overexpression data vs. control data. Requires utils.py and Java library to run GSEA 

(here gsea2-2.2.2.jar was used).

