# ATLASx Analysis tools 

The data and scripts contained in this repository allow the user to reproduce the figures in the ATLASx manuscript from the original data files.
Reference to the article (doi): https://doi.org/10.1101/2021.02.17.431583


## 1) Installation

The installation can be completed in less than 10 minutes, including installation of dependencies and fetching the data from the git repository.

### Requirements

- python 2.7 or higher (Tested and recommended: python 3.9.6)
- pip

The code is adapted for python 3, but it can also be executed in python 2.s

Runtimes are indicated for each script and were determined normal desktop computer (macOS).

### Download repository

`$ git clone https://github.com/EPFL-LCSB/ATLASxAnalyses`

To install the required dependencies: 

`$ cd ATLASxAnalysis`
`$ make`

### Note
Data files are stored using git large file storage (lfs). The make file will install git lfs automatically. However, if lfs was not installed previously, the repository has to be updated after installation:

`$ git pull`

This is needed to retrieve the data files from the repository after installation.


## 2) Reproduce Network Analysis

`$ cd NetworkAnalysis/Source`

Plot component distribution for database scopes

`$ python3 get_component_distribution.py #Runtime: 97s`

By default, the database scopes will be plotted. For a resolution by data source, add data_sources as an argument to the above command:

`$ python3 get_component_distribution.py data_sources #Runtime: 11s`

### CSV file conversion
CVS files of networks are quite practical for visualisation, e.g. in the Gephi software.
To convert the gpickle files for database scopes to CSV, run the following:

`$ python3 print_csv_from_gpickle.py #Runtime: 52s`

As above, for single data sources use:

`$ python3 print_csv_from_gpickle.py data_sources #Runtime: 50s`

The output is automatically written to a new folder called ATLASxAnalyses_output created in the same directory as the repository.

The data files used in the repository is the same as the one used in the manuscript, which has been downloaded on 8 November 2020.
For updated network files, please contact the authors of the paper directly.



## 3) Reproduce MetaCyc pathway coverage plot

The data and scripts contained in this repository allow the user to reproduce the
Figure 3: "Pathway search comparison to dataset of pathways extracted from MetaCyc"
of the ATLASx manuscript from the original data files and the Supplementary Figures S3 and S4:
"Pathway search comparison to dataset of pathways extracted from MetaCyc for
BNICE.ch-curated reactions" and "Distribution of the length (i.e., number of reaction steps)
of reconstructed MetaCyc pathways."

`$ cd MetaCycPWanalysis`

Perform the pathway search for MetaCyc, MetaCyc BNICE-curated, ATLASx and ATLASx BNICE-curated
(running time: 5 min for edges coverage, approx. 2.5 days on 10 cores for pathway search)

`$ python3 check_coverage_and_rank_pathways.py`

Running this script will produce 5 output files within the output folder:
- networkEdgesCoverage.csv
- pw_ranking_chemATLAS_hp.csv
- pw_ranking_chemATLAS.csv
- pw_ranking_MetaCyc_hp.csv
- pw_ranking_MetaCyc.csv

As the script takes long time to be executed we recommend to use the existing files
that are stored in "output_generated_for_article" folder for plotting (copy them to output folder):

To generate the plots run:

`$ python3 plot_edges_coverage.py`

`$ python3 plot_pathway_rank.py`

The 2 scripts will generate the  figures used in the publication within the "plots" folder:
MetaCycVsATLASx.png - figure 3A of the manuscript
MetaCycVsATLASxBNICE.png - Supplementary Figure 3A
pathways_rank_chemATLAS_hp.png  - figure 3B of the manuscript
pathways_rank_chemATLAS.png - Supplementary Figure 3B
pathways_rank_MetaCyc_hp.png  - figure 3B of the manuscript
pathways_rank_MetaCyc.png - Supplementary Figure 3B
plotMetacycPathwayLengthDistributionAll_all.png - main part of the Supplementary Figure 4
plotMetacycPathwayLengthDistributionAll_crop.png - cropped in part of the Supplementary Figure 4
