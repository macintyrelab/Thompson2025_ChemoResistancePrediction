# Predicting response to cytotoxic chemotherapy

This repo holds the code and data needed to reproduce the analysis for *Predicting resistance to chemotherapy using chromosomal instability signatures* (Thompson, Madrid, Hernando, et al. Nat Genetics, 2025).

## Data

All data for reproducing analyses and figures can be downloaded from [figshare](https://figshare.com/s/09bbc25f748991a248a6). Once downloaded, please add all of them in a folder named as `Input_Data` at the root of the repo to ensure smooth execution of the code.

## Linking code to manuscript figures (main, extended and supplement)

All the code is in the `Analysis` folder, which contains 8 directories which we recommend to run in the following order to reproduces each figure in the manuscript:

-   [Implementing_Platins_Response_Predictor](https://github.com/macintyrelab/ChemoResistantPredictionCIN/tree/main/Analysis/Implementing_Platins_Response_Predictor): This repository contains code to apply the biomarker for predicting platin resistance at pan-cancer level, and reproduce Figure 1b

-   [Optimizing_Taxanes_Response_Predictor](https://github.com/macintyrelab/ChemoResistantPredictionCIN/tree/main/Analysis/Optimizing_Taxanes_Response_Predictor): This repository contains code to optimize the biomarker for predicting taxane resistance using cell lines from the DepMap project. This code reproduces Figure 1d, Supplementary Figure 1, Supplementary Table 2

-   [Identifying_an_Anthracycline_Response_Predictor](https://github.com/macintyrelab/ChemoResistantPredictionCIN/tree/main/Analysis/Identifying_an_Anthracycline_Response_Predictor): This repository contains code to develop and optimize the biomarker for predicting anthracycline resistance:

    -   The `dox_mechanism.Rmd` script reproduces analysis of Supplementary Note 1, Supplementary Figures 78-79

    -   The `transcriptomic_tcga_data.R` script is needed to get transcriptomic data from TCGA

    -   The `dox_optimal_threshold.R` script is used to optimize the threshold for biomarker classification, and reproduces Figure 1f

-   [Predicting_Chemotherapy_Response_in_Ovarian_Cancer](https://github.com/macintyrelab/ChemoResistantPredictionCIN/tree/main/Analysis/Predicting_Chemotherapy_Response_in_Ovarian_Cancer): This repository contains code to reproduce the pilot study performed in the OV04 cohort:

    -   The `OV04_REMARK_Diagram.R` script reproduces Extended Data Figure 2

    -   The `generate_clin_history_plots.R` script reproduces plots summarizing clinical history of OV04 patients (Supplementary Figures 24-73). This also reproduces swimmer plots of Supplementary Figures 2-4.

    -   The `OV04_Survival_Analysis.R` script is used to perform survivial and complementary analyses for the pilot study, and reproduces Figure 2, Extended Data Figure 3, Supplementary Table 6

-   [Predicting_Chemotherapy_Response_Pancancer](https://github.com/macintyrelab/ChemoResistantPredictionCIN/tree/main/Analysis/Predicting_Chemotherapy_Response_Pancancer): This repository contains code to reproduce the performance assessment of our biomarkers by emulating trials:

    -   The `Curating_HMF_Clinical_Data.R` and `Curating_TCGA_Clinical_Data.R` scripts are used to curate clinical data from public datasets

    -   The `Generating_supplementary_tables.R` script is used to create Supplementary Table 5, which shows clinical traits of the different tumour-type HMF and TCGA cohorts.

    -   The `Phase_III_general_filtering_diagram.R` script reproduces Supplementary Figure 80, which is a flowchart summarizing inclusion and exclusion criteria

    -   The `Pancancer_survival_analysis.R` script is used to perform survivial analyses across primary and metastatic samples from different tumour types. This script reproduces Figure 3-4, Extended Data Figures 6-7, Extended Data Tables 1-2, Supplementary Note 2, Supplementary Figures 5-21, Supplementary Figures 74-77, Supplementary Figure 4.

    -   The `HRDetect_MyriadmyChoice_analysis.R` script reproduces Supplementary Note 3

-   [Predicting_Response_Using_TSO500](https://github.com/macintyrelab/ChemoResistantPredictionCIN/tree/main/Analysis/Predicting_Response_Using_TSO500): This repository contains code to reproduce the clinical utility of our biomarkers by using TSO500 data (Figure 5a-b, Supplementary Figure 22, Extended Data Table 3)

-   [Predicting_Response_Using_Liquid_Biopsy](https://github.com/macintyrelab/ChemoResistantPredictionCIN/tree/main/Analysis/Predicting_Response_Using_Liquid_Biopsy): This repository contains code to reproduce the clinical utility of our biomarkers by using liquid biopsies (Figure 5c-d, Supplementary Figure 23, Extended Data Table 3)

-   [Helper_Scripts](https://github.com/macintyrelab/ChemoResistantPredictionCIN/tree/main/Analysis/Helper_Scripts): This repository contains all functions needed for running all scripts

## Contact
If you experience any issues or have questions about the code, please open a Github issue with a minimum reproducible example. For questions around collaborations or sensitive patient data, please contact us directly at Geoff Macintyre gmacintyre@cnio.es or Anna M. Piskorz anna.piskorz@cruk.cam.ac.uk.

## License
The contents of this repository are copyright (c) 2025, Tailor Bio Ltd and Spanish National Cancer Research Centre (CNIO).

The contents of this repository are published and distributed under the GAP Available Source License v1.0 (ASL).

The contents of this repository are distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the ASL for more details.

The methods implemented in the code are the subject of a patent on using copy number signatures to predict response to doxorubicin treatment in ovarian cancer (PCT/EP2021/065058), and a patent on a method for identifying pan-cancer copy number signatures (PCT/EP2022/077473).

Any commercial use of this code is prohibited.
