

<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--

-->


<!-- PROJECT SHIELDS -->
<!--

-->
[![Contributors](https://github.com/othneildrew/Best-README-Template/graphs/contributors)[contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



Table of Contents
* [Installation](#Installation)
  * [Getting the code](#getting-the-code)
  * [Creating the enviroment with required dependencies](#Creating-the-enviroment-with-required-dependencies)
* [Getting QC report](#Getting-QC-report)
  * [Configure run](#Configure-run)
  * [Filtration](#Filtration)
* [Output](#Output)

# Ribose-Map_QC
To get quality control of files generated by Ribose-Map
Also supported by percentage and composition of rNMPs from ribose-seq output 

## Installation

### Getting the code
The development version from [GitHub](https://github.com/) with:

```sh
git clone https://github.com/DKundnani/Ribose-Map-QC.git
```

### Creating the enviroment with required dependencies

```sh
conda env create --file rNMP_match_analysis/mm_removal.yml
```

## Getting QC report
### Configure run
```bash
vim MMremoval_configure_run.sh
#Change the variables for you run as per mentions in the bash configure file
```

### Filtration
```bash
bash MMremoval_configure_run.sh

```
featgroup<-grepl( "RNASE",rownames(logTCGA40)) #optional, a set of features to separated
corr_pair<-pairwisecorr <- pairwise_corr(df=logTCGA40,featurelist=rownames(logTCGA40),featuregroup=featgroup)
```

<img src="man/figures/README-multiple-pairwise-correlation-2.png" width="100%" />
