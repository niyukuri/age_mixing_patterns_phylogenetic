# Inferring age mixing patterns in HIV transmission from phylogenetic trees

A simulation study to understand age mixing patterns in HIV transission networks inferred from phylogenetic trees. From the phylogenetic transmission clusters, we estimated a transmission network using phylogenetic linkage between pairs of sequences. We then inferred measurements that depict age-mixing patterns in transmission i.e. proportions of men/women of different age groups paired with women/men of another age group, and mean and standard deviation of an average age difference between women/men with their respective pairs. The uncertainty around these measurements was assessed as a function of sampling coverage in different sequence missingness scenarios was also assessed.


## CONTENTS

This repo contains the information necessary to reproduce the simulation study:

* [Code and data files](#code-and-data-files)
   * Code files for simulation and post-simulaiton analysis
   * Data files 
   * Results files
* [System and software requirements](#system-and-software-requirements)
* [Copyright and licensing information](#copyright-and-licensing-information)
* [Contact information](#contact-information)

## CODE AND DATA FILES 

Three types of files (not including this readme file) are stored on this repository: code files, data files and figure files.


### Code files

All code is written in R. R is a statistical programming language and software package that is distributed under a GNU General Public License. R documentation and software is available for free download through the R Project for Statistical Computing website at http://www.r-project.org. The software is available as an executable file for a wide variety of operating systems and computer architectures, and a compilable binary is also available should you need it.

  ***sim_code*** -- contains R scripts to generate the data for age mixing patterns in HIV transmission networks
  
  ***post_sim_code*** -- contains R markdown scripts to analyse and visualise the simulation outputs, and an R script to compare metrics in between different data missingness scenarios


### Data files
  
  ***sim_outputs*** -- contains csv files one per each job submission at CHPC which contains statistics on age mixing patterns from records of the transmission networks and phylogenetic analysis.
  
  
### Results files

  ***MCAR*** -- results of age mixing pattenrs when we consider scenario of missing completly at random
  
  ***MAR_a*** -- results of age mixing pattenrs when we consider scenario of missing at random with at most 70% of sample being women
  
  ***MAR_b*** -- results of age mixing pattenrs when we consider scenario of missing at random with at most 30% of sample being women
  
  ***MAR_c*** -- results of age mixing pattenrs when we consider scenario of missing at random with at most 50% of sample being women
  
  ***MCAR_MAR_comparison*** -- results of comparison of results between MCAR and MAR scenarios

  

## SYSTEM AND SOFTWARE REQUIREMENTS

### Operating system


  We run the simulation study at the Cape Town Centre for High Performance Computing (CHPC) and tested it on personal computers (OS X Version 10.11.6 and Linux Ubuntu Version 16.04) and the golett cluster of the Flemish Supercomputer Centre (VSC).

### Required software

  **R version 3.4.4** <www.r-project.org> For statistical computing. To Install R, do the following:
  
  1. Open an internet browser and go to www.r-project.org.
  2.  Click the "download R" link in the middle of the page under "Getting Started."
  3. Select a CRAN location (a mirror site) and click the corresponding link.
  4. Click on the "Download R for ***your OS***" link at the top of the page.
  
  

  **Seq-Gen version 1.3.4** <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4> Simulates viral evolution across a transmission network and produces a sequence alignment. To install Seq-Gen, do the following:
  
  1. Visit the following Github repository to download the latest version of Seq-Gen: <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4>
  2. Click on the "Source Code" zip file to download
  3. Click on the zip file to unzip the folder
  4. Navigate to the source folder to confirm there is a file called "Makefile"
  5. Now you will need to compile the program using the Terminal on your computer
  6. Via the Terminal, change your working directory to the source folder by typing after the prompt: `cd "file/path/here/Seq-Gen-1.3.4 2/source"`
  7. Once your working directory has been set to the source folder, type after the prompt: `make`
  8. Now open the source folder and verify that a new file is present called "seq-gen"
  9. Copy that file and paste it into your R working directory
  
  If installing on HPC facility, you may follow the instructions from 1 up to 5. And you will load the tool via the the PBS file, for example `module add /apps/chpc/scripts/modules/bio/app/Seq-Gen/1.3.4`.
  

  **FastTree version 2.1.10** <http://www.microbesonline.org/fasttree/#Install> Reconstructs a phylogenetic tree from a large alignment dataset. To install FastTree, do the following:
  
  1. Visit the website for downloading instructions: <http://www.microbesonline.org/fasttree/#Install>
  2. If you have a Linux operating system, you can directly download the executable files that are linked on that website. Those downloaded files can then be placed in your R working directory
  3. If you are using an OS X operating system, open the link "FastTree.c" in a new browser window
  4. Right-click on the program and click "Save as"
  5. Save anywhere on your computer
  6. Open the Terminal on your computer and change your working directory to the folder that contains "FastTree.c". After the prompt type:  `cd "file/path/here"`
  7. After the directory has been changed, after the prompt type: `gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm`
  8. Now check to see if a new executable file has been created in that folder
  9. Copy that file and paste it into your R working directory

If installing on HPC facility, you may follow the instructions from 1 up to 7. And you will load the tool via the the PBS file, for example `module add /apps/chpc/scripts/modules/bio/app/FastTree/2.1.10`.


  **SimpactCyan version 0.21** and RSimpactCyan. SimpactCyan is the core program that allows fast simulation of HIV transmission across a sexual network. RSimpactCyan is the R package that enables initiation and running of models built by SimpactCyan. Installation instructions for both are at: <https://github.com/j0r1/RSimpactCyan/blob/master/INSTALLATION.md>


  **ClusterPicker version 1.2.3** <http://hiv.bio.ed.ac.uk/software.html> Cluster Picker identifies clusters in newick-formatted phyogenetic trees containing thousands of sequences. Cut-offs for within cluster genetic distance and bootstrap support are selected by the user.

  To use ClusterPicker, do the following:
  
  1. Install Java 1.6.0 or higherU
  2. Visit the website for downloading instructions for ClusterPicker: <http://hiv.bio.ed.ac.uk/software.html>
  3. Donwload the ClusterPicker command line version

    
  A long list of auxiliary R packages is required to run the post-simulation analysis for the simulation study.

  install.packages("devtools")
  
  install.packages("pacman")
  
  library(devtools)

  install_github("j0r1/readcsvcolumns/pkg")

  install_github("wdelva/RSimpactHelp", dependencies = TRUE)

  p_load(RSimpactCyan, RSimpactHelper, Rcpp, ape, expoTree, data.table, readr, phangorn, lme4, nlme, dplyr, adephylo, treedater, geiger, picante, igraph, phyloTop, phytools, Rsamtools, robustbase, intergraph, lubridate, tidyr)
  

To run the simulation on your desktop, you need to add the executable tools "Seq-Gen", and "FastTree", and "ClusterPicker" command line version to your working directory, as well as the root viral gene sequence (hiv.seq.C.pol.j.fasta). Run the the
`wrapper.age.mix.R` file for simulations. That file sources other needed files: `age.mix.MCAR.MAR.comput.R`, `age.mixing.MAR.fun.R`, `age.mixing.MCAR.fun.R`, `advanced.transmission.network.builder.R`, `needed.functions.RSimpactHelp.R`.

If you are running the simulation on HPC facility make sure you change the working directory accordingly, and ipload in your working directory also `ClusterPicker_1.2.3.jar` and execute the PBS file `run_large_AD_sacema.pbs`, you can rename is as you wnat and make sure also you modify the content according to your working environment on your HPC. 

Note that you need to verify if the working directory is set in all R files for both simulation plateforms (desktop or HPC), and also in the pbs file if you are on HPC.


## COPYRIGHT AND LICENSING INFORMATION

All files are copyright protected and are made available under the GPL 3.0 License <https://www.gnu.org/licenses/gpl-3.0.en.html>. This means that this work is suitable for commercial use, that licensees can modify the work, that they must release the source alongside with Derivative Work, and that Derivative Work must be released under the same terms.


## CONTACT INFORMATION

David Niyukuri
Email: <niyukuri@sun.ac.za>



