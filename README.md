# GenIn

GenIn is an open-source tool based on neural nets to rapidly and automatically assign genotypes to clade 2.3.4.4b H5 viruses collected in Europe since October 2020, starting from complete nucleotide genome sequences. Genotypes are assigned using the methods described [here](https://doi.org/10.1093%2Fve%2Fveae027).

The updated list of all genotypes recognized by GenIn is available in the [latest version of file "genotype2composition.txt" available in the "data" folder](https://github.com/izsvenezie-virology/genin/blob/main/data/version1/genotype2composition.txt).


## Installation

GenIn is available for Linux and requires R, nnet library, perl, mafft as dependacies. The easiest way to install it is to clone current github reposititory with the following command
```
git clone https://github.com/izsvenezie-virology/genin
```
and to create a conda environment using the yml file provided
```
conda env create -p genin_env -f environment.yml
```

## Usage

### Input
The tool requires a FASTA file as input, with each sequence formatted as described below.

GenIn analyzes multiple A(H5Nx x=1,2,3,4,5,8) Influenza virus sequences simultaneously. It can handle partial and complete genome sequences of multiple samples, but it needs all the eight segments to assess the genotype correctly. You must provide a single file containing all the nucleotide sequences in FASTA format. Sequences must adhere to the [IUPAC code](https://www.bioinformatics.org/sms/iupac.html). GenIn relies on the FASTA header to assign the sequence to a specific segment and sample. For this reason, the header must contain both a sample ID (consistent among sequences of the same sample) and one of the following segment names: `PB2`, `PB1`, `PA`, `HA`, `NP`, `NA`, `MP`, `NS`, separated by an underscore (i.e. ">sequenza_name_PA"). Genin assumes to analyze only H5Nx clade 2.3.4.4b sequences, for this reason the HA segment will be ignored. The segments belonging to the same sample must have the exact same sample ID (remember that it is case-sensitive).

### Basic usage
With the conda environment active, you can use GenIn to predict the genotypes of your samples running:
```
genin.pl --fasta YOUR_FASTA.FA --output OUTPUT.FOLDER
```
where
* YOUR_FASTA.FA is the fasta file with all the samples you want to predict genotype
* OUTPUT.FOLDER is directory where computation will be done and results will be saved in “result.txt” file

In addition, to speed up the computation, you can add the following option:
```
--Ncpu NUMBER_OF_CPU
```
to set up the number of cpus to be used by GenIn.

### Update database
GenIn is constantly fueled with new identified genotypes, thus you should always use the latest version of our database; to fullful such requirements, simply re-download this current github reposititory in the same directory with the command
```
git clone https://github.com/izsvenezie-virology/genin
```

## Outputs
In the output folder, GenIn saves a file called "results.txt" (csv format), which can be properly opened and viewed in Excel. Such file contains eighteen columns, a line as header with the information about the content of the column and a line for each analysed sample. The first column contains the Sample Name; the second contains the predicted genotype. If GenIn identifies the genotype correctly, the genotype is reported by one or two capital letters (e.g. AB), otherwise, an asterisk (\*) after the genotype is reported and details regarding the uncertainty of the assignment are printed in the last column called “Note”. The third column contains the probability associated with the given genotype.
    
The other columns (from the 4th to the 17th) contain a number associated with a reference genome for each segment (you can find it in [10.1093/ve/veae027](https://doi.org/10.1093%2Fve%2Fveae027)), and the next column contains the probability associated with the prediction of the single segment. Segment HA is ignored.

## Cite GenIn
We are currently writing the paper. 
Until the publication please cite the GitHub repository:

[https://github.com/izsvenezie-virology/genin](https://github.com/izsvenezie-virology/genin)

## License
GenIn is licensed under the GNU Affero v3 license (see [LICENSE](LICENSE)).

# Fundings

This work was supported by KAPPA-FLU HORIZON-CL6-2022-FARM2FORK-02-03 (grant agreement No 101084171) and by the NextGeneration EU-MUR PNRR Extended Partnership initiative on Emerging Infectious Diseases (Project no. PE00000007, INF-ACT).

>Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Health and Digital Executive Agency (HEDEA). 
>Neither the European Union nor the granting authority can be held responsible for them
