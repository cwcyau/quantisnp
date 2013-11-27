QuantiSNP v2.3 Beta


Last updated: 03 April 2011


Changes:



1. Introduction


QuantiSNP identifies putative copy number alterations from Illumina Infinium I/II SNP genotyping data. In version 2, QuantiSNP uses an improved data model with improved robustness and adaptability for use with Affymetrix 500K and SNP 6.0 data.


2. Installation


The following files are supplied:
	- Binary quantisnp (Linux) or quantisnp.exe (Windows).
	- Configuration files levels.dat and params.dat (levels-hd.dat is available for Infinium HD users).
	- Local GC content files for Build 35 and 36 of the human genome (b35.zip/b36.zip).
	- MATLAB Run-Time Libraries Installer (MCRInstallerXXX.bin - Linux, MCRInstallerXXX.exe - Windows).

The MATLAB Run-Time Libraries must be installed before QuantiSNP can be used.


3. Usage


A machine with at least 2 Gb memory is required and 4 Gb is recommended for processing SNP data from the Illumina 1M and Affymetrix 6.0 SNP arrays.

3.1 Input data


For single file processing, the input files must be plain tab-delimited text files with the following columns:
- Probe ID / SNP Name
- Chromosome
- Position
- Log R Ratio
- B Allele Frequency

Suitable text files can be obtained directly from Illumina's BeadStudio or GenomeStudio software or manually extracted from other file formats.

Affymetrix users can make use of the PennCNV Affymetrix extraction routines to generate suitable files (see http://www.openbioinformatics.org/penncnv/).


For batch file processing using BeadStudio reports, the report files must contain the following columns:

- Sample ID
- SNP Name
- Chromosome
- Position
- Log R Ratio
- B Allele Freq

In addition, a tab-delimited gender information file can be created with the following columns:

- Sample ID
- Gender (Male/Female)



3.2 Command line options


QuantiSNP maybe invoked with the command line arguments specified in Table \ref{tab:CommandLineOptions} from a Windows MS-DOS box or a Linux terminal.

	--outdir {directory name} (Required) Directory in which output files are stored (must exist).

	--levels {filename} (Required) Path to a configuration file containing list of copy number states and associated mean levels for the Log R Ratio. Default: levels.dat is assumed to be in the same directory as the executable.

	--config {filename} (Required) Path to configuration file containing list of hyperparameter settings. Default: params.dat is assumed to be in the same directory as the executable.

	--gcdir {directory name} (Optional) Directory in which local GC content files are stored. If not specified, then local GC-based correction of the Log R Ratio is not performed.

	--lsetting {number} (Optional) The characteristic length used to calculate transition probabilities. Default: 2,000,000.

	--emiters {number} (Optional) The number of iterations used for the EM algorithm during learning. Default: 10.
				
	--plot (Optional) Generates a series of plots (gzipped Postscript format) of putative copy number alterations found.

	--genotype (Optional) Generates a gzipped text file containing list of Generalised Genotypes.

	--chr {1-N} (Optional) Processes specified chromosomes only, e.g. --chr [1, 4:23] would run QuantiSNP using the data from chromosomes 1 and 4 to X only.

	--chrX {number} (Optional) Specifies the X chromosome, e.g. --chrX 23. Default: chrX is 23. This setting can be altered for analysis of non-human species.

	--doXcorrect Optional & Specifies whether to do correction of the Log R Ratio for the X chromosome. Default: No X correction.

	--isaffy (Optional) Specify for the processing of Affymetrix data. Default: Illumina data processing assumed.

[ for single file processing ]

	--sampleid {name} (Required) Sample ID - this is used to generate the name of the output file.

	--gender {male/female} (Required) Specifies gender of the sample. Adjusts processing for the X chromosome for males. If not specified, then automatic gender calling is used to predict gender.

	--input-file {filename} (Required) Path to text file containing input data.

[ for BeadStudio file processing ]

	--logfile {name} (Required) - writes a report of containing the samples in the BeadStudio report that were processed

	--genderfile {name} (Optional) - file containing sample IDs and gender information. If no gender file is supplied or a sample contains no gender information in the gender file then automatic gender calling is used.

	--beadstudio-files {name} (Required) - BeadStudio report file

3.3 Output


QuantiSNP generates up to four output files for each sample:

- A list of putative copy number alterations in {samplename}.cnv.
- A list of putative loss-of-heterozygosity regions in {samplename}.loh 
- A list of genotypes in {samplename}.gn.
- Plots of putative CNVs in {samplename}.ps.gz.
- Quality control parameters are stored in {samplename}.qc.


3.3.1 CNV file

The CNV output file is a plain text tab-delimited file containing a list of putative CNVs. The following columns are contained in the file:

	- Sample Name.
	- Chromosome.
	- Start Position (bp).
	- End Position (bp).
	- Start Probe Name. Probe name of first probe in CNV region.
	- End Probe name. Probe name of last probe in CNV region.
	- Length (bp).
	- Number of Probes.
	- Copy number.
	- Maximum Log Bayes Factor. Log Bayes Factor of most probable copy number state.
	- Log Bayes Factor (six columns). Log Bayes Factors for all copy number states.



3.3.2 LOH file

The LOH output file is a plain text tab-delimited file containing a list of putative copy-neutral LOH events. The following columns are contained in the file:

	- Sample Name.
	- Chromosome.
	- Start Position (bp).
	- End Position (bp).
	- Start Probe Name. Probe name of first probe in CNV region.
	- End Probe name. Probe name of last probe in CNV region.
	- Length (bp).
	- Number of Probes.
	- Copy number.
	- Maximum Log Bayes Factor. Log Bayes Factor of most probable copy number state.
	- Log Bayes Factor (six columns). Log Bayes Factors for all copy number states.


3.3.3 Genotype file

The genotype output file is a plain text tab-delimited file containing a list of putative generalised and diploid genotype calls. The following columns are contained in the file:

	- Probe Name.
	- Chromosome.
	- Position (bp).
	- Log R Ratio (corrected for local GC content)
	- B allele frequency.
	- Copy number (of most probable copy number state).
	- Maximum Log Bayes Factor (of most probable copy number state).
	- Generalised genotype call.
	- Generalised genotype call probability.
	- Diploid genotype call.
	- Diploid genotype call probability.



3.3.4 QC file

The quality control output file is a plain text tab-delimited file containing a list of model fit and quality control parameters. The following data are contained in the file for each chromosome:

	- Sample ID.
	- Outlier rate. Estimated probability of outliers in the data.
	- Std. Dev. LRR. A measure of the spread of Log R Ratio values.
	- Std. Dev. BAF. A measure of the spread of distribution of B allele frequencies for heterozygote genotypes.
	- Gender. Predicted gender of the sample if using automatic gender calling.

4. Analysis


For each copy number alteration identified a score is assigned to it. This measure is  called the Log Bayes Factor and is a quantity that represents the supports for the existence of a copy number alteration in the specified location given by the available SNP data.

A stringent threshold for the Log Bayes Factor, of at least 30, is recommended to obtain low false positive rates (< 1%). Thresholds between 10-30 increases power to detect smaller copy number alterations at the expense of an increase in false positive calls (up to 10%). Copy number alterations identified with a Log Bayes Factor of less than 10 are generally insignificant and it is recommended that they be filtered out.


5. Quality Control


Low quality data can be detected by examining the parameters in the quality control (QC) files. There are a number of key metrics:

Outlier rate. A high outlier rate suggests the noise model assumed by QuantiSNP does not match the actual noise characteristics of the data.

Std. Dev. (Log R Ratio/B allele frequency). A measure of the spread of the Log R Ratio and the heterozygote clusters of the B allele frequency model provide a measure of the noise in the data. Typical values associated with high quality data are 0.1-0.25 (LRR) and 0.025-0.04 (BAF). If we plot the two metrics for a large number of samples it is possible to identify groups of bad quality samples by clustering.



6. Additional notes


6.1 Local GC content based correction


Low quality Illumina SNP data often contains wave-like artefacts in the Log R Ratio. QuantiSNP uses local GC content to remove these wave-like artefacts in order to reduce false CNV calls.

6.2 Chromosome X processing


An excess number of CNV calls may be made on Chromosome X with data processed using an older version of Illumina's cluster file generation process. For the Infinium HD product line, the cluster file generation for the X chromosome has been improved. For further information:

http://www.illumina.com/downloads/XChrClustering_TN.pdf

QuantiSNP can centred the Log R Ratio values for the X chromosome based on the specific gender of the sample. Use the switch --doXcorrect to centre the Log R Ratio at zero for human females or the deletion level for human males.

6.3 Parameter settings


The parameter files levels.dat and params.dat contain important default model parameters. Modifications to these values could improve the functionality of QuantiSNP. For Infinium HD users a file levels-hd.dat is provided. These parameter files assume that the data is normalised using Illumina's own BeadStudio software, if you have used your own custom normalisation method or you are analysing Affymetrix data it is necessary to define your own levels.

6.3.1 Copy number levels

Each copy number state is associated with a mean value of the Log R Ratio. These are contained in the files:

levels.dat - older Illumina non-HD SNP arrays (HumanHapXXX)
levels-hd.dat - Illumina 6XX-Quad, Illumina 1M, Illumina-Omni
levels-affy.dat - Affymetrix SNP 6.0

6.3.2 Hyperparameters

The statistical model underlying QuantiSNP involves a number of preset hyperparameters. These are contained in the file params.dat.


