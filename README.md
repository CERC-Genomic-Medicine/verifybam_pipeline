# Pipeline to verify BAM/CRAM files from WES/WGS

## Requirements

1. Download and compile [VerifyBamId2](https://github.com/Griffan/VerifyBamID).
2. Verify that you are running Python 3 and the following packages are installed: `argparse`, `pandas`, `matplotlib`, `fpdf`.
3. Verify that you have `samtools` installed.

## Running

1. Clone this repository to the directory where you will run the pipeline:
   ```
   git clone https://github.com/CERC-Genomic-Medicine/verifybam_pipeline.git
   ```

2. Modify `nextflow.config` configuration file.
     * `params.bams` -- path to your BAM/CRAM file(s). You can use `glob` expressions to selecect multiple files.
     * `params.targets_bed` -- path to the WES target BED file. Used only when analyzing WES BAM/CRAM.
     * `params.fasta` -- path to the genome reference FASTA file.
     * `params.genome_build` -- set to "b37" or "b38".
     * `params.verifyBamId_dir` -- path to the `VerifyBamID` directory.
     * `params.reported_sex` -- path to the file with the reported sex. It is a tab-delimied file with two columns named: "NAME", "SEX". The "NAME" column stores the sample name/id. The sample name/id must match the BAM/CRAM filename without suffixes. For example, if the BAM file name is "Sample123.run2.phase1.bam", then the sample will be referred to as "Sample123". The "SEX" column stores self-reported sex values encoded as: F - female, M - male, NA - missing.
     * `executor.$slurm.queueSize` -- maximal number of SLURM jobs to submit at once.
  
3. Run pipeline.
   
   WGS:
   ```
   module load nextflow
   module load samtools
   nextflow run Verify.nf -w ~/scratch/work_directory
   ```
   WES:
   ```
   module load nextflow
   module load samtools
   nextflow run Verify_wes.nf -w ~/scratch/work_directory
   ```
   Important: when working on Compute Canada HPC, set working directory to ~/scratch/\<new directory name\>. This will speed up IO and also save space on your `project` partition. After the execution, if there were no errors and you are happy with the results, you can remove this working directory.

## Results
Results will saved in the `results` directory:
  * `results/report.pdf` -- PDF file with all plots.
  * `results/summary.txt` -- Space-delimited table with all the results, which was used to generate PDF.
  * `results/plots/` -- All plots in JPEG format from PDF.
  * `results/verifyBamId_logs/` -- All log files from the verifyBamId tool.
  * `results/dp/` -- Files with averagte DP values by chromosome for each sample. These files were used to generate `summary.txt`.

Below you can find description of the `summary.txt` columns:
 | Column | Description |
 | --- | --- |
 | NAME | Sample name, which is equal to the corresponding BAM/CRAM name without a suffix. |
 | HGDP_SNPS_OVERLAP | number of SNPs used by VerifyBamId for estimating contamination using the HGDP reference panel. This columns is not present when analyzing WES data. |
 | HGDP_SNPS | number of SNPs within the HGDP reference panel. This columns is not present when analyzing WES data. |
 | HGDP_AVG_DP | average read depth at SNPs used by VerifyBamId for estimating contamination using the HGDP reference panel. This columns is not present when analyzing WES data. |
 | HGDP_FREEMIX | VerifyBamId estimated contamination using the HGDP reference panel. This columns is not present when analyzing WES data. |
 | HGDP_PC1_CONTAMINATING | Contaminating sample's PC1 projected to the HGDP PCA space. This columns is not present when analyzing WES data. |
 | HGDP_PC1_INTENDED | Intended sample's PC1 projected to the HGDP PCA space. This columns is not present when analyzing WES data. |
 | HGDP_PC2_CONTAMINATING | Contaminating sample's PC2 projected to the HGDP PCA space. This columns is not present when analyzing WES data. |
 | HGDP_PC2_INTENDED | Intended sample's PC2 projected to the HGDP PCA space. This columns is not present when analyzing WES data. |
 | HGDP_PC3_CONTAMINATING | Contaminating sample's PC3 projected to the HGDP PCA space. This columns is not present when analyzing WES data. |
 | HGDP_PC3_INTENDED | Intended sample's PC3 projected to the HGDP PCA space. This columns is not present when analyzing WES data. |
 | HGDP_PC4_CONTAMINATING | Contaminating sample's PC4 projected to the HGDP PCA space. This columns is not present when analyzing WES data. |
 | HGDP_PC4_INTENDED | Intended sample's PC4 projected to the HGDP PCA space. This columns is not present when analyzing WES data. |
 | 1000G_SNPS_OVERLAP | number of SNPs used by VerifyBamId for estimating contamination using the 1000G reference panel. This columns is not present when analyzing WES data. |
 | 1000G_SNPS | number of SNPs within the 1000G reference panel. |
 | 1000G_AVG_DP | average read depth at SNPs used by VerifyBamId for estimating contamination using the 1000G reference panel. |
 | 1000G_FREEMIX | VerifyBamId estimated contamination using the 1000G reference panel. |
 | 1000G_PC1_CONTAMINATING | Contaminating sample's PC1 projected to the 1000G PCA space. |
 | 1000G_PC1_INTENDED | Intended sample's PC1 projected to the 1000G PCA space. |
 | 1000G_PC2_CONTAMINATING | Contaminating sample's PC2 projected to the 1000G PCA space. |
 | 1000G_PC2_INTENDED | Intended sample's PC2 projected to the 1000G PCA space. |
 | 1000G_PC3_CONTAMINATING | Contaminating sample's PC3 projected to the 1000G PCA space. |
 | 1000G_PC3_INTENDED | Intended sample's PC3 projected to the 1000G PCA space. |
 | 1000G_PC4_CONTAMINATING | Contaminating sample's PC4 projected to the 1000G PCA space. |
 | 1000G_PC4_INTENDED | Intended sample's PC4 projected to the 1000G PCA space. |
 | ALL_CHROMS | Total number of chromosomes. |
 | ALL_N_BASES | Total number of bases across all chromosomes. |
 | ALL_PRCNT_COVERED | Percent of covered bases across all chromosomes. |
 | ALL_PRCNT_DP1 | Percent of bases covered at DP>=1 across all chromosomes. |
 | ALL_PRCNT_DP2 | Percent of bases covered at DP>=2 across all chromosomes. |
 | ALL_PRCNT_DP3 | Percent of bases covered at DP>=3 across all chromosomes. |
 | ALL_PRCNT_DP4 | Percent of bases covered at DP>=4 across all chromosomes. |
 | ALL_PRCNT_DP5 | Percent of bases covered at DP>=5 across all chromosomes. |
 | ALL_PRCNT_DP10 | Percent of bases covered at DP>=10 across all chromosomes. |
 | ALL_PRCNT_DP20 | Percent of bases covered at DP>=20 across all chromosomes. |
 | ALL_PRCNT_DP30 | Percent of bases covered at DP>=30 across all chromosomes. |
 | ALL_PRCNT_DP40 | Percent of bases covered at DP>=40 across all chromosomes. |
 | ALL_PRCNT_DP50 | Percent of bases covered at DP>=50 across all chromosomes. |
 | ALL_PRCNT_DP60 | Percent of bases covered at DP>=60 across all chromosomes. |
 | ALL_PRCNT_DP70 | Percent of bases covered at DP>=70 across all chromosomes. |
 | ALL_PRCNT_DP80 | Percent of bases covered at DP>=80 across all chromosomes. |
 | ALL_PRCNT_DP90 | Percent of bases covered at DP>=90 across all chromosomes. |
 | ALL_PRCNT_DP100 | Percent of bases covered at DP>=100 across all chromosomes. |
 | ALL_PRCNT_DP110 | Percent of bases covered at DP>=110 across all chromosomes. |
 | ALL_AVG_DP | Average depth across all chromosomes. |
 | AUTO_CHROMS | Total number of autosomal chromosomes. |
 | AUTO_N_BASES | Total number of bases across autosomal chromosomes. |
 | AUTO_PRCNT_COVERED | Percent of covered bases across autosomal chromosomes. |
 | AUTO_PRCNT_DP1 | Percent of bases covered at DP>=1 across autosomal chromosomes. |
 | AUTO_PRCNT_DP2 | Percent of bases covered at DP>=2 across autosomal chromosomes. |
 | AUTO_PRCNT_DP3 | Percent of bases covered at DP>=3 across autosomal chromosomes. |
 | AUTO_PRCNT_DP4 | Percent of bases covered at DP>=4 across autosomal chromosomes. |
 | AUTO_PRCNT_DP5 | Percent of bases covered at DP>=5 across autosomal chromosomes. |
 | AUTO_PRCNT_DP10 | Percent of bases covered at DP>=10 across autosomal chromosomes. |
 | AUTO_PRCNT_DP20 | Percent of bases covered at DP>=20 across autosomal chromosomes. |
 | AUTO_PRCNT_DP30 | Percent of bases covered at DP>=30 across autosomal chromosomes. |
 | AUTO_PRCNT_DP40 | Percent of bases covered at DP>=40 across autosomal chromosomes. |
 | AUTO_PRCNT_DP50 | Percent of bases covered at DP>=50 across autosomal chromosomes. |
 | AUTO_PRCNT_DP60 | Percent of bases covered at DP>=60 across autosomal chromosomes. |
 | AUTO_PRCNT_DP70 | Percent of bases covered at DP>=70 across autosomal chromosomes. |
 | AUTO_PRCNT_DP80 | Percent of bases covered at DP>=80 across autosomal chromosomes. |
 | AUTO_PRCNT_DP90 | Percent of bases covered at DP>=90 across autosomal chromosomes. |
 | AUTO_PRCNT_DP100 | Percent of bases covered at DP>=100 across autosomal chromosomes. |
 | AUTO_PRCNT_DP110 | Percent of bases covered at DP>=110 across autosomal chromosomes. |
 | AUTO_AVG_DP | Average depth across autosomal chromosomes. |
 | X | 1 - X chromose is present, 0 - otherwise. |
 | X_AVG_DP | average read depth for X chromosome. |
 | X_NORM_AVG_DP | average read depth for X chromosome normalized by the average read depth across autosomal chromosomes. |
 | Y | 1 - Y chromosome is present, 0 - othewise. |
 | Y_AVG_DP | average read depth for Y chromosome. |
 | Y_NORM_AVG_DP | average read depth for Y chromosome normalized by the average read depth across autosomal chromosomes. |
