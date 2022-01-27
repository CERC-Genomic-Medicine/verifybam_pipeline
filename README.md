# Pipeline to verify BAM/CRAM files from WES/WGS

## Requirements

1. Download and compile [VerifyBamId2](https://github.com/Griffan/VerifyBamID).
2. Verify that you are running Python 3 and the following packages are installed: argparse, pandas, matplotlib, fpdf.
3. Verify thay you have samtools installed.

## Running

1. Clone this repository to the directory where you will run the pipeline:
   ```
   git clone https://github.com/CERC-Genomic-Medicine/verifybam_pipeline.git
   ```

2. Modify `nextflow.config` configuration file.
     * `params.bams` -- path to your BAM/CRAM file(s). You can use `glob` expressions to selecect multiple files.
     * `params.fasta` -- path to the genome reference FASTA file.
     * `params.genome_build` -- set to "b37" or "b38".
     * `params.verifyBamId_dir` -- path to the VerifyBamID directory.
     * `params.reported_sex` -- path to the file with the reported sex. It is a tab-delimied file with two columns named: "NAME", "SEX". The "NAME" column stores the sample name/id. The sample name/id must match the BAM/CRAM filename without suffixes. For example, if the BAM file name is "Sample123.run2.phase1.bam", then the sample will be referred to as "Sample123". The "SEX" columns stores self-reported sex values encoded as: F - female, M - male, NA - missing.
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
