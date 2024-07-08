#!/usr/bin/env nextflow
/*
* AUTHOR: Daniel Taliun <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2024
*/


process verify_bams {
	cache "lenient"
	
	executor 'slurm'
	scratch '$SLURM_TMPDIR'

	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return "retry" }
	maxRetries 3
	cpus 1
   	memory "3GB"
   	time "2h"

	input:
	tuple val(name), file(bam), file(bam_index)
	path fasta

	output:
	path("${name}.1000G.Ancestry")
	path("${name}.1000G.selfSM")
	path("${name}.1000G.log")

	publishDir "results/verifyBamId_logs/HGDP", pattern: "*.HGDP.log", mode: "copy"
	publishDir "results/verifyBamId_logs/1000G", pattern: "*.1000G.log", mode: "copy"

	"""

${params.verifyBamId_dir}/bin/VerifyBamID --SVDPrefix ${params.verifyBamId_dir}/resource/exome/1000g.phase3.10k.${params.genome_build}.exome.vcf.gz.dat --Reference ${fasta} --BamFile ${bam} --NumThread 1 --DisableSanityCheck  --NumPC 4 --Output ${name}.1000G 2>&1 | tee ${name}.1000G.log


	"""
}

process compute_bams_dp {
	cache "lenient"
	scratch true

        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return "retry" }
        maxRetries 3
        cpus 1

        input:
        tuple val(name), path(bam), path(bam_index)
        each file(fasta)
	each file(targets)

        output:
        file "${name}.dp.txt"
	file "*.by_chrom.txt"

	publishDir "results/dp/", pattern: "*.by_chrom.txt", mode: "copy"

        """
	samtools depth -a -s -q20 -Q20 --reference ${fasta} -b ${targets} ${bam} | aggregate_dp.py -o ${name}.dp
        """
}

process concat_1000g_ancestry {

	executor "local"
	memory '1 GB'
	cpus 1

	input:
	path ancestry_files

	output:
	path("1000g_ancestry.txt")

	"""
	printf "NAME" > 1000g_ancestry.txt
	for i in {1..4}; do
		printf "\t1000G_PC%d_CONTAMINATING\t1000G_PC%d_INTENDED" \${i} \${i}
	done >> 1000g_ancestry.txt
	printf "\n" >> 1000g_ancestry.txt
	for f in `find . -name "*.Ancestry" -printf "%f\n" | sort`; do
		name=\${f%%.*}
		printf "%s" "\${name}"
		for i in {1..4}; do
			printf "\t%s" `grep "^\${i}" \${f} | cut -f2`
			printf "\t%s" `grep "^\${i}" \${f} | cut -f3`
		done
		printf "\n"
	done >> 1000g_ancestry.txt
	"""
}

process concat_1000g_log {
	
	executor "local"
        memory '1 GB'
        cpus 1

        input:
        path log_files

        output:
        path("1000g_log.txt")

	"""
	printf "NAME\t1000G_SNPS_OVERLAP\n" > 1000g_log.txt
	for f in `find . -name "*.1000G.log" -printf "%f\n" | sort`; do
		name=\${f%%.*}
		printf "%s" "\${name}"
		snps_overlap=\$(grep -Po "(?<=\\[SimplePileup\\] Total Number Markers: )[0-9]+"  \${f})
		printf "\t%s" "\${snps_overlap}"
		printf "\n"
	done >> 1000g_log.txt
	"""

}

process concat_1000g_contamination {

	executor "local"
	memory '1 GB'
	cpus 1

	input:
	path contamination_files

	output:
	path("1000g_contamination.txt")

	"""
	printf "NAME\t1000G_SNPS\t1000G_AVG_DP\t1000G_FREEMIX\n" > 1000g_contamination.txt
	for f in `find . -name "*.selfSM" -printf "%f\n" | sort`; do
                name=\${f%%.*}
                printf "%s" "\${name}"
                printf "\t%s" `cut -f4 \${f} | tail -n1`
                printf "\t%s" `cut -f6 \${f} | tail -n1`
                printf "\t%s" `cut -f7 \${f} | tail -n1`
                printf "\n"
	done >> 1000g_contamination.txt
	"""
}



process concat_bams_dp {

	executor "local"
	memory '1 GB'
	cpus 1

	input:
	path  bam_dp_files

	output:
	path("dp.txt")

	"""
	filename=`find . -name "*.dp.txt" | head -n1`
	printf "NAME\t%s\n" "`head -n1 \${filename}`" > dp.txt
	for f in `find . -name "*.dp.txt" -printf "%f\n" | sort`; do
		name=\${f%%.*}
		printf "%s\t%s\n" "\${name}" "`tail -n1 \${f}`" >> dp.txt
	done
	"""
}


process joined {

	executor "local"
	memory '1 GB'
	cpus 1

	input:
	path dp_file
	path contamination_1000g_file
	path log_1000g_file
	path ancestry_1000g_file

	output:
	path("summary.txt")

	publishDir "results/", pattern: "*.txt", mode: "copy"

	"""
	join --check-order --header ${contamination_1000g_file} ${log_1000g_file}  | join --check-order --header - ${ancestry_1000g_file} > 1000g.txt
	join --check-order --header 1000g.txt ${dp_file} > summary.txt
	"""
}


process plot {

	executor "local"
	memory '4 GB'
	cpus 1

	input:
	path summary
	path reported_sex

	output:
	path("*.jpeg")
	path("report.pdf")

	publishDir "results/plots/", pattern: "*.jpeg", mode: "copy"
	publishDir "results/", pattern: "*.pdf", mode: "copy"

    """
    generate_plots.py -s ${summary} -pca1 ${params.verifyBamId_dir}/resource/exome/1000g.phase3.10k.${params.genome_build}.exome.vcf.gz.dat.V -pop1 ${workflow.projectDir}/Populations/1000g_populations.txt -sex ${reported_sex}
    make_pdf.py -o report.pdf
    """
}

workflow {
	bams = Channel.fromFilePairs("${params.bams}", size: -1) { file -> file.getName().replaceAll("(.bai|.crai)\$", "") }.map {it -> [it[0].replaceAll("(.bam|.cram)\$", ""), it[1][0], it[1][1]] }

	verify_estimates = verify_bams(bams, Channel.fromPath(params.fasta).collect())
	depths = compute_bams_dp(bams, Channel.fromPath(params.fasta).collect(), Channel.fromPath(params.targets_bed))

	all = joined(
		concat_bams_dp(depths[0].collect()),
		concat_1000g_contamination(verify_estimates[1].collect()),
		concat_1000g_log(verify_estimates[2].collect()),
 		concat_1000g_ancestry(verify_estimates[0].collect())
	) 
	plot(all, Channel.fromPath(params.reported_sex))
}

