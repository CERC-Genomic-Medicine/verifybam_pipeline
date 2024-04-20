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
   	memory "4GB"
   	time "4h"

	input:
	tuple val(name), file(bam), file(bam_index)
	path fasta

	output:
	path("${name}.HGDP.Ancestry")
	path("${name}.HGDP.selfSM")
	path("${name}.1000G.Ancestry")
	path("${name}.1000G.selfSM")
	path("*.log")

	publishDir "results/verifyBamId_logs/HGDP", pattern: "*.HGDP.log", mode: "copy"
	publishDir "results/verifyBamId_logs/1000G", pattern: "*.1000G.log", mode: "copy"


	"""
	${params.verifyBamId_dir}/bin/VerifyBamID --SVDPrefix ${params.verifyBamId_dir}/resource/hgdp.100k.${params.genome_build}.vcf.gz.dat --Reference ${fasta} --BamFile ${bam} --NumThread 1 --NumPC 4 --Output ${name}.HGDP 2>&1 | tee ${name}.HGDP.log
	${params.verifyBamId_dir}/bin/VerifyBamID --SVDPrefix ${params.verifyBamId_dir}/resource/1000g.phase3.100k.${params.genome_build}.vcf.gz.dat --Reference ${fasta} --BamFile ${bam} --NumThread 1 --NumPC 4 --Output ${name}.1000G 2>&1 | tee ${name}.1000G.log
	"""
}

process compute_bams_dp {
	cache "lenient"

        executor 'slurm'
        scratch '$SLURM_TMPDIR'

        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return "retry" }
        maxRetries 3
        cpus 1
	memory "64GB"
	time "12h"

        input:
        tuple val(name), path(bam), path(bam_index)
	path fasta

        output:
        path("${name}.dp.txt")
	path("*.by_chrom.txt")

	publishDir "results/dp/", pattern: "*.by_chrom.txt", mode: "copy"

        """
	samtools depth --reference ${fasta} -a -s -q10 -Q10 ${bam} | aggregate_dp.py -o ${name}.dp
        """
}


process concat_hgdp_ancestry {

	executor "local"
	memory '4 GB'
	cpus 1

	input:
	path ancestry_files

	output:
	path("hgdp_ancestry.txt")

	"""
	printf "NAME" > hgdp_ancestry.txt
	for i in {1..4}; do
		printf "\tHGDP_PC%d_CONTAMINATING\tHGDP_PC%d_INTENDED" \${i} \${i}
	done >> hgdp_ancestry.txt
	printf "\n" >> hgdp_ancestry.txt
	for f in `find . -name "*.Ancestry" -printf "%f\n" | sort`; do
		name=\${f%%.*}
		printf "%s" "\${name}"
		for i in {1..4}; do
			printf "\t%s" `grep "^\${i}" \${f} | cut -f2`
			printf "\t%s" `grep "^\${i}" \${f} | cut -f3`
		done
		printf "\n"
	done >> hgdp_ancestry.txt
	"""
}


process concat_hgdp_contamination {

	executor "local"
	memory '4 GB'
	cpus 1

	input:
	path contamination_files

	output:
	path("hgdp_contamination.txt")

	"""
	printf "NAME\tHGDP_SNPS\tHGDP_AVG_DP\tHGDP_FREEMIX\n" > hgdp_contamination.txt
	for f in `find . -name "*.selfSM" -printf "%f\n" | sort`; do
		name=\${f%%.*}
		printf "%s" "\${name}"
		printf "\t%s" `cut -f4 \${f} | tail -n1`
		printf "\t%s" `cut -f6 \${f} | tail -n1`
		printf "\t%s" `cut -f7 \${f} | tail -n1`
		printf "\n"
	done >> hgdp_contamination.txt
	"""
}


process concat_1000g_ancestry {

	executor "local"
	memory '4 GB'
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


process concat_1000g_contamination {

	executor "local"
	memory '4 GB'
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
	memory '4 GB'
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


process join {

	executor "local"
	memory '4 GB'
	cpus 1

	input:
	path dp_file
	path contamination_hgdp_file
	path ancestry_hgdp_file 
	path contamination_1000g_file
	path ancestry_1000g_file

	output:
	path("summary.txt")

	publishDir "results/", pattern: "*.txt", mode: "copy"

	"""
	join --check-order --header ${contamination_hgdp_file} ${ancestry_hgdp_file} > hgdp.txt
	join --check-order --header ${contamination_1000g_file} ${ancestry_1000g_file} > 1000g.txt
	join --check-order --header hgdp.txt 1000g.txt > contamination_ancestry.txt
	join --check-order --header contamination_ancestry.txt ${dp_file} > summary.txt
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
	generate_plots.py -s ${summary} -pca1 ${params.verifyBamId_dir}/resource/1000g.phase3.100k.${params.genome_build}.vcf.gz.dat.V -pop1 ${workflow.projectDir}/Populations/1000g_populations.txt -pca2 ${params.verifyBamId_dir}/resource/hgdp.100k.${params.genome_build}.vcf.gz.dat.V -pop2 ${workflow.projectDir}/Populations/HGDP_populations.txt -sex ${reported_sex}
	make_pdf.py -o report.pdf
	"""
}

workflow {
	bams = Channel.fromFilePairs("${params.bams}", size: -1) { file -> file.getName().replaceAll("(.bai|.crai)\$", "") }.map {it -> [it[0].replaceAll("(.bam|.cram)\$", ""), it[1][0], it[1][1]] }

	verify_estimates = verify_bams(bams, Channel.fromPath(params.fasta).collect())
	depths = compute_bams_dp(bams, Channel.fromPath(params.fasta).collect())

	all = join(
		concat_bams_dp(depths[0].collect()),
		concat_hgdp_contamination(verify_estimates[1].collect()),
		concat_hgdp_ancestry(verify_estimates[0].collect()),
		concat_1000g_contamination(verify_estimates[3].collect()),
 		concat_1000g_ancestry(verify_estimates[2].collect())
	) 
	plot(all, Channel.fromPath(params.reported_sex))
}

