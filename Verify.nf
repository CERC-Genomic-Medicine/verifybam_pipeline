

process verify_bams {
	cache "lenient"

	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return "retry" }
	maxRetries 3
	cpus 1

	input:
	tuple val(name), file(bam), file(bam_index) from Channel.fromPath(params.bams).map{ bam -> [ bam.getSimpleName(), bam, bam + (bam.getExtension() == "bam" ? ".bai" : ".crai") ] }
	each file(fasta) from Channel.fromPath(params.fasta)

	output:
	file "${name}.HGDP.Ancestry" into bam_hgdp_ancestry
	file "${name}.HGDP.selfSM" into bam_hgdp_contamination
	file "${name}.1000G.Ancestry" into bam_1000g_ancestry
	file "${name}.1000G.selfSM" into bam_1000g_contamination
	file "*.log" into verify_bam_logs

	publishDir "results/verifyBamId_logs/HGDP", pattern: "*.HGDP.log", mode: "copy"
	publishDir "results/verifyBamId_logs/1000G", pattern: "*.1000G.log", mode: "copy"


	"""
	${params.verifyBamId_dir}/bin/VerifyBamID --SVDPrefix ${params.verifyBamId_dir}/resource/hgdp.100k.${params.genome_build}.vcf.gz.dat --Reference ${fasta} --BamFile ${bam} --NumThread 1 --NumPC 4 --Output ${name}.HGDP 2>&1 | tee ${name}.HGDP.log
	${params.verifyBamId_dir}/bin/VerifyBamID --SVDPrefix ${params.verifyBamId_dir}/resource/1000g.phase3.100k.${params.genome_build}.vcf.gz.dat --Reference ${fasta} --BamFile ${bam} --NumThread 1 --NumPC 4 --Output ${name}.1000G 2>&1 | tee ${name}.1000G.log
	"""
}


process concat_hgdp_ancestry {

	executor "local"
	cpus 1

	input:
	file ancestry_files from bam_hgdp_ancestry.collect()

	output:
	file("hgdp_ancestry.txt") into hgdp_ancestry

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
	cpus 1

	input:
	file contamination_files from bam_hgdp_contamination.collect()

	output:
	file("hgdp_contamination.txt") into hgdp_contamination

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
	cpus 1

	input:
	file ancestry_files from bam_1000g_ancestry.collect()

	output:
	file("1000g_ancestry.txt") into _1000g_ancestry

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
	cpus 1

	input:
	file contamination_files from bam_1000g_contamination.collect()

	output:
	file("1000g_contamination.txt") into _1000g_contamination

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


process compute_bams_dp {
	cache "lenient"

        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return "retry" }
        maxRetries 3
        cpus 1

        input:
        tuple val(name), file(bam), file(bam_index) from Channel.fromPath(params.bams).map{ bam -> [ bam.getSimpleName(), bam, bam + (bam.getExtension() == "bam" ? ".bai" : ".crai") ] }

        output:
        file "${name}.dp.txt" into bam_dp
	file "*.by_chrom.txt" into bam_dp_by_chrom

	publishDir "results/dp/", pattern: "*.by_chrom.txt", mode: "copy"

        """
	samtools depth -a -s -q10 -Q10 ${bam} | aggregate_dp.py -o ${name}.dp
        """
}


process concat_bams_dp {

	executor "local"
	memory '4 GB'
	cpus 1

	input:
	file bam_dp_files from bam_dp.collect()

	output:
	file("dp.txt") into dp

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
	file dp_file from dp
	file contamination_hgdp_file from hgdp_contamination
	file ancestry_hgdp_file from hgdp_ancestry
	file contamination_1000g_file from _1000g_contamination
	file ancestry_1000g_file from _1000g_ancestry

	output:
	file "summary.txt" into summarized

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
	file summary from summarized
	file reported_sex from Channel.fromPath(params.reported_sex)

	output:
	file "*.jpeg" into plots
	file "report.pdf" into report

	publishDir "results/plots/", pattern: "*.jpeg", mode: "copy"
	publishDir "results/", pattern: "*.pdf", mode: "copy"

	"""
	generate_plots.py -s ${summary} -pca1 ${params.verifyBamId_dir}/resource/1000g.phase3.100k.${params.genome_build}.vcf.gz.dat.V -pop1 ${workflow.projectDir}/Populations/1000g_populations.txt -pca2 ${params.verifyBamId_dir}/resource/hgdp.100k.${params.genome_build}.vcf.gz.dat.V -pop2 ${workflow.projectDir}/Populations/HGDP_populations.txt -sex ${reported_sex}
	make_pdf.py -o report.pdf
	"""
}

