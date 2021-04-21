rule options:
	params:
		threads = 4
options = rules.options.params


# Define file names
rule files:
	params:
		sequence_dataset = "input_files/gisaid_hcov-19.fasta",
		new_genomes = "input_files/genome_sequences.fasta",
		aligned = "input_files/aligned.fasta",
		metadata_gisaid = "input_files/metadata_nextstrain.tsv",
		corelab_metadata = "input_files/metadata_impacc.csv",
		sample_metadata = "input_files/impacc-virology-clin-sample.csv",
		patient_metadata = "input_files/impacc-virology-clin-individ.csv",
		batch_layout = "input_files/batch_layout.csv",
		reference = "config/reference.gb",
		refgenome_size = "29903",
		max_missing = "30"


files = rules.files.params

rule lineages:
	message:
		"""
		Getting list of representative genomes per pangolin lineage:
		- Pick one representative genome per lineage
		- Export list with sequence names
		"""
	input:
		genomes = files.sequence_dataset,
		metadata = files.metadata_gisaid
	params:
		download_file = "yes",
		howmany = 1
	output:
		list = "output_files/sequences/rep_genomes.txt",
		base_dataset = "output_files/sequences/base_dataset.fasta",
		base_metadata = "output_files/metadata/base_metadata.tsv"
	shell:
		"""
		python3 scripts/get_lineage_reps.py \
			--download {params.download_file} \
			--howmany {params.howmany} \
			--output {output.list}
		
		echo Wuhan/Hu-1/2019 >> {output.list}
		echo Wuhan/WH01/2019 >> {output.list}
		
		python3 scripts/seqtree_handler.py \
			--input {input.genomes} \
			--format fasta \
			--action keep \
			--list {output.list} \
			--output {output.base_dataset}
		
		python3 scripts/seqtree_handler.py \
			--input {input.metadata} \
			--format tsv \
			--action keep \
			--index strain \
			--list {output.list} \
			--output {output.base_metadata}

		mv pangolin_lineages.csv output_files/pangolin_lineages.csv
		"""


rule filter_coverage:
	message:
		"""
		Filtering sequence files to:
		- Identify and remove poor quality genomes
		- Generate initial quality assurance matrix
		"""
	input:
		genomes = files.new_genomes
	params:
		size = files.refgenome_size,
		index = "sample_id",
		missing = files.max_missing
	output:
		matrix = "output_files/qa/qa_matrix1.tsv",
		rename = "output_files/sequences/rename.tsv",
		new_sequences = "output_files/sequences/renamed_genomes.fasta"
	shell:
		"""
		python3 scripts/filter_by_coverage.py \
			--genomes {input.genomes} \
			--index {params.index} \
			--refgenome-size {params.size} \
			--max-missing {params.missing} \
			--output {output.matrix} \
			--output2 {output.rename}
			
		
		python3 scripts/seqtree_handler.py \
			--input {input.genomes} \
			--format fasta \
			--action rename \
			--list {output.rename} \
			--output {output.new_sequences}	
		"""


rule inspect_metadata:
	message:
		"""
		Inspecting metadata file to:
		- Detect genomes with missing metadata
		- Expand quality assurance matrix
		"""
	input:
		metadata1 = files.corelab_metadata,
		metadata2 = files.sample_metadata,
		metadata3 = files.patient_metadata,
		batch = files.batch_layout,
		matrix1 = rules.filter_coverage.output.matrix
	params:
		index = "sample_id"
	output:
		sample_metadata = "output_files/assured_data/sample_metadata.tsv",
		matrix = "output_files/qa/qa_matrix2.tsv",
		rename = "output_files/sequences/rename2.tsv"
	shell:
		"""
		python3 scripts/inspect_metadata.py \
			--metadata1 {input.metadata1} \
			--metadata2 {input.metadata2} \
			--metadata3 {input.metadata3} \
			--batch {input.batch} \
			--index {params.index} \
			--matrix {input.matrix1} \
			--output1 {output.sample_metadata} \
			--output2 {output.matrix} \
			--output3 {output.rename}
		
		"""


rule multifasta:
	message:
		"""
		Combining sequence files as a multifasta file
		"""
	input:
		base_dataset = rules.lineages.output.base_dataset,
		new_genomes = "output_files/sequences/renamed_genomes.fasta",
		qamatrix = rules.inspect_metadata.output.matrix,
		rename_file = rules.inspect_metadata.output.rename
	params:
		format = "fasta"
	output:
		list_seqs = "output_files/sequences/full_genomes_list.txt",
		filtered_seqs = "output_files/sequences/filtered_seqs.fasta",
		ren_sequences = "output_files/sequences/renamed_seqs.fasta",
		combined_seqs = "output_files/sequences/sequences.fasta"
	shell:
		"""
		grep -v FAIL {input.qamatrix} | cut -d$'\t' -f 1 | sed -e 1d > {output.list_seqs}

		python3 scripts/seqtree_handler.py \
			--input {input.new_genomes} \
			--format {params.format} \
			--action keep \
			--list {output.list_seqs} \
			--output {output.filtered_seqs}
		
		python3 scripts/seqtree_handler.py \
			--input {output.filtered_seqs} \
			--format {params.format} \
			--action rename \
			--list {input.rename_file} \
			--output {output.ren_sequences}	

		cat {input.base_dataset} {output.ren_sequences} > {output.combined_seqs}
		"""


rule combine_metadata:
	message:
		"""
		Combining metadata files
		"""
	input:
		base_metadata = rules.lineages.output.base_metadata,
		sample_metadata = rules.inspect_metadata.output.sample_metadata
	output:
		combined_metadata = "output_files/assured_data/metadata.tsv"
	shell:
		"""
		python3 scripts/metadata_merger.py \
			--metadata1 {input.base_metadata} \
			--metadata2 {input.sample_metadata} \
			--output {output.combined_metadata}
		"""



## Aligning the sequences using MAFFT
rule align:
	message:
		"""
		Aligning sequences to {input.reference}
		    - gaps relative to reference are considered real
		"""
	input:
		sequences = rules.multifasta.output.combined_seqs,
		aligned = files.aligned,
		reference = files.reference
	params:
		threads = options.threads
	output:
		alignment = "output_files/sequences/aligned.fasta"
	shell:
		"""
		augur align \
			--sequences {input.sequences} \
			--existing-alignment {files.aligned} \
			--reference-sequence {input.reference} \
			--nthreads {params.threads} \
			--output {output.alignment} \
			--remove-reference \
			--debug
		"""


### Masking alignment sites
rule mask:
	message:
		"""
		Mask bases in alignment
		  - masking {params.mask_from_beginning} from beginning
		  - masking {params.mask_from_end} from end
		  - masking other sites: {params.mask_sites}
		"""
	input:
		alignment = rules.align.output.alignment
	params:
		mask_from_beginning = 55,
		mask_from_end = 100,
		mask_sites = "150 153 635 1707 1895 2091 2094 2198 2604 3145 3564 3639 3778 4050 5011 5257 5736 5743 5744 6167 6255 6869 8022 8026 8790 8827 8828 9039 10129 10239 11074 11083 11535 13402 13408 13476 13571 14277 15435 15922 16290 16887 19298 19299 19484 19548 20056 20123 20465 21550 21551 21575 22335 22516 22521 22661 22802 24389 24390 24622 24933 25202 25381 26549 27760 27761 27784 28253 28985 29037 29039 29425 29553 29827 29830"
	output:
		alignment = "output_files/tree/masked_alignment.fasta"
	shell:
		"""
		python3 scripts/mask-alignment.py \
			--alignment {input.alignment} \
			--mask-from-beginning {params.mask_from_beginning} \
			--mask-from-end {params.mask_from_end} \
			--mask-sites {params.mask_sites} \
			--output {output.alignment}
		"""


### Inferring Maximum Likelihood tree using the default software (IQTree)

rule tree:
	message: "Building tree"
	input:
		alignment = rules.mask.output.alignment
	params:
		threads = options.threads
	output:
		tree = "output_files/tree/tree_raw.nwk"
	shell:
		"""
		
		augur tree \
			--alignment {input.alignment} \
			--nthreads {params.threads} \
			--output {output.tree}
		"""


### Running TreeTime to estimate time for ancestral genomes

rule root2tip:
	message:
		"""
		Perform root-to-tip analysis:
		  - detect molecular clock outliers
		  - generate root-to-tip regression plot
		"""
	input:
		tree = rules.tree.output.tree,
		alignment = rules.align.output.alignment,
		metadata = rules.combine_metadata.output.combined_metadata
	params:
		clock_rate = 0.0008,
		clock_std_dev = 0.0004,
		root = "Wuhan/Hu-1/2019 Wuhan/WH01/2019",
		outdir = "./output_files/root2tip"
	output:
		outliers = "output_files/root2tip/outliers_list.txt"
	shell:
		"""
		treetime \
			--tree {input.tree} \
			--aln {input.alignment} \
			--dates {input.metadata} \
			--clock-filter 4 \
			--reroot {params.root} \
			--gtr JC69 \
			--clock-rate {params.clock_rate} \
			--clock-std-dev {params.clock_std_dev} \
			--max-iter 2 \
			--coalescent skyline \
			--plot-rtt root2tip_plot.pdf \
			--tip-labels \
			--verbose 1 \
			--outdir {params.outdir}
		
		grep "\-\-" {params.outdir}/dates.tsv | cut -d$'\t' -f 1 > {output.outliers}
		"""


rule assurance:
	message:
		"""
		Perform last step of quality assurance:
		  - generate final QA matrix
		  - move sequence and metadata files to final directory
		"""
	input:
		qamatrix2 = rules.inspect_metadata.output.matrix,
		outliers = rules.root2tip.output.outliers,
		lineage_seqs = rules.lineages.output.list,
		genomes = rules.multifasta.output.filtered_seqs,
		rename_file = rules.inspect_metadata.output.rename
	params:
		index = "sample_id",
		format = "fasta"
	output:
		qamatrix = "output_files/assured_data/qa_matrix.tsv",
		quality_seqs = "output_files/assured/quality_genomes.fasta",
		final_seqs = "output_files/assured_data/sequences.fasta",
	shell:
		"""
		python3 scripts/get_QAmatrix.py \
			--matrix {input.qamatrix2} \
			--genomes {input.genomes} \
			--outliers {input.outliers} \
			--index {params.index} \
			--output1 {output.qamatrix} \
			--output2 {output.quality_seqs}

		python3 scripts/seqtree_handler.py \
			--input {output.quality_seqs} \
			--format {params.format} \
			--action rename \
			--list {input.rename_file} \
			--output {output.final_seqs}	
		"""


### Clearing the working directory (only executed when needed)

rule clean:
	message: "Removing directories: {params}"
	params:
		"output_files"

	shell:
		"""
		rm -rfv {params}
		"""


