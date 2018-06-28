## paintor wdl

task preprocess {
	Array[String] interval_array
	String interval_string
	File gds_file
	Array[File] sample_ids
	Array[File] assoc_files
	File? meta_file
	File? annotation_file
	String? anno_cols
	Int? mac
	String pval_col
	String effect_col
	String outpref

	Int memory
	Int disk

	String interval = sub(interval_string, "\t", ".")

	command {
		echo ${interval}
		R --vanilla --args ${sep=":" interval_array} ${gds_file} ${sep="," sample_ids} ${sep="," assoc_files} ${meta_file} ${default="NA" annotation_file} ${default="NA" anno_cols} ${default="40" mac} ${pval_col} ${effect_col} ${outpref} < /fineMap/paintor/preprocess.R
	}

	runtime {
		docker: "manninglab/finemap:paintor.v3.0"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		Array[File] ld_files = glob("${outpref}.${interval}.LD.*")
		File ld_avg = "${outpref}.${interval}.LD.all"
		File variant_list = "${outpref}.${interval}.markers.csv"
		File annotation_out = "${outpref}.${interval}.annotations"
		File assoc_out = "${outpref}.${interval}"
		File zcol_names = "zcol.txt"
		String meta_zcol = read_lines("meta_zcol.txt")
		File ld_names = "ld.txt"
		File anno_names = "anno.txt"
		File out_message = "out_message.txt"
	}
}

task runPaintor {
	String interval_string
	String interval = sub(interval_string, "\t", ".")
	Array[File] ld_files
	File annotation_out
	File assoc_out
	File zcol_names
	File ld_names
	File anno_names
	Int max_causal
	String outpref

	Int memory
	Int disk

	Array[String] zcol = read_lines(zcol_names)
	Array[String] ld = read_lines(ld_names)
	Array[String] anno = read_lines(anno_names)

	command {
		mv -t ./ ${annotation_out} ${assoc_out} ${sep = " " ld_files} && \
		echo ${outpref}.${interval} >> input.txt && \
		PAINTOR -input input.txt \
		-in . \
		-out . \
		-Gname "${outpref}.${interval}.enrichment.estimate" \
		-Lname "${outpref}.${interval}.log.bayesfactor" \
		-Zhead ${sep = "," zcol} \
		-LDname ${sep = "," ld} \
		-annotations ${sep="," anno} \
		-max_causal ${max_causal} \
		-set_seed 1
	}

	runtime {
		docker: "manninglab/finemap:paintor.v3.0"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}
	
	output {
		File results = "${outpref}.${interval}.results"
		File enrichment = "${outpref}.${interval}.enrichment.estimate"
		File bayes = "${outpref}.${interval}.enrichment.estimate"		
	}
}

task summaryLD {
	String outpref
	String interval_string
	File paintor_results
	File annotation_out
	File anno_names
	File zcol_names
	Array[File] ld_files
	Float? pval_thresh
	String meta_ld_file
	String meta_zcol

	Int memory
	Int disk

	String interval = sub(interval_string, "\t", ".")
	Array[String] anno = read_lines(anno_names)
	Array[String] zcol = read_lines(zcol_names)

	command {
		python /fineMap/paintor/CANVIS.py \
		-l ${paintor_results} \
		-v ${sep=" " zcol} ${meta_zcol} \
		-a ${annotation_out} \
		-c ${sep=" " anno} \
		-r ${sep=" " ld_files} ${meta_ld_file} \
		-o ${outpref}.${interval} \
		-t 99 \
		-T ${default="1" pval_thresh} \
		-L y && \
		inkscape ${outpref}.${interval}.svg -e ${outpref}.${interval}.png
	}

	runtime {
		docker: "manninglab/finemap:paintor.canvis"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}
	
	output {
		File html = "${outpref}.${interval}.html"
		File svg = "${outpref}.${interval}.svg"
		File pdf = "${outpref}.${interval}.png"
	}

}

task summaryNoLD {
	String outpref
	String interval_string
	File paintor_results
	File annotation_out
	File anno_names
	File zcol_names
	Float? pval_thresh
	String meta_zcol

	Int memory
	Int disk

	String interval = sub(interval_string, "\t", ".")
	Array[String] anno = read_lines(anno_names)
	Array[String] zcol = read_lines(zcol_names)

	command {
		python /fineMap/paintor/CANVIS.py \
		-l ${paintor_results} \
		-v ${sep=" " zcol} ${meta_zcol} \
		-a ${annotation_out} \
		-c ${sep=" " anno} \
		-o ${outpref}.${interval} \
		-t 99 \
		-T ${default="1" pval_thresh} && \
		inkscape ${outpref}.${interval}.svg -e ${outpref}.${interval}.png
	}

	runtime {
		docker: "manninglab/finemap:paintor.canvis"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}
	
	output {
		File html = "${outpref}.${interval}.html"
		File svg = "${outpref}.${interval}.svg"
		File pdf = "${outpref}.${interval}.png"
	}

}

task catResults {
	String outpref
	Array[File] all_results

	command {
		R --vanilla --args ${sep="," all_results} ${outpref} < /fineMap/paintor/catResults.R
	}

	runtime {
		docker: "manninglab/finemap:paintor.v3.0"
		disks: "local-disk 50 SSD"
		memory: "5G"
	}

	output {
		File cat_results = "${outpref}.top.variants.csv"
	}

}


workflow group_assoc_wf {
	# File this_gds_list
	# File this_sample_id_list
	# File this_assoc_list
	Array[File] these_gds
	Array[File] these_sample_ids
	Array[File] these_assoc
	File this_meta_file
	String this_outpref

	File this_interval_file
	File? this_annotation_file
	String? these_anno_cols
	Int? this_mac
	String this_pval_col
	String this_effect_col
	Int this_max_causal
	Boolean plot_ld
	Float? this_pval_thresh

	Int pre_memory
	Int paintor_memory
	Int summary_memory
	Int this_disk

	# read file lis
	# Array[File] these_gds = read_lines(this_gds_list)
	# Array[File] these_sample_ids = read_lines(this_sample_id_list)
	# Array[File] these_assoc = read_lines(this_assoc_list)

	Array[Array[String]] these_intervals = read_tsv(this_interval_file)
	Array[String] these_interval_lines = read_lines(this_interval_file)

	Array[Pair[Array[String],String]] these_interval_pairs = zip(these_intervals, these_interval_lines)

	scatter( this_interval_pair in these_interval_pairs ) {

		Int this_chr = sub(this_interval_pair.left[0], "chr", "")
		Int this_chr_int = this_chr - 1

		call preprocess {
			input: interval_array = this_interval_pair.left, interval_string = this_interval_pair.right, gds_file = these_gds[this_chr_int], sample_ids = these_sample_ids, assoc_files = these_assoc, meta_file = this_meta_file, annotation_file = this_annotation_file, anno_cols = these_anno_cols, mac = this_mac, pval_col = this_pval_col, effect_col = this_effect_col, outpref = this_outpref, memory = pre_memory, disk = this_disk
		}

		call runPaintor {
			input: interval_string = this_interval_pair.right, ld_files = preprocess.ld_files, annotation_out = preprocess.annotation_out, assoc_out = preprocess.assoc_out, zcol_names = preprocess.zcol_names, ld_names = preprocess.ld_names, anno_names = preprocess.anno_names, max_causal = this_max_causal, outpref = this_outpref, memory = paintor_memory, disk = this_disk
		}

		if (plot_ld) {
			call summaryLD {
				input: outpref = this_outpref, interval_string = this_interval_pair.right, paintor_results = runPaintor.results, annotation_out = preprocess.annotation_out, anno_names = preprocess.anno_names, zcol_names = preprocess.zcol_names, meta_zcol = preprocess.meta_zcol, ld_files = preprocess.ld_files, pval_thresh = this_pval_thresh, meta_ld_file = preprocess.ld_avg, memory = summary_memory, disk = this_disk
			}	
		}
		if (!plot_ld){
			call summaryNoLD {
				input: outpref = this_outpref, interval_string = this_interval_pair.right, paintor_results = runPaintor.results, annotation_out = preprocess.annotation_out, anno_names = preprocess.anno_names, zcol_names = preprocess.zcol_names, meta_zcol = preprocess.meta_zcol, pval_thresh = this_pval_thresh, memory = summary_memory, disk = this_disk
			}	
		}
		
	}

	call catResults {
		input: all_results = runPaintor.results, outpref = this_outpref
	}
}