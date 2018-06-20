## paintor wdl

task preprocess {
	Array[String] interval_array
	String interval_string
	File gds_file
	Array[File] sample_ids
	Array[File] assoc_files
	File annotation_file
	Array[String] anno_cols
	Int? mac
	String pval_col
	String effect_col
	Float? pval_thresh

	Int memory
	Int disk

	String interval = sub(interval_string, "\t", ".")

	command {
		echo ${interval}
		R --vanilla --args ${sep=":" interval_array} ${gds_file} ${sep="," sample_ids} ${sep="," assoc_files} ${annotation_file} ${sep="," anno_cols} ${default="10" mac} ${pval_col} ${effect_col} ${default="0.01" pval_thresh} < /fineMap/paintor/preprocess.R
	}

	runtime {
		docker: "manninglab/finemap:paintor.v3.0"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		Array[File] ld_files = glob("${interval}.LD.*")
		File ld_avg = "${interval}.all.LD"
		Array[File] ld_passed = glob("pval.passed.${interval}.LD.*")
		File variant_list = "${interval}.markers.csv"
		File variant_list_passed = "pval.passed.${interval}.markers.csv"
		File annotation_out = "${interval}.annotations"
		File annotation_out_passed = "pval.passed.${interval}.annotations"
		File assoc_out = "${interval}"
		File assoc_out_passed = "pval.passed.${interval}"
		File zcol_names = "zcol.txt"
		File ld_names = "ld.txt"
		File anno_names = "anno.txt"
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

	Int memory
	Int disk

	Array[String] zcol = read_lines(zcol_names)
	Array[String] ld = read_lines(ld_names)
	Array[String] anno = read_lines(anno_names)

	command {
		mv -t ./ ${annotation_out} ${assoc_out} ${sep = " " ld_files} && \
		echo ${interval} >> input.txt && \
		PAINTOR -input input.txt \
		-in . \
		-out . \
		-Gname "${interval}.enrichment.estimate" \
		-Lname "${interval}.log.bayesfactor" \
		-Zhead ${sep = "," zcol} \
		-LDname ${sep = "," ld} \
		-annotations ${sep = "," anno} \
		-max_causal ${max_causal} \
		-set_seed 1
	}

	runtime {
		docker: "manninglab/finemap:paintor.v3.0"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}
	
	output {
		File results = "${interval}.results"
		File enrichment = "${interval}.enrichment.estimate"
		File bayes = "${interval}.enrichment.estimate"		
	}
}

task summary {
	String interval_string
	File paintor_results
	File annotation_out
	File anno_names
	File ld_avg

	Int memory
	Int disk

	Array[String] anno = read_lines(anno_names)
	String interval = sub(interval_string, "\t", ".")

	command {
		python /fineMap/paintor/CANVIS.py \
		-l ${paintor_results} \
		-z meta_p \
		-a ${annotation_out} \
		-s ${sep=" " anno} \
		-r ${ld_avg} \
		-o ${interval} \
		-t 99 \
		-p \
		-L y
	}

	runtime {
		docker: "manninglab/finemap:paintor.canvis"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}
	
	output {
		File html = "${interval}.html"
		File svg = "${interval}.svg"
	}

}

task catResults {
	Array[File] all_results

	command {
		R --vanilla --args ${sep="," all_results} < /fineMap/paintor/catResults.R
	}

	runtime {
		docker: "manninglab/finemap:paintor.v3.0"
		disks: "local-disk 50 SSD"
		memory: "5G"
	}

	output {
		File cat_results = "top.variants.csv"
	}

}


workflow group_assoc_wf {
	
	File this_interval_file
	Array[File] these_gds_files
	Array[File] these_sample_ids
	Array[File] these_assoc_files
	File this_annotation_file
	Array[String] these_anno_cols
	Int? this_mac
	String this_pval_col
	String this_effect_col
	Int this_max_causal
	Float? this_pval_thresh

	Int pre_memory
	Int paintor_memory
	Int summary_memory
	Int this_disk

	Array[Array[String]] these_intervals = read_tsv(this_interval_file)
	Array[String] these_interval_lines = read_lines(this_interval_file)

	Array[Pair[Array[String],String]] these_interval_pairs = zip(these_intervals, these_interval_lines)

	scatter( this_interval_pair in these_interval_pairs ) {

		Int this_chr = sub(this_interval_pair.left[0], "chr", "")
		Int this_chr_int = this_chr - 1

		call preprocess {
			input: interval_array = this_interval_pair.left, interval_string = this_interval_pair.right, gds_file = these_gds_files[this_chr_int], sample_ids = these_sample_ids, assoc_files = these_assoc_files, annotation_file = this_annotation_file, anno_cols = these_anno_cols, mac = this_mac, pval_col = this_pval_col, effect_col = this_effect_col, pval_thresh = this_pval_thresh, memory = pre_memory, disk = this_disk
		}

		call runPaintor {
			input: interval_string = this_interval_pair.right, ld_files = preprocess.ld_files, annotation_out = preprocess.annotation_out, assoc_out = preprocess.assoc_out, zcol_names = preprocess.zcol_names, ld_names = preprocess.ld_names, anno_names = preprocess.anno_names, max_causal = this_max_causal, memory = paintor_memory, disk = this_disk
		}

		call summary {
			input: interval_string = this_interval_pair.right, paintor_results = runPaintor.results, annotation_out = preprocess.annotation_out, anno_names = preprocess.anno_names, ld_avg = preprocess.ld_avg, memory = summary_memory, disk = this_disk
		}
	}

	call catResults {
		input: all_results = runPaintor.results
	}
}