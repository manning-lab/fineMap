## paintor wdl
### preprocess 
# identify shared variants = 6.11.18
# LD calculation
# prepare assoc results = association file script (add alt ref)
# prepare annotation file = annotations.r
# output zcol names
# output ld names
# output annotation names


task preprocess {
	Array[String] interval
	File gds_file
	Array[File] sample_ids
	Array[File] assoc_files
	File annotation_file
	Array[String] anno_cols
	Int? mac
	String pval_col
	String effect_col

	Int memory
	Int disk

	command {
		R --vanilla --args ${sep = ":" interval} ${gds_file} ${sep="," sample_ids} ${sep="," assoc_files} ${annotation_file} ${sep="," anno_cols} ${default="10" mac} ${pval_col} ${effect_col} < /fineMap/paintor/preprocess.R
	}

	runtime {
		docker: "tmajarian/paintor:0.3"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		Array[File] ld_files = glob("Locus1.LD.*")
		File variant_list = "Locus1.markers.csv"
		File annotation_out = "Locus1.annotations"
		File assoc_out = "Locus1"
		File zcol_names = "zcol.txt"
		File ld_names = "ld.txt"
		File anno_names = "anno.txt"
	}
}

task runPaintor {
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
		echo "Locus1" >> input.txt && \
		PAINTOR -input input.txt \
		-in . \
		-out . \
		-Zhead ${sep = ", " zcol} \
		-LDname ${sep = ", " ld} \
		-annotations ${sep = ", " anno} \
		-mcmc \
		-max_causal ${max_causal}
	}

	runtime {
		docker: "tmajarian/paintor:0.3"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}
	
	output {
		File results = "Locus1.results"
	}
}

task summary {
	File paintor_results
	Array[File] ld_files
	File annotation_out
	File assoc_out
	File zcol_names
	File ld_names
	File anno_names

	Int memory
	Int disk

	Array[String] zcol = read_lines(zcol_names)
	Array[String] ld = read_lines(ld_names)
	Array[String] anno = read_lines(anno_names)
	# R --vanilla --args ${paintor_results} ${sep="," ld_files} ${annotation_out} ${assoc_out} ${sep="," zcol} ${sep="," ld} ${sep="," anno} < /fineMap/paintor/summary.R

	# CANVIS.py -l output/Locus1.results -z ZSCORE.EU -a RunDirectory/Locus1.annotations -s 6_Weak_transcription 10_Active_enhancer_2 1_Active_TSS 2_Weak_TSS 9_Active_enhancer_1

	command {
		touch Locus1.plots.png
	}

	runtime {
		docker: "tmajarian/paintor:0.3"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}
	
	output {
		File plots = "Locus1.plots.png"
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

	Int pre_memory
	Int paintor_memory
	Int summary_memory
	Int this_disk

	Array[Array[String]] these_intervals = read_tsv(this_interval_file)

	scatter( this_interval in these_intervals ) {

		Int this_chr = sub(this_interval[0], "chr", "")
		Int this_chr_int = this_chr - 1

		call preprocess {
			input: interval = this_interval, gds_file = these_gds_files[this_chr_int], sample_ids = these_sample_ids, assoc_files = these_assoc_files, annotation_file = this_annotation_file, anno_cols = these_anno_cols, mac = this_mac, pval_col = this_pval_col, effect_col = this_effect_col, memory = pre_memory, disk = this_disk
		}

		call runPaintor {
			input: ld_files = preprocess.ld_files, annotation_out = preprocess.annotation_out, assoc_out = preprocess.assoc_out, zcol_names = preprocess.zcol_names, ld_names = preprocess.ld_names, anno_names = preprocess.anno_names, max_causal = this_max_causal, memory = paintor_memory, disk = this_disk
		}

		call summary {
			input: paintor_results = runPaintor.results, ld_files = preprocess.ld_files, annotation_out = preprocess.annotation_out, assoc_out = preprocess.assoc_out, zcol_names = preprocess.zcol_names, ld_names = preprocess.ld_names, anno_names = preprocess.anno_names, memory = summary_memory, disk = this_disk
		}
	}
}