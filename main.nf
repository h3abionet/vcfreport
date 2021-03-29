nextflow.enable.dsl=2

def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}

process count_markers {
    tag "count_markers_${dataset}_${chrm}"

    input:
      tuple val(dataset), val(chrm), file(dataset_vcf_set)
    output:
      tuple val(dataset), file(dataset_vcf_set), file(count)
    script:
      count = "${chrm}.txt"
      dataset_vcf = dataset_vcf_set[0]
      """
      echo -e "${chrm}: \$(bcftools index --nrecords ${dataset_vcf})" > ${count}
      """
}

process combine_counts {
    tag "combine_counts_${dataset}"
    input:
        tuple val(dataset), val(vcf_files), val(counts_files)
    output:
        tuple val(dataset), file(count_summary)
    script:
        count_summary = "${dataset}_summary.txt"
        """
        cat ${counts_files.join(" ")} | sort -V >> ${count_summary}
        TOTAL=\$(cut -f 2 -d' ' ${count_summary}| bc)
        echo -e "Total: ${TOTAL}" >> ${count_summary}
        """
}

def unppack = {
     (chromosome, files) -> chromosome
}

workflow {
    datasets = []

    params.globDatasets.each { dataset_name, dataset_glob ->
        println dataset_name
        pattern = dataset_glob + '{,.tbi}'
        vcfs_channel2 = Channel
            .fromFilePairs(pattern)
            .map { [dataset_name, it[0], it[1]] }
        count_markers(vcfs_channel2)
        }

    // Count per chrm
    //count_markers(datasets_cha)

    // Combine
    counts_cha = count_markers.out.groupTuple(by:0)
    combine_counts(counts_cha).view()
}
