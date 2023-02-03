#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include {check_files; get_sites_only_1; combine_vcfs_1; bcftools_stats; parse_bcftools_stats; combine_csv } from './modules/qc'

workflow report{
    take:

    main:
        datasets = []
        params.datasets.each { dataset, description, dataset_vcf, dataset_sample ->
            datas = []
            dataset_vcfs = file(dataset_vcf)
            if (dataset_vcfs instanceof List){
                dataset_vcfs.each{ vcf ->
                    datas = [ dataset, description.split().join('-') ]
                    check_files([vcf, dataset_sample])
                    datas << vcf
                    datas << dataset_sample
                    datasets << datas 
                }
            }
            else{
                datas = [ dataset, description.split().join('-') ]
                check_files([dataset_vcf, dataset_sample])
                datas << dataset_vcf
                datas << dataset_sample
                datasets << datas 
            }
            
        }
        datasets_sites = []
        params.datasets_sites.each { dataset, description, dataset_vcf, dataset_sample ->
            datas = []
            dataset_vcfs = file(dataset_vcf)
            if (dataset_vcfs instanceof List){
                dataset_vcfs.each{ vcf ->
                    datas = [ dataset, description.split().join('-') ]
                    check_files([vcf, dataset_sample])
                    datas << vcf
                    datas << dataset_sample
                    datasets_sites << datas 
                }
            }
            else{
                datas = [ dataset, description.split().join('-') ]
                check_files([dataset_vcf, dataset_sample])
                datas << dataset_vcf
                datas << dataset_sample
                datasets_sites << datas 
            }
            
        }
        datasets_ch = Channel.from(datasets)
        datasets_sites_ch = Channel.from(datasets_sites)
        
        get_sites_only_1(datasets_ch.map{ dataset, description, vcf, sample -> [ dataset, description, file(vcf) ] })
        sites_datas = get_sites_only_1.out.groupTuple(by:[0,1])
            .mix(datasets_sites_ch.groupTuple(by:[0,1])
            .map{ dataset, desc, vcfs, samples -> [ dataset, desc, vcfs ]})
        combine_vcfs_1(sites_datas)
        bcftools_stats(combine_vcfs_1.out)
        parse_bcftools_stats(bcftools_stats.out)
        all_stats = parse_bcftools_stats.out.map{ dataset, step, stats -> [ 'H3A_report', dataset, step, stats ] }
            .groupTuple(by:0)
            .map{ group, datasets, steps, stats -> [ group, stats, 'csv' ] }
        combine_csv(all_stats)

    emit:
        datasets_ch
}


workflow{
    report()
}