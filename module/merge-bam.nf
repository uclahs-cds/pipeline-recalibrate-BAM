include { gatherbams } from './gather-bams.nf'
include { mergesamfiles } from './merge-sam-files.nf'

workflow merge_bams {
    take:
    bams_to_merge

    main:
    if (params.parallelize_by_chromosome) {
        gatherbams(bams_to_merge)
        output_ch_merge_bams = gatherbams.out.gathered_bams
    } else {
        mergesamfiles(bams_to_merge)
        output_ch_merge_bams = mergesamfiles.out.merged_bams
    }

    emit:
    output_ch_merge_bams = output_ch_merge_bams
}
