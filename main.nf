#!/usr/bin/env nextflow

params.input = false
params.help = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    engine = new groovy.text.SimpleTemplateEngine()

    bindings = ["nb_threads":"$params.nb_threads",
                "para_diff": "$params.para_diff",
                "iso_diff": "$params.iso_diff",
                "perp_diff_min": "$params.perp_diff_min",
                "perp_diff_max": "$params.perp_diff_max",
                "lambda1": "$params.lambda1",
                "lambda2": "$params.lambda2",
                "output_dir":"$params.output_dir",
                "b_thr":"$params.b_thr"]

    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "FreeWater and FW corrected DTI metrics pipeline"
log.info "==============================================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "[Inputs]"
log.info "Input: $params.input"
log.info "Output directory: $params.output_dir"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[FreeWater fitting]"
log.info "Parallel diff: $params.para_diff"
log.info "Iso diff: $params.iso_diff"
log.info "Perpendicular diff min: $params.perp_diff_min"
log.info "Perpendicular diff max: $params.perp_diff_max"
log.info "Lambda 1: $params.lambda1"
log.info "Lambda 2: $params.lambda2"
log.info "b-threshold: $params.b_thr"
log.info ""

log.info "Number of processes per tasks"
log.info "============================="
log.info "FreeWater fitting: $params.nb_threads"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if (params.input){
    log.info "Input: $params.input"
    input = file(params.input)
    in_data = Channel
        .fromFilePairs("$input/**/*{brain_mask.nii.gz,bval,bvec,dwi.nii.gz}",
                       size: 4,
                       maxDepth:1,
                       flat: true) {it.parent.name}
}

(all_data_for_kernels, data_for_fw, grads_mask_for_metrics) = in_data
    .map{sid, brain_mask, bval, bvec, dwi -> [tuple(sid, brain_mask, bval, bvec, dwi),
                                                tuple(sid, brain_mask, bval, bvec, dwi),
                                                tuple(sid, brain_mask, bval, bvec)]}
    .separate(3)

all_data_for_kernels.first().set{unique_data_for_kernels}

process Compute_Kernel {
  cpus 1
  publishDir = "${params.output_dir}/Compute_Kernel"

  input:
    set sid, file(brain_mask), file(bval), file(bvec), file(dwi) from unique_data_for_kernels

  output:
    file("kernels/") into kernel_for_fw

  script:
    """
    scil_compute_freewater.py $dwi $bval $bvec\
      --mask $brain_mask\
      --para_diff $params.para_diff\
      --perp_diff_min $params.perp_diff_min\
      --perp_diff_max $params.perp_diff_max\
      --iso_diff $params.iso_diff\
      --lambda1 $params.lambda1\
      --processes $params.nb_threads\
      --lambda2 $params.lambda2\
      --b_thr $params.b_thr\
      --save_kernels kernels/ \
      --compute_only
    """
}

data_for_fw
  .combine(kernel_for_fw)
  .set{data_with_kernel_for_fw}

process Compute_FreeWater {
    cpus { params.nb_threads * task.attempt }
    memory { 6.GB * task.attempt }

    input:
    set sid, file(brain_mask), file(bval), file(bvec), file(dwi), file(kernels) from data_with_kernel_for_fw

    output:
    set sid, "${sid}__dwi_fw_corrected.nii.gz" into fw_corrected_dwi
    file "${sid}__FIT_dir.nii.gz"
    file "${sid}__FIT_FiberVolume.nii.gz"
    file "${sid}__FIT_FW.nii.gz"
    file "${sid}__FIT_nrmse.nii.gz"

    script:
    """
    scil_compute_freewater.py $dwi $bval $bvec\
        --mask $brain_mask\
        --para_diff $params.para_diff\
        --perp_diff_min $params.perp_diff_min\
        --perp_diff_max $params.perp_diff_max\
        --iso_diff $params.iso_diff\
        --lambda1 $params.lambda1\
        --processes $params.nb_threads\
        --lambda2 $params.lambda2\
        --load_kernels $kernels

    mv results/dwi_fw_corrected.nii.gz ${sid}__dwi_fw_corrected.nii.gz
    mv results/FIT_dir.nii.gz ${sid}__FIT_dir.nii.gz
    mv results/FIT_FiberVolume.nii.gz ${sid}__FIT_FiberVolume.nii.gz
    mv results/FIT_FW.nii.gz ${sid}__FIT_FW.nii.gz
    mv results/FIT_nrmse.nii.gz ${sid}__FIT_nrmse.nii.gz
    rm -rf results
    """
}

grads_mask_for_metrics
    .join(fw_corrected_dwi)
    .set{data_for_dti_metrics}

process FW_Corrected_Metrics {
    cpus 3

    input:
      set sid, file(brain_mask), file(bval), file(bvec), file(fw_corrected_dwi) from data_for_dti_metrics

    output:
    file "${sid}__fw_corr_ad.nii.gz"
    file "${sid}__fw_corr_evecs.nii.gz"
    file "${sid}__fw_corr_evecs_v1.nii.gz"
    file "${sid}__fw_corr_evecs_v2.nii.gz"
    file "${sid}__fw_corr_evecs_v3.nii.gz"
    file "${sid}__fw_corr_evals.nii.gz"
    file "${sid}__fw_corr_evals_e1.nii.gz"
    file "${sid}__fw_corr_evals_e2.nii.gz"
    file "${sid}__fw_corr_evals_e3.nii.gz"
    file "${sid}__fw_corr_fa.nii.gz"
    file "${sid}__fw_corr_ga.nii.gz"
    file "${sid}__fw_corr_rgb.nii.gz"
    file "${sid}__fw_corr_md.nii.gz"
    file "${sid}__fw_corr_mode.nii.gz"
    file "${sid}__fw_corr_norm.nii.gz"
    file "${sid}__fw_corr_rd.nii.gz"
    file "${sid}__fw_corr_tensor.nii.gz"
    file "${sid}__fw_corr_nonphysical.nii.gz"
    file "${sid}__fw_corr_pulsation_std_dwi.nii.gz"
    file "${sid}__fw_corr_residual.nii.gz"
    file "${sid}__fw_corr_residual_iqr_residuals.npy"
    file "${sid}__fw_corr_residual_mean_residuals.npy"
    file "${sid}__fw_corr_residual_q1_residuals.npy"
    file "${sid}__fw_corr_residual_q3_residuals.npy"
    file "${sid}__fw_corr_residual_residuals_stats.png"
    file "${sid}__fw_corr_residual_std_residuals.npy"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_image_math.py convert $brain_mask $brain_mask -f --data_type int16
    scil_compute_dti_metrics.py $fw_corrected_dwi $bval $bvec --mask $brain_mask\
        --ad ${sid}__fw_corr_ad.nii.gz --evecs ${sid}__fw_corr_evecs.nii.gz\
        --evals ${sid}__fw_corr_evals.nii.gz --fa ${sid}__fw_corr_fa.nii.gz\
        --ga ${sid}__fw_corr_ga.nii.gz --rgb ${sid}__fw_corr_rgb.nii.gz\
        --md ${sid}__fw_corr_md.nii.gz --mode ${sid}__fw_corr_mode.nii.gz\
        --norm ${sid}__fw_corr_norm.nii.gz --rd ${sid}__fw_corr_rd.nii.gz\
        --tensor ${sid}__fw_corr_tensor.nii.gz\
        --non-physical ${sid}__fw_corr_nonphysical.nii.gz\
        --pulsation ${sid}__fw_corr_pulsation.nii.gz\
        --residual ${sid}__fw_corr_residual.nii.gz\
        -f
    """
}
