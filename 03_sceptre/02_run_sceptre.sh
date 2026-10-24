#!/bin/bash
# Run sceptre-pipeline (https://github.com/timothy-barry/sceptre-pipeline) on one SCEPTRE object.
# gRNAs are assigned with the mixture model; resource settings are those used for each run.
#
# Usage: bash 02_run_sceptre.sh <analysis>
#   cis:   pilot_full_day0 | full_day2 | full_day4 | full_day2_day4 | pilot_full
#          pilot_day0_ds | full_day0_ds | full_day2_ds | full_day4_ds
#   trans: pilot_full_day0_trans | pilot_full_trans | full_day2_day4_trans
#   individual gRNAs (objects from 04_results/02_individual_grnas.R): indv_cis | indv_trans

PROJECT_DIR="${PROJECT_DIR:-path/to/230126_BMD-STINGseq}"

# Resource settings shared by the larger runs
large=(--set_analysis_parameters_memory 8GB --prepare_assign_grnas_memory 8GB --combine_assign_grnas_memory 8GB
       --combine_assign_grnas_time 30m --prepare_association_analysis_time 30m --combine_association_analysis_time 30m
       --prepare_association_analysis_memory 8GB)
massive=(--grna_pod_size 25 --assign_grnas_time_per_grna 15m --run_association_analysis_time_per_pair 15m
         --set_analysis_parameters_time 60m --set_analysis_parameters_memory 8GB --prepare_assign_grnas_memory 8GB
         --assign_grnas_memory 32GB --combine_assign_grnas_memory 8GB --combine_assign_grnas_time 100m
         --prepare_association_analysis_time 60m --combine_association_analysis_time 60m
         --prepare_association_analysis_memory 12GB --run_association_analysis_memory 32GB
         --combine_association_analysis_memory 24GB --run_qc_time 100m --discovery_pairs trans)
small=(--grna_pod_size 25 --assign_grnas_time_per_grna 8m --run_association_analysis_time_per_pair 8m --set_analysis_parameters_time 18m)
pilot_full=(--grna_pod_size 25 --assign_grnas_time_per_grna 10m --run_association_analysis_time_per_pair 20m
            --set_analysis_parameters_time 18m --combine_association_analysis_time 20m --combine_assign_grnas_time 20m)

NXF_VERSION=""
case "$1" in
  pilot_full_day0)  RUN=240831_BMDSS_pilot_full_day0_cis_pca; JAVA="-Xms1G -Xmx240G"
                    ARGS=(--grna_pod_size 5 --assign_grnas_time_per_grna 12m --run_association_analysis_time_per_pair 8m --set_analysis_parameters_time 25m
                          --assign_grnas_memory 20GB --run_association_analysis_memory 20GB --combine_association_analysis_memory 8GB "${large[@]}") ;;
  full_day2)        RUN=240831_BMDSS_full_day2_cis_pca; JAVA="-Xms500M -Xmx100G"
                    ARGS=(--grna_pod_size 5 --assign_grnas_time_per_grna 10m --run_association_analysis_time_per_pair 8m --set_analysis_parameters_time 18m) ;;
  full_day4)        RUN=240831_BMDSS_full_day4_cis_pca; JAVA="-Xms500M -Xmx100G"
                    ARGS=(--grna_pod_size 5 --assign_grnas_time_per_grna 8m --run_association_analysis_time_per_pair 8m --set_analysis_parameters_time 18m) ;;
  full_day2_day4)   RUN=250210_BMDSS_full_day2_day4_cis_pca; JAVA="-Xms500M -Xmx10G"; ARGS=("${small[@]}") ;;
  pilot_full)       RUN=240831_BMDSS_pilot_full_cis_pca25; JAVA="-Xms500M -Xmx8G"; ARGS=("${pilot_full[@]}") ;;
  pilot_day0_ds|full_day0_ds|full_day2_ds|full_day4_ds)
                    RUN=250114_BMDSS_$1; JAVA="-Xms500M -Xmx10G"; ARGS=("${small[@]}") ;;
  pilot_full_day0_trans)
                    RUN=240831_BMDSS_pilot_full_day0_trans_pca; JAVA="-Xms500M -Xmx240G"
                    ARGS=(--grna_pod_size 5 --assign_grnas_time_per_grna 12m --run_association_analysis_time_per_pair 8m --set_analysis_parameters_time 25m
                          --assign_grnas_memory 20GB --run_association_analysis_memory 20GB --combine_association_analysis_memory 8GB "${large[@]}"
                          --discovery_pairs trans) ;;
  pilot_full_trans) RUN=250131_BMDSS_pilot_full_trans_pca_massive; JAVA="-Xms500M -Xmx16G"; ARGS=("${massive[@]}") ;;
  full_day2_day4_trans)
                    RUN=250210_BMDSS_full_day2_day4_trans_pca_massive; JAVA="-Xms500M -Xmx16G"; ARGS=("${massive[@]}") ;;
  indv_cis)         RUN=260529_BMDSS_pilot_full_cis_pca_IndvGuides; JAVA="-Xms500M -Xmx8G"; NXF_VERSION=23.10.1; ARGS=("${pilot_full[@]}") ;;
  indv_trans)       RUN=260603_BMDSS_pilot_full_trans_pca_IndvGuides; JAVA="-Xms500M -Xmx16G"; NXF_VERSION=23.10.1
                    ARGS=(--grna_pod_size 25 --assign_grnas_time_per_grna 12m --run_association_analysis_time_per_pair 8m --set_analysis_parameters_time 25m
                          --assign_grnas_memory 8GB --run_association_analysis_memory 24GB --combine_association_analysis_memory 24GB "${large[@]}"
                          --run_qc_time 2h --run_qc_memory 16GB --discovery_pairs trans) ;;
  *) echo "Unknown analysis: $1"; exit 1 ;;
esac

[[ -n "${NXF_VERSION}" ]] && export NXF_VER=${NXF_VERSION}
export NXF_OPTS="${JAVA}"

output_directory="${PROJECT_DIR}/sceptre/${RUN}/"
data_directory="${output_directory}ondisc/"

nextflow pull timothy-barry/sceptre-pipeline
nextflow run timothy-barry/sceptre-pipeline -r main -resume \
 --sceptre_object_fp ${data_directory}sceptre_object.rds \
 --response_odm_fp ${data_directory}response.odm \
 --grna_odm_fp ${data_directory}grna.odm \
 --output_directory ${output_directory} \
 --grna_assignment_method mixture \
 --pair_pod_size 100 \
 "${ARGS[@]}"
