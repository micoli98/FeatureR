## De novo extraction of features for de novo extraction of signature

featureR_denovo <- function(df, output_path, seed, showPlots=F) {
  ############## Segmentation production
  segments <- extract_segments(df, verbose =T)
  write.table(segments, file.path(output_path, "segmentation.tsv"))
  ############## Extract SV features
  print("############## Starting SV features extraction ##############")
  ## Complex and single SVs
  cosim_feat <- cosim(df)

  ## Intra-chromosomal and inter-chromosomal SVs
  trater_feat <- trater(df)

  ## Complex and sparse SVs
  breakpoints <- bk_distance(df)
  compact_sparse_feat <- extract_compact_sparse(breakpoints, SEED)

  #Deletion and duplication lenghts
  deldups_long <- deldup(df)
  dels_feat <- extract_classes(deldups_long$deletions, "deletion", 6, SEED, plot=showPlots)
  dups_feat <- extract_classes(deldups_long$duplications, "duplication", 5, SEED, plot=showPlots)

  # SVs types: insertions, unbalanced translocation, reciprocal translocation, reciprocal inversions, LINE insertions
  sv_types <-SV_categories(df)

  # Magnitudes: quantiles 0.25 and 0.75 of the logR
  magnitudes_feat <- logR_quantiles(df, segments)

  # Chromotripsis using Shatterseek
  chrmtrps_feat <- is.chromotripsis(df, segments)

  ############## Extract CNV features
  print("############## Starting CNV features extraction ##############")
  ss <- segments %>%
    mutate(length = end - start) %>%
    filter(length>0) %>%
    filter(sample %in% samples$sample) %>%
    dplyr::select(sample, length)
  ss_feat <- extract_classes(ss, "segsize", 4, SEED, plot=showPlots)

  bp5MB_counts <- bkp5MB_extract(df, segments)
  bp5MB_feat <- bkp5MB_discretize(bp5MB_counts, SEED, plot=showPlots)

  oscil_counts<- oscil_extract(df, segments)
  oscil_feat <- extract_classes(oscil_counts, "oscillation", 4 ,SEED, plot=showPlots)

  bpArm_counts <- bpArm_extract(df, segments)
  bpArm_feat <- extract_classes(bpArm_counts, "arm_breakpoints", 6, SEED, plot=showPlots)

  changep_counts <- changp_extract(df, segments)
  changep_feat <- extract_classes(changep_counts, "changepoint", 4, SEED, plot=showPlots)

  ######## Aggregation in one features dataframe and one models list
  features <- Reduce(function(df1, df2) merge(df1, df2, by = "sample", all = TRUE),
                     list(cosim_feat_feat, trater_feat, compact_sparse_feat[[2]], dels_feat[[2]], dups_feat[[2]], sv_types,
                          magnitudes_feat, chrmtrps_feat, ss_feat[[2]], bp5MB_feat[[2]], oscil_feat[[2]],
                          bpArm_feat[[2]], changep_feat[[2]]))
  write.table(features, paste0(output_path, "features.tsv"), sep='\t', col.names = TRUE)

  models <- list(compact_sparse = compact_sparse_feat[[1]], deletions = dels_feat[[1]], duplications = dups_feat[[1]],
                 segsize = ss_feat[[1]], bp5MB = bp5MB_feat[[1]], oscill = oscil_feat[[1]], bpArm = bpArm_feat[[1]],
                 changepoint = changep_feat[[1]])
  saveRDS(models, paste0(output_path, "models.rds"))
  print(paste0("Features.tsv, segmentation.tsv and models.rds have been saved in ", output_path))
  return(c(features, models))
}

## Extraction of features from existing models

featureR_fromModels <- function(df, output_path, models) {
  ############## Segmentation production
  segments <- extract_segments(df, verbose =T)
  write.table(segments, file.path(output_path, "segmentation.tsv"))

  ############## Extract SV features
  print("############## Starting SV features extraction ##############")
  ## Complex and single SVs
  cosim_feat <- cosim(df)

  ## Intra-chromosomal and inter-chromosomal SVs
  trater_feat <- trater(df)

  ## Complex and sparse SVs
  breakpoints <- bk_distance(df)
  compact_sparse_feat <- MM_fromModel(breakpoints, "compact_sparse",models[["compact_sparse"]])

  #Deletion and duplication lenghts
  deldups_long <- deldup(df)
  dels_feat <- classes_fromModel(deldups_long$deletions, "deletion", models[["deletions"]])
  dups_feat <- classes_fromModel(deldups_long$duplications, "duplication", models[["duplications"]])

  # SVs types: insertions, unbalanced translocation, reciprocal translocation, reciprocal inversions, LINE insertions
  sv_types <-SV_categories(df)

  # Magnitudes: quantiles 0.25 and 0.75 of the logR
  magnitudes_feat <- logR_quantiles(df, segments)

  # Chromotripsis using Shatterseek
  chrmtrps_feat <- is.chromotripsis(df, segments)

  ############## Extract CNV features
  print("############## Starting CNV features extraction ##############")
  ss <- segments %>%
    mutate(length = end - start) %>%
    filter(length>0) %>%
    filter(sample %in% samples$sample) %>%
    dplyr::select(sample, length)
  ss_feat <- classes_fromModel(ss, "segsize", models[["segsize"]])

  bp5MB_counts <- bkp5MB_extract(df, segments)
  bp5MB_feat <- MM_fromModel(bp5MB_counts, "bp5MB", models[["bp5MB"]])

  oscil_counts<- oscil_extract(df, segments)
  oscil_feat <- classes_fromModel(oscil_counts, "oscillation", models[["oscill"]])

  bpArm_counts <- bpArm_extract(df, segments)
  bpArm_feat <- classes_fromModel(bpArm_counts, "arm_breakpoints", models[["bpArm"]])

  changep_counts <- changp_extract(df, segments)
  changep_feat <- classes_fromModel(changep_counts, "changepoint", models[["changepoint"]])

  ######## Aggregation in one features dataframe and one models list
  features <- Reduce(function(df1, df2) merge(df1, df2, by = "sample", all = TRUE),
                     list(cosim_feat, trater_feat, compact_sparse_feat, dels_feat, dups_feat, sv_types,
                          magnitudes_feat, chrmtrps_feat, ss_feat, bp5MB_feat, oscil_feat,
                          bpArm_feat, changep_feat))
  write.table(features, paste0(output_path, "features.tsv"), sep='\t', col.names = TRUE)

  print(paste0("Features.tsv and segmentation.tsv have been saved in ", output_path))
  return(features)
}
