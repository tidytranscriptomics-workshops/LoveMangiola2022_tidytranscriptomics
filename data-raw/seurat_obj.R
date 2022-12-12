seurat_obj <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/cancer_only_analyses/lymphoid/cancer_lymphoid_cell_type_curated.rds")

seurat_obj = seurat_obj |> RunPCA() |> select(-contains("UMAP")) |> RunUMAP(dims=1:20)
seurat_obj = seurat_obj |> select(cell, file, 3, 8, 9,S.Score, G2M.Score , Phase , curated_cell_type , contains("UMAP"))
seurat_obj = seurat_obj %>% filter(.cell %in% (seurat_obj %>% sample_frac(0.5) %>%  pull(.cell) %>% c(seurat_obj %>% filter(grepl("Delta", curated_cell_type)) %>% pull(.cell)) %>% unique))
seurat_obj = seurat_obj[VariableFeatures(seurat_obj),]
seurat_obj[["SCT"]]@scale.data = seurat_obj[["SCT"]]@scale.data[c("CD3D", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B"),]
seurat_obj[["SCT"]]@data = seurat_obj[["SCT"]]@data[c("CD3D", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B"),]
DefaultAssay(seurat_obj) = "SCT"
seurat_obj[["integrated"]] = NULL

# Update treatment 
slot(seurat_obj$SCT@SCTModel.list[[1]], 'median_umi') = median(seurat_obj$SCT@SCTModel.list[[1]]@cell.attributes$umi)

seurat_obj = seurat_obj |> 
	tidyseurat::extract(file, "sample", "../data/.*/([a-zA-Z0-9_-]+)/outs.+", remove = FALSE) |> 
	left_join(
		readr::read_csv("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/OMBC - Sheet1.csv") |>
			select(sample, type)
	)  |>
	mutate(treatment = if_else(type=="MBC", TRUE, FALSE)) |> 
	select(-sample, -type) |> 
	mutate(treatment = factor(treatment))

# Filter ribosome rich
seurat_obj = seurat_obj |>  filter(curated_cell_type != "CD4+_ribosome_rich")

# Save
save(seurat_obj , file="data/seurat_obj.rda", compress = "xz")
