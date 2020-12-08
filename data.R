# load from CytoGLMM vignette
load_data = function() {
  load("se_aghaeepour2017immune.Rdata")
  exprs = assay(se_aghaeepour2017immune)
  sample_info = rowData(se_aghaeepour2017immune)
  sample_info_names = names(sample_info)
  df = cbind(as.data.frame(exprs), as.data.frame(sample_info))
  df %<>% as_tibble
  protein_names = colData(se_aghaeepour2017immune) %>% 
    as.data.frame %>% 
    dplyr::filter(type == "function") %>%
    .$protein_name
  gate_protein_names = colData(se_aghaeepour2017immune) %>% 
    as.data.frame %>% 
    dplyr::filter(type == "phenotype") %>%
    .$protein_name
  df %<>% select_at(c("donor", "term", "celltype", protein_names))
  df %<>% dplyr::rename(condition = term)
  df$celltype %<>% as.factor()
  df$donor %<>% as.factor()
  df
}