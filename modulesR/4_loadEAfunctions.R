##Enrichment analysis functions

source("modulesR/0_loadLibraries.R")
loadpkg("org.Hs.eg.db")
loadpkg("biomaRt")
loadpkg("dplyr")
loadpkg("clusterProfiler")

orgdb = "org.Hs.eg.db"
biomart_dataset = "hsapiens_gene_ensembl"
keggname = "hsa"
mart = biomaRt::useMart(biomart = "ensembl", dataset = biomart_dataset)
#entrez = biomaRt::getBM(attributes = c("refseq_mrna", "entrezgene"), mart = mart)
#entrez$entrezgene = as.character(entrez$entrezgene)
entrezsymbol = biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), mart = mart)
entrezsymbol$entrezgene = as.character(entrezsymbol$entrezgene)


summarize_cp = function(res, comparison) {
  summaries = data.frame()
  for (ont in names(res)) {
    ontsum = as.data.frame(res[[ont]])
    ontsum$ont = ont
    summaries = rbind(summaries, ontsum)
  }
  summaries$comparison = comparison
  return(summaries)
}

enrich_cp = function(res, comparison, type="over") {
  res = res %>% data.frame()  %>% left_join(entrezsymbol, by = "hgnc_symbol") %>% filter(!is.na(entrezgene))
  # universe = brcaData@GISTIC@AllByGene$Gene.Symbol
  if(type=="all"){
    res <- res %>% filter(abs(logFC) > lfc & adj.P.Val < pval) # lfc and pval threshold defined above in the volcano plot
    genes = res$entrezgene
    
    mf = clusterProfiler::enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH",
                                   qvalueCutoff = 1, pvalueCutoff = 1)
    cc = clusterProfiler::enrichGO(genes,  OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH",
                                   qvalueCutoff = 1, pvalueCutoff = 1)
    bp = clusterProfiler::enrichGO(genes,  OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH",
                                   qvalueCutoff = 1, pvalueCutoff = 1)
    kg = clusterProfiler::enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1,
                                     qvalueCutoff = 1, pAdjustMethod = "BH")
    all = list(mf = mf, cc = cc, bp = bp, kg = kg)
    all[["summary"]] = summarize_cp(all, comparison)
    return(all)
  }
  if(type=="over"){
    res.over <- res %>% filter(logFC > lfc  & adj.P.Val < pval)
    genes = res.over$entrezgene
    mf = clusterProfiler::enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH",
                                   qvalueCutoff = 1, pvalueCutoff = 1)
    cc = clusterProfiler::enrichGO(genes,  OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH",
                                   qvalueCutoff = 1, pvalueCutoff = 1)
    bp = clusterProfiler::enrichGO(genes,  OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH",
                                   qvalueCutoff = 1, pvalueCutoff = 1)
    kg = clusterProfiler::enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1,
                                     qvalueCutoff = 1, pAdjustMethod = "BH")
    all = list(mf = mf, cc = cc, bp = bp, kg = kg)
    all[["summary"]] = summarize_cp(all, comparison)
    return(all)
  }
  
  if(type=="under"){
    res.under <- res %>% filter(logFC < -lfc & adj.P.Val < pval)
    genes = res.under$entrezgene
    mf = clusterProfiler::enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH",
                                   qvalueCutoff = 1, pvalueCutoff = 1)
    cc = clusterProfiler::enrichGO(genes,  OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH",
                                   qvalueCutoff = 1, pvalueCutoff = 1)
    bp = clusterProfiler::enrichGO(genes,  OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH",
                                   qvalueCutoff = 1, pvalueCutoff = 1)
    kg = clusterProfiler::enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1,
                                     qvalueCutoff = 1, pAdjustMethod = "BH")
    all = list(mf = mf, cc = cc, bp = bp, kg = kg)
    all[["summary"]] = summarize_cp(all, comparison)
    return(all)
  }
}

convert_enriched_ids = function(res, entrezsymbol) {
  res = res %>% mutate(geneID = strsplit(as.character(geneID), "/")) %>% tidyr::unnest(geneID) %>% 
    left_join(entrezsymbol, by = c(geneID = "entrezgene")) %>% group_by(ID, 
                                                                        Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, ont, 
                                                                        comparison) %>% summarise(geneID = paste(geneID, collapse = "/"), symbol = paste(hgnc_symbol, 
                                                                                                                                                         collapse = "/"))
  return(res)
}

