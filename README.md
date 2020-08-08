# snip
useful snippets


# remove a subset of genes based on regex 
```{r}
# starting from an indexed list (eg. by cell type) 
gene_union = unlist(genes_scmod) %>% unique()
gene_rm2 = str_detect(string = gene_union, pattern = "RP[0-9]-" )
gene_union = gene_union[!gene_rm2]
gene_rm = gene_union[grep(pattern = "RP11|MT-", gene_union)]
gene_union = gene_union[gene_union %ni% gene_rm]
```
