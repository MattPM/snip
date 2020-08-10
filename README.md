# snip
code snippets used frequently for reference. 

## tidyverse 

Tidy syntax based 

## map values indexed by a list 
```{r}
# read Meta data table
meta_table = read_delim(file = "git_ignore/meta_table.txt", delim = "\t")
ids_map = names(meta_table)

# some meta data table e.g. seurat or bioconductor pheno data 
md = s@meta.data

# add metadata based on map values list wise 
meta = list()
for (i in 1:length(ids_map)) {
  md = s@meta.data
  md[[ids_map[i]]] = plyr::mapvalues(md$sample, from = meta_table$sample, to = meta_table[[ids_map[i]]] )
  meta[[i]] = md %>% select(ids_map[i])
} 

# metadata to add to object 
meta = do.call(cbind, meta) 
rownames(meta) = NULL
```

## get specific columns matching string 

```{r}
Get specific columns matching strings 

x = x %>%
select(-ENSEMBL63_GENE_ID) %>% 
select(matches('D000|D001|ENSEMBL63_GENE_NAME')) %>% 
column_to_rownames("ENSEMBL63_GENE_NAME") 

```

## ggplot 

ggplot things 

## rotate axis 
```{r}
 + theme(axis.text.x=element_text(angle = -90, hjust = 0))
```



## strings 

For doing annoying things with strings and regex 

## separate the BTM names by the first space to get the short name 

```
# btm_names = the names of the list BTM e.g. names(btm) str(btm) -- list 

vapply(strsplit(btm_names," "), `[`, 1, FUN.VALUE=character(1))
```

## remove a subset of genes based on regex 
```{r}
# starting from an indexed list (eg. by cell type) 
gene_union = unlist(genes_scmod) %>% unique()
gene_rm2 = str_detect(string = gene_union, pattern = "RP[0-9]-" )
gene_union = gene_union[!gene_rm2]
gene_rm = gene_union[grep(pattern = "RP11|MT-", gene_union)]
gene_union = gene_union[gene_union %ni% gene_rm]

```

## put quotes around each element of a vector 
```{r}
vector = c('BUB1B, CCL22, CD58, CD59, DBI')
cat(gsub("\\b", '"', vector, perl=T)) 

```


## single cell 


## var genes in scran from seurat remove mt rp

```
bsce = Convert(from = bc, to = "sce")
poiss_trend = scran::trendVar(bsce, use.spikes = FALSE)
var_fit = scran::decomposeVar(x = bsce, fit = poiss_trend)
var_gene = var_fit %>% as.data.frame() %>% rownames_to_column("gene") %>% filter(bio > 0.001) %$% gene
gene_rm = str_detect(string = var_gene, pattern = "RP[0-9]-" )
var_gene = var_gene[!gene_rm]
gene_rm2 = var_gene[grep(pattern = "RP11|MT-", var_gene)]
var_gene = var_gene[var_gene %ni% gene_rm2]
```


## Git 

```
cd _path_to_repo_
git init 

git pull https://github.com/MattPM/repo

git add .
git commit -m "add dsb norm funcion"

git config --global user.email my_email@__.com
git commit --amend --reset-author

git remote add origin https://github.com/MattPM/repo

git pull origin master 
git push origin master 

git status 
# On branch master
# nothing to commit, working tree clean

```
