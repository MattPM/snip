# snip
code snippets used frequently for reference. 

## tidyverse 

Tidy syntax based 

### map values indexed by a list 
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

### dplyr::group_by() %>% dplyr::summarize() create multiple summary vars at once
```
# just separate the new vars being added within the summarize call after the pipe! 
# df is a long dataframe with celltype, timepoint, sampleid, gene, norm_count 
summary_df = df %>% 
  group_by(celltype, timepoint, sampleid, gene) %>% 
  summarize(
    gene_sd = sd(norm_count),
    gene_mad = mad(norm_count), 
    gene_mean = mean(norm_count), 
    ) 
```
### select specific columns matching string 
```{r}
Get specific columns matching strings 

x = x %>%
select(-ENSEMBL63_GENE_ID) %>% 
select(matches('D000|D001|ENSEMBL63_GENE_NAME')) %>% 
column_to_rownames("ENSEMBL63_GENE_NAME") 

```

## ggplot 

ggplot things 

### rotate axis 
```{r}
 + theme(axis.text.x=element_text(angle = -90, hjust = 0))
```
### add a marginal plot (like a histogram) 
```
# main plot
p1 = ggplot(df, aes(x = x, y = y)) + geom_point(shape = 21, size = 2.5, color = "grey" ,stroke = 0.1)

# marginal x axis for top
xtop = axis_canvas(p1, axis = "x") + geom_density(data = df, aes(x = x, fill = timepoint))

# marginal x axis for bottom 
xbottom = axis_canvas(p1, axis = "x") +
  geom_point(data = df, aes(x = x, y = whatever), shape = 16, alpha = 0.2) + 
  geom_smooth(data = df, aes(x = x, y = whatever),method = "loess", color = "blue3") 

# marginal boxplot for y axis 
ybox = axis_canvas(p1, axis = "y", coord_flip = FALSE) +
  geom_boxplot(data = df, aes(y = component_2, fill = group_separator), show.legend = TRUE)

# combine 
p2 = insert_xaxis_grob(p1, xtop, grid::unit(.2, "null"), position = "top")
p3 = insert_yaxis_grob(p2, ybox, grid::unit(.2, "null"), position = "right")
p4 =  insert_xaxis_grob(p3, xbottom, grid::unit(.6, "null"), position = "bottom")
p6 = ggdraw(p4)
p6
# save 
ggsave(p6, filename = paste0(figpath, "whatever.pdf"), width = 8, height = 10)

#
```


## strings 

For doing annoying things with strings and regex 

### separate something (BTM list names) by the first space to get the short name 

```
# btm_names = the names of the list BTM e.g. names(btm) str(btm) -- list 

new_btm_names = vapply(strsplit(btm_names," "), `[`, 1, FUN.VALUE=character(1))
```
### remove a subset of genes based on regex 
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

### var genes in scran from seurat remove mt rp 
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
## Lists 

### reorder a list based on vector
```
mylist = mylist[order(match(names(mylist), vector_with_desired_order))]
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

######### 
# to ignore entire file path in project directory recursively in the gitignore file add: 
path/I/Want/To/Ignore*

# to ignore a single file in each sub directory with the same name (ignore all generated data and figures)
generated_data/
data/
figures/

```
## links 
https://bioconductor.org/about/release-announcements/  
https://awesomeopensource.com/project/EmilHvitfeldt/r-color-palettes  
http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
