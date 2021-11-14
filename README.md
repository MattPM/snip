# snip

Code snippets for convenience.

## genomics stuff 

### convert to entrez 

```{r}
ent =
  tryCatch(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys = res$CC10$gene, columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL"),
    error = function(e) return(NA)
  )
dup = duplicated(ent$SYMBOL)
ent = ent[!dup, ]$ENTREZID
ent = ent[!is.na(ent)]
```

## tidyverse 

### calculate average of a list of gene stes in a list of gene matrices 

```
# av = list of RNAseq pbulk matrices, sig = gene list 

d_sig = lapply(av, function(x){
  lapply(sig, function(y){
    Matrix::colMeans(x[intersect(y, rownames(x)),  ])
    }) %>% 
    do.call(rbind, .)
  })

  
```



Tidy syntax based 

### map values indexed by a list 

```{r}
# map a vector of values from a metatable met directly to aother dataframe with data d  
# with the shared variable subject (in d) and `Subject ID` (in met)
# d = some data frame 
# met = some metadata columns to append onto d matching by subject ids 
ids_map = names(met)[2:length(names(met))] # 1 is `Subject ID`, dont' map that 

# add metadata based on map values list wise 
for (i in 1:length(ids_map)) {
  d[[ids_map[i]]] = plyr::mapvalues(ddf$subject, from = met$`Subject ID`, to = met[[ids_map[i]]])
} 
# rbind .. 
```

### refer to a function argument (environment-variable) as a data-variable to use tidy non standard eval within a function
#### this is also a good example of collapsing by unique values e.g. to create a design matrix 
[see here](https://dplyr.tidyverse.org/articles/programming.html) we need to **Embrace** the argument in double brackets {{}} 
```{r}
collapsedown = function(dat, grouping_var){
  gvar = rlang::sym(grouping_var)
  dat %>% 
    dplyr::group_by({{gvar}}) %>% 
    dplyr::summarise_each(list(~unique(.)))  
}
yy = collapsedown(dat = cell_metadata[ ,vars_all], grouping_var = 'sample')

```

### dplyr::group_by() %>% dplyr::summarize() create multiple summary vars at once
remember one can just separate the new vars being added within the summarize call after the pipe. 
```
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

### filter rows matching a string 
just use grepl inside the filter argument 
```{r}
# combine a list of dataframes and filter based on the group variable matching a string 
df2 = list %>%
  bind_rows() %>% 
  filter(grepl("this_string|or_this_string", x = group)) 
```

### lead and lag to calculate a fold change
assuming log transformed data so subtracting to calculate a fold change. 
Calculate the log fold change of the variable "average" over 2 timepoints.  
The filter removes the unused rows because they are now summarized by the fold change. 
```{r}
new_df = 
  old_df %>% 
  arrange(category, sampleid) %>% # optional
  mutate(fold_change = lead(average) - average) %>% 
  filter(timepoint == "second_timepoint") 
 
```

### full join with frame-specific variable rename 
joining similar named dataframes by a single variable and retaining frame specific columns by replacing ".x" ".y" string with something that is more obvious 

```{r}
data = full_join(df1, df2, by = "var1") 
oldnames = colnames(data)
newnames = str_replace(string = oldnames, pattern = "\\.x", replacement = "_DF1")
newnames = str_replace(string = newnames, pattern = "\\.y", replacement = "_DF2")
colnames(data) = newnames
# use case 
dsub = data %>% select("gene", "logFC_DF1"    "P.Value_DF2"    "logFC_DF2" )
```
## plotting

### Make the legend key shape a different size than the shape on the plot 

```{r}
guides(color = guide_legend(override.aes = list(size=2))) 

```

### adjust a single color to manually add transparency 
```
adjustcolor( "red", alpha.f = 0.2)
```

### stripped down clean aesthetic mtheme

```
theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        legend.position = 'none')) +
```


### reorder categorical axis by manual order 
new_order = c("x3", "x1", "x2")
df$x = factor(df$x, levels = new_order)

### reorder categorical y axis by average y axis value of tidy data
av_y = df %>% group_by(x) %>% summarize(meany = mean(y)) %>% arrange(meany) %$% y
df$y = factor(df$y, levels = av_y)

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
```
### gsea bubble individual 
```
p = ggplot(df, aes(x = module, y = padj )) +
  geom_point(aes(size = NES),  color = "#195e83") +
  theme_bw() + 
  theme(axis.text.y = element_text(face = "bold", color = "black", size = 5)) + 
  theme(axis.text.x = element_text(face = "bold", color = "black")) + 
  theme(axis.title.x = element_text(face = "bold", color = "black")) + 
  ylab(" -log10 enrichment p value ") + 
  xlab("") + 
  geom_hline(yintercept = 1.3, linetype = "dashed", size = 0.25) + 
  coord_flip()
```
### label points by group on a scatterplot a la seurat 
```{r}
# taken from seurat V2 DimPlot
df = c(x, y, celltype) 

# for example
ggplot(df, aes(x = x, y = y) + geom_point()

# define cluster center 
centers = data.plot %>% 
  dplyr::group_by(ident) %>% 
    summarize(x = median(x = x), y = median(x = y))

# add that data layer to the plot 
p = p + geom_point(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0) + 
            geom_text(data = centers, mapping = aes(label = celltype))
            
            
            
centers = df2 %>% dplyr::group_by(cluster_cohort) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# add that data layer to the plot 

p = ggplot(df2, aes(x = UMAP_1, y = UMAP_2, color = cluster_cohort)) +
  geom_point(size = 0.2, alpha = 0.3, show.legend = FALSE) 
p = p + geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0, show.legend = FALSE) + 
  ggrepel::geom_text_repel(data = centers, mapping = aes(label = cluster_cohort), show.legend = FALSE) + 
  scale_color_manual(values = distinct_60)
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
### put quotes around each element of a vector 
```{r}
vector = c('BUB1B, CCL22, CD58, CD59, DBI')
cat(gsub("\\b", '"', vector, perl=T)) 

```

## error handling 

## Generic trycatch 
```
tryCatch(function(object = x, args = args), error = function(e) return(NA))

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

### simple git workflow 
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
### generating the r markdown included in git repo as a pdf 

```
# https://stackoverflow.com/questions/7694887/is-there-a-command-line-utility-for-rendering-github-flavored-markdown
# install grip 
# cd to dir where markdown is 
grip README.md
# go to localhost/whatever and print -> save as pdf 


```


## Linux and HPC 

### Show directory structure as a tree 

```

# without the files just directories 
ls -R | grep ":$" | sed -e 's/:$//' -e 's/[^-][^\/]*\//--/g' -e 's/^/   /' -e 's/-/|/'

# with the directories 
find . | sed -e "s/[^-][^\/]*\//  |/g" -e "s/|\([^ ]\)/|-\1/"

```

### R 3.5 linux workflow 
```
# mem request eg 
qrsh -l mem_free=240G,h_vmem=30G -pe threaded 8 

step 3 change to your working directory where you want to start the scripts form
cd /path/to/your_directory

# R 3.5 
module load r/3.5.1_20190921
```
## python virtual env 
```
# eg for umap in r must have python installed.
virtualenv_create("r-reticulate")
virtualenv_install("r-reticulate", "umap-learn")
use_virtualenv("r-reticulate")
library(umap)
config = umap.defaults
config$n_neighbors = 35
config$min_dist = 0.6

# run umap
ump = umap(mymatrix,config = config)
```
## system stuff

### speed up time machine backups 
```
sudo sysctl debug.lowpri_throttle_enabled=1
```

### stay awake 
```
# last command is time in seconds - this for 27hours, e.g. 86400 is stay awa
caffeinate -d -i -m -s -t 1000000
```

## links 
https://bioconductor.org/about/release-announcements/  
https://awesomeopensource.com/project/EmilHvitfeldt/r-color-palettes  
https://bioconductor.org/packages/3.8/workflows/html/simpleSingleCell.html
http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
