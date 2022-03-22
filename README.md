# Code Snippets

Convenience snippets

## package dev 
### CRAN submission
- push all function changes to github  
- change the version number in DESCRIPTION  
- update readme.md, news.md  
- push these changes  
- devtools::check(cran = TRUE)   
- devtools::check_win_release()   
- update cran-comments.md with the output of cran checks  
- push this change and any remaining changes to github  
- devtools::submit_cran()  
- NOW CHECK EMAIL CONFIRM SUBMISSION via email (maintiner stable email)  
- once accepted by CRAN, tag a release with the same version number as the updated DESCRIPTION  
- tag release instructions: https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository  

## genomics related  

### correlation matrix p adjust

Given some data of variables (cols) with observations (rows), use Hmisc::rcorr() to get a symmetric correlation matrix of each pairwise comparison between columns.  object returned from rcorr contains a correlation matrix of rho values, p values and n valuee. The function below will calculate the adjusted p value. 

```{r}
p.adjust.cormat = function(hmisc.cor, method = 'fdr'){ 
  stopifnot(isTRUE(isSymmetric(hmisc.cor$P)))
  p.adj =  p.adjust(hmisc.cor$P[lower.tri(hmisc.cor$P)], method = method)
  p.adj.mx <- matrix(rep(0,ncol(hmisc.cor$P)*ncol(hmisc.cor$P)), nrow = ncol(hmisc.cor$P))
  p.adj.mx[lower.tri(p.adj.mx)] <- p.adj
  p.adj.mx[upper.tri(p.adj.mx)] <- p.adj
  diag(p.adj.mx) = 1
  colnames(p.adj.mx) = rownames(p.adj.mx) = colnames(hmisc.cor$P)
  return(p.adj.mx)
}
```

### manually calculated bonferonni of pairwise correlations 
given a matrix with number of columns Y 
the number of possible pairwise comparisons is (Y * (Y - 1) ) / 2  
The simple conservative bonferonni comparison multiplies the number of p values by the number of tests. 

```{r}

pmat = Hmisc::rcorr(x)$P
pvals = lower.tri(pmat)
ntest = ( ncol(x) * ncol(x) - 1 ) / 2 ) 
bonferonni.corrected = pvals * ntest

```

### convert entrez IDs to gene symbols visa versa 

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

## data wrangling

### Non standard evaluation for programming with dplyr verbs

[see here](https://dplyr.tidyverse.org/articles/programming.html) need to **Embrace** the argument in double brackets {{}}. 
```{r}
collapsedown = function(dat, grouping_var){
  gvar = rlang::sym(grouping_var)
  dat %>% 
    dplyr::group_by({{gvar}}) %>% 
    dplyr::summarise_each(list(~unique(.)))  
}
yy = collapsedown(dat = cell_metadata[ ,vars_all], grouping_var = 'sample')

```

### create multiple grouped aggregation and summary stats simultaneously
After group by summarize, just separate the new vars being added with commas within the summarize call. 
```
 # e.g. df is a long dataframe with celltype, timepoint, sampleid, gene, norm_count 
summary_df = df %>% 
  group_by(celltype, timepoint, sampleid, gene) %>% 
  summarize(
    gene_sd = sd(norm_count),
    gene_mad = mad(norm_count), 
    gene_mean = mean(norm_count), 
    ) 
```

### select specific columns matching string 
use matches 

```{r}
x = x %>%
  dplyr::select(matches('somestring')) 

```

### filter rows matching a string 

grepl inside the filter argument 
```{r}
x = x %>% 
  dplyr::filter(grepl("string1|string2", x = group)) 
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

## data visualization


```
  celltypes             cu              
"BC_Mem"              "lightslateblue"
"BC_Naive"            "#2B3D26"       
"CD103_Tcell"         "#E25822"       
"CD14_Mono"           "#654522"       
"CD16_Mono"           "#8DB600"       
"CD38_Bcell"          "#882D17"       
"CD4_CD161_Mem_Tcell" "#DCD300"       
"CD4_CD25_Tcell"      "#B3446C"       
"CD4_CD56_Tcell"      "maroon1"       
"CD4_CD57_Tcell"      "#604E97"       
"CD4_Efct_Mem_Tcell"  "#F99379"       
"CD4Naive_Tcell"      "#0067A5"       
"CD8_CD161_Tcell"     "darkseagreen1" 
"CD8_Mem_Tcell"       "#008856"       
"CD8_Naive_Tcell"     "#848482"       
"CD8_NKT"             "#C2B280"       
"HSC"                 "#BE0032"       
"IgA_CD14_Mono"       "#A1CAF1"       
"MAIT_Like"           "#F38400"       
"mDC"                 "#875692"       
"NK"                  "#F3C300"       
"pDC"                 "#222222"       
"BC_Mem"              "midnightblue"  

```


### Make sure the 0 point of a diverging color palette is white in a heatmap

A diverging palette should always be used to display standardized vars or data above and below 0. The mid point of the heatmap should ALWAYS be the mid point in the diverging palette or else you can confuse the reader. 

```{r}
# got this simple solution from the post below
# https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r

# any diverging palette 
cu = RColorBrewer::brewer.pal(name = 'PRGn',n = 12)
range <- max(abs(mtx)); 
pheatmap(mtx, color = cu, breaks = seq(-range, range, length.out = 12))

```
to find a diverging palette: 
`RColorBrewer::display.brewer.all()`


### Make the legend key shape a different size than the shape on the plot 

```{r}
guides(color = guide_legend(override.aes = list(size=2))) 

```

### adjust a single color to manually add transparency 
```
grDevices::adjustcolor( "red", alpha.f = 0.2)
```

### clean theme 

```
boxbox = list(
  theme_bw(),
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank())
        )
```

### rotate axis 
```{r}
 + theme(axis.text.x=element_text(angle = -90, hjust = 0))
```

## strings 

For doing annoying things with strings and regex 

### separate something by the first space get the short name 
e.g. btm names...
```
new_names = vapply(strsplit(oldnames," "), `[`, 1, FUN.VALUE=character(1))
```

### remove a subset of vector elements based on regex 
```{r}
gene_rm = str_detect(string = gene_union, pattern = "RP[0-9]-" )
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
```

### remove accidental committed file from git history 
```
#https://rtyley.github.io/bfg-repo-cleaner/#usage
java -jar  Downloads/bfg-1.14.0.jar --delete-files id_{NAME_OF_HUGE_HTML_FILE.html}  PATH/TO/LOCAL_GIT_REPO/myrepo 
```

### to ignore entire file path in project directory recursively in the gitignore file add: 
path/I/Want/To/Ignore*

### to ignore a single file in each sub directory with the same name (ignore all generated data and figures)
generated_data/
data/
figures/

### generating the r markdown included in git repo as a pdf 

```
# https://stackoverflow.com/questions/7694887/is-there-a-command-line-utility-for-rendering-github-flavored-markdown
# install grip 
# cd to dir where markdown is 
grip README.md
# go to localhost/whatever and print -> save as pdf 

```


## Linux and HPC related 

### Show directory structure as a tree 

```

# without the files just directories 
ls -R | grep ":$" | sed -e 's/:$//' -e 's/[^-][^\/]*\//--/g' -e 's/^/   /' -e 's/-/|/'

# with the directories 
find . | sed -e "s/[^-][^\/]*\//  |/g" -e "s/|\([^ ]\)/|-\1/"

```

## Reticulate 
### python virtual env 
```# eg for umap in r must have python installed.
virtualenv_create("r-reticulate")
virtualenv_install("r-reticulate", "umap-learn")
use_virtualenv("r-reticulate")
library(umap)
```
## system related

### speed up time machine backups 
```
sudo sysctl debug.lowpri_throttle_enabled=1
```
### stay awake 
```
# last command is time in seconds - this for 27hours, e.g. 86400 is stay awa
caffeinate -d -i -m -s -t 1000000
```
