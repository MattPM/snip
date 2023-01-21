# Code Snippets

Convenience snippets

## R package development  
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

### update packages with a packagedown site 
- push all function changes to github  
- change the version number in DESCRIPTION  
- update readme.md, news.md  
- push these changes  
- pkgdown::build_site()
- push to github
- tag release instructions: https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository  

## stats related 

extract variable T statistics on a matrix, with a grouping factor as the first column. Useful for sanirt check on direction of variable importance. 
Assuming data structure of dataframe `d` looks like: 
```
group Var1 Var2
"low" 5   3
"low" 5   5
"high"  3   7
"high"  3   8

```

code 
```
# get t stats of individual features 
do.ttest = function(x, y) {
  ttest = t.test(x ~ y)
  out <- c(tstat = ttest$statistic, pval = ttest$p.value)
  out
}
tval = apply(d[ , -c(1)], MARGIN = 2, FUN = do.ttest, y = d$group) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('feature')
```

Quick visualization of group differences for select vars in the matrix. 

```
d.long = d %>%
  gather(variable, value, Var1:VarN)

p = ggplot(data = d.long, aes(x = value, y = group,  fill = group, color = group)) +
  geom_boxplot(size = 0.2, outlier.size = 0.1)  +
  facet_grid(rows = vars(variable), scales = 'free', space = 'free') +
  theme(strip.text.y = element_text(angle = 0))

# save long fig size 

```

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


## data visualization

### Make the 0 point of a diverging color palette white

A diverging palette should always be used to display standardized vars or data above and below 0. The mid point of the heatmap should ALWAYS be the mid point in the diverging palette or else you can confuse the reader. 

```{r}
# got this simple solution from the post below
# https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r

# any diverging palette 
cu = RColorBrewer::brewer.pal(name = 'PRGn',n = 12)
range <- max(abs(mtx)); 
pheatmap(mtx, color = cu, breaks = seq(-range, range, length.out = 12))

# to find a diverging palette: 
RColorBrewer::display.brewer.all()

```
### color adjust a single color to manually add transparency 
```
grDevices::adjustcolor( "red", alpha.f = 0.2)

```

## ggplot related 


### ggplot change legend key shape from the main geom on the plot

```{r}
guides(color = guide_legend(override.aes = list(size=2))) 

```

### gplot Change all the legend element sizes 
```{r}
theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)) #change legend text font size
```

### ggplot clean theme 
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

### ggplot rotate axis 
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
gene.rm = str_detect(string = gene.vector, pattern = "RP[0-9]-" )
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

### git workflow 
```
cd _path_to_repo_
git init 

git pull https://github.com/MattPM/repo

git add .
git commit -m "add funcion"

# if error related to user auth etc. try this
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
