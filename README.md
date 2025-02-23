# Code Snippets

Miscellaneous functions and documentation.  


## Website workflow

### github pages rendering
After making changes use `rmarkdown::render_site()`. Then commit + push. 

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
- [instructions on tagging a release](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)
- reference: [R package release](https://r-pkgs.org/release.html)


### update packages with a packagedown site 
- push all function changes to   
- change the version number in DESCRIPTION  
- update readme.md, news.md  
- push these changes  
- pkgdown::build_site()
- push to github
- tag release instructions: https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository  

## stats related 

### quick linear model of predictors association with outcome 

Useful for quickly testing the direction of association between features and response variable, for example, get the direction of association with the outcome of the features with high variable importance scores from a ML model. 

structure of `dat`

| group.factor    | gene1  | gene2 |
| -------- | ------- |-------- |
| high  | 0.5   |0.6|
| low | 0.4   |0.3|
| high    | 0.3    |0.2 |

```{r}

# classification example with 2 classes - ensure factor levels are ordered correctly 
dat$group = factor($group, levels = c('low', 'high'))

get.linear.coef = function(x, y) {
  test = lm(x ~ y)
  out <- c(beta = coef(test)[2])
}

# optionally define feature subset in place of -1 indexing below, which removes the group var from the features. 

feature.coef = apply(dat[ ,-1], 2,  get.linear.coef, y = dat$group)
dplot = data.frame(beta = feature.coef)
plot(dplot) 

```

### quick t test of predictors association with outcome 

Same as above but use a Wilcoxon rank test or t test.  
Here i also collect the test p value...   
Again assuming the structure of the dat `dat` is the same as above.  
 
```
# get t stats of individual features 
do.ttest = function(x, y) {
  ttest = t.test(x ~ y)
  out <- c(tstat = ttest$statistic, pval = ttest$p.value)
  return(out)
}
tval = apply(dat[ , -1], 2, do.ttest, y = dat$group) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column('feature')
```

Quick visualization of group differences for select vars in the matrix. 

```
d.long = d %>%
  gather(variable, value, Var1:VarN)

p = ggplot(data = d.long, aes(x = value, y = group,  fill = group, color = group)) +
  geom_boxplot(size = 0.2, outlier.size = 0.1)  +
  facet_grid(rows = vars(variable), scales = 'free', space = 'free') +
  theme(strip.text.y = element_text(angle = 0))

```

### simple scale function

```
# simple function for scaling because scale return format is annoying
scale.simple = function(x) { 
  (x - base::mean(x)) / stats::sd(x)
  }

```

## genomics related  


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

### Seurat v2 conversion to v3 or 4 + 

```{r}

# R 4.2 Seurat v 3 or 4 
suppressMessages(library(Seurat))

# your path datapath = my_path
# dir.create(datapath)

# read main object Seurat 2.4 
# s = readRDS(file = 'my_path/'))

# save cell IDs 
cells = rownames(s@meta.data)

# extract data structure from Seurat V2 object
md = s@meta.data[cells, ]
rna.raw = s@raw.data[ ,cells]
rna.norm = s@data[ ,cells]
adt.raw = s@assay$CITE@raw.data[ ,cells]
adt.norm = s@assay$CITE@data[ ,cells]

# make Seurat v4 object 
ss = CreateSeuratObject(counts = rna.raw, meta.data = md)
ss@assays$RNA@data = rna.norm
adt.assay = CreateAssayObject(counts = adt.raw)
ss[['CITE']] <- adt.assay
ss@assays$CITE@data <- adt.norm
saveRDS(ss, file = paste0(datapath, 'new.rds'))

```

## data wrangling

## specify a conflicted function preference
prefer to use dplyr filter over another package. 
```
conflicted::conflict_prefer("filter", "dplyr")
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

### Non standard evaluation for programming with dplyr verbs
Note dont do this anymore better to not program with tidyverse functions.... 

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
  theme(panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank() 
        )
```

### ggplot make the axis tics more natural and uniform 

```
# e.g. y axis 
scale_y_continuous(breaks= scales::pretty_breaks())

```

### ggplot rotate axis 
```{r}
 theme(axis.text.x=element_text(angle = -90, hjust = 0))
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

### put quotes around each element of a vector - useful for making gene lists
```{r}
vector = c('BUB1B, CCL22, CD58, CD59, DBI')
cat(gsub("\\b", '"', vector, perl=T)) 

```

## error handling and control flow with trycatch 

When iterating through and testing want the function to continue if it encounters an error. 

### using trycatch 
generic 
```
tryCatch(function(object = x, args = args), error = function(e) return(NA))
```
Machine learning pipeline example with additional options 


```
  # Loop over each algorithm in the list
  for (i in seq_along(caret.algorithm)) {
    method.use <- caret.algorithm[i]
    grid.use   <- grid.list[[i]]
    
    print(paste0("training ", method.use))
    model.fits[[i]] <- tryCatch({
      train(
        group ~ .,
        data = data.combined,
        method = method.use,
        tuneGrid = grid.use,
        trControl = ctrl,
        metric = "ROC"
      )
    }, error = function(e) {
      message("Error for method ", method.use, ": ", e$message)
      NA
    })
    message("done")
  }

```



## Lists 

### combine a list into a named dataframe 

Binding named lists into rows of a data.frame with conversions into different data containers like data.frame or tibble

```{r}

function(x) {
  # tibble::as_tibble(x <- data.table::setDF(
  # data.frame(x <- data.table::setDF(
    data.table::rbindlist(x, use.names = TRUE, fill = TRUE, idcol = "id")
))
}

```

### flatten a list of lists into the top level list base R only  
Optionally append the name of the first level list with the names of second level list being combined.  

For example the object, `pathways` with the following structure:

$CD4Tcell_naive  
"HALLMARK_TNFA_SIGNALING_VIA_NFKB"  

$CD8Tcell_dim  
"HALLMARK_TNFA_SIGNALING_VIA_NFKB"  
"HALLMARK_KRAS_SIGNALING_UP"  

Gets combined into a single list:  
"CD4Tcell_naive--HALLMARK_TNFA_SIGNALING_VIA_NFKB"  
"CD8Tcell_dim--HALLMARK_TNFA_SIGNALING_VIA_NFKB"  
"CD8Tcell_dim--HALLMARK_KRAS_SIGNALING_UP"  

```{r}
pathways = do.call(
  c, # combine
  lapply(names(pathways), function(x) {
    setNames(pathways[[x]], paste0(x, "--", names(pathways[[x]]))) 
  })
)
```



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

### remove accidental COMITTED file from git history 

```
# cd to the root dir of the project 
git reset HEAD ~git reset HEAD~
```

### remove accidental PUSHED file from git history  
If you already comitted the file it is a bit more challenging, but this tool has worked. 
```
# download this tool: https://rtyley.github.io/bfg-repo-cleaner/#usage
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
