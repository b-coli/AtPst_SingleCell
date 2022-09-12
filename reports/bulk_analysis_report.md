Bulk Data Analysis
================
2022-09-08

6 whole-tissue and 2 protoplasted (bulked) samples were processed by:

-   Trimming adapters and filtering low-quality sequence (trimgalore)
-   Mapping to the Arabidopsis genome using STAR (same STAR version and
    reference used for single-cell data)

To analyze differential expression as a function of protoplasting and P.
syringae treatment, I first loaded the data into a digital expression
matrix. Here, I plot expression distributions for each sample.

``` r
bulk_dge_raw %>% 
  tidyr::gather("Sample_Name", "Expression") %>% 
  left_join(sample_metadata, by = "Sample_Name") %>%
  ggplot(aes(x = Sample_Name, y = log1p(Expression))) + 
  geom_boxplot(aes(fill = interaction(Treatment, Sample_Type))) + 
  scale_fill_discrete(name = "Sample Type") +
  theme_bw()
```

![](/home/bjcole/AtPst_SingleCell/reports/bulk_analysis_report_files/figure-gfm/plot_raw_distribution-1.png)<!-- -->

Expression seems fairly even across samples, but could do with a little
normalization/filtering

To do this, I imported all data into a new DGEList object (edgeR). I
filtered out low-expressed or low-variance genes, and defined a design
matrix to identify genes induced by protoplasting, independent of pst3k
treatment.

``` r
design <- bulk_erobj$design
row.names(design) <- row.names(bulk_erobj$samples)

kableExtra::kable(design, format = "latex") %>% kableExtra::kable_styling(full_width = FALSE)
```

I then normalized the expression matrix. Here, I plot the normalized
expression distribution per sample.

``` r
bulk_dge_norm %>% 
  as_tibble(rownames = "Locus") %>%
  tidyr::gather("Sample_Name", "Expression", -Locus) %>%
  left_join(sample_metadata) %>%
  ggplot(aes(x = Sample_Name, y = log1p(Expression))) + 
  geom_boxplot(aes(fill = interaction(Treatment, Sample_Type))) + 
  scale_fill_discrete(name = "Sample Type") +
  theme_bw()
```

    ## Joining, by = "Sample_Name"

![](/home/bjcole/AtPst_SingleCell/reports/bulk_analysis_report_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

I then used edgeR to build generalized linear models (GLMs) to the data,
given the design matrix specified above. I used an adjusted p-value
threshold of 0.01 and a log-fold change cutoff of 0.5 to define genes
differentially expressed by protoplasting (controlling for DC3000
treatment and interaction between protoplasting and DC3000 treatment). I
then repeated the same procedure for genes differentially expressed by
DC3000 treatment, using a log-fold change cutoff of 2 and an adjusted
p-value of 0.01. Here, I plot the expression distributions of genes
based on the effect protoplasting has on them…

``` r
expr_data <- bulk_dge_norm %>%
  as_tibble(rownames = "Locus") %>%
  tidyr::gather("Sample_Name", "Expression", -Locus) %>%
  left_join(sample_metadata, by = "Sample_Name") %>%
  left_join(bulk_de_table, by = "Locus")

ggplot(expr_data, aes(x = Sample_Name, y = log1p(Expression))) + 
  geom_boxplot(aes(fill = de_type_proto)) + 
  scale_fill_discrete(name = "Protoplast DE Gene") +
  theme_bw()
```

![](/home/bjcole/AtPst_SingleCell/reports/bulk_analysis_report_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

…and also by mean expression between protoplasted and bulk samples

``` r
expr_data %>%
  group_by(Locus, Sample_Type, Treatment, de_type_proto) %>%
  summarize(mean_expr = mean(Expression)) %>%
  spread(Sample_Type, mean_expr) %>%
  ggplot(aes(x = log1p(Leaf), y = log1p(Protoplast))) + 
  geom_point(aes(color = de_type_proto), size = 0.5) +
  facet_wrap("Treatment")
```

    ## `summarise()` has grouped output by 'Locus', 'Sample_Type', 'Treatment'. You
    ## can override using the `.groups` argument.

![](/home/bjcole/AtPst_SingleCell/reports/bulk_analysis_report_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Here are the same plots for the DC3000 treatment:

``` r
ggplot(expr_data, aes(x = Sample_Name, y = log1p(Expression))) + 
  geom_boxplot(aes(fill = de_type_treat)) + 
  scale_fill_discrete(name = "DC3000 DE Gene") +
  theme_bw()
```

![](/home/bjcole/AtPst_SingleCell/reports/bulk_analysis_report_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
expr_data %>%
  group_by(Locus, Sample_Type, Treatment, de_type_treat) %>%
  summarize(mean_expr = mean(Expression)) %>%
  spread(Treatment, mean_expr) %>%
  ggplot(aes(x = log1p(Control), y = log1p(DC3000))) + 
  geom_point(aes(color = de_type_treat), size = 0.5) +
  facet_wrap("Sample_Type")
```

    ## `summarise()` has grouped output by 'Locus', 'Sample_Type', 'Treatment'. You
    ## can override using the `.groups` argument.

![](/home/bjcole/AtPst_SingleCell/reports/bulk_analysis_report_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->
