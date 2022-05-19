# Amplification is the Primary Mode of Gene-by-Sex Interaction in Complex Human Traits
Carrie Zhu, Matthew J. Ming, Jared M. Cole, Mark Kirkpatrick, Arbel Harpak

Provided below are instructions and details for scripts used to generate the results and figures in ["Amplification is the Primary Mode of Gene-by-Sex Interaction in Complex Human Traits"](https://www.biorxiv.org/content/10.1101/2022.05.06.490973v1).

## Outline:
1. Download GWAS summary statistics files from [UTBox](https://utexas.box.com/s/ef25198jq6owpcq5j2wq6najovlq75b8)
2. Install software: [plink 1.9](https://www.cog-genomics.org/plink/), [plink 2.0](https://www.cog-genomics.org/plink/2.0/), [ldsc](https://github.com/bulik/ldsc), [Ensembl VEP](https://useast.ensembl.org/info/docs/tools/vep/script/index.html), [mashr](https://github.com/stephenslab/mashr)
3. Follow file path outline shown in directory_outline
4. Update configuration file with own file paths: config.R
5. For each section in Documentation, it is best to follow the code in order

## Software
*plink v1.9 beta*
*plink v2.0 alpha*
*LD SCore v1.0.1*
Ensembl command line *variant effect predictor (VEP) v106*
- We used the command line VEP tool to annotate SNPs, following documentation listed on the website
- We downloaded cache files for human genome assembly GRCh37 using INSTALL.pl
- perl module Set::IntervalTree also needs to be installed to use the --nearest flag for VEP  
*mashr* package in R

## Documentation
### General Flags
- ```-p``` or ```--pheno``` flag indicates the phenotype code
- ```-n``` or ```--name``` flag indicates a formated phenotype name, often used for a title of a plot

### Phenotype files
Phenotype files are obtained from UK Biobank and renamed pheno_(phenotype code).txt. A list of phenotype codes and formatted names for labels are provided in pheno_names.txt

### Single snp analysis
##### Miami plots from GWAS summary statistics estimated in males and females only
Download sex-specific summary statistics
Code Example: ```./manhattan.R -p arm_fatfree_mass_L -n "Arm fat-free mass (L)"```
<br> <br/>
##### SNP annotation for list of SNPs after clumping and thresholding, removing SNPs with p-value>5e-8, pairwise LD threshold r<sup>2</sup>>0.1, or within 250kb  
Download 1000G phase 3 genotype data, all_phase3 files, which were created using the following code:  
```plink2 --pfile all_phase3 --chr 1-22 --max-alleles 2 --keep eur_ids.txt --rm-dup exclude-all --king-cutoff 0.0442 --make-bed --out all_phase3```  
eur_ids.txt contained subpopulation codes CEU and GBR on separate lines to be kept in the sample  

Code Example: ```./snp_annotation.sh -p height```  

### LD Score regression
##### Estimate heritability and genetic correlation from sex-specific GWAS summary statistics
Download the pre-computed LD scores from the [ldsc tutorial](https://github.com/bulik/ldsc) (Bulik-Sullivan et al. 2015).  
Code Example: ```./ldsc_basic.sh -p height``` 

#### Create plot for Figure 1
Download ldsc_results.txt, which contains sex-specific heritability estimates and male-female genetic correlations, estimated in the previous step.  
Code: ```./r2_by_h2.R```

### Multivariate adaptive shrinkage (mashr)
The ```-m``` or ```--mode``` flag may be used in the three scripts below. The flag indicates whether or not to produce results specifically for the polygenic score pipeline (**Text S1**), which uses a smaller sample size for the cross-validation procedure. Use the flag and input 'pgs' only if performing the PGS pipeline, otherwise do not use the flag. 

Read in the data to create a matrix of effect estimates and standard errors by sex from sex-specific GWAS summary statistics.  
Code Example: ```./mash_setup.R -p height```

Incorporate pre-specified hypothesis covariance matrices and estimate mixture proportions.  
Code Example: ```./mash_100.R -p height```

Compute posterior estimates: posterior mean, standard deviation, weight, and lfsr.  
Code Example: ```./mash_posterior.R -p height```

#### Plots using mixture weights
Create a overall and compact heatmaps of mixture weights. Plots for **Fig. S4**.   
Code Example: ```./mash_heatmap.R -p height -n Height```

Plot for **Fig. 4A**.  Download pheno_names.txt and sex_ids.txt.  
Code: ```./phenovar_by_phenomean.R```

Plot for **Fig. 4B**. Download mash_weights.txt, which summarizes weights from all traits.  
Code: ```./phenovar_by_amplification.R```

Plot for **Fig. S7**  
Code: ```nontrivial.R```

#### Test different p-value threshold (Methods and Fig. S3)
The ```-s``` or ```--same``` flag may be used in the three scripts below. To keep the same random sample size for input to *mash*, use the flag, with the parameter '_same'. Otherwise, do not use the flag. 

```./mash_setup.R``` needs to be run first for the same trait.  
Code Example: ```mash_p_threshold.R -p height```  

Plots for **Fig. S3B**.  
Code Example: ```mash_pvalue_plot.R -p height -n Height```  

Plots for **Fig. S3A**. Download noeffect_weight.txt and noeffect_weight_same.txt, which has the weight on the no effect matrix for each phenotype and p-value threshold used in this particular analysis.  
Code Example: ```mash_pvalue_null_plot.R```  

### mash simulations
environ_matrix.R
    # need maf_sample_20k.txt in QC dir, output to GWAS dir
environ_small.R
environ_large.R
environ_mash.R
environ_heatmap.R

### PGS 
PGS_testset_1.R
PGS_GWAS_2.R
PGS_CT_SCORE_4.R
PGS_predict_5_linear.R

PGS_plot_final.R
    # need pgs_linear_results_five.txt, output pgs_combined_r2.txt
pheno_pgs.R
    # directions to move best pgs file to gwas/pheno directory    #!
pheno_pgs_overall.$
    # need sexspecific_pheno_pgs_lm.txt from previous

### testosterone as underlier
G_testosterone.R
G_testosterone_pgs.R
G_corr_testosterone.R
G_corr_testosterone_pgs.R
G_corr_testosterone_age.R

### shared amplification
gen_env_bootstrap.R 
    # need ldsc_results.txt in LDSC directory

### selection

