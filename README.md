# TargetedCN
A R package to infer copy number states from targeted capture data.  Used to call CN changes from targeted sequencing data in [7].

## Using the package

The package can be installed in the standard way using either devtools or install.packages.  Once installed refer to the package documentation for each function for details.  In brief, to use this package you need:
 - BAM files for tumours and a pool of normals captured with the same bait set from which coverage is calculated.
 - A GRanges object (from something like a bed file) detailing the capture regions to associate with each gene.
 - A table of SNPs with counts of reads supporting the reference and alternate allele.

## Brief description of methods

Due to the difficulty in precisely identifying copy number (CN) states from targeted capture data, we instead aim to classify each captured gene as either CN neutral (i.e., having CN equal to the sample ploidy), amplified (i.e., CN gain), and deleted (i.e., CN loss).  For those genes with many captured single nucleotide polymorphisms (SNPs) that were heterozgyous in the matched normal, we further split this classification into genes with homozygous and heterozgyous loss.  To classify genes into these categories we combined data from the coverage and from the biallelic frequencies at SNPs heterozgyous in the normal (henceforth BAFs).

For each gene, we counted the number of 5' fragment ends that fell within the targeted capture region in both the tumour sample of interest and a panel of normals captured using the same capture panel.  Each sample was normalised by the total number of mapping fragments for that sample.  A negative binomial generalised linear model with log-link function was fit for each gene with a factor to indicate tumour/normal status and an offset to account for library size difference.  We optionally apply further normalisation via an additional offset in the model to account for GC bias [1] and composition bias (using a variation on the RNA-seq TMM method [2]).  

We treat the normal panel as replicates and used maximum likelihood estimation to calculate global and gene-specific negative binomial over-dispersion parameters using the empirical Bayes method of [3].  Finally, a likelihood ratio test was performed to compare the model with and without a factor indicating tumour/normal membership for each sample and multiple hypothesis testing was applied using the method of Benjamini and Hochberg [4].  Genes with a corrected p-value of less than 0.01 and a normalised log fold-change in fragment counts less than/greater than 0 were tentatively marked as being deleted/amplified.

Genes with a genuine difference in CN state will also contain a change in the BAFs.  Therefore, we look for additional evidence of the CN gain/loss identified from differences in counts between tumour/normal by looking at the BAF data.  The one exception to this is regions of homozygous deletion where all sequenced reads originate from normal contamination and so show no difference in BAF from CN neutral regions (this makes identifying homozgyously deleted regions more difficult).  

Any location identified as containing a difference from the reference genome in either the tumour or normal at the location of a SNP was identified as a potentially informative SNP.  The only SNPs that contain information about the CN state are those that are heterozygous in the normal.  At each of the potentially informative SNPs we calculated the BAF as the number of reads supporting the alternative allele divided by the number of reads supporting either the alternative or the normal allele.  Any location with reads supporting a third allele were discarded.  Where a matched normal was available we calculated the BAF in the matched normal.  For each remaining SNP we performed a binomial test followed by multiple hypothesis correction under a null hypothesis of 5% deviation from homozygousity (to exclude SNPs with a high error rate or sub-clonal variants).  SNPs with a corrected p-value of less than 0.01 were considered to be heterozygous in the normal and the BAF was calculated from the reads in the tumour.

As non-reference reads are less likely to be captured and to map than those containing the reference sequence, we calculated a correction factor for reference strand bias as the sum of reads in matched normals supporting the alternative allele at SNPs marked as informative divided by those supporting the reference allele. We statistically phase the informative SNPs using IMPUTE [5].  To identify genes for which the allelic ratio is unbalanced (which will be the case for genes with a CN gain/heterozygous loss) we fit a beta-binomial model (following the method of [6]).  The probability of success and over-dispersion parameter were fit by maximising the gene-wide and sample-wide joint likelihood respectively.  Genes with allelic imbalance due to changes in CN were then identified using a likelihood ratio test comparing the beta-binomial model with the maximum likelihood values of the probability of success to a beta-binomial model where the probability of success was set to 0.465 (0.5  times the reference strand bias correction factor).  Any gene with a p-value less than 0.01 after multiple hypothesis correction was marked as  having allelic imbalance.

Generally the BAF data has lower variability and fewer biases than the coverage data and so the most reliable CN calls are for those genes that have a reasonable number of heterozgyous SNPs.  In practice we consider only those genes that have at least 20 informative SNPs when a high sensitivity is required.

[1] - Benjamini Y, Speed TP. Summarizing and correcting the GC content bias in high-throughput sequencing. Nucleic Acids Research. 2012;40(10):e72. doi:10.1093/nar/gks001.

[2] - Robinson MD, Oshlack A 2010. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol 11: R25 10.1186/gb-2010-11-3-r25

[3] -   Chen, Y, Lun, ATL, and Smyth, GK (2014). Differential expression analysis of complex RNA-seq experiments using edgeR. In: Statistical Analysis of Next Generation Sequence Data, Somnath Datta and Daniel S Nettleton (eds), Springer, New York

[4] - Benjamini, Yoav; Hochberg, Yosef (1995). "Controlling the false discovery rate: a practical and powerful approach to multiple testing" (PDF). Journal of the Royal Statistical Society, Series B. 57 (1): 289–300. MR 1325392.

[5] - B. N. Howie, P. Donnelly, and J. Marchini (2009) A flexible and accurate genotype imputation method for the next generation of genome-wide association studies. PLoS Genetics 5(6): e1000529

[6] - Martincorena I, Roshan A, Gerstung M, et al. High burden and pervasive positive selection of somatic mutations in normal human skin. Science (New York, NY). 2015;348(6237):880-886. doi:10.1126/science.aaa6806

[7] - Tarpey PS, Behjati S, Young MD, et al. The driver landscape of sporadic chordoma. Nat Commun. 2017 Oct 12;8(1):890. doi: 10.1038/s41467-017-01026-0.
