# yeasts-dosage

### Selection of strains procedure
- 301 strains from *Peter et. al  2018* => regular diploids, heterozygous, available 
- VCF filtering and annotation: variants splitting, normalization, removing of sites  with low genotyping rate (below 0.1) or with only reference or missing genotypes, annotation with VEP 105 => add script
- Calculation of Hamming (IBS) distances with PLINK v1.90b6.18 64-bit (16 Jun 2020) on the filtered VCF file and other distances to compare [PLINK v2.00a5.6LM 64-bit Intel (29 Oct 2023) --sample-diff, --make-king-table ] 
  => commands in the `run-plink.sh`
 - Further analyses with the `king.R` script => strains selection based on PLINK 1.9 IBS DISTANCES, also some plots of distances:
   * **Cut-off method**   
     Remove strains that have the above distance lower than the 10th distance percentyle to any other strain == 0.01884688.## 106 strains remained => for them I recalculated the lof metrics. 
     Their IDS (one per line) are in the `selected-ids.txt` file. PTV stats summarized in the `summary_per_gene_selected_strains_only.tsv`, `summary_per_strain_selected_strains_only.tsv`, and `summary_per_variant_selected_strains_only.tsv` files.
    * **Selection based on dendogram**  We will use this method!!!   
       Analysis based on the `1011-dist.mdist` and `1011-dist.mdist.id` (matrix from PLINK 1.9)  
         i) generation of dendograms (hclust, three different linking methods);                
         ii) cutting of the upgma ("average") dendogram into 98 groups (with distance from the final nodes smaller than 0.02);   
         iii) picking one strain from each group (saved in `selected-ids-dendogram.txt` **THIS IS THE FINAL STRAIN SET**).   
      Recalculated PTV metrics for these strains are in the `summary_per_gene_dend.tsv` , `summary_per_variant_dend.tsv`, `summary_per_strain_dend.tsv` files.  
      Selected strains shown (red labels) also on the nj tree `nj.pdf` (to ensure that the simple hierarchical clustering works fine for this data).   

 * PTV summaries => `script.R`      


All files and scripts are in the `/mnt/qnap/projects/Yeast_ZGE/distances/` directory.
