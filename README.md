# From Local Stressors to Regional Sources: Eutrophication and Connectivity Shape a Coastal Metacommunity


#### author list hidden for double-anonymised review 

## Abstract 

<img align="right" src="2_figure/figure2forshow.png" width=450> 

Coastal urbanization exposes marine ecosystems to nutrient pollution by increasing wastewater and urban runoff. Yet, how eutrophication shapes benthic metacommunities over the complex hydrological network inherent to coastal seas - remains poorly understood. Here, we integrated 6 years of biodiversity data from 88 standardized biodiversity samplers (ARMS) with high-resolution water quality records across a subtropical eutrophic urban coast. Chlorophyll-a concentration—a proxy for nutrient enrichment—was strongly and negatively correlated with normalized benthic richness (GLM, p = 1.58e-06), explaining 23.3% of compositional variation along the primary ordination axis (envfit, p = 0.001). This environmental filter systematically restructured communities: relative richness of benthic primary producers (Bacillariophyta, Rhodophyta) increased with total OTU richness (R² = 0.21), while that of Annelida declined sharply (R² = 0.60), driving a concurrent rise in evenness (R² = 0.72). Network analysis revealed that spatial connectivity further shaped metacommunity structure. Communities clustered into four distinct modules with 93.2% concordance to predefined hydrological regions. A highly connected 'biodiversity corridor' in the south and east, characterized by low nutrient levels and high network centrality, harbored all putative source communities. In contrast, isolated sinks in the west and Tolo Harbor exhibited low richness and connectivity. Our findings demonstrate that eutrophication and spatial connectivity jointly structure coastal metacommunities, identifying source corridors critical for conservation in an era of rapid coastal change.





## Table of Contents

### Supporting Materials 
  1. [Raw sequence](https://doi.org/10.6084/m9.figshare.29481053) 
  2. [Data](3_data)
  3. [Figures](2_figure/260304_mergeFigure.pdf)
  4. [Tables](2_figure/260304_mergeTable.pdf)
  5. [Supplementary Materials](2_figure/260305_supplementaryMaterials_FIN.pdf)

### Sequence processing pipeline 
1. [Import & cutadap](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.1_importAndCutAdapt.sh): import raw sequence data (.fastq) into Qiime artefacts (.qza) and remove PCR adaptors.
2. [Denoise-paired](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.2_denoiseAndPair.sh): remove sequences likely induced by error and merge the reverse/forward reads.
3. [Decontam](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.3_decontam.r): a process to look into the negative control and remove sequences that might have come from sample contamination.
4. [Amino Acid translation](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.4_aaTranslate.r): translate DNA sequence into amino acid and remove sequences with one of the following conditions: 1) any STOP codon, 2) >3 deletion, 3) any frameshift, 4) any insertion.
5. [Cluster all sequences](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.5_clusterReads.sh) by 97% similarity into operational taxonomic units (OTUs) for downstream data analysis.
6. [Taxonomic assignment](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.6_taxAssign.sh) with BLAST against two different libraries: 1) McIlroy et al. 2024 & 2) Medori2 (GB260).

### Data Analysis 
1. Environmental data
   - [Heatmap](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.1_eData_heatmap.r) (Figure 1d, Table 1)
   - [MPA east vs west](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.2_eastVSwest.r) (Table S2)
2. Species richness by ARMS 
   - [Merge richness from all three fractions](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.3_combinFractionbyARMS.r) (Table S1)
   - [Environmental data ~ species richness](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.4_eDATAvsRichness.r) (Table 2) 
3. Community composition
   - [PCoA](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.5_PCoA.r) (Figure 2)
   - [Permutational Multivariate Analysis of Variance (adonis2)](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.6_adonis2.r) 
   - [Diverging Bar Chart & Chi-Square analysis](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.8_sidewayBar.r) (Figure 3)

     
