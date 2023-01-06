-----

## Title

Differential Expression

## Methods:

The two main steps in performing differential expression analysis are to
estimate the relative abundance of transcripts, and to apply statistical
models to test for differential expression between treatment groups.
Estimating relative abundance is basically determining how many NGS
reads match a given gene within a genome. In this module we have used
Salmon (Patro et al. 2017) to estimate relative abundance, tximport
(Soneson, Love, and Robinson 2015) to import the Salmon abundance
estimates, and DESeq2 (Love, Huber, and Anders 2014) to perform
statistical tests to identify differentially expressed genes. Salmon is
a method for quantifying transcript abundance from RNA-seq reads that is
accurate and fast First we built a salmon index from the de-novo
transcriptome. We did the Salmon abundance estimation all the samples,
stored it in quant directory and aligned them to AipIndex using Salmon.
To use Salmon input in tximport, we have to build a table mapping
transcripts of gene. Next we have written R script to import Salmon
alignments into DESeq2 and performed differential expression analysis.An
important task here in analysis of RNA sequencing (RNA-seq) data is the
aim of finding genes that are differentially expressed across groups of
samples. DESeq2 performs this differential analysis. Finally we have
filtered the genes and made a method section that has the differentially
expressed genes in a kable-formatted
table.

Results:

| ko        | pathway      | path                                             |         log |      padj | factor                        |
| :-------- | :----------- | :----------------------------------------------- | ----------: | --------: | :---------------------------- |
| ko:K00333 | path:ko01100 | Metabolic pathways                               | \-0.9083018 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K00333 | path:ko00190 | Oxidative phosphorylation                        | \-0.9083018 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K00522 | path:ko04978 | Mineral absorption                               | \-0.7371782 | 0.3515705 | Menthol\_Menthol\_vs\_Control |
| ko:K00522 | path:ko04217 | Necroptosis                                      | \-0.7371782 | 0.3515705 | Menthol\_Menthol\_vs\_Control |
| ko:K00522 | path:ko04216 | Ferroptosis                                      | \-0.7371782 | 0.3515705 | Menthol\_Menthol\_vs\_Control |
| ko:K00819 | path:ko00330 | Arginine and proline metabolism                  |   0.7144250 | 0.2530751 | Menthol\_Menthol\_vs\_Control |
| ko:K00819 | path:ko01110 | Biosynthesis of secondary metabolites            |   0.7144250 | 0.2530751 | Menthol\_Menthol\_vs\_Control |
| ko:K00819 | path:ko01130 | Biosynthesis of antibiotics                      |   0.7144250 | 0.2530751 | Menthol\_Menthol\_vs\_Control |
| ko:K00819 | path:ko01100 | Metabolic pathways                               |   0.7144250 | 0.2530751 | Menthol\_Menthol\_vs\_Control |
| ko:K00967 | path:ko01100 | Metabolic pathways                               |   0.3491718 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K00967 | path:ko00564 | Glycerophospholipid metabolism                   |   0.3491718 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K00967 | path:ko00440 | Phosphonate and phosphinate metabolism           |   0.3491718 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K01251 | path:ko01100 | Metabolic pathways                               |   1.2977136 | 0.1471194 | Menthol\_Menthol\_vs\_Control |
| ko:K01251 | path:ko00270 | Cysteine and methionine metabolism               |   1.2977136 | 0.1471194 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko01130 | Biosynthesis of antibiotics                      |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko04922 | Glucagon signaling pathway                       |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko01110 | Biosynthesis of secondary metabolites            |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko01120 | Microbial metabolism in diverse environments     |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko01100 | Metabolic pathways                               |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko00010 | Glycolysis / Gluconeogenesis                     |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko00620 | Pyruvate metabolism                              |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko04964 | Proximal tubule bicarbonate reclamation          |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko03320 | PPAR signaling pathway                           |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko04920 | Adipocytokine signaling pathway                  |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko04151 | PI3K-Akt signaling pathway                       |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko04152 | AMPK signaling pathway                           |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko04931 | Insulin resistance                               |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko00020 | Citrate cycle (TCA cycle)                        |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko04068 | FoxO signaling pathway                           |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01596 | path:ko04910 | Insulin signaling pathway                        |   0.7606325 | 0.0955962 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko01130 | Biosynthesis of antibiotics                      | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko01120 | Microbial metabolism in diverse environments     | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko01110 | Biosynthesis of secondary metabolites            | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko00680 | Methane metabolism                               | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko01230 | Biosynthesis of amino acids                      | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko00710 | Carbon fixation in photosynthetic organisms      | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko00051 | Fructose and mannose metabolism                  | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko00010 | Glycolysis / Gluconeogenesis                     | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko01100 | Metabolic pathways                               | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko01200 | Carbon metabolism                                | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01623 | path:ko00030 | Pentose phosphate pathway                        | \-1.0362436 | 0.1870460 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko01230 | Biosynthesis of amino acids                      |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko00680 | Methane metabolism                               |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko01110 | Biosynthesis of secondary metabolites            |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko03018 | RNA degradation                                  |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko01200 | Carbon metabolism                                |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko01100 | Metabolic pathways                               |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko00010 | Glycolysis / Gluconeogenesis                     |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko01130 | Biosynthesis of antibiotics                      |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko04066 | HIF-1 signaling pathway                          |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K01689 | path:ko01120 | Microbial metabolism in diverse environments     |   0.8458592 | 0.3312406 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04722 | Neurotrophin signaling pathway                   |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04745 | Phototransduction - fly                          |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04015 | Rap1 signaling pathway                           |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04713 | Circadian entrainment                            |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05133 | Pertussis                                        |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04915 | Estrogen signaling pathway                       |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04925 | Aldosterone synthesis and secretion              |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04626 | Plant-pathogen interaction                       |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04022 | cGMP-PKG signaling pathway                       |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04261 | Adrenergic signaling in cardiomyocytes           |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05200 | Pathways in cancer                               |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04720 | Long-term potentiation                           |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04921 | Oxytocin signaling pathway                       |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04114 | Oocyte meiosis                                   |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05170 | Human immunodeficiency virus 1 infection         |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04728 | Dopaminergic synapse                             |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04014 | Ras signaling pathway                            |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04218 | Cellular senescence                              |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04016 | MAPK signaling pathway - plant                   |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05167 | Kaposi sarcoma-associated herpesvirus infection  |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04744 | Phototransduction                                |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04922 | Glucagon signaling pathway                       |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04912 | GnRH signaling pathway                           |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05031 | Amphetamine addiction                            |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04270 | Vascular smooth muscle contraction               |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05034 | Alcoholism                                       |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05418 | Fluid shear stress and atherosclerosis           |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05152 | Tuberculosis                                     |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04371 | Apelin signaling pathway                         |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04625 | C-type lectin receptor signaling pathway         |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05010 | Alzheimer disease                                |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04971 | Gastric acid secretion                           |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05163 | Human cytomegalovirus infection                  |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04070 | Phosphatidylinositol signaling system            |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04924 | Renin secretion                                  |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04740 | Olfactory transduction                           |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04910 | Insulin signaling pathway                        |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04024 | cAMP signaling pathway                           |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04916 | Melanogenesis                                    |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04020 | Calcium signaling pathway                        |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04970 | Salivary secretion                               |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko04750 | Inflammatory mediator regulation of TRP channels |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02183 | path:ko05214 | Glioma                                           |   0.4003852 | 0.4258071 | Menthol\_Menthol\_vs\_Control |
| ko:K02980 | path:ko03010 | Ribosome                                         | \-1.1641770 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K02993 | path:ko03010 | Ribosome                                         | \-0.7752986 | 0.2530751 | Menthol\_Menthol\_vs\_Control |
| ko:K03236 | path:ko03013 | RNA transport                                    | \-0.7910526 | 0.4862660 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko04728 | Dopaminergic synapse                             | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko04111 | Cell cycle - yeast                               | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko04151 | PI3K-Akt signaling pathway                       | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko04391 | Hippo signaling pathway - fly                    | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko04530 | Tight junction                                   | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko04071 | Sphingolipid signaling pathway                   | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko05160 | Hepatitis C                                      | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko05142 | Chagas disease (American trypanosomiasis)        | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko04152 | AMPK signaling pathway                           | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko04390 | Hippo signaling pathway                          | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko03015 | mRNA surveillance pathway                        | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko04261 | Adrenergic signaling in cardiomyocytes           | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K04354 | path:ko05165 | Human papillomavirus infection                   | \-2.0689098 | 0.0272952 | Menthol\_Menthol\_vs\_Control |
| ko:K05759 | path:ko04015 | Rap1 signaling pathway                           | \-0.7621169 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K05759 | path:ko04013 | MAPK signaling pathway - fly                     | \-0.7621169 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K05759 | path:ko04810 | Regulation of actin cytoskeleton                 | \-0.7621169 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K05759 | path:ko05131 | Shigellosis                                      | \-0.7621169 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K05759 | path:ko05132 | Salmonella infection                             | \-0.7621169 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K11251 | path:ko05034 | Alcoholism                                       |   1.0680166 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K11251 | path:ko05322 | Systemic lupus erythematosus                     |   1.0680166 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K11251 | path:ko04217 | Necroptosis                                      |   1.0680166 | 0.4836718 | Menthol\_Menthol\_vs\_Control |
| ko:K16186 | path:ko04140 | Autophagy - animal                               |   1.6274960 | 0.1461510 | Menthol\_Menthol\_vs\_Control |
| ko:K16186 | path:ko04150 | mTOR signaling pathway                           |   1.6274960 | 0.1461510 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko04110 | Cell cycle                                       | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko04013 | MAPK signaling pathway - fly                     | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko04151 | PI3K-Akt signaling pathway                       | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko04212 | Longevity regulating pathway - worm              | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko04390 | Hippo signaling pathway                          | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko05203 | Viral carcinogenesis                             | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko05160 | Hepatitis C                                      | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko04114 | Oocyte meiosis                                   | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko05161 | Hepatitis B                                      | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K16197 | path:ko04391 | Hippo signaling pathway - fly                    | \-1.3002193 | 0.1939680 | Menthol\_Menthol\_vs\_Control |
| ko:K17920 | path:ko04144 | Endocytosis                                      |   1.1188034 | 0.0000011 | Menthol\_Menthol\_vs\_Control |

## References

<div id="refs" class="references">

<div id="ref-Love">

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12): 550–50.
<https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Patro">

Patro, Rob, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and Carl
Kingsford. 2017. “Salmon Provides Fast and Bias-Aware Quantification of
Transcript Expression.” *Nature Methods* 14 (4): 417–19.
<https://doi.org/10.1038/nmeth.4197>.

</div>

<div id="ref-Soneson">

Soneson, Charlotte, Michael I. Love, and Mark D. Robinson. 2015.
“Differential Analyses for RNA-Seq: Transcript-Level Estimates Improve
Gene-Level Inferences.” *F1000Research* 4 (December): 1521–1.
<https://www.ncbi.nlm.nih.gov/pubmed/26925227>.

</div>

</div>
