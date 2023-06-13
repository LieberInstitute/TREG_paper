library("here")
library("sessioninfo")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")
library("writexl")

## Plot output directories
dir_plots <- here::here(
    "plots",
    "04_gene_ontology"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir_rdata <- here::here(
    "processed-data",
    "04_gene_ontology"
)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Load ranked candidate TREGs
load(here("processed-data", "01_find_tregs", "rank_invar.Rdata"), verbose = TRUE)

## Get ENTREZID
sigGene <- lapply(c(20, 50), function(k) { 
    bitr(head(rank_invar_df$ensembl_id, n = k), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
})
names(sigGene) <- c(20, 50)
lengths(sigGene)

geneUniverse <- unique(bitr(rownames(gene_propZero), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID)
length(geneUniverse)
# [1] 10981

## The following GO code was adapted from
## https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/05_GO_enrichment/01_GO_enrichment.R

## Run GO and KEGG enrichment analysis
go <- compareCluster(sigGene,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = "org.Hs.eg.db",
    ont = "ALL",
    pAdjustMethod = "BH",
    readable = TRUE
)

kegg <- compareCluster(sigGene,
    fun = "enrichKEGG",
    universe = geneUniverse,
    organism = "hsa"
)

## Save Rdata with gene ontology enrichment results
save(go, kegg, file = file.path(dir_rdata,
    "gene_ontology_enrichment_objects.Rdata"
))

## Plot BP, CC and MF
plot_go <- function(ont, title_p, path, filename) {
    dotplot_1 <- ggplot(filter(go, ONTOLOGY == ont), aes(Cluster, Description)) +
        theme_bw() +
        geom_point(aes(color = p.adjust, size = Count)) +
        scale_color_gradientn(
            colours = c("#f7ca64", "#46bac2", "#7e62a3"),
            trans = "log10",
            guide = guide_colorbar(reverse = TRUE, order = 1),
            name = "FDR"
        ) +
        scale_size_continuous(range = c(2, 10)) +
        xlab("Cluster") +
        ylab("") +
        ggtitle(title_p)

    ggsave(filename = filename, path = path, dotplot_1, height = 6, width = 5)
}

plot_go(ont = "BP", title_p = "Biological Process", filename = "GOenrichment_BP.pdf", path = dir_plots)
plot_go(ont = "CC", title_p = "Cellular Component", filename = "GOenrichment_CC.pdf", path = dir_plots)
plot_go(ont = "MF", title_p = "Molecular Function", filename = "GOenrichment_MF.pdf", path = dir_plots)

## Plot KEGG
dotplot_1 <- ggplot(kegg, aes(Cluster, Description)) +
    theme_bw() +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_gradientn(
        colours = c("#f7ca64", "#46bac2", "#7e62a3"),
        trans = "log10",
        guide = guide_colorbar(reverse = TRUE, order = 1),
        name = "FDR"
    ) +
    scale_size_continuous(range = c(2, 10)) +
    xlab("Cluster") +
    ylab("") +
    ggtitle("KEGG")

ggsave(filename = "KEGGenrichment.pdf", path = dir_plots, dotplot_1, height = 6, width = 5)

## Print some result
as.data.frame(go)[, -10]
#    Cluster ONTOLOGY         ID                        Description GeneRatio  BgRatio       pvalue     p.adjust       qvalue Count
# 1       20       BP GO:0006325             chromatin organization      7/18 455/9460 1.139519e-05 7.691753e-03 6.381306e-03     7
# 2       20       MF GO:0003712 transcription coregulator activity      5/16 374/9777 2.456529e-04 1.066311e-02 7.672322e-03     5
# 3       20       MF GO:0003713 transcription coactivator activity      4/16 203/9777 2.699521e-04 1.066311e-02 7.672322e-03     4
# 4       50       BP GO:0006325             chromatin organization     12/47 455/9460 1.488362e-06 1.705663e-03 1.655999e-03    12
# 5       50       BP GO:0051568          histone H3-K4 methylation      5/47  50/9460 4.361377e-06 2.499069e-03 2.426303e-03     5
# 6       50       BP GO:0016570               histone modification     10/47 357/9460 7.678988e-06 2.599430e-03 2.523742e-03    10
# 7       50       BP GO:0016571                histone methylation      6/47 100/9460 9.073055e-06 2.599430e-03 2.523742e-03     6
# 8       50       BP GO:0006479                protein methylation      6/47 132/9460 4.422845e-05 7.168094e-03 6.959379e-03     6
# 9       50       BP GO:0008213                 protein alkylation      6/47 132/9460 4.422845e-05 7.168094e-03 6.959379e-03     6
# 10      50       BP GO:0034968         histone lysine methylation      5/47  81/9460 4.701925e-05 7.168094e-03 6.959379e-03     5
# 11      50       BP GO:0018205       peptidyl-lysine modification      8/47 272/9460 5.003905e-05 7.168094e-03 6.959379e-03     8
# 12      50       BP GO:0018022        peptidyl-lysine methylation      5/47  90/9460 7.800975e-05 9.933242e-03 9.644013e-03     5
# 13      50       CC GO:0035097  histone methyltransferase complex      4/46  51/9829 8.934249e-05 1.366940e-02 1.269604e-02     4
# 14      50       CC GO:0034708          methyltransferase complex      4/46  75/9829 4.003629e-04 2.738343e-02 2.543352e-02     4
# 15      50       CC GO:0070603   SWI/SNF superfamily-type complex      4/46  81/9829 5.369299e-04 2.738343e-02 2.543352e-02     4
# 16      50       CC GO:0016607                      nuclear speck      7/46 340/9829 9.320704e-04 3.510826e-02 3.260829e-02     7
# 17      50       CC GO:0044666                     MLL3/4 complex      2/46  11/9829 1.147329e-03 3.510826e-02 3.260829e-02     2
# 18      50       MF GO:0003712 transcription coregulator activity     12/46 374/9777 9.874611e-08 1.757681e-05 1.507177e-05    12
# 19      50       MF GO:0003713 transcription coactivator activity      7/46 203/9777 4.039527e-05 3.595179e-03 3.082797e-03     7
# 20      50       MF GO:0003682                  chromatin binding      8/46 383/9777 3.638758e-04 2.158996e-02 1.851298e-02     8

sort(table(as.data.frame(go)$geneID), decreasing = TRUE)
#                                                 KANSL1/KMT2C/RLF/ASH1L/NCOA6                                           KANSL1/KMT2C/RLF/ASH1L/NCOA6/SUZ12 
#                                                                            3                                                                            3 
#                                                     KANSL1/KMT2C/NCOA6/SUZ12                                                      ARID1B/RERE/KMT2C/EP300 
#                                                                            2                                                                            1 
#                                    ARID1B/RERE/KMT2C/EP300/DDX17/CCAR1/NCOA6                                                     ARID1B/SF3B1/TRRAP/SUZ12 
#                                                                            1                                                                            1 
#                            JMJD1C/ARID1B/RERE/EP300/ASH1L/NCOA6/SUZ12/NAP1L4                                  JMJD1C/ARID1B/RERE/KANSL1/SF3B1/KMT2C/EP300 
#                                                                            1                                                                            1 
#    JMJD1C/ARID1B/RERE/KANSL1/SF3B1/KMT2C/EP300/TRRAP/ASH1L/TLK2/SUZ12/NAP1L4                                               JMJD1C/ARID1B/RERE/KMT2C/EP300 
#                                                                            1                                                                            1 
# JMJD1C/ARID1B/RERE/KMT2C/EP300/MED13L/DDX17/TRRAP/CCAR1/N4BP2L2/NCOA6/BCLAF1                  JMJD1C/KANSL1/SF3B1/KMT2C/EP300/TRRAP/RLF/ASH1L/NCOA6/SUZ12 
#                                                                            1                                                                            1 
#                               KANSL1/SF3B1/KMT2C/EP300/TRRAP/RLF/ASH1L/NCOA6                                                                  KMT2C/NCOA6 
#                                                                            1                                                                            1 
#                                  SF3B1/THOC2/DDX17/TCF12/SRSF4/PRPF4B/BCLAF1 
#                                                                            1

## I do see ARID1B in several ones here, but not AKT3.

as.data.frame(kegg)[, -9]
#    Cluster       ID                                              Description GeneRatio  BgRatio       pvalue   p.adjust      qvalue Count
# 1       20 hsa05211                                     Renal cell carcinoma      3/10  58/4501 0.0002285934 0.01409436 0.006390966     3
# 2       20 hsa04630                               JAK-STAT signaling pathway      3/10  65/4501 0.0003209489 0.01409436 0.006390966     3
# 3       20 hsa05215                                          Prostate cancer      3/10  66/4501 0.0003358373 0.01409436 0.006390966     3
# 4       20 hsa05164                                              Influenza A      3/10  80/4501 0.0005931352 0.01409436 0.006390966     3
# 5       20 hsa04935           Growth hormone synthesis, secretion and action      3/10  84/4501 0.0006846501 0.01409436 0.006390966     3
# 6       20 hsa05152                                             Tuberculosis      3/10  86/4501 0.0007336206 0.01409436 0.006390966     3
# 7       20 hsa04068                                   FoxO signaling pathway      3/10  87/4501 0.0007589272 0.01409436 0.006390966     3
# 8       20 hsa04919                        Thyroid hormone signaling pathway      3/10  92/4501 0.0008938768 0.01452550 0.006586461     3
# 9       20 hsa05161                                              Hepatitis B      3/10 102/4501 0.0012078379 0.01744655 0.007910985     3
# 10      20 hsa05225                                 Hepatocellular carcinoma      3/10 109/4501 0.0014646539 0.01904050 0.008633749     3
# 11      20 hsa05167          Kaposi sarcoma-associated herpesvirus infection      3/10 115/4501 0.0017104341 0.02021422 0.009165963     3
# 12      20 hsa04024                                   cAMP signaling pathway      3/10 132/4501 0.0025438388 0.02717131 0.012320596     3
# 13      20 hsa04370                                   VEGF signaling pathway      2/10  38/4501 0.0029931309 0.02717131 0.012320596     2
# 14      20 hsa04664                          Fc epsilon RI signaling pathway      2/10  43/4501 0.0038219132 0.02717131 0.012320596     2
# 15      20 hsa04330                                  Notch signaling pathway      2/10  44/4501 0.0039991595 0.02717131 0.012320596     2
# 16      20 hsa04917                              Prolactin signaling pathway      2/10  44/4501 0.0039991595 0.02717131 0.012320596     2
# 17      20 hsa05213                                       Endometrial cancer      2/10  44/4501 0.0039991595 0.02717131 0.012320596     2
# 18      20 hsa05218                                                 Melanoma      2/10  44/4501 0.0039991595 0.02717131 0.012320596     2
# 19      20 hsa05221                                   Acute myeloid leukemia      2/10  44/4501 0.0039991595 0.02717131 0.012320596     2
# 20      20 hsa04929                                           GnRH secretion      2/10  45/4501 0.0041802022 0.02717131 0.012320596     2
# 21      20 hsa05230                      Central carbon metabolism in cancer      2/10  48/4501 0.0047459666 0.02920304 0.013241864     2
# 22      20 hsa04662                        B cell receptor signaling pathway      2/10  49/4501 0.0049420527 0.02920304 0.013241864     2
# 23      20 hsa04625                 C-type lectin receptor signaling pathway      2/10  55/4501 0.0061963113 0.03084564 0.013986689     2
# 24      20 hsa04720                                   Long-term potentiation      2/10  55/4501 0.0061963113 0.03084564 0.013986689     2
# 25      20 hsa05223                               Non-small cell lung cancer      2/10  55/4501 0.0061963113 0.03084564 0.013986689     2
# 26      20 hsa04916                                            Melanogenesis      2/10  56/4501 0.0064181805 0.03084564 0.013986689     2
# 27      20 hsa05235   PD-L1 expression and PD-1 checkpoint pathway in cancer      2/10  56/4501 0.0064181805 0.03084564 0.013986689     2
# 28      20 hsa05214                                                   Glioma      2/10  57/4501 0.0066436771 0.03084564 0.013986689     2
# 29      20 hsa05212                                        Pancreatic cancer      2/10  59/4501 0.0071054971 0.03185223 0.014443116     2
# 30      20 hsa05220                                 Chronic myeloid leukemia      2/10  61/4501 0.0075816607 0.03279280 0.014869611     2
# 31      20 hsa05165                           Human papillomavirus infection      3/10 200/4501 0.0082264373 0.03279280 0.014869611     3
# 32      20 hsa04914                  Progesterone-mediated oocyte maturation      2/10  64/4501 0.0083225599 0.03279280 0.014869611     2
# 33      20 hsa04660                        T cell receptor signaling pathway      2/10  65/4501 0.0085765793 0.03279280 0.014869611     2
# 34      20 hsa05210                                        Colorectal cancer      2/10  65/4501 0.0085765793 0.03279280 0.014869611     2
# 35      20 hsa01521                EGFR tyrosine kinase inhibitor resistance      2/10  67/4501 0.0090951159 0.03284347 0.014892587     2
# 36      20 hsa04613                  Neutrophil extracellular trap formation      2/10  67/4501 0.0090951159 0.03284347 0.014892587     2
# 37      20 hsa04922                               Glucagon signaling pathway      2/10  68/4501 0.0093596060 0.03288510 0.014911463     2
# 38      20 hsa04012                                   ErbB signaling pathway      2/10  69/4501 0.0096275593 0.03293639 0.014934718     2
# 39      20 hsa01522                                     Endocrine resistance      2/10  71/4501 0.0101738015 0.03391267 0.015377406     2
# 40      20 hsa04066                                  HIF-1 signaling pathway      2/10  74/4501 0.0110188029 0.03410582 0.015464987     2
# 41      20 hsa04550 Signaling pathways regulating pluripotency of stem cells      2/10  74/4501 0.0110188029 0.03410582 0.015464987     2
# 42      20 hsa04666                         Fc gamma R-mediated phagocytosis      2/10  74/4501 0.0110188029 0.03410582 0.015464987     2
# 43      20 hsa05224                                            Breast cancer      2/10  77/4501 0.0118942501 0.03514210 0.015934880     2
# 44      20 hsa05226                                           Gastric cancer      2/10  77/4501 0.0118942501 0.03514210 0.015934880     2
# 45      20 hsa05231                             Choline metabolism in cancer      2/10  78/4501 0.0121927696 0.03522356 0.015971815     2
# 46      20 hsa04915                               Estrogen signaling pathway      2/10  79/4501 0.0124946191 0.03531088 0.016011411     2
# 47      20 hsa04210                                                Apoptosis      2/10  83/4501 0.0137350538 0.03799057 0.017226495     2
# 48      20 hsa04926                                Relaxin signaling pathway      2/10  84/4501 0.0140533563 0.03806117 0.017258508     2
# 49      20 hsa04071                           Sphingolipid signaling pathway      2/10  87/4501 0.0150277212 0.03986946 0.018078462     2
# 50      20 hsa04371                                 Apelin signaling pathway      2/10  89/4501 0.0156933831 0.04080280 0.018501673     2
# 51      20 hsa05160                                              Hepatitis C      2/10  91/4501 0.0163717933 0.04173202 0.018923022     2
# 52      20 hsa04722                           Neurotrophin signaling pathway      2/10  92/4501 0.0167157470 0.04178937 0.018949025     2
# 53      20 hsa04910                                Insulin signaling pathway      2/10 100/4501 0.0195796687 0.04718624 0.021396189     2
# 54      20 hsa04062                              Chemokine signaling pathway      2/10 101/4501 0.0199515066 0.04718624 0.021396189     2
# 55      20 hsa04022                               cGMP-PKG signaling pathway      2/10 102/4501 0.0203263799 0.04718624 0.021396189     2
# 56      20 hsa04218                                      Cellular senescence      2/10 102/4501 0.0203263799 0.04718624 0.021396189     2
# 57      20 hsa04072                        Phospholipase D signaling pathway      2/10 106/4501 0.0218559809 0.04984697 0.022602676     2
# 58      50 hsa04935           Growth hormone synthesis, secretion and action      4/20  84/4501 0.0004355245 0.04715475 0.027251455     4
# 59      50 hsa04919                        Thyroid hormone signaling pathway      4/20  92/4501 0.0006164020 0.04715475 0.027251455     4
# 60      50 hsa04929                                           GnRH secretion      3/20  45/4501 0.0009454515 0.04821803 0.027865939     3


sort(table(as.data.frame(kegg)$geneID), decreasing = TRUE)


#      10000/5894       10000/5894/2033            10000/2033             5894/2033      10000/57492/5894 10000/5894/2033/23389  10000/5894/2033/2776 
#              39                    12                     2                     2                     1                     1                     1 
# 10000/5894/2776             6310/2033 
#               1                     1 

subset(rank_invar_df, Symbol == "AKT3")
#   Symbol      ensembl_id rank_invar
# 5   AKT3 ENSG00000117020        873

bitr("ENSG00000117020", fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
#           ENSEMBL ENTREZID
# 1 ENSG00000117020    10000

bitr("10000", fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
#  ENTREZID SYMBOL
# 1    10000   AKT3
bitr("5894", fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
#   ENTREZID SYMBOL
# 1     5894   RAF1
bitr("2033", fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
#   ENTREZID SYMBOL
# 1     2033  EP300

## Below we see the full results Note that the defaults for
## enrichGO() and enrichKEGG() are:
## * pvalueCutoff = 0.05
## * qvalueCutoff = 0.2
options(width = 300)
as.data.frame(go)
#    Cluster ONTOLOGY         ID                        Description GeneRatio  BgRatio       pvalue     p.adjust       qvalue                                                                       geneID Count
# 1       20       BP GO:0006325             chromatin organization      7/18 455/9460 1.139519e-05 7.691753e-03 6.381306e-03                                  JMJD1C/ARID1B/RERE/KANSL1/SF3B1/KMT2C/EP300     7
# 2       20       MF GO:0003712 transcription coregulator activity      5/16 374/9777 2.456529e-04 1.066311e-02 7.672322e-03                                               JMJD1C/ARID1B/RERE/KMT2C/EP300     5
# 3       20       MF GO:0003713 transcription coactivator activity      4/16 203/9777 2.699521e-04 1.066311e-02 7.672322e-03                                                      ARID1B/RERE/KMT2C/EP300     4
# 4       50       BP GO:0006325             chromatin organization     12/47 455/9460 1.488362e-06 1.705663e-03 1.655999e-03    JMJD1C/ARID1B/RERE/KANSL1/SF3B1/KMT2C/EP300/TRRAP/ASH1L/TLK2/SUZ12/NAP1L4    12
# 5       50       BP GO:0051568          histone H3-K4 methylation      5/47  50/9460 4.361377e-06 2.499069e-03 2.426303e-03                                                 KANSL1/KMT2C/RLF/ASH1L/NCOA6     5
# 6       50       BP GO:0016570               histone modification     10/47 357/9460 7.678988e-06 2.599430e-03 2.523742e-03                  JMJD1C/KANSL1/SF3B1/KMT2C/EP300/TRRAP/RLF/ASH1L/NCOA6/SUZ12    10
# 7       50       BP GO:0016571                histone methylation      6/47 100/9460 9.073055e-06 2.599430e-03 2.523742e-03                                           KANSL1/KMT2C/RLF/ASH1L/NCOA6/SUZ12     6
# 8       50       BP GO:0006479                protein methylation      6/47 132/9460 4.422845e-05 7.168094e-03 6.959379e-03                                           KANSL1/KMT2C/RLF/ASH1L/NCOA6/SUZ12     6
# 9       50       BP GO:0008213                 protein alkylation      6/47 132/9460 4.422845e-05 7.168094e-03 6.959379e-03                                           KANSL1/KMT2C/RLF/ASH1L/NCOA6/SUZ12     6
# 10      50       BP GO:0034968         histone lysine methylation      5/47  81/9460 4.701925e-05 7.168094e-03 6.959379e-03                                                 KANSL1/KMT2C/RLF/ASH1L/NCOA6     5
# 11      50       BP GO:0018205       peptidyl-lysine modification      8/47 272/9460 5.003905e-05 7.168094e-03 6.959379e-03                               KANSL1/SF3B1/KMT2C/EP300/TRRAP/RLF/ASH1L/NCOA6     8
# 12      50       BP GO:0018022        peptidyl-lysine methylation      5/47  90/9460 7.800975e-05 9.933242e-03 9.644013e-03                                                 KANSL1/KMT2C/RLF/ASH1L/NCOA6     5
# 13      50       CC GO:0035097  histone methyltransferase complex      4/46  51/9829 8.934249e-05 1.366940e-02 1.269604e-02                                                     KANSL1/KMT2C/NCOA6/SUZ12     4
# 14      50       CC GO:0034708          methyltransferase complex      4/46  75/9829 4.003629e-04 2.738343e-02 2.543352e-02                                                     KANSL1/KMT2C/NCOA6/SUZ12     4
# 15      50       CC GO:0070603   SWI/SNF superfamily-type complex      4/46  81/9829 5.369299e-04 2.738343e-02 2.543352e-02                                                     ARID1B/SF3B1/TRRAP/SUZ12     4
# 16      50       CC GO:0016607                      nuclear speck      7/46 340/9829 9.320704e-04 3.510826e-02 3.260829e-02                                  SF3B1/THOC2/DDX17/TCF12/SRSF4/PRPF4B/BCLAF1     7
# 17      50       CC GO:0044666                     MLL3/4 complex      2/46  11/9829 1.147329e-03 3.510826e-02 3.260829e-02                                                                  KMT2C/NCOA6     2
# 18      50       MF GO:0003712 transcription coregulator activity     12/46 374/9777 9.874611e-08 1.757681e-05 1.507177e-05 JMJD1C/ARID1B/RERE/KMT2C/EP300/MED13L/DDX17/TRRAP/CCAR1/N4BP2L2/NCOA6/BCLAF1    12
# 19      50       MF GO:0003713 transcription coactivator activity      7/46 203/9777 4.039527e-05 3.595179e-03 3.082797e-03                                    ARID1B/RERE/KMT2C/EP300/DDX17/CCAR1/NCOA6     7
# 20      50       MF GO:0003682                  chromatin binding      8/46 383/9777 3.638758e-04 2.158996e-02 1.851298e-02                            JMJD1C/ARID1B/RERE/EP300/ASH1L/NCOA6/SUZ12/NAP1L4     8

as.data.frame(kegg)
#    Cluster       ID                                              Description GeneRatio  BgRatio       pvalue   p.adjust      qvalue                geneID Count
# 1       20 hsa05211                                     Renal cell carcinoma      3/10  58/4501 0.0002285934 0.01409436 0.006390966       10000/5894/2033     3
# 2       20 hsa04630                               JAK-STAT signaling pathway      3/10  65/4501 0.0003209489 0.01409436 0.006390966       10000/5894/2033     3
# 3       20 hsa05215                                          Prostate cancer      3/10  66/4501 0.0003358373 0.01409436 0.006390966       10000/5894/2033     3
# 4       20 hsa05164                                              Influenza A      3/10  80/4501 0.0005931352 0.01409436 0.006390966       10000/5894/2033     3
# 5       20 hsa04935           Growth hormone synthesis, secretion and action      3/10  84/4501 0.0006846501 0.01409436 0.006390966       10000/5894/2033     3
# 6       20 hsa05152                                             Tuberculosis      3/10  86/4501 0.0007336206 0.01409436 0.006390966       10000/5894/2033     3
# 7       20 hsa04068                                   FoxO signaling pathway      3/10  87/4501 0.0007589272 0.01409436 0.006390966       10000/5894/2033     3
# 8       20 hsa04919                        Thyroid hormone signaling pathway      3/10  92/4501 0.0008938768 0.01452550 0.006586461       10000/5894/2033     3
# 9       20 hsa05161                                              Hepatitis B      3/10 102/4501 0.0012078379 0.01744655 0.007910985       10000/5894/2033     3
# 10      20 hsa05225                                 Hepatocellular carcinoma      3/10 109/4501 0.0014646539 0.01904050 0.008633749      10000/57492/5894     3
# 11      20 hsa05167          Kaposi sarcoma-associated herpesvirus infection      3/10 115/4501 0.0017104341 0.02021422 0.009165963       10000/5894/2033     3
# 12      20 hsa04024                                   cAMP signaling pathway      3/10 132/4501 0.0025438388 0.02717131 0.012320596       10000/5894/2033     3
# 13      20 hsa04370                                   VEGF signaling pathway      2/10  38/4501 0.0029931309 0.02717131 0.012320596            10000/5894     2
# 14      20 hsa04664                          Fc epsilon RI signaling pathway      2/10  43/4501 0.0038219132 0.02717131 0.012320596            10000/5894     2
# 15      20 hsa04330                                  Notch signaling pathway      2/10  44/4501 0.0039991595 0.02717131 0.012320596             6310/2033     2
# 16      20 hsa04917                              Prolactin signaling pathway      2/10  44/4501 0.0039991595 0.02717131 0.012320596            10000/5894     2
# 17      20 hsa05213                                       Endometrial cancer      2/10  44/4501 0.0039991595 0.02717131 0.012320596            10000/5894     2
# 18      20 hsa05218                                                 Melanoma      2/10  44/4501 0.0039991595 0.02717131 0.012320596            10000/5894     2
# 19      20 hsa05221                                   Acute myeloid leukemia      2/10  44/4501 0.0039991595 0.02717131 0.012320596            10000/5894     2
# 20      20 hsa04929                                           GnRH secretion      2/10  45/4501 0.0041802022 0.02717131 0.012320596            10000/5894     2
# 21      20 hsa05230                      Central carbon metabolism in cancer      2/10  48/4501 0.0047459666 0.02920304 0.013241864            10000/5894     2
# 22      20 hsa04662                        B cell receptor signaling pathway      2/10  49/4501 0.0049420527 0.02920304 0.013241864            10000/5894     2
# 23      20 hsa04625                 C-type lectin receptor signaling pathway      2/10  55/4501 0.0061963113 0.03084564 0.013986689            10000/5894     2
# 24      20 hsa04720                                   Long-term potentiation      2/10  55/4501 0.0061963113 0.03084564 0.013986689             5894/2033     2
# 25      20 hsa05223                               Non-small cell lung cancer      2/10  55/4501 0.0061963113 0.03084564 0.013986689            10000/5894     2
# 26      20 hsa04916                                            Melanogenesis      2/10  56/4501 0.0064181805 0.03084564 0.013986689             5894/2033     2
# 27      20 hsa05235   PD-L1 expression and PD-1 checkpoint pathway in cancer      2/10  56/4501 0.0064181805 0.03084564 0.013986689            10000/5894     2
# 28      20 hsa05214                                                   Glioma      2/10  57/4501 0.0066436771 0.03084564 0.013986689            10000/5894     2
# 29      20 hsa05212                                        Pancreatic cancer      2/10  59/4501 0.0071054971 0.03185223 0.014443116            10000/5894     2
# 30      20 hsa05220                                 Chronic myeloid leukemia      2/10  61/4501 0.0075816607 0.03279280 0.014869611            10000/5894     2
# 31      20 hsa05165                           Human papillomavirus infection      3/10 200/4501 0.0082264373 0.03279280 0.014869611       10000/5894/2033     3
# 32      20 hsa04914                  Progesterone-mediated oocyte maturation      2/10  64/4501 0.0083225599 0.03279280 0.014869611            10000/5894     2
# 33      20 hsa04660                        T cell receptor signaling pathway      2/10  65/4501 0.0085765793 0.03279280 0.014869611            10000/5894     2
# 34      20 hsa05210                                        Colorectal cancer      2/10  65/4501 0.0085765793 0.03279280 0.014869611            10000/5894     2
# 35      20 hsa01521                EGFR tyrosine kinase inhibitor resistance      2/10  67/4501 0.0090951159 0.03284347 0.014892587            10000/5894     2
# 36      20 hsa04613                  Neutrophil extracellular trap formation      2/10  67/4501 0.0090951159 0.03284347 0.014892587            10000/5894     2
# 37      20 hsa04922                               Glucagon signaling pathway      2/10  68/4501 0.0093596060 0.03288510 0.014911463            10000/2033     2
# 38      20 hsa04012                                   ErbB signaling pathway      2/10  69/4501 0.0096275593 0.03293639 0.014934718            10000/5894     2
# 39      20 hsa01522                                     Endocrine resistance      2/10  71/4501 0.0101738015 0.03391267 0.015377406            10000/5894     2
# 40      20 hsa04066                                  HIF-1 signaling pathway      2/10  74/4501 0.0110188029 0.03410582 0.015464987            10000/2033     2
# 41      20 hsa04550 Signaling pathways regulating pluripotency of stem cells      2/10  74/4501 0.0110188029 0.03410582 0.015464987            10000/5894     2
# 42      20 hsa04666                         Fc gamma R-mediated phagocytosis      2/10  74/4501 0.0110188029 0.03410582 0.015464987            10000/5894     2
# 43      20 hsa05224                                            Breast cancer      2/10  77/4501 0.0118942501 0.03514210 0.015934880            10000/5894     2
# 44      20 hsa05226                                           Gastric cancer      2/10  77/4501 0.0118942501 0.03514210 0.015934880            10000/5894     2
# 45      20 hsa05231                             Choline metabolism in cancer      2/10  78/4501 0.0121927696 0.03522356 0.015971815            10000/5894     2
# 46      20 hsa04915                               Estrogen signaling pathway      2/10  79/4501 0.0124946191 0.03531088 0.016011411            10000/5894     2
# 47      20 hsa04210                                                Apoptosis      2/10  83/4501 0.0137350538 0.03799057 0.017226495            10000/5894     2
# 48      20 hsa04926                                Relaxin signaling pathway      2/10  84/4501 0.0140533563 0.03806117 0.017258508            10000/5894     2
# 49      20 hsa04071                           Sphingolipid signaling pathway      2/10  87/4501 0.0150277212 0.03986946 0.018078462            10000/5894     2
# 50      20 hsa04371                                 Apelin signaling pathway      2/10  89/4501 0.0156933831 0.04080280 0.018501673            10000/5894     2
# 51      20 hsa05160                                              Hepatitis C      2/10  91/4501 0.0163717933 0.04173202 0.018923022            10000/5894     2
# 52      20 hsa04722                           Neurotrophin signaling pathway      2/10  92/4501 0.0167157470 0.04178937 0.018949025            10000/5894     2
# 53      20 hsa04910                                Insulin signaling pathway      2/10 100/4501 0.0195796687 0.04718624 0.021396189            10000/5894     2
# 54      20 hsa04062                              Chemokine signaling pathway      2/10 101/4501 0.0199515066 0.04718624 0.021396189            10000/5894     2
# 55      20 hsa04022                               cGMP-PKG signaling pathway      2/10 102/4501 0.0203263799 0.04718624 0.021396189            10000/5894     2
# 56      20 hsa04218                                      Cellular senescence      2/10 102/4501 0.0203263799 0.04718624 0.021396189            10000/5894     2
# 57      20 hsa04072                        Phospholipase D signaling pathway      2/10 106/4501 0.0218559809 0.04984697 0.022602676            10000/5894     2
# 58      50 hsa04935           Growth hormone synthesis, secretion and action      4/20  84/4501 0.0004355245 0.04715475 0.027251455  10000/5894/2033/2776     4
# 59      50 hsa04919                        Thyroid hormone signaling pathway      4/20  92/4501 0.0006164020 0.04715475 0.027251455 10000/5894/2033/23389     4
# 60      50 hsa04929                                           GnRH secretion      3/20  45/4501 0.0009454515 0.04821803 0.027865939       10000/5894/2776     3


write_xlsx(as.data.frame(go), path = file.path(dir_rdata, "TableSxx_GO.xlsx"))
write_xlsx(as.data.frame(kegg), path = file.path(dir_rdata, "TableSxx_KEGG.xlsx"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.0 (2023-04-21)
#  os       macOS Ventura 13.4
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2023-06-13
#  rstudio  2023.06.0+421 Mountain Hydrangea (desktop)
#  pandoc   2.17.1.1 @ /opt/homebrew/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package          * version   date (UTC) lib source
#  AnnotationDbi    * 1.62.1    2023-05-02 [1] Bioconductor
#  ape                5.7-1     2023-03-13 [1] CRAN (R 4.3.0)
#  aplot              0.1.10    2023-03-08 [1] CRAN (R 4.3.0)
#  Biobase          * 2.60.0    2023-04-25 [1] Bioconductor
#  BiocGenerics     * 0.46.0    2023-04-25 [1] Bioconductor
#  BiocParallel       1.34.2    2023-05-28 [1] Bioconductor
#  Biostrings         2.68.1    2023-05-16 [1] Bioconductor
#  bit                4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
#  bit64              4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
#  bitops             1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
#  blob               1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
#  brio               1.1.3     2021-11-30 [1] CRAN (R 4.3.0)
#  cachem             1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
#  callr              3.7.3     2022-11-02 [1] CRAN (R 4.3.0)
#  cli                3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
#  clusterProfiler  * 4.8.1     2023-05-03 [1] Bioconductor
#  codetools          0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
#  colorout           1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace         2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
#  cowplot            1.1.1     2020-12-30 [1] CRAN (R 4.3.0)
#  crayon             1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
#  data.table         1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
#  DBI                1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
#  devtools         * 2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
#  digest             0.6.31    2022-12-11 [1] CRAN (R 4.3.0)
#  DOSE               3.26.1    2023-05-03 [1] Bioconductor
#  downloader         0.4       2015-07-09 [1] CRAN (R 4.3.0)
#  dplyr              1.1.2     2023-04-20 [1] CRAN (R 4.3.0)
#  ellipsis           0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
#  enrichplot         1.20.0    2023-04-25 [1] Bioconductor
#  fansi              1.0.4     2023-01-22 [1] CRAN (R 4.3.0)
#  farver             2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
#  fastmap            1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
#  fastmatch          1.1-3     2021-07-23 [1] CRAN (R 4.3.0)
#  fgsea              1.26.0    2023-04-25 [1] Bioconductor
#  fs                 1.6.2     2023-04-25 [1] CRAN (R 4.3.0)
#  generics           0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
#  GenomeInfoDb       1.36.0    2023-04-25 [1] Bioconductor
#  GenomeInfoDbData   1.2.10    2023-05-06 [1] Bioconductor
#  ggforce            0.4.1     2022-10-04 [1] CRAN (R 4.3.0)
#  ggfun              0.0.9     2022-11-21 [1] CRAN (R 4.3.0)
#  ggplot2          * 3.4.2     2023-04-03 [1] CRAN (R 4.3.0)
#  ggplotify          0.1.0     2021-09-02 [1] CRAN (R 4.3.0)
#  ggraph             2.1.0     2022-10-09 [1] CRAN (R 4.3.0)
#  ggrepel            0.9.3     2023-02-03 [1] CRAN (R 4.3.0)
#  ggtree             3.8.0     2023-04-25 [1] Bioconductor
#  glue               1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
#  GO.db              3.17.0    2023-05-06 [1] Bioconductor
#  GOSemSim           2.26.0    2023-04-25 [1] Bioconductor
#  graphlayouts       1.0.0     2023-05-01 [1] CRAN (R 4.3.0)
#  gridExtra          2.3       2017-09-09 [1] CRAN (R 4.3.0)
#  gridGraphics       0.5-1     2020-12-13 [1] CRAN (R 4.3.0)
#  gson               0.1.0     2023-03-07 [1] CRAN (R 4.3.0)
#  gtable             0.3.3     2023-03-21 [1] CRAN (R 4.3.0)
#  HDO.db             0.99.1    2023-05-06 [1] Bioconductor
#  here             * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
#  hms                1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
#  htmltools          0.5.5     2023-03-23 [1] CRAN (R 4.3.0)
#  htmlwidgets        1.6.2     2023-03-17 [1] CRAN (R 4.3.0)
#  httpuv             1.6.11    2023-05-11 [1] CRAN (R 4.3.0)
#  httr               1.4.6     2023-05-08 [1] CRAN (R 4.3.0)
#  igraph             1.4.3     2023-05-22 [1] CRAN (R 4.3.0)
#  IRanges          * 2.34.0    2023-04-25 [1] Bioconductor
#  jsonlite           1.8.5     2023-06-05 [1] CRAN (R 4.3.0)
#  KEGGREST           1.40.0    2023-04-25 [1] Bioconductor
#  later              1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
#  lattice            0.21-8    2023-04-05 [1] CRAN (R 4.3.0)
#  lazyeval           0.2.2     2019-03-15 [1] CRAN (R 4.3.0)
#  lifecycle          1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
#  lubridate          1.9.2     2023-02-10 [1] CRAN (R 4.3.0)
#  magrittr           2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
#  MASS               7.3-60    2023-05-04 [1] CRAN (R 4.3.0)
#  Matrix             1.5-4.1   2023-05-18 [1] CRAN (R 4.3.0)
#  memoise            2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
#  mime               0.12      2021-09-28 [1] CRAN (R 4.3.0)
#  miniUI             0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
#  munsell            0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
#  nlme               3.1-162   2023-01-31 [1] CRAN (R 4.3.0)
#  org.Hs.eg.db     * 3.17.0    2023-05-06 [1] Bioconductor
#  patchwork          1.1.2     2022-08-19 [1] CRAN (R 4.3.0)
#  pillar             1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
#  pkgbuild           1.4.0     2022-11-27 [1] CRAN (R 4.3.0)
#  pkgconfig          2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
#  pkgload            1.3.2     2022-11-16 [1] CRAN (R 4.3.0)
#  plyr               1.8.8     2022-11-11 [1] CRAN (R 4.3.0)
#  png                0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
#  polyclip           1.10-4    2022-10-20 [1] CRAN (R 4.3.0)
#  prettyunits        1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
#  processx           3.8.1     2023-04-18 [1] CRAN (R 4.3.0)
#  profvis            0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
#  promises           1.2.0.1   2021-02-11 [1] CRAN (R 4.3.0)
#  prompt             1.0.1     2023-05-06 [1] Github (gaborcsardi/prompt@7ef0f2e)
#  ps                 1.7.5     2023-04-18 [1] CRAN (R 4.3.0)
#  purrr              1.0.1     2023-01-10 [1] CRAN (R 4.3.0)
#  qvalue             2.32.0    2023-04-25 [1] Bioconductor
#  R6                 2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
#  RColorBrewer       1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp               1.0.10    2023-01-22 [1] CRAN (R 4.3.0)
#  RCurl              1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
#  remotes            2.4.2     2021-11-30 [1] CRAN (R 4.3.0)
#  reshape2           1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
#  rlang              1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
#  rprojroot          2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
#  RSQLite            2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
#  rsthemes           0.4.0     2023-05-06 [1] Github (gadenbuie/rsthemes@34a55a4)
#  rstudioapi         0.14      2022-08-22 [1] CRAN (R 4.3.0)
#  S4Vectors        * 0.38.1    2023-05-02 [1] Bioconductor
#  scales             1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
#  scatterpie         0.2.1     2023-06-07 [1] CRAN (R 4.3.0)
#  sessioninfo      * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
#  shadowtext         0.1.2     2022-04-22 [1] CRAN (R 4.3.0)
#  shiny              1.7.4     2022-12-15 [1] CRAN (R 4.3.0)
#  stringi            1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
#  stringr            1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
#  suncalc            0.5.1     2022-09-29 [1] CRAN (R 4.3.0)
#  testthat         * 3.1.8     2023-05-04 [1] CRAN (R 4.3.0)
#  tibble             3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
#  tidygraph          1.2.3     2023-02-01 [1] CRAN (R 4.3.0)
#  tidyr              1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
#  tidyselect         1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
#  tidytree           0.4.2     2022-12-18 [1] CRAN (R 4.3.0)
#  timechange         0.2.0     2023-01-11 [1] CRAN (R 4.3.0)
#  treeio             1.24.1    2023-05-31 [1] Bioconductor
#  tweenr             2.0.2     2022-09-06 [1] CRAN (R 4.3.0)
#  urlchecker         1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
#  usethis          * 2.2.0     2023-06-06 [1] CRAN (R 4.3.0)
#  utf8               1.2.3     2023-01-31 [1] CRAN (R 4.3.0)
#  vctrs              0.6.2     2023-04-19 [1] CRAN (R 4.3.0)
#  viridis            0.6.3     2023-05-03 [1] CRAN (R 4.3.0)
#  viridisLite        0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
#  withr              2.5.0     2022-03-03 [1] CRAN (R 4.3.0)
#  writexl          * 1.4.2     2023-01-06 [1] CRAN (R 4.3.0)
#  xtable             1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
#  XVector            0.40.0    2023-04-25 [1] Bioconductor
#  yulab.utils        0.0.6     2022-12-20 [1] CRAN (R 4.3.0)
#  zlibbioc           1.46.0    2023-04-25 [1] Bioconductor
# 
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
