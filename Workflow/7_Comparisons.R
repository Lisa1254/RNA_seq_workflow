##
## Comparisons of Multiple Datasets ----
##

#Data required for this script:
##   > significant genes determined in scripts 4 & 5
load("Output/sig_deg_asd_control.Rdata")
load("Output/sig_drim_asd_control.Rdata")
##   > other gene lists for comparison (can use SPARK ASD-risk genes as specified below)
##   > Gene to symbol & description mappings (Script 2)
load("Output/gene_symbols.Rdata")

##
# Libraries & Functions----
##

#homologene is used to map genes to known orthologs in other species for comparison
library(homologene)

#Function used to 
source()

##
# Define Key Variables ----
##

#Output directory for saving files
out_dir = paste0(getwd(),"/Output/")

#Taxonomic ID number for species being compared
human_taxID <- 9606
mouse_taxID <- 10090
#Zebrafish ID, as another common species is below
#zebrafish_taxID <- 7955


#For comparison, will be using list of genes with demonstrated strong association with ASD-risk, according to Simons Foundation Powering Autism Research for Knowledge
#https://www.sfari.org/grant/2021-genomics-of-asd-pathways-to-genetic-therapies-request-for-applications/
human_SPARK <- c("ADNP", "AHDC1", "ANKRD11", "ARID1B", "ASXL3", "BCKDK", "BCL11A", "CASK", "CHAMP1", "CHD2", "CHD8", "CNOT3", "CREBBP", "CSNK2A1", "CTCF", "CTNNB1", "DDX3X", "DEAF1", "DYNC1H1", "DYRK1A", "EHMT1", "EIF3F", "EP300", "FOXG1", "FOXP1", "GRIN2B", "HIVEP2", "HNRNPH2", "HNRNPU", "IQSEC2", "KMT2A", "MED13L", "MYT1L", "NSD1", "PACS1", "POGZ", "PPP2R5D", "SCN1A", "SCN2A", "SCN8A", "SETBP1", "SETD5", "SLC6A1", "SMARCA2", "SRCAP", "STXBP1", "SYNGAP1", "TRIO", "USP9X")


##
# Compare known homologs ----
##

#Consider streamlining the next series of steps into a function

#Express significant genes by their symbols
sig_deg_symbol <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_asd.control)),"external_gene_name"]
sig_drim_symbol <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_drim_asd.control)),"external_gene_name"]

#Check for duplicated gene names
sum(duplicated(sig_deg_symbol))
sum(duplicated(sig_drim_symbol))
#Both have zero duplicated names

#If any are duplicated, can either locate and manually adjust symbols, or
sig_deg_symbol <- unique(sig_deg_symbol)
sig_drim_symbol <- unique(sig_drim_symbol)

#Use homologene to get homologs
deg_homologs <- homologene(sig_deg_symbol, inTax = mouse_taxID, outTax = human_taxID,  db = homologene::homologeneData2)

drim_homologs <- homologene(sig_drim_symbol, inTax = mouse_taxID, outTax = human_taxID,  db = homologene::homologeneData2)

#Check if any of genes in this set match the SPARK list
deg_homologs[which(deg_homologs$`9606` %in% human_SPARK),]
#None
drim_homologs[which(drim_homologs$`9606` %in% human_SPARK),]
#Cask, Iqsec2, and Myt1l match

homol_overlap_drim <- drim_homologs[which(drim_homologs$`9606` %in% human_SPARK),c(1,2)]
colnames(homol_overlap_drim) <- c("Mouse", "Human")
homol_overlap_drim$Description <- gene_symbols[match(homol_overlap_drim$Mouse, gene_symbols$external_gene_name), 3]


##
# Comparison Plots ----
##

# Would like to add Venn diagram style comparison

# This is a paste from my other script of creating a heatmap style of overlap in significant genes in several sets
# This is more useful with more than 2 comparison sets
# Consider modifying the code to keep as function or process so that it can accommodate desired number of gene sets for comparison and/or adding another set of relevant genes to the import for comparison

any_asdSet_sym <- annot_sig_dge$external_gene_name
any_vanSet_sym <- annot_sig_dge_vanco$external_gene_name

top50_risk_drSym <- annot_ens_dr$external_gene_name

gene_sets <- c("any_asdSet_sym", "any_vanSet_sym", "top50_risk_drSym", "resc_genes_update")

all_genes_any_set <- unique(c(any_asdSet_sym, any_vanSet_sym, top50_risk_drSym, resc_genes_update))

SigGene_table_sym <- matrix(nrow=length(all_genes_any_set), ncol=length(gene_sets))
rownames(SigGene_table_sym) <- all_genes_any_set
colnames(SigGene_table_sym) <- c("Any ASD-NT", "VANCO", "top 50 Risk", "Prev. Rescue Study")

for (gl in 1:length(gene_sets)) {
  temp_list <- get(gene_sets[gl])
  for (gene in 1:nrow(SigGene_table_sym)) {
    if (rownames(SigGene_table_sym)[gene] %in% temp_list) {
      SigGene_table_sym[gene,gl] <- 1
    } else {
      SigGene_table_sym[gene,gl] <- 0
    }
  }
}

SigGene_table_sym2 <- SigGene_table_sym[which(rowSums(SigGene_table_sym)>1),]

gene_annot <- gene_symbol_all[match(rownames(SigGene_table_sym2),gene_symbol_all$external_gene_name),"description"]
gene_annot <- strsplit(gene_annot, " \\[")
gene_annot <- unlist(gene_annot)
gene_annot <- gene_annot[seq(1,length(gene_annot),2)]

ha_gene <- rowAnnotation(
  Gene = anno_text(paste(rownames(SigGene_table_sym2), gene_annot, sep = ": "),
                   rot = 0,
                   just = 'left',
                   gp = gpar(fontsize = 8))
)

hmap <- Heatmap(SigGene_table_sym2,
                name = 'Genes in Common',
                col = c('0' = 'white', '1' = 'forestgreen'),
                rect_gp = gpar(col = 'grey85'),
                border_gp = gpar(col = "black"),
                
                show_row_dend = FALSE)
draw(hmap + ha_gene)





##