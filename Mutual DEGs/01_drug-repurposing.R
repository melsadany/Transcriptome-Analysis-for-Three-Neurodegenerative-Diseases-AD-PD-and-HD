################################################################################
#                             drug repurposing                                 #
################################################################################
# set
rm(list = ls())
gc()
pacman::p_load(tidyverse, ggh4x, doMC)
theme_set(theme_minimal() +
            theme(axis.line = element_line(colour = "black", linewidth = 0.5), 
                  axis.text = element_text(face = "bold"), 
                  axis.text.y = element_text(size = 8),
                  axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
                  axis.ticks = element_blank(), 
                  axis.title = element_text(face = "bold"), 
                  strip.text = element_text(face = "bold"),
                  strip.switch.pad.grid = unit(0, "points"),
                  strip.placement = "outside", 
                  panel.spacing = unit(1, "points"), 
                  panel.grid = element_blank(),
                  legend.position = "bottom", 
                  # legend.title = element_blank(), 
                  legend.text = element_text(face = "bold", size = 5),
                  plot.title = element_text(size = 10),
                  plot.subtitle = element_text(size = 8),
                  title = element_text(face = "bold", size = 10),
                  plot.caption = element_text(hjust = 0)
            ))
my.guides <- guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5))
redblack.col <- c("#800000", "black")
six.colors <- c("#800000", "#cc7277", "#4f6162", "#e65236", "#56483a", "#73937e")
pload <- function(fname,envir=.GlobalEnv){
  con <- pipe(paste("/opt/homebrew/bin/pixz -d <",fname),"rb")
  load(con,envir=envir); close(con)
}
psave <- function(...,file){  
  con = pipe(paste("/opt/homebrew/bin/pixz -2 -q 80 -f 3 > ",file,".pxz",sep=""),"wb") 
  save(...,file=con,envir=.GlobalEnv); close(con) 
} 
################################################################################
project.dir <- "/Volumes/Mac/ZC/GP-paper/Transcriptome-Analysis-for-Three-Neurodegenerative-Diseases-AD-PD-and-HD/"
setwd(project.dir)
################################################################################
# read the CMAP signature, and keep FDA-approed drugs only
fda <- read_delim("../../FDA_approved_Products.txt") %>%
  mutate(DrugName = sub(" ", "-", tolower(DrugName)),
         ActiveIngredient = sub(" ", "-", tolower(ActiveIngredient)))
# full CMAP matrix with rows as drug/molecule names and genes as colnames
pload("../../cp_mean_coeff_mat_tsv.Rdata.pxz")
fda.mdrug <- mdrug[c(tolower(rownames(mdrug)) %in% tolower(fda$DrugName)|
                       tolower(rownames(mdrug)) %in% tolower(fda$ActiveIngredient))
                   ,]
# save fda drug matrix
psave(fda.mdrug, file = "../../cp_mean_coeff_mat_tsv_FDA-only.Rdata.pxz")
gc()
################################################################################
# read list of mutual DEGS in NDDs
mut.genes <- read_csv("Mutual DEGs/MutualsValuesAll_V2.csv", 
                      col_names = c("gene", 
                                    "AD_lfc", "AD_q",
                                    "PD_lfc", "PD_q",
                                    "HD_lfc", "HD_q"), 
                      skip = 2) %>%
  mutate(AD = ifelse(AD_lfc < 0, "D", "U"),
         PD = ifelse(PD_lfc < 0, "D", "U"),
         HD = ifelse(HD_lfc < 0, "D", "U")) %>%
  filter(AD == PD & AD == HD)
################################################################################
# filter cmap drug matrix to only keep genes in the mutual DEGs list
fda.mdrug.filtered <- fda.mdrug[,colnames(fda.mdrug) %in% mut.genes$gene]
################################################################################
# write the function to calculate jaccard index, assuming you're entering just a direction for disease genes
jaccard.index <- function(subject, drug) {
  # I'm assuming genes are in the first column of each of the input vectors/dataframes
  df <- right_join(subject, drug)
  colnames(df) <- c("gene", "subject", "drug")
  df <- df %>%
    mutate(subject = ifelse(is.na(subject), 0, subject)) %>%
    mutate(UU = ifelse(subject=="U" & drug>0, 1, 0)) %>%
    mutate(DD = ifelse(subject=="D" & drug<0, 1, 0)) %>%
    mutate(UD = ifelse(subject=="U" & drug<0, 1, 0)) %>%
    mutate(DU = ifelse(subject=="D" & drug>0, 1, 0))
  UU <- df %>% filter(subject =="U" | drug > 0)
  DD <- df %>% filter(subject =="D" | drug < 0)
  UD <- df %>% filter(subject =="U" | drug < 0)
  DU <- df %>% filter(subject =="D" | drug > 0)
  jaccard <- ((sum(df$UU)/nrow(UU))+(sum(df$DD)/nrow(DD))-(sum(df$UD)/nrow(UD))-(sum(df$DU)/nrow(DU)))/2
  df2 <- data.frame(drug = jaccard,
                    UU = sum(df$UU),
                    DD = sum(df$DD),
                    UD = sum(df$UD),
                    DU = sum(df$DU))
  # colnames(df2)[1] <- colnames(drug)[2]
  return(df2)
}
################################################################################
# get correlation/jaccard index by drug
registerDoMC(cores = 3)
d.scores <- foreach(i = 1:nrow(fda.mdrug.filtered), .combine = rbind) %dopar% {
  d <- rownames(fda.mdrug.filtered)[i]
  j <- jaccard.index(subject = mut.genes%>%select(gene, all=AD), 
                     drug = fda.mdrug.filtered[i,]%>%as.data.frame()%>%rownames_to_column("gene"))
  cbind(drug = d,j %>% rename(score =1))
  
}
# plot distribution of computed scores
d.scores %>%
  ggplot(aes(x=score)) +
  geom_histogram(bins = 50) +
  labs(x="Jaccard Index")
ggsave("Mutual DEGs/Visualization/distribution-of-Jaccard-Index-with-FDA-approved-drugs.svg",
       width = 3, height = 2, units = "in", bg = "white", dpi = 360)

# plot number of genws with changes, and total score as fill color
d.scores %>%
  ggplot(aes(x=UD, y=DU, color = score)) +
  geom_point(size=1) +
  scale_color_gradient2(low = redblack.col[2], high = redblack.col[1], 
                        name = "Jaccard Index") +
  guides(color = guide_colorbar(barwidth = 6, barheight = 0.5)) +
  labs(x = "number of UP genes in NDDs and get DOWN by the drug",
       y = "number of DOWN genes in NDDs and get UP by the drug",
       caption = "a positive value in Jaccard Index means the drug has a similar transcriptomic profile to genes' lfc")
ggsave("Mutual DEGs/Visualization/UD-VS-DU-from-Jaccard-Index.svg",
       width = 6, height = 5, units = "in", bg = "white", dpi = 360)

# save the calculated Jaccard Index scores for all drugs
write_csv(d.scores, "Mutual DEGs/computed-jaccard-index-scores.csv")
################################################################################
# heatmap for 5 genes of interest
genes.of.int <- c("NFKB1", "NFKBIA", "RELA", "TRIM4", "SMAD4")
drugs.of.int <- d.scores %>%
  filter(score < -0.1)
fda.mdrug.filtered[drugs.of.int$drug,] %>%
  as.data.frame() %>%
  select(any_of(genes.of.int)) %>%
  rownames_to_column("drug") %>%
  pivot_longer(cols = any_of(genes.of.int), names_to = "gene", values_to = "drug_exp") %>%
  ggplot(aes(x=gene, y=drug, fill=drug_exp))+
  geom_tile()+
  scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1])+
  my.guides+
  theme(axis.text.x.bottom = element_text(angle = 0, hjust = 0.5))
ggsave("Mutual DEGs/Visualization/gens-exp-in-drugs-of-int.svg",
       width = 5, height = 11, units = "in", dpi = 360, bg = "white")
# mut.genes %>%
#   filter(gene %in% genes.of.int) %>%
#   select(gene, ends_with("lfc")) %>%
#   pivot_longer(cols = c("AD_lfc", "PD_lfc", "HD_lfc"), names_to = "NDD", values_to = "NDD_exp") %>%
#   mutate(NDD = sub("_lfc", "", NDD)) %>%
#   left_join(mdrug.filtered.2, relationship = "many-to-many") %>% drop_na()
# ################################################################################


################################################################################


