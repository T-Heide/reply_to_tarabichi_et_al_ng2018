################################################################################
# FILENAME:    '03-create_figures.R'
################################################################################

# Options and variables:
dataDir   <- "./results/datasets"
simResDir <- "./results/simulations"
figDir    <- "./results/figures"

# Libs: ########################################################################

library(cowplot)
library(dplyr)
library(reshape2)

# Functions: ###################################################################

# Modification function for spread matrices:
modMatrix <- function(m) {
  rownames(m)   <- m[,1]   # Set resonable rownames.
  m             <- m[,-1]
  m             <- t(as.matrix(m))
  m[is.nan(m)]  <- NA
  return(m)
}

charLog2 <- function(x) {
  sapply(x, function(y){
    y <- log2(as.numeric(as.character(y)))
    bquote(2^.(y))
  })
}

plotHistogram <- function(datFileBase) {
  
  # Load simulated sequencing data:
  datA <- read.delim(paste0(datFileBase, "-simulated_sequencing.tsv"))
  pDat_hist <- filter(datA, VAF > 0.05)
  pDat_hist$CLONE  <- as.character(pDat_hist$CLONE)
  
  # Load cell count data:
  datB <- read.delim(paste0(datFileBase, "-cell_number.tsv"))
  pDat_lines <- data.frame(lty=2, x=datB$clone2 / datB$total * 0.5)
  
  # Plot data:  
  p <- ggplot(pDat_hist, aes(x=VAF, fill=CLONE, group=CLONE)) + 
    geom_histogram(breaks=seq(from=0, to=1, by=0.01)) + 
    geom_vline(data=pDat_lines, aes(xintercept=x), linetype=2) +
    xlim(0, NA) +
    scale_fill_manual(name="Clone",
                      breaks=c("1","0","2"),
                      labels=c("Background","Mutant (ancestral)",
                               "Mutant (private)"),
                      values=c("1"="#009E73", "0"="#D55E00", "2"="#0072B2"))
    
  return(p)
}

################################################################################

# Load the data:
modFits    <- readRDS(file.path(dataDir, "1f_model_fits.rds"))
modFitsExt <- readRDS(file.path(dataDir, "1f_model_fits_ext.rds"))
cellCounts <- readRDS(file.path(dataDir, "cell_counts.rds"))


# Create figure dir:
dir.create(figDir, showWarnings=FALSE, recursive=TRUE)

# 1) Create heatmaps of 1f-model tests #########################################

# Spread results as matrix containing fractions of neutral calls:
matrixRes    <- modMatrix(dcast(modFits, sc_mutation_rate ~ sc_deltaS, 
                                value.var="non_neutral", fun.aggregate=mean,
                                na.rm=TRUE))

matrixResExt <- modMatrix(dcast(modFitsExt, sc_mutation_rate ~ sc_deltaS, 
                                value.var="non_neutral", fun.aggregate=mean,
                                na.rm=TRUE))


# Heatmap of dS and subclone mutation rate time combinations:
pheat <- melt(matrixRes) %>%
  ggplot(aes(x=factor(Var2), y=Var1 - 1.0, fill=value)) + 
  geom_tile(color="gray30", linetype=0) +
  scale_x_discrete(breaks=2^seq(1, 10, by=1), labels=charLog2) +
  scale_y_continuous(breaks=seq(0, 1, by=0.2)) +
  xlab(bquote('subclone mutation rate ('*mu['sc']*')')) +
  ylab(bquote('subclone selective advantage ('*Delta*S['sc']*')')) +
  ggtitle(label="Stochastic simulations: 1/f test") +
  scale_fill_gradientn(colors=c(colorRampPalette(c("blue","white"))(27),
                                colorRampPalette(c("white","red"))(100-27),
                                rep("white", 10)),
                       limits=c(0,1.1), 
                       breaks=seq(0, 1, by=0.25),
                       na.value = 'gray75') +
  labs(fill=expression(R^2<0.98))

ggsave(file.path(figDir, "heatmap-frac_neutral-dS-mu.pdf"), 
       pheat, width=6, height=5)


pheatExt <- melt(round(matrixResExt*2, 1)/2) %>%
  filter(Var2 %in% as.character(2^(1:9))) %>%
  filter(Var1 %in% as.character(seq(1, 2, by=0.05))) %>%
  ggplot(aes(x=factor(Var2), y=Var1-1.0, fill=value)) + 
    geom_tile(color="gray30", linetype=1) +
    scale_x_discrete(breaks=2^seq(1, 10, by=1), labels=charLog2) + 
    scale_y_continuous(breaks=seq(0, 1, by=0.2)) + 
    xlab(bquote('subclone mutation rate ('*mu['sc']*')')) + 
    ylab(bquote('subclone selective advantage ('*Delta*S['sc']*')')) + 
    ggtitle(label="Stochastic simulations: 1/f test [0.025,0.45]") + 
    labs(fill="Sensitivity") +
    scale_fill_gradientn(colors=colorRampPalette(c("white","darkgreen"))(100),
                         limits=c(0,1)) + 
    geom_text(aes(label=sprintf("%.2f", value)))


ggsave(file.path(figDir, "heatmap_frac_neutral_dS-mu-full_range.pdf"), 
       pheatExt, width=6.2, height=5)


# 2) Plot subclone frequency for all dS values: ################################

# Plot the relationship of selective advantage and subclone fraction:
data_text_layer <- data.frame(x=c(0.5,0.99), 
                              y=c(1.05,0.025), 
                              label=c("Fixation","LOD (100X)"))

plot_cell_fracs <- cellCounts %>%
  filter(as.character(sc_deltaS) %in% seq(from=1, to=2, by=0.1)) %>%
  ggplot(aes(x=sc_deltaS-1, y=subcloneFrac, group=as.character(sc_deltaS))) +
  geom_violin(width=0.05, fill="black", alpha=0.3, scale="width") +
  geom_boxplot(width=0.02, outlier.shape=NA, fill="black", linetype=0) +
  stat_summary(fun.y=median, geom="point", colour="gray95", size=0.75) +
  xlab(bquote('subclone selective advantage (adv)')) + 
  ylab(bquote('subclone cell fraction')) +
  geom_hline(yintercept=0.1, linetype=2, color="#7F7F7F", alpha=0.8) +
  geom_hline(yintercept=1, linetype=2, color="#7F7F7F", alpha=0.8) +
  scale_x_continuous(breaks=seq(from=0, to=1, by=0.2), limits=c(-0.05,1.1)) +
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.2))+
  background_grid(major="xy", size.major=0.5) +
  geom_text(data=data_text_layer, aes(x=x, y=y, label=label), inherit.aes=0)

ggsave(file.path(figDir, "correlation-dS-frac_clone_alt.pdf"), 
       plot_cell_fracs, width=4.5, height=3.5)


# 3) Create per simulation histograms of selected combinations: ################

# Parameter combinations to plot as histogram:
selParams  <- list(c("16","1.5","1006"), 
                   c("16","1.25","1012"), 
                   c("128","1.5","1015"), 
                   c("1024","1.5","1013"))

null <- lapply(selParams, function(params) {
  # Determine base path of simulation files for current parameters:
  patF  <- "simulation-mmr_%s-mbr_%s-seed_%s-cst_%s"
  baseN <- sprintf(patF, params[1], params[2],  params[3], "256")
  baseP <- file.path(simResDir, 
                     sprintf("mmr_%s", params[1]), 
                     sprintf("mbr_%s", params[2]), 
                     baseN)
  
  # Load data and plot as histogram:
  histoP  <- plotHistogram(baseP)
  figName <- sprintf("histogram_simulation_mmr%s-mbr%s-seed%s.pdf", 
                     params[1], params[2],  params[3])
  ggsave(file.path(figDir, figName), histoP, width=5, height=2.5)
})
