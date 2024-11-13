setwd("~/lustre_sl37/cardinal/dynamic_eQTLs/gaspacho_fork/GASPACHO")
library("Matrix")
for (file in list.files("R/", pattern="*.R")) {
  source(paste0("R/", file))
}
library(reticulate)
np <- import("numpy")
df = read.csv("Bcells_res/test.csv")[1:5000,]
gplvm = readRDS("Bcells_res/test_gplvm.updatedDCSE.rds")
metadata = read.csv("Bcells_res/bcells_seacells.meta.csv")[1:5000,]
library(stringr)
df$donor = sapply(str_split(df$X, "-"), `[`, 1)
df$donor_id = as.numeric(as.factor(df$donor))
gt_donor_df = unique(subset(df, select=c(gt, donor_id)))
ordered_gt_donor_df = gt_donor_df[order(gt_donor_df$donor_id),]

# replicate an eQTL with SELL gene
bfs = getBF_withnan(yj = df$SELL_exp,
                    gplvm = gplvm, 
                    G = as.matrix(rbind(ordered_gt_donor_df$gt, ordered_gt_donor_df$gt)), 
                    did = as.numeric(as.factor(df$donor)))
bfs
y=df$SELL_exp
g0 = df$gt
af = mean(df$gt, na.rm = TRUE)/2
deltaG = 0.001
betas = getBeta(
  y, 
  gplvm, 
  g0, 
  af, 
  deltaG = 0.001
)
plotdf = as.data.frame(cbind(gplvm$Param$Xi[[2]][, 1:4], df$SELL_exp, betas))
ggplot(plotdf, aes(x = V1, y = V2, color = betas)) +
  geom_point(size=0.1) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme_minimal()

# test a random input as expression 
random_vector <- rnorm(5000, mean = 0, sd = 1)
bfs = getBF_withnan(yj = random_vector,
                    gplvm = gplvm, 
                    G = as.matrix(rbind(ordered_gt_donor_df$gt, ordered_gt_donor_df$gt)), 
                    did = as.numeric(as.factor(df$donor)))
bfs
y=random_vector
g0 = df$gt
af = mean(df$gt, na.rm = TRUE)/2
deltaG = 0.001
betas = getBeta(
  y, 
  gplvm, 
  g0, 
  af, 
  deltaG = 0.001
)
plotdf = as.data.frame(cbind(gplvm$Param$Xi[[2]][, 1:4], df$SELL_exp, betas))
ggplot(plotdf, aes(x = V1, y = V2, color = betas)) +
  geom_point(size=0.1) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme_minimal()
