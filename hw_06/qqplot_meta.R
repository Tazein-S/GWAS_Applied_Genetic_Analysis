library(ggplot2)
library(data.table)

pvals <- fread("METAANALYSIS1.TBL", select="P-value")[[1]]

# thin points - only keep every 10th point in the bulk, all extreme ones
n <- length(pvals)
qq_df <- data.frame(
  exp = sort(-log10(ppoints(n)), decreasing=TRUE),
  obs = sort(-log10(pvals), decreasing=TRUE)
)

# thin the non-extreme points to speed up plotting
thin <- c(which(qq_df$obs > 3), seq(1, nrow(qq_df), by=10))
qq_df <- qq_df[unique(thin), ]

p <- ggplot(qq_df, aes(x=exp, y=obs)) +
  geom_point(color="cyan4", size=0.5) +
  geom_abline(slope=1, intercept=0) +
  labs(x="-log10(expected P)", y="-log10(observed P)", title="Meta QQ Plot: Uganda, DCC, DDS, AADM") +
  theme_bw()

ggsave("qq_metaanalysis.jpeg", plot=p, dpi=720)
