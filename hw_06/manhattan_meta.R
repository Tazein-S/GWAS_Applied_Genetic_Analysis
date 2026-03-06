library(ggplot2)
library(data.table)

dat <- fread("thin.TBL", col.names=c("CHR", "POS", "P.value"))

dat[, CHR := as.numeric(CHR)]
dat <- dat[!is.na(P.value) & P.value > 0]
dat <- dat[order(CHR, POS)]

dat[, cumpos := POS + (CHR - 1) * 3e8]
chrom_mids <- dat[, .(mid=mean(cumpos)), by=CHR]

ggplot(dat, aes(x=cumpos, y=-log10(P.value), color=factor(CHR%%2))) +
  geom_point(size=0.3) +
  scale_color_manual(values=c("cyan4", "navy"), guide="none") +
  scale_x_continuous(breaks=chrom_mids$mid, labels=chrom_mids$CHR) +
  geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed") +
  labs(x="Chromosome", y="-log10(P-value)", title="Meta-analysis Manhattan Plot") +
  theme_bw()

ggsave("manhattan_meta.jpeg", dpi=720, width=12, height=6)
