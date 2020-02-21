source("./general_functions.R")
OUT_DIR = Sys.getenv('OUT_DIR')
library(ggplot2)

df = read.csv(paste0(OUT_DIR, "/astle_compare/astle_var_ukbb500k.csv"), stringsAsFactors=F)
source("a")
df = df[!is.element(df$trait, c("BASO#", "BASO%", "RDW")),]
mini = min(df$Minor.Allele.Frequency)
maxi = max(df$Minor.Allele.Frequency)
increm = (maxi-mini)/10
df$X8 <- cut(df$Minor.Allele.Frequency, breaks = seq(mini, maxi, increm))
p = ggplot(df, aes(y =r2, x = factor(X8))) +
      geom_boxplot( outlier.shape = NA,outlier.colour = NA)+
      geom_point(color="red")+
      labs(x = "MAF", y = "R2" )+
      scale_x_discrete(breaks=seq(0,360,15))


png("~/a.png")
print(p)
dev.off()
