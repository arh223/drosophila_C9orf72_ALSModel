# man.dat is data.frame
# names(man.dat): chr, pos, pval
# hit.idx: optional index of SNP hits (e.g., SNPs with p-values below some threshold); these points will be colored red
# main.title: optional plot title

manhattan.plot <- function(manhat.df, main.title="Manhattan Plot", hit.idx=NULL)
{
  names(manhat.df) <- c("chr", "pos", "pval")
  chr.list <- unique(manhat.df$chr)
  col.vec <- rep(c(1,8), length=length(chr.list))
  full.col.vec <- col.vec[match(manhat.df$chr, chr.list)]
  
  if (!is.null(hit.idx))
  {
    full.col.vec[hit.idx] <- 2
  }
  
  chr.tab <- table(manhat.df$chr)
  mid.val <- round(chr.tab/2)
  start.idx <- sapply(1:length(chr.tab), function(x) sum(chr.tab[1:x]))
  tick.idx <- start.idx - mid.val
  tick.lab <- names(chr.tab)
  
  
  plot(-log10(manhat.df$pval), ylim=c(0,12.5), xaxt='n', col=full.col.vec, pch=16, ylab="-log10(p-value)", main=main.title, xlab="chromosome")
  axis(1, at=tick.idx, labels=c("1", "2", "3", "4", "5", "6"), las=2, cex.axis=0.8)
  abline(h=5, col="red", lty=2)
}



