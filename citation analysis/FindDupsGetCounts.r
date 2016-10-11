#=====================================================
#Web of Science was used to locate the citations
#of individual papers that used or cited the
#following software packages:
#ADMB, Stan, TMB, JAGS, Classic BUGS, WinBUGS,
#OpenBUGS
#Since multiple methods were used to find citing
#papers, duplicates must be removed before returning
#the cites by year.
#=====================================================
XX <- read.csv("Citation analysis.csv", stringsAsFactors=F)
software <- c("Classic BUGS", "WinBUGS", "OpenBUGS",
              "JAGS", "ADMB", "TMB", "Stan")
nsoftware <- length(software)
years <- seq(min(XX$PY), 2015, 1)
nyears <- length(years)
results <- matrix(nrow=nyears, ncol=nsoftware,
                  dimnames=list(years,software))
#get data in the right format
for (i in 1:nsoftware) {
   #go through the data for each software package in turn
   Xdata <- XX[XX$Software.package==software[i],]
   #find the duplicated doi values
   doi.dups <- duplicated(Xdata$DI, incomparables=c(""))
   #find the duplicated journal+vol+pagestart
   journal.dups <- duplicated(paste(Xdata$SO, Xdata$VL, Xdata$BP),
                              incomparables=c(""))
   #find the duplicated WOS codes
   WOS.dups <- duplicated(Xdata$UT, incomparables=c(""))
   #extract the years (PY) from the data where there
   #are no duplicates
   #results[,i] <- hist(Xdata[!(doi.dups | journal.dups | WOS.dups),"PY"],
   #     breaks=seq(years[1]-0.5,years[nyears]+0.5,1), main="",
   #     xlab="Year", ylab="Citations", col="gray50")$counts
   results[,i] <- tabulate(factor(Xdata[!(doi.dups | journal.dups | WOS.dups),"PY"], levels = years)) #alt way of getting results
   ##   mtext(side=3,line=-1,software[i], cex=1.3)
}
write.csv(file="CountsByYear.csv",results)

## #create plot of citations
## pdf("Draft.pdf",width=12,height=12)
## mat <- matrix(1:nsoftware, nrow=nsoftware, ncol=1)
## par(mar=c(0,0,0,0), oma=c(3,3,1,1))
## ylims <- 110+apply(results, MARGIN=2, FUN=max)
## layout(mat=mat, widths=1, heights=ylims)
## for (i in 1:nsoftware) {
##    if (i==nsoftware) {
##       names.arg <- years
##    } else {
##       names.arg <- NA
##    }
##    barplot(results[,i], col="gray", ann=F, axes=F, names.arg=names.arg,
##            xaxs="i", yaxs="i", ylim=c(0,ylims[i]))
##    box()
##    axis(side=2, las=1, at=seq(0,600,100))
##    mtext(side=3, line =-1.7, outer=F, software[i], cex=1)
## }
## mtext(side=1, line=3, outer=T, "Year", cex=1.3)
## mtext(side=2, line=3, outer=T, "Number of citations", cex=1.3)
## dev.off()
