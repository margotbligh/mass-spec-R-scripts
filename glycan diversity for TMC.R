# Figure glycan diversity for TMC 05May2021

require(xcms)
require(dplyr)

rawdata<-readMSData("C:/Users/admin/ownCloud/Deep-water kelp/Presentations/Thesis committee meeting 11May2021/MS31_20210205_procA_Pos531_2_33b_11.mzML", 
                    mode="onDisk")

rtrange=c(800,1400)
dp2<-chromatogram(rawdata, mz= c(562.2, 562.4), rt=rtrange)
dp3<-chromatogram(rawdata, mz= c(724.25, 724.45), rt=rtrange)
dp4<-chromatogram(rawdata, mz= c(886.3, 886.5), rt=rtrange)
dp5<-chromatogram(rawdata, mz=c(1048.35, 1048.55), rt=rtrange)

ggplot(dp2, aes(x=))

png(filename="C:/Users/admin/ownCloud/Deep-water kelp/Presentations/Thesis committee meeting 11May2021/Fig2_oligosaccharides_POS531-2_Ex33b_procA 05May2021.png",
    height = 10, width=8.9, units = "cm", res = 200)
par(mfrow=c(4,1), mai=c(0.2,0.7,0,0.1), omi=c(0.3,0,0,0))
plot(dp2, xaxt='n', main=NA, las=2, ylab='', yaxt='n')
axis(1, seq(800,1400,100), labels = NA)
axis(side=2, at=c(0, 5e+5, 1e+6),
     labels=expression(0, 5%*%10^5, 1%*%10^6), las=2)
title(ylab="Intensity (a.u.)", line=4, cex.lab=1)


plot(dp3, xaxt='n', main=NA, las=2, ylab='', yaxt='n', ylim=c(0, 4e+5))
axis(1, seq(800,1400,100), labels = NA)
axis(side=2, at=c(0, 2e+5, 4e+5),
     labels=expression(0, 2%*%10^5, 4%*%10^5), las=2)
title(ylab="Intensity (a.u.)", line=4, cex.lab=1)


plot(dp4, xaxt='n', main=NA, las=2, ylab='', yaxt='n', ylim=c(0, 8e+4))
axis(1, seq(800,1400,100), labels = NA)
axis(side=2, at=c(0, 4e+4, 8e+4),
     labels=expression(0, 4%*%10^4, 8%*%10^4), las=2)
title(ylab="Intensity (a.u.)", line=4, cex.lab=1)


plot(dp5, main=NA, ylab='', xaxt='n', xlab='', 
     yaxt='n', ylim=c(0, 8e+3))
axis(1, seq(800,1400,100))
axis(side=2, at=c(0, 4e+3, 8e+3),
     labels=expression(0, 4%*%10^3, 8%*%10^3), las=2)
title(ylab="Intensity (a.u.)", line=4, cex.lab=1)
mtext("Retention time (sec)", side=1, outer=T, at=0.56, cex = 0.7, line = 0.8)
dev.off()
