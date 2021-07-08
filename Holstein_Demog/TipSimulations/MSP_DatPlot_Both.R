# 26th September
# Reading in and printing results for Python Holstein Simulations

# Update 27th Feb 2019
# Plotting both high- and low-N0 estimates

# Update 1st Apr 2018
# Plotting separate high and low N0 figures

# 11th Feb 2020
# Setting input samples to 204.

# Function to turn log-log fit into actual function
# For prediction
xyll <- function(x,int,coeff){
	y <- exp(int)*(x^coeff)
	return(y)
}

# Mean tip length
# Function for plotting with error bars;
# Also adding log-log line of best fit
# And prediction of tip-age, given number of samples
xy.error.bars <- function(x,y,ybot,ytop,xin,xt,yt){
	dfit <- lm(log(y)~log(x))
	dpred <- predict(dfit)
	par(mar=c(5,7,4,2)+0.1)
	plot(x, y, pch=16, ylim=c(min(y-ybot),max(y+ytop)),log="xy",xlab="Haploid Sample Size",ylab="",xaxt="n",yaxt="n")
	axis(1,at=c(10,20,50,100,250,500,1000))
	axis(2,at=c(10,50,100,500,1000,5000,10000),las=2)
	ymx <- max(y+ytop)
	ymn <- min(y-ybot)
	ylab <- exp((log(max(ymx)) - log(min(ymn)))/2 + log(min(ymn)))
	text(x=3.75,y=ylab,labels=paste("Mean","Tip","Age",sep="\n"),xpd=NA)
	arrows(x, y-ybot, x, y+ytop, code=3, angle=90, length=0.1)
	lines(exp(dpred)~x)
	xout <- round(xyll(xin,dfit$coefficients[1],dfit$coefficients[2]))
	segments(xin,min(y-ybot)/2,xin,xout,lty=2)
	segments(min(x)/2,xout,xin,xout,lty=2)
	text(xt,yt,sprintf("%d samples covers %d generations",xin,xout))
}

dat <- read.table("MSP_Results.dat")
dat2 <- read.table("MSP_Results_Low.dat")
pdf(file="TipAgeHOLHighN0.pdf",width=8,height=8)
xy.error.bars(dat$V1,dat$V2,abs(dat$V3-dat$V2),(dat$V4-dat$V2),2*102,225,5000)
dev.off()
pdf(file="TipAgeHOLLowN0.pdf",width=8,height=8)
xy.error.bars(dat2$V1,dat2$V2,abs(dat2$V3-dat2$V2),(dat2$V4-dat2$V2),2*102,225,3000)
dev.off()
