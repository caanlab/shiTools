# PLOT -----

## get ribbon from effects:Effect ----
## remove random intercepts ----
## get ggplot colors ----

# source("~/GoogleDrive/# Documents/2021/20210901_VivitrolCorticalThickness/Data_Prep.R")
# source("~/GoogleDrive/# Documents/2021/20210901_VivitrolCorticalThickness/Early/Data_Stat_IPW.R") # just run PREP

Tb <- data.frame(
	Brain = X$Brain,
	Crav = X$Crav,
	Session = X$Session,
	invwt = X$invwt,
	sqrtinvwt = sqrt(X$invwt),
	Subject = X$Subject
)

xMdl <- lmer(Brain ~ Session*Crav + (1|Subject), data = Tb, weights = Tb$invwt)
xFx <- Effect(c("Crav","Session"), mod = xMdl, xlevels = list(Crav = seq(from=0,to=9,by=.05))) %>% as.data.frame() # get ribbon and stuff

x <- coef(xMdl)$Subject # get random intercepts, one per subject
x <- cbind(row.names(x) %>% as.data.frame(), x$`(Intercept)` %>% as.data.frame()) # combine subjects (from x rownames) and intercepts
names(x) <- c("Subject","Const")
Tb <- left_join(Tb,x,by="Subject")
Tb$BrainAdj <- Tb$Brain - Tb$Const + mean(x$Const) # remove random intercepts, add back overall mean

xxGp <- ggplot() +
	geom_point(dat = Tb, aes(x=Crav,y=BrainAdj,color=Session,shape=Session,size=sqrtinvwt),alpha=.7,position=position_jitter(w=0.1)) +
	scale_shape_manual(values = c(16,16,16)) +
	scale_color_manual(values = scales::hue_pal()(3)) + # gets default colors
	scale_fill_manual(values = scales::hue_pal()(3)) +
	geom_line(data = xFx, aes(x=Crav,y=fit,color=Session)) +
	geom_ribbon(data = xFx, aes(x=Crav,ymin=lower,ymax=upper,fill=Session), alpha=.3)
print(xxGp)

## add significance to bar plots ----
## break y-axis ----

xDat <- LmeEmm %>% summary() %>% dplyr::select(Session,emmean,SE)
xDat$Session <- c("pre-treatment","on-treatment","post-treatment") %>% {factor(.,levels=.)}
xDat$emmean <- xDat$emmean-2.6 # prepare to start from zero

TEXT_12 <- "p<0.05 (whole-brain)"
TEXT_13 <- shiRFormatPval(summary(LmeEmmCon)$p.value[1])
TEXT_23 <- shiRFormatPval(summary(LmeEmmCon)$p.value[2])

YLIM <- c(0,0.48)
YGAP <- YLIM[2]-YLIM[1]
YDOWNSHIFT <- 0.015
YGAPFACTOR <- 1.1/13
XBAREXTEND <- 0.1

plot_FontSize_Note <- 18

xxGp <- ggplot() +
	geom_bar(dat=xDat, aes(x=Session, y=emmean, fill=Session),stat="identity",alpha=1,width=0.8) +
	geom_errorbar(dat=xDat, aes(x=Session, y=emmean, ymin=emmean-SE, ymax=emmean+SE, xmin=NA, xmax=NA), width=.2) +
	scale_fill_manual(values = c("#F8766D", "#00BA42", "#7997FF")) +
	coord_cartesian(ylim = YLIM) +
	scale_x_discrete(limits=xDat$Session) +
	scale_y_continuous(
		breaks = seq(from=0,to=0.5,by=0.1),
		labels = c("0.0"="0.0","0.1"="2.7","0.2"="2.8","0.3"="2.9","0.4"="3.0","0.5"="3.1") # add 2.6 back except zero for y-axis breaking in photoshop
	) +
	annotate("segment", x=1-XBAREXTEND,xend=2+XBAREXTEND, y=YLIM[2]-0*(YGAP*YGAPFACTOR)-YDOWNSHIFT,yend=YLIM[2]-0*(YGAP*YGAPFACTOR)-YDOWNSHIFT) + # add lines and stars for pairwise comparisons
	annotate("segment", x=1-XBAREXTEND,xend=3+XBAREXTEND, y=YLIM[2]-1*(YGAP*YGAPFACTOR)-YDOWNSHIFT,yend=YLIM[2]-1*(YGAP*YGAPFACTOR)-YDOWNSHIFT) +
	annotate("segment", x=2-XBAREXTEND,xend=3+XBAREXTEND, y=YLIM[2]-2*(YGAP*YGAPFACTOR)-YDOWNSHIFT,yend=YLIM[2]-2*(YGAP*YGAPFACTOR)-YDOWNSHIFT) +
	annotate("text", x=(1+2)/2, y=YLIM[2]-0*(YGAP*YGAPFACTOR)-YDOWNSHIFT, label=TEXT_12, parse=FALSE, size=plot_FontSize_Note/.pt, color="black", vjust=-.5) +
	annotate("text", x=(1+3)/2, y=YLIM[2]-1*(YGAP*YGAPFACTOR)-YDOWNSHIFT, label=TEXT_13, parse=FALSE, size=plot_FontSize_Note/.pt, color="black", vjust=-.5) +
	annotate("text", x=(2+3)/2, y=YLIM[2]-2*(YGAP*YGAPFACTOR)-YDOWNSHIFT, label=TEXT_23, parse=FALSE, size=plot_FontSize_Note/.pt, color="black", vjust=-.5)

print(xxGp)

## remove legend ----
gg <- ggplot() + guides(colour = "none", shape = "none")


# LM, LME, and GEE ----

## anova notes ----

data.frame(
	FUNC = c("summary()","anova()","car::Anova(type=3)"),
	LM = c("[sen. to fact]","[type I]","[sen. to fact]"),
	LMER = c("[sen. to fact]","","[sen. to fact]"),
	GEEGLM = c("[sen. to fact]","[type I]","[not supp]")
) %>% print()

## get tangent slope ----

model <- lm(ants_diff ~ poly(x, 3), data = X)

x0 <- seq(from=0,to=1.7,by=0.01) # Calculate the predicted value of the tangent line at x = x0
y0 <- predict(model, newdata = data.frame(x = x0))

slope <- predict(model, newdata = data.frame(x = x0), deriv = 1, se.fit = TRUE, interval = "confidence") # CI
slope_pi <- predict(model, newdata = data.frame(x = x0), deriv = 1, se.fit = TRUE, interval = "prediction") # PI

## get tangent slope ----

dat = data.frame(x=rnorm(20),y=rnorm(20),m=rnorm(20))
mod = lm(y ~ (x+I(x^2)+I(x^3))*m, data=dat)

t=-10
dm1=c();
dm2=c();
for (i in 1:50) {
	t = t + 0.4;
	dm1[i] = summary(emtrends(mod,~m,"x",at=list(x=t,m=-1)))$x.trend
	dm2[i] = summary(emtrends(mod,~m,"x",at=list(x=t,m=1)))$x.trend
}
plot(dm1)
plot(dm2)

## get marginal trends from interaction between continuous variables ----

dat=data.frame(x=rnorm(20),y=rnorm(20),m1=rnorm(20),m2=rnorm(20),m3=rnorm(20))
mod=lm(y~x*m1*m2*m3,data=dat)

emtrends(mod,~m1+m2+m3,"x",at=list(m1=c(-1,1),m2=c(-1,1),m3=c(-1,1)))


# BOOTSTRAP ----

## get bca p-value ----

k <- 0
Bca.Pval <- c()
for (i in 1:length(Region)) {
	for (j in 1:length(MMVariable)) {
		k <- k+1
		STOP = F
		pMin <- 0
		pMax <- 1
		pNow <- 0.05
		pInc <- 0.00001 # stop when p-value or absolute p-value change is smaller than this
		Sig0 <- 0
		while (!STOP) {
			Sig <- prod( boot::boot.ci(BOOT,index=k,type="bca",conf=1-pNow)$bca[4:5] ) # BOOT comes from boot::boot
			pDelta <- pNow - pNow0
			STOP <- (abs(pDelta) < pInc) && (prod(Sig,Sig0)<0 || pNow < pInc)
			pNow0 <- pNow
			Sig0 <- Sig
			if (Sig<=0) { # if non-sig, next try new p-value halfway from current to max
				pMin <- pNow
			} else { # if sig, next try new p-value halfway from min to current
				pMax <- pNow
			}
			pNow <- (pMin + pMax)/2 # halfway
			# sprintf("k = %2d   pNow = %5f   pDelta = %+5f   Sig = %s   pMin = %5f   pMax = %5f\n",k,pNow0,pDelta,shiRIf(Sig>0,"+","-"),pMin,pMax) %>% cat
		}
		Bca.Pval[k] <- pNow
	}
}


# STRING IN DPLYR ----

Grp = "a.group.variable"
df %<>% dplyr::mutate(!!sym(Grp):="All")
df %<>% dplyr::select(-!!sym(Grp))