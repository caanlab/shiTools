library(readxl)
library(effectsize)
library(plyr)
library(magrittr)
library(crayon)
library(ppcor)
library(boot)
library(rms)
library(fastDummies)
library(emmeans)
library(geepack)
library(reshape2)
library(effects)
library(ggeffects)
library(lme4)
library(lmerTest)
library(Amelia)
library(car)
library(scales)
library(tidyverse)

# library(merTools)

####

shiRBootBcaConfInt <- function(
		theta,theta_boot,theta_jack,alpha=0.05
) {
	
	# theta				: scalar value of parameter estimated from original data
	# theta_boot	: vector of parameter estimates from bootstrap
	# theta_jack	: vector of parameter estimates from jackknife
	# alpha				: alpha level (default=0.05, i.e. 95% BCa CI)
	#
	# citation: 
	# Efron, B. (1987). Better bootstrap confidence intervals. Journal of the American statistical Association, 82(397), 171-185.
	
	theta_jack_diff <- mean(theta_jack) - theta_jack
	
	a <- sum(theta_jack_diff^3) / (6 * sum(theta_jack_diff^2) ^ 1.5)
	
	z0 <- qnorm(mean(theta_boot <= theta))
	
	z_alpha1 <- qnorm(alpha/2)
	z_alpha2 <- qnorm(1-alpha/2)
	
	alpha1 <- pnorm( z0 + (z0 + z_alpha1) / (1 - a * (z0 + z_alpha1)) )
	alpha2 <- pnorm( z0 + (z0 + z_alpha2) / (1 - a * (z0 + z_alpha2)) )
	
	return( quantile(theta_boot, probs=c(alpha1,alpha2)) )
	
}

####

shiRDescriptive <- function(
		df, GrpVar=NULL
) {
	FUNC_LIST <- list(
		"N"    = function(x)sum(!is.na(x)),
		"MEAN" = function(x)mean(x,na.rm=TRUE),
		"SD"   = function(x)sd(x,na.rm=TRUE),
		"SE"   = function(x)(sd(x,na.rm=TRUE))/sqrt(sum(!is.na(x))),
		"MED"  = function(x)quantile(x,probs=0.50,na.rm=TRUE),
		"Q25"  = function(x)quantile(x,probs=0.25,na.rm=TRUE),
		"Q75"  = function(x)quantile(x,probs=0.70,na.rm=TRUE),
		"MIN"  = function(x)min(x,na.rm=TRUE),
		"MAX"  = function(x)max(x,na.rm=TRUE),
		"SUM"  = function(x)sum(x,na.rm=TRUE)
	)
	if (length(GrpVar)==0) {
		df %<>% lapply( shiRFact2Num ) %>% tibble::as.tibble()
		df %<>% dplyr::select(-c(names(df)[sapply(df,function(x)all(is.na(x)))]))
		Descr <- data.frame(Var=names(df))
		for (i in 1:length(FUNC_LIST)) {
			Descr[[names(FUNC_LIST)[i]]]<- sapply(df,FUNC_LIST[[i]])
		}
	} else {
		BY_LIST <- df %>% dplyr::select(!!sym(GrpVar))
		df %<>% dplyr::select(-!!sym(GrpVar)) %>% lapply( shiRFact2Num ) %>% tibble::as.tibble()
		df %<>% dplyr::select(-c(names(df)[sapply(df,function(x)all(is.na(x)))]))
		Descr <- aggregate(x=df, by=BY_LIST, FUN=FUNC_LIST[[1]]) %>% reshape2::melt(id.vars=GrpVar, value.name=names(FUNC_LIST)[1])
		for (i in 2:length(FUNC_LIST)) {
			Descr %<>% cbind( aggregate(x=df, by=BY_LIST, FUN=FUNC_LIST[[i]]) %>% reshape2::melt(id.vars=GrpVar, value.name=names(FUNC_LIST)[i]) %>% dplyr::select(-c(GrpVar,"variable")) )
		}
		Descr %<>% dplyr::rename("Var" = "variable")
	}
	return(Descr)
}

####

shiRDescriptivePrint <- function(
		df, # data frame
		GrpVar = NULL, # which group variable(s)
		Ignore = NULL, # which variables to be ignored (e.g. subject IDs)
		DoTest = TRUE, # whether to run group comparisons (always false unless exactly 2 groups and NumTest and CatTest are correctly specified)
		NumTest = "wilcox.test", # or t.test
		CatTest = "fisher.test", # or chisq.test
		digits = 2 # for rounding purpose, must be non-negative integer, may be Inf for no rounding
) {
	
	sym=dplyr::sym
	
	if (digits == Inf) {
		ROUND <- function(num)sprintf("%g",num)
	} else {
		ROUND <- function(num)sprintf(paste0("%.",digits,"f"),num)
	}
	
	GROUP = ".shiRDescriptivePrint.Grp"
	
	if (length(GrpVar)==0) {
		df %<>% dplyr::mutate(!!sym(GROUP):="All")
		df[[GROUP]] %<>% factor()
	} else if (length(GrpVar)==1) {
		df %<>% dplyr::mutate(!!sym(GROUP):=!!sym(GrpVar[1]))
	} else {
		df %<>% dplyr::mutate(!!sym(GROUP):=!!sym(GrpVar[1]))
		for (i in 2:length(GrpVar)) {
			df[[GROUP]] %<>% interaction(df[[GrpVar[i]]])
		}
	}
	
	df[[GROUP]] <- droplevels(df[[GROUP]])
	NGRP <- nlevels(df[[GROUP]])
	GRPLEV <- levels(df[[GROUP]])
	
	DOTEST <- DoTest && NGRP == 2
	
	VarList <- names(df) %>% setdiff(c(GrpVar,GROUP,Ignore))
	
	for (i in 1:length(VarList)) {
		if (!is.numeric(df[[VarList[i]]])) {
			df[[VarList[i]]] %<>% factor()
		}
	}
	ISFACTOR <- sapply(df %>% dplyr::select(VarList),is.factor)
	
	# header
	
	N <- dplyr::count(df,!!sym(GROUP))$n
	TABLE <- "Variable"
	for (grp in 1:NGRP) {
		xTAB <- GRPLEV[grp] %>% paste0(" (n=",N[grp],")")
		TABLE %<>% c(xTAB)
	}
	
	# each variable
	for (i in 1:length(VarList)) {
		xdf <- df %>% dplyr::filter(!is.na(!!sym(VarList[i])))
		# variable name
		xTAB <- c(VarList[i], rep("",NGRP))
		TABLE %<>% rbind(xTAB)
		# no. of subjects
		xN <- dplyr::count(xdf, !!sym(GROUP))$n
		if (!identical(N,xN)) {	TABLE %<>% rbind(c("No. of participants", xN)) }
		# stats for factor
		if (ISFACTOR[i]) {
			# column 0
			xTAB <- paste0(levels(df[[VarList[i]]]),", No. (%)")
			# column for each group
			for (grp in 1:NGRP) {
				xTAB1 <- dplyr::count(xdf %>% dplyr::filter(!!sym(GROUP)==GRPLEV[grp]), !!sym(VarList[i]))$n
				xTAB1 %<>% paste0(" (",(xTAB1/xN[grp]*100) %>% ROUND,"%)")
				xTAB %<>% cbind(xTAB1)
			}
			# stats for numerical
		} else {
			# column 0
			xTAB <- c("Mean (SD)","Median (IQR)")
			# column for each group
			for (grp in 1:NGRP) {
				xx <- xdf[[VarList[i]]][xdf[[GROUP]]==GRPLEV[grp]]
				xTAB1 <- c(mean(xx) %>% ROUND, median(xx) %>% ROUND)
				xTAB2 <- c(sd(xx) %>% ROUND,paste0(quantile(xx,.25) %>% ROUND,"–",quantile(xx,.75) %>% ROUND))
				xTAB1 %<>% paste0(" (", xTAB2, ")")
				xTAB %<>% cbind(xTAB1)
			}
		}
		# add to main table
		TABLE %<>% rbind(xTAB)
	}
	
	rownames(TABLE) <- NULL
	colnames(TABLE) <- NULL
	
	if ((tolower(NumTest) != "t.test" && tolower(NumTest) != "wilcox.test") || (tolower(CatTest) != "chisq.test" && tolower(CatTest) != "fisher.test")) {
		DOTEST = FALSE
	}
	if (!DOTEST) {
		TABLE0 <- TABLE
		TABLE0[!sapply(TABLE0[,2],function(x)nchar(x)==0),1] %<>% {paste0("  ",.)}
		TABLE0[1,1] <- "Variable"
		print.table(TABLE0)
		return(TABLE)
	}
	
	if (tolower(NumTest) == "t.test") {
		NUMTEST = function(x,y)t.test(x,y)
	} else if (tolower(NumTest) == "wilcox.test") {
		NUMTEST = function(x,y)wilcox.test(x,y)
	}
	
	if (tolower(CatTest) == "chisq.test") {
		CATTEST = function(x,y)chisq.test(x,y)
	} else if (tolower(CatTest) == "fisher.test") {
		CATTEST = function(x,y)fisher.test(x,y)
	}
	
	xTAB <- "p-value"
	xCNT <- 0
	for (i in 2:nrow(TABLE)) {
		if (TABLE[i,2] != "") {xTAB %<>% c("")} else {
			xCNT <- xCNT+1
			xVAR <- TABLE[i,1]
			if (ISFACTOR[xCNT]) {
				xSTAT <- CATTEST(df[[GROUP]],df[[xVAR]])
			} else {
				xSTAT <- NUMTEST(df[[xVAR]][df[[GROUP]]==GRPLEV[1]],df[[xVAR]][df[[GROUP]]==GRPLEV[2]])
			}
			xTAB %<>% c(paste0("p",shiRFormatPval(xSTAT$p.value)))
		}
	}
	
	TABLE %<>% cbind(xTAB)
	
	rownames(TABLE) <- NULL
	colnames(TABLE) <- NULL
	
	TABLE0 <- TABLE
	TABLE0[!sapply(TABLE0[,2],function(x)nchar(x)==0),1] %<>% {paste0("  ",.)}
	TABLE0[1,1] <- "Variable"
	print.table(TABLE0)
	return(TABLE)
	
}

####

shiRFact2Num <- function(
		xFact
) {
	if (is.factor(xFact))  {
		if (lapply(levels(xFact),is.character) %>% unlist %>% any) {
			xFact <- as.numeric(xFact)
		} else {
			xFact <- as.numeric(as.character(xFact))
		}
	}
	return(xFact)
}

####

shiRFormatPval <- function(
    pval,
    sep=""
) {
  pstring <- rep("",length(pval))
  for (i in 1:length(pval)) {
    pstring[i] <- shiRIf(
      pval[i]>=0.995,paste(">","0.99",sep=sep),shiRIf(
        pval[i]>=0.0995,paste("=",format(round(pval[i], 2), nsmall = 2),sep=sep),shiRIf(
          pval[i]>=0.001,paste("=",format(round(pval[i], 3), nsmall = 3),sep=sep),paste("<","0.001",sep=sep)
        )))
  }
  return(pstring)
}

####

shiRFormula <- function(
    y,
    x
) {
  pt1 <- paste(y,"~")
  pt2 <- paste(x,collapse="+")
  return(as.formula(paste(pt1,pt2)))
}

####

shiRLogisticReg <- function(
    formula,
    data = environment(formula),
    ...
) {
  
  xMdl1 <- rms::lrm(formula, data = data, ...) # get R2, Coef
  xMdl2 <- glm(formula, data = data, family = "binomial", ...) # get LR X2 and p-value
  xCI <- confint(xMdl2) %>% suppressMessages
  xAnv <- car::Anova(xMdl2,type=3)
  
  OUT <- list("stats"=NULL,"coefficients"=NULL,"mdl_glm"=NULL,"mdl_lrm"=NULL)
  
  OUT$stats <- xMdl1$stats[c("Obs","Model L.R.","d.f.","P","R2","C")]
  names(OUT$stats) <- c("N","LR_Chi2","df","pval","Nagelkerke_R2","ROC_AUC")
  
  OUT$coefficients <- data.frame(
    B = xMdl1$coefficients,
    B_SE = xMdl1$var %>% {.[row(.)==col(.)]} %>% sqrt(),
    B_CI_2.5 = xCI[,1],
    B_CI_97.5 = xCI[,2],
    ExpB = xMdl1$coefficients %>% exp(),
    ExpB_CI_2.5 = xCI[,1] %>% exp(),
    ExpB_CI_97.5 = xCI[,2] %>% exp(),
    LR_Chi2 = c(NaN,xAnv[,1]),
    Df = c(NaN,xAnv[,2]),
    pval = c(NaN,xAnv[,3])
  )
  
  OUT$mdl_glm <- xMdl2
  OUT$mdl_lrm <- xMdl1
  
  return(OUT)
  
}

####

shiRLogisticReg_PurpSelect <- function(
    data,
    y,
    x,
    xForceIncl = NULL,
    cutoff_CandSelect_p = 0.25,
    cutoff_CandRetain_p = 0.1,
    cutoff_CandRetain_pct = 15,
    cutoff_NonCandIncl_p = 0.15,
    ...
) {
  
  # deal with input
  
  if (cutoff_CandRetain_pct < 1) {
    warning(paste("cutoff_CandRetain_pct = ", cutoff_CandRetain_pct, "* 100 =", cutoff_CandRetain_pct*100))
  }
  
  for (i in 1:length(xForceIncl)) {
    if (length(xForceIncl)==0) {break}
    if (!(xForceIncl[i] %in% x)) {
      x <- c(x,xForceIncl[i])
    }
  }
  
  cat(crayon::yellow("\ncutoff_CandSelect_p  ","="),cutoff_CandSelect_p)
  cat(crayon::yellow("\ncutoff_CandRetain_p  ","="),cutoff_CandRetain_p)
  cat(crayon::yellow("\ncutoff_CandRetain_pct","="),cutoff_CandRetain_pct)
  cat(crayon::yellow("\ncutoff_CandSelect_p  ","="),cutoff_NonCandIncl_p)
  cat("\n\n")
  
  param <- list(
    "y"=y,
    "x"=x,
    "xForceIncl"=xForceIncl,
    "cutoff_CandSelect_p"=cutoff_CandSelect_p,
    "cutoff_CandRetain_p"=cutoff_CandRetain_p,
    "cutoff_CandRetain_pct"=cutoff_CandRetain_pct,
    "cutoff_NonCandIncl_p"=cutoff_NonCandIncl_p
  )
  
  # dummy code
  
  for (i in 1:length(x)) {
    if (is.factor(data[[x[i]]]) && nlevels(data[[x[i]]])>2) {
      data_dumm <- fastDummies::dummy_cols(data, select_columns = x[i])
      x_dumm <- variable.names(data_dumm)[!(variable.names(data_dumm) %in% variable.names(data))]
      x <- x[!(x %in% x[i])]
      x <- c(x,x_dumm)
    }
  }
  
  # STEP 1: UNIVARIABLE
  
  Step1 <- data.frame()
  
  for (i in 1:length(x)) {
    fml <- shiRFormula(y, x[i])
    OUT <- shiRLogisticReg(formula = fml, data = data, ...)
    stopifnot(nrow(OUT$coefficients) == 2)
    Step1 <- rbind(Step1,OUT$coefficients[2,])
  }
  
  Step1 <- cbind(data.frame(Var=x), Step1)
  stopifnot(!any(is.na(Step1$pval)))
  Step1 <- Step1[order(Step1$pval),]
  Step1$Candidate <- (Step1$pval < cutoff_CandSelect_p) | (Step1$Var %in% xForceIncl)
  Step1.compact <- Step1[,c("Var","B","B_SE","LR_Chi2","Df","pval","Candidate")]
  
  x_cand <- Step1$Var[Step1$Candidate]
  x_noncand <- Step1$Var[!Step1$Candidate]
  
  Step1.log <- NULL
  
  for (i in 1:length(Step1$Var)) {
    if (Step1$Candidate[i]) {
      Step1.log <- rbind(Step1.log,paste("Include: ",Step1$Var[i],sep=""))
    } else {
      Step1.log <- rbind(Step1.log,paste("Exclude: ",Step1$Var[i],sep=""))
    }
  }
  
  cat(crayon::yellow("Step1:\n  "))
  Step1.log %>% cat(sep="\n  ")
  
  if (length(x_cand)==0) {
    return(list(
      "param"=param,
      "Step1"=Step1,
      "Step1.compact"=Step1.compact,
      "Step1.log"=Step1.log))
  }
  
  
  # STEP 2: MULTIVARIABLE
  
  Step2 <- data.frame()
  
  iter <- 0
  fml <- shiRFormula(y, x_cand)
  OUT <- shiRLogisticReg(formula = fml, data = data, ...)
  
  xVar <- data.frame(Var = x_cand)
  xCoef <- OUT$coefficients %>% {.[2:nrow(.),]}
  xPctCoefCh <- data.frame(PctCoefCh = rep(0,nrow(xCoef)))
  xIter <- data.frame(Iter = rep(paste(iter," ",sep=""),nrow(xCoef)))
  
  xTAB <- cbind(xIter,xVar,xCoef,xPctCoefCh)
  xTAB <- xTAB[order(-xTAB$pval),]
  Step2 <- rbind(Step2,xTAB)
  
  xTAB0 <- xTAB
  
  nvar <- 1
  xMaxPval <- xTAB0$pval[nvar]
  xx_pick <- xTAB0$Var[nvar]
  xx_safe <- setdiff(xTAB0$Var,xx_pick)
  
  Step2.log <- NULL
  
  while ( xMaxPval >= cutoff_CandRetain_p && length(xx_safe) > 0 && nvar <= nrow(xTAB0) && iter<1000 ) {
    
    if (xx_pick %in% xForceIncl) {
      nvar <- nvar + 1
      if (nvar<nrow(xTAB0)) {
        xMaxPval <- xTAB0$pval[nvar]
        xx_pick <- xTAB0$Var[nvar]
        xx_safe <- setdiff(xTAB0$Var,xx_pick)
      }
      next
    }
    
    
    iter <- iter + 1
    fml <- shiRFormula(y, xx_safe)
    OUT <- shiRLogisticReg(formula = fml, data = data, ...)
    
    xVar <- data.frame(Var = xx_safe)
    xCoef <- OUT$coefficients %>% {.[2:nrow(.),]}
    
    xPctCoefCh <- data.frame( PctCoefCh = rep(0,nrow(xCoef)) )
    rownames(xPctCoefCh) <-  rownames(xCoef)
    xPctCoefCh$PctCoefCh <- (xCoef[rownames(xPctCoefCh),"B"]/xTAB0[rownames(xPctCoefCh),"B"]-1)*100
    xIter <- data.frame(Iter = rep(paste(iter," ",sep=""),nrow(xCoef)))
    
    
    xTAB <- cbind(xIter,xVar,xCoef,xPctCoefCh)
    xTAB <- xTAB[order(-xTAB$pval),]
    
    if ( any(abs(xPctCoefCh)>=cutoff_CandRetain_pct) ) {
      
      Step2.log <- rbind(Step2.log,paste("Iter ",iter,": include ",xx_pick,sep=""))
      
      nvar <- nvar + 1
      
      if (nvar<nrow(xTAB0)) {
        xMaxPval <- xTAB0$pval[nvar]
        xx_pick <- xTAB0$Var[nvar]
        xx_safe <- setdiff(xTAB0$Var,xx_pick)
      }
      
    } else {
      
      Step2.log <- rbind(Step2.log,paste("Iter ",iter,": exclude ",xx_pick,sep=""))
      
      xTAB0 <- xTAB
      
      nvar <- 1
      xMaxPval <- xTAB0$pval[nvar]
      xx_pick <- xTAB0$Var[nvar]
      xx_safe <- setdiff(xTAB0$Var,xx_pick)
      
    }
    
    Step2 <- rbind(Step2,xTAB0)
    
  }
  
  Step2.compact <- Step2[,c("Iter","Var","B","B_SE","LR_Chi2","Df","pval","PctCoefCh")]
  
  cat(crayon::yellow("\nStep2:\n  "))
  Step2.log %>% cat(sep="\n  ")
  
  if (length(x_noncand)==0) {
    return(list(
      "param"=param,
      "Step1"=Step1,
      "Step1.compact"=Step1.compact,
      "Step1.log"=Step1.log,
      "Step2"=Step2,
      "Step2.compact"=Step2.compact,
      "Step2.log"=Step2.log
    ))
  }
  
  # STEP 3: NONCANDIDATE
  
  Step3 <- data.frame()
  
  xx_safe <- xTAB0$Var
  iter <- 0
  
  fml <- shiRFormula(y, xx_safe)
  OUT <- shiRLogisticReg(formula = fml, data = data, ...)
  xVar <- data.frame(Var = xx_safe)
  xCoef <- OUT$coefficients %>% {.[2:nrow(.),]}
  xPctCoefCh <- data.frame( PctCoefCh = rep(0,length(xx_safe)) )
  xIter <- data.frame(Iter = rep(paste(iter," ",sep=""),nrow(xCoef)))
  xTAB0 <- cbind(xIter,xVar,xCoef,xPctCoefCh)
  
  Step3 <- rbind(Step3,xTAB0)
  
  Step3.log <- NULL
  
  for (i in 1:length(x_noncand)) {
    
    iter <- iter + 1
    
    fml <- shiRFormula(y, c(xx_safe,x_noncand[i]))
    OUT <- shiRLogisticReg(formula = fml, data = data, ...)
    
    xVar <- data.frame(Var = c(xx_safe,x_noncand[i]))
    xCoef <- OUT$coefficients %>% {.[2:nrow(.),]}
    
    xPctCoefCh <- data.frame( PctCoefCh = rep(0,nrow(xCoef)) )
    rownames(xPctCoefCh) <-  rownames(xCoef)
    xPctCoefCh[rownames(xTAB0),1] <- xTAB0$B/xCoef[rownames(xTAB0),"B"]
    xIter <- data.frame(Iter = rep(paste(iter," ",sep=""),nrow(xCoef)))
    
    xTAB <- cbind(xIter,xVar,xCoef,xPctCoefCh)
    xTAB <- xTAB[order(xTAB$pval),]
    
    if ( any(abs(xPctCoefCh)>=cutoff_CandRetain_pct) || xTAB[x_noncand[i],"pval"]<cutoff_NonCandIncl_p ) {
      
      Step3.log <- rbind(Step3.log,paste("Iter ",iter,": include ",x_noncand[i],sep=""))
      
      xx_safe <- c(xx_safe,x_noncand[i])
      xTAB0 <- xTAB
      
    } else {
      
      Step3.log <- rbind(Step3.log,paste("Iter ",iter,": exclude ",x_noncand[i],sep=""))
      
    }
    
    Step3 <- rbind(Step3,xTAB)  
    
  }
  
  Step3 <- rbind(Step3,xTAB0)
  
  Step3.compact <- Step3[,c("Iter","Var","B","B_SE","LR_Chi2","Df","pval","PctCoefCh")]
  
  cat(crayon::yellow("\nStep3:\n  "))
  Step3.log %>% cat(sep="\n  ")
  
  return(list(
    "param"=param,
    "Step1"=Step1,
    "Step1.compact"=Step1.compact,
    "Step1.log"=Step1.log,
    "Step2"=Step2,
    "Step2.compact"=Step2.compact,
    "Step2.log"=Step2.log,
    "Step3"=Step3,
    "Step3.compact"=Step3.compact,
    "Step3.log"=Step3.log
  ))
  
}

####

shiRCat <- function(
    ...
) {
  cat("\n")
  cat(rep(c(crayon::bgCyan(" "),crayon::bgMagenta(" "),crayon::bgYellow(" "),crayon::bgRed(" "),crayon::bgGreen(" "),crayon::bgBlue(" ")),10),sep="")
  cat("\n")
  cat(crayon::bgWhite(crayon::black(" ")))
  cat(crayon::bgWhite(crayon::black(...)))
  cat(crayon::bgWhite(crayon::black(" ")))
  cat("\n")
}

####

shiRIf <- function(
    TF,ifT,ifF
) {
  if (TF) {
    return(ifT)
  } else {
    return(ifF)
  }
}

####

shiRResid <- function(
    data=NULL,
    covname,
    ignore=NULL,
    ...
) {
  covname <- intersect(covname,names(data))
  dataNew <- data[,NULL]
  for (i in names(data)) {
    if ((!is.numeric(data[[i]])) || (i %in% ignore)) {
      dataNew[[i]] <- data[[i]]
      next
    }
    if (i %in% covname) {
      next
    }
    dataTemp <- data.frame(y=data[,i],data[,covname])
    xxFit <- glm(y~., data=dataTemp, ...) %>% predict(., dataTemp, na.action=na.pass, ...)
    xresid <-  xxFit %>% {dataTemp$y - . + mean(., na.rm=TRUE)}
    dataNew[[i]] <- xresid
  }
  return(dataNew)
}