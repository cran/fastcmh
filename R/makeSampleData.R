#' Create sample data for FastCMH
#'
#' This function creates sample data for use with the FastCMH method.
#'
#' @param folder The folder in which the data will be saved. Default is current directory "./".
#' @param xfilename The name of the data file. Default is "data.txt".
#' @param yfilename The name of the label file. Default is "label.txt".
#' @param covfilename The name of the file containing the covariate categories. This file actually just contains K numbers, where K is the number of covariates. Default is "cov.txt".
#' @param K The number of covariates (a positive integer). Default is K=2.
#' @param L The number of features (length of each sequence). Default is L=1000.
#' @param n The number of samples (cases and controls combined). Default is n=200, i.e. 100 cases and 100 controls.
#' @param noiseP The background noise in the data (as a probability of 0/1 being flipped). Default is noiseP=0.3.
#' @param corruptP The probability of data corruption: each bit has probability correuptP of being flipped. Default is corruptP=0.05.
#' @param rho The strength of the confounding in the confounded interval (as a probability). Default is rho=0.8 (i.e. a very strong signal).
#' @param tau1 The location of the significant interval (starting point). Default value is tau1=100.
#' @param taulength1 The length of the significant interval. Default value is taulength1=4, so default significant interval is [100, 103].
#' @param tau2 The location of the confounded significant interval (starting point). Default value is tau2=200.
#' @param taulength2 The length of the confounded significant interval. Default value is taulength2=4, so default significant interval is [200, 203].
#' @param seednum The seed used for generating the data. Default value is sednum=2.
#' @param truetaufilename The file where the location of the true significant intervals are saved (as opposed to the detected significant intervals). Default is "truetau.txt".
#' @param showOutput Flag to decide whether or not to show output, where files are created, their names, etc. Default is FALSE
#'
#' @section See Also:
#' 'FastCMH'
#'
#'
#'
#' @examples
#' 
#' #make a small sample data set, using the default parameters
#' makeSampleData("~/sampledata/", showOutput=TRUE)
#' 
#' #make a very small sample data set
#' makeSampleData("~/sampledata/small/", n=20, L=10, tau1=2, taulength1=2, tau2=6, taulength2=2)
#' 
#' 
#' @export
makeSampleData <- function(folder="./", xfilename="data.txt", yfilename="label.txt", covfilename="cov.txt", K=2, L=1000, n=200, noiseP=0.3, corruptP=0.05, rho=0.8, tau1=100, taulength1=4, tau2=200, taulength2=4, seednum=2, truetaufilename="truetau.txt", showOutput=FALSE){



    #check that bindata is installed
    if (!requireNamespace("bindata", quietly = TRUE)) {
      stop("Package 'bindata' needed for this function to work. Please install it.", call. = FALSE)
    }

    #----------------------------------------------------------------#

    #bindata is installed; can start to make sample data
    if (seednum > 0){
        set.seed(seednum)
    }

    if (showOutput){
        cat("#--------------------------------------------------------#\n")
        cat("FASTCMH: CREATING SAMPLE DATA\n")
        cat("#--------------------------------------------------------#\n")
        if (seednum > 0){
                cat("setting seed: ", seednum, "\n")
        } else{
            cat("not setting seed\n")
        }
    }
    #now need to normalise path for windows
    folder <- fixSeparator(folder)

    #make sure folder has "/" at end:
    folder <- checkdir(folder)


    #need to preappend folder:
    xfilename = paste0(folder, xfilename)
    yfilename = paste0(folder, yfilename)
    covfilename = paste0(folder, covfilename)
    truetaufilename = paste0(folder, truetaufilename)

    eachCov <- trunc(n/2)

    truetau <- rbind( c(tau1, tau1+taulength1-1))
                        


    ####second method
    ##marginal probs
    p1 <- 0.5
    p2 <- 0.5
    p3 <- 0.5
    marg <- c(p1, p2, p3)
#    epsval <- 0.05
    epsval <- corruptP

    sigma <- cbind(c(1,0,rho/2),c(0,1,rho/2),c(rho/2,rho/2,1))
    #using bindata package
    x <- bindata::rmvbin(n,margprob=marg,bincorr=sigma)
    ylabels <- x[, 3]
    covlabels <- x[, 2]
    isSignif1 <- x[, 1]
    epsvec <- rbinom(n, size=1, prob=epsval)
    isSignif2 <- 1*xor(covlabels, epsvec)

    #create the matrix..
    mat <- fillWithNoise(n, L, noiseP)
    mat <- insertSignifSubseq(mat, tau1, taulength1, isSignif1)
    mat <- insertSignifSubseq(mat, tau2, taulength2, isSignif2)

    if (showOutput){
        trueSeqStr <- paste0("[", tau1, ", ", tau1+taulength1-1, "]")
        confoundedSeqStr <- paste0("[", tau2, ", ", tau2+taulength2-1, "]")
        cat("data has TRUE siginificant subsequence:       ", trueSeqStr, "\n", sep="")
        cat("data has CONFOUNDED siginificant subsequence: ", confoundedSeqStr, "\n", sep="")
    }
    
    ##no need to sort - fine in this format
    #now do ordering
    covorder <- order(covlabels, decreasing=TRUE)
    covlabels <- covlabels[covorder]
    ylabels <- ylabels[covorder]
    mat <- mat[covorder, ]

    numCov1 <- sum(covlabels)
    numCov0 <- length(covlabels) - numCov1
    covlabelsShort <- c(numCov0, numCov1)
    
    if (showOutput){
        cat("\n")
        cat("saving sample data in folder: ", folder, sep="", "\n")
    }
    
    write.table(t(mat), file=xfilename, sep=" ", row.names=FALSE, col.names=FALSE)

    if (showOutput){
        cat("writing x (data) to: ", xfilename, "\n")
    }
    
    write.table(ylabels, file=yfilename, sep=" ", row.names=FALSE, col.names=FALSE)

    if (showOutput){
        cat("writing y (phenotype labels) to: ", yfilename, "\n")
    }

    #writing the short version...
    write.table(covlabelsShort, file=covfilename, sep=" ", row.names=FALSE, col.names=FALSE)

    if (showOutput){
        cat("writing c (covariate labels) to: ", covfilename, "\n")
    }
            

    write.table(truetau, file=truetaufilename, sep=" ", row.names=FALSE, col.names=FALSE)
#    cat("writing true tau's (but not needed for FastCMH) to: ", truetaufilename, "\n")

    if (showOutput){
        cat("writing true tau's to: ", truetaufilename, "\n")
        cat("#--------------------------------------------------------#\n")
    }

}



fixSeparator <- function(path){
    mysep <- .Platform$file.sep
    if (mysep=="/"){
        path <- gsub("\\", mysep, path, fixed=TRUE)
    }
    if (mysep=="\\"){
        path <- gsub("/", mysep, path, fixed=TRUE)
    }
    return(path)
}


#insert a significant subsequence into the data
insertSignifSubseq <- function(mat, tau, taulength, isSignif){
    taustart <- tau
    tauend <- tau+taulength-1
    for (i in 1:length(isSignif)){
        y <- isSignif[i]
        if(y==0){
            ##do nothing; make zeros
            mat[i, taustart:tauend] <- rep(0, taulength)
        } else {
            ##add significant interval
            mat[i, taustart:tauend] <- genSignifInt(taulength)
        }
    }
    return(mat)
}


#generate a significant interval of length taulength
#an interval of zeros with one 1, e.g.
#[0 0 0 1 0]
genSignifInt <- function(taulength){
    index <- sample.int(taulength, size=1)
    vec <- rep(0, taulength)
    vec[index] <- 1
    return (vec)
}


#probability of a significant interval
getCaseP <- function(p, n){
    return (  1 - (1.0-p)^(1.0/n)  )
}



#fill a data matrix with background noise
fillWithNoise <- function(rows, cols, p){
    #rows is samples, cols is number of snps
    mat <- matrix(0, ncol=cols, nrow=rows)
    for (i in 1:rows){
        mat[i, ] <- rbinom(cols, size=1, prob=p)
    }
    return(mat)
}



#check if a matrix is a covariance matrix
isCovariance <- function(mat){
    eig1 <- eigen(mat)
    val <- eig1$values
    print(val)
    allpos <- all(val > 0)
    return(allpos)
}


#make sure it is a directory (last character is "/")
#also, if folder does not exist, create it.
checkdir <- function(str, showOutput=FALSE){
    returnString <- str
    #make platform independent
    mysep <- .Platform$file.sep

    if( substr(str, nchar(str), nchar(str)) == mysep ){
        #do nothing
    } else{
        returnString <- paste0(returnString, mysep)
    }

    #now check if folder exists; if not, create it
    if (dir.exists(returnString)==FALSE){
        if (showOutput){
            cat("folder ", returnString, " does not exist. Creating folder now...\n", sep="")
        }
        dir.create(returnString, showWarnings=FALSE, recursive=TRUE)
    }
    return(returnString)
}



#makeSampleData <- function(folder="./", xfilename="data.txt", yfilename="label.txt", covfilename="cov.txt", K=2, L=1000, n=400, noiseP=0.1, rho=0.8, tau1=100, taulength1=4, tau2=200, taulength2=4, seednum=2, truetaufilename="truetau.txt"){


#make some sample data - only for a quick test of filtering function
#two significant intervals:
#[9,10] 1e-5
#[40,42] 1e-6
makeSampleFilteringData <- function(folder="./", sigIntFilename="sampleSigInts.csv"){
    #sort out filename
    folder <- checkdir(folder)
    filename = paste0(folder, sigIntFilename)

    set.seed(2)
    # columns l, tau, a, x, P.value
    # length, tau, count a, count x (a <= x), p-value double
    l <- c(3, 2, 1, 3, 3, 4, 3)
    tau <- c(10, 9, 10, 11, 40, 38, 41)
    a <- c(2, 3, 3, 2, 3, 2, 3)
    x <- c(3, 4, 4, 4, 3, 4, 4)
    P.values <- c(1e-4, 1e-5, 1e-3, 1e-3, 1e-6, 1e-3, 1e-4)
    #df <- data.frame(l, tau, a, x, P.values)
    #df should not have a and x
    df <- data.frame("l"=l, "tau"=tau,  "Pvalues"=P.values)

    #rearrange rows
    n <- nrow(df)
    perm <- sample(x=n, size=n, replace=FALSE)
    df <- df[perm, ]


#    print(df)
    write.table(df, file=filename, quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)
    cat("saved sample data to ", filename, "\n")
}
