#' @useDynLib fastcmh
#' @importFrom Rcpp sourceCpp
#' @import bindata
#' @importFrom stats rbinom
#' @importFrom utils str write.table


fastcmhDoc <- function(){
    cat("For FastCMH documentation, type \"?FastCMH\" \n")
}
#the above function exists to avoid roxygen2 messing up NAMESPACE


#' Run the FastCMH algorithm 
#'
#' This function runs the FastCMH algorithm on a particular data set. 
#'
#' @param folder The folder in which the data is saved. If the any of 'data', 'label' and 'pvalue' arguments are not specified, then filenames must have following a naming convention inside the folder: data file is "data.txt" (i.e. full path is "folder/data.txt"), phenotype label file is "label.txt", and covariate label file is "cov.txt". More details on the structure of these files is given below, or the user can use the 'makeSampleData' function to see an example of the correct data formats. If folder="/data/", the data in fastcmh/inst/extdata is used.
#'
#' @param data The filename for the data file. Default is NULL. The data file must be an L x n txt file containing only 0s and 1s, which are space-separated in each row, while each row is on a separate newline.
#'
#' @param label The filename for the phenotype label file. Default is NULL. The label file should consist of a single column (i.e. each row is on a separate line) of 0s and 1s.
#'
#' @param cov The filename for the covariate label file. Default is NULL. The cov filecontains a single column of positive integers. The first row, containing value n_1, specifies that the first n_1 columns have covariate value 1; the second row, containing n_2, specifies that the next n_2 rows have covariate value 2, etc.
#'
#' @param alpha The value of the FWER; must be a number between 0 and 1. Default is 0.05.
#'
#' @param Lmax The maximum length of significant intervals which is considered. Must be a non-negative integer. For example, Lmax=10 searches for significant intervals up to length 10. Setting Lmax=0 will search for significant intervals up to any length (with algorithm pruning appropriately). Default is Lmax=0.
#'
#' @param showProcessing A flag which will turn printing to screen on/off. Default is FALSE (which is "off")
#'
#' @param saveAllPvals A flag which controls whether or not all the intervals (less than minimum attainable pvalue) will be returned. Default is FALSE (which is "no, do notreturn all intervals")
#'
#' @param doFDR A flag which controls whether or not Gilbert's Tarone FDR procedure (while accounting for positive regression dependence) is performed. Default is FALSE (which is "no, do not do FDR")
#'
#' @param useDependenceFDR A flag which controls whether or not Gilbert's Tarone FDR procedure uses the dependent formulation by Benjamini and Yekutieli (2001), which further adjusts alpha by dividing by the harmonic mean. This flag is only used if doFDR==TRUE. Default is FALSE.
#'
#' @param saveToFile A flag which controls whether or not the results are saved to file. By default, saveToFile=FALSE, and the data frame is returned in R. See the examples below.
#'
#' @param saveFilename A string which gives the filename to which the output is saved (needs to have saveToFile=TRUE) as an RData file. Default is "fastcmhresults.RData".
#'
#' @param saveFolder A string which gives the path to which the output will be saved (needs to have saveToFile=TRUE). Default is "./".
#'
#' 
#' @section Details:
#' This function runs the FastCMH algorithm on a particular data set in order to discover intervals that are statistically significantly associated with a particular label, while accounting for categorical covariates. 
#' The user must either supply the folder, which contains files named "data.txt", "label.txt" and "cov.txt", or the non-default filenames must be specified individually. See the descriptions of arguments "data", "label" and "cov" to see the format of the input files, or make a small sample data file using the makeSampleData() function.
#' By default, filtered results are provided. The user also has the option of using an FDR procedure rather than the standard FWER-preserving procedure. 
#'
#' 
#' @section Value:
#' FastCMH will return a list if saveToFile=FALSE (default setting), otherwise it will save the list in an .RData file. The fields of the list are:
#' 
#' sig: a dataframe listing the significant intervals, after filterting. Columns 'start', 'end' and 'pvalue' indicate the start and end points of the interval (inclusive), and the p-value for that interval
#' 
#' unfiltered: a dataframe listing all the significant intervals before filtering. The filtering compares the overlapping intervals and returns the interval with the smallest p-value in each cluster of overlapping intervals. Dataframe has has structure as sig
#'
#' fdr: (if doFDR==TRUE) significant intervals using Gilbert's FDR-Tarone procedure, after filtering. Dataframe has same structure as sig
#'
#' unfilteredFdr: (if doFDR==TRUE) a dataframe listing all the significant intervals before filtering. See description of 'unfiltered'
#' @section Author(s):
#' 
#' allTestablle: (if saveAllPvals==TRUE) a dataframe listing all the testable intervals, many of which will not be significant. Dataframe has same structure as sig
#' 
#' histObs: Together with histFreq gives a histogram of maximum attainable CMH statistics
#' 
#' histFreq: Histogram of maximum attainable CMH statistics (only reliable in the testable range)
#' 
#' summary: a character string summarising the results. Use 'cat(...$summary)' to print the results with the correct indentation/new lines.
#' 
#' timing: a list containing (i) details, a character string summarising the runtime values for the experiment - use cat(...$details) for correct indentation, etc. (ii) exec, the total execution time. (iii) init, the time to initialise the objects. (iv) fileIO, the time to read the input files. (v) compSigThresh, the time to compute the significance threshold. (vi) compSigInt, the time to compute the significant intervals. (vii) peakMemUsage, a value indicating the total amount of memory used. However, in this version, the value will always be 0.
#' 
#' 
#' 
#' Felipe Llinares Lopez, Dean Bodenham
#' 
#'
#' @section See Also:
#' 'makeSampleData'
#'
#'
#' @section References:
#'
#' Gilbert, P. B. (2005) A modified false discovery rate multipl-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(1), 143-158.
#'
#' Benjamini, Y., Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of statistics, 29(4), 1165-1188.
#'
#'
#' @examples
#' #Example with default naming convention used for data, label and cov files
#' # Note: using "/data/" as the argument for folder 
#' #       accesses the data/ directory in the fastcmh package folder
#' mylist <- FastCMH("/data/") 
#'
#' #Example where the progress will be shown
#' mylist <- FastCMH(folder="/data/", showProcessing=TRUE) 
#'
#' #Example where many parameters are specified
#' mylist <- FastCMH(folder="/data/", data="data2.txt", alpha=0.01, Lmax=7)
#'
#' #Example where Gilbert's Tarone-FDR procedure is used
#' mylist <- FastCMH("/data/", doFDR=TRUE) 
#'
#' #Example where FDR procedure takes some dependence structures into account
#' mylist <- FastCMH("/data/", doFDR=TRUE, useDependenceFDR=TRUE) 
#'
#' #Example where the data frame is saved to file
#' FastCMH("/data/", saveToFile=TRUE, saveFolder="./", saveFilename="output.RData") 
#'
#'
#' @export
FastCMH <- function(folder=NULL, data=NULL, label=NULL, cov=NULL, alpha=0.05, Lmax=0, showProcessing=FALSE, saveAllPvals=FALSE, doFDR=FALSE, useDependenceFDR=FALSE, saveToFile=FALSE, saveFilename="fastcmhresults.RData", saveFolder=NULL){

    #check if any of the file names are NULL; not specified
    if ( (is.null(data)) || (is.null(label)) || (is.null(cov)) )
        notAllFileNamesProvided <- TRUE
    else
        notAllFileNamesProvided <- FALSE

    #if folder not specified AND not all filenames are provided, stop and throw error.
    if ( (is.null(folder)) && (notAllFileNamesProvided)){
    filenameErrorString <- paste0("Folder is not provided and not all files specified. FastCMH aborted. \nPlease either specify the folder, and then use 'data.txt', 'label.txt', and 'cov.txt' as the filenames for the data, phenotype labels and covariate labels, respectively, or specify all filenames, including paths.\n")

      stop(filenameErrorString, call. = FALSE)
    }


    #get default filenames, in case some are missing
    if (!is.null(folder)){

        #special case:
        #now make it access data folder in R package
        if (folder=="/data/"){
            #for some reason, automatically goes to fastcmh/inst/
            folder <- file.path(system.file(package="fastcmh"), "extdata")
            folder <- paste0(folder, .Platform$file.sep)
        }

        #normalise path for Windows
#        folder <- normalizePath(folder)
        folder <- fixSeparator(folder)
        #add / to end of folder, if not present
        folder <- checkdir(folder)

        #get defaults
        df <- loadDefaultFileNames(folder) 

        #if data filename is null, use default - later will check existence    
        if ( is.null(data) ) {
            data <- df$xfilename
        } else {
            data <- paste0(folder, data)
        }

        #if label filename is null, use default - later will check existence    
        if ( is.null(label) ) {
            label <- df$yfilename
        } else {
            label <- paste0(folder, label)
        }

        #if cov filename is null, use default - later will check existence    
        if ( is.null(cov) ) {
            cov <- df$cfilename
        } else {
            cov <- paste0(folder, cov)
        }
    }

    #check that data file exists
    if (!file.exists(data)){
        dataErrorString <- paste0("data file ", data, " does not exist. FastCMH aborted\n")
        stop(dataErrorString, call.=FALSE)
    } else {
        #expand path in case of tilde ~
        data <- path.expand(data)
    }

    #check that label file exists
    if (!file.exists(label)){
        labelErrorString <- paste0("label file ", label, " does not exist. FastCMH aborted\n")
        stop(labelErrorString, call. = FALSE)
    } else {
        #expand path in case of tilde ~
        label <- path.expand(label)
    }

    #check that cov file exists
    if (!file.exists(cov)){
        covErrorString <- paste0("cov file ", cov, " does not exist. FastCMH aborted\n")
        stop(covErrorString, call. = FALSE)
    } else {
        #expand path in case of tilde ~
        cov <- path.expand(cov)
    }


    #check alpha
    if (!is.numeric(alpha)){
        stop("alpha is not numeric. FastCMH aborted. Please ensure alpha is numeric and in interval(0, 1).\n")
    }
    if (  (is.numeric(alpha)) && ( (alpha <= 0) || (alpha >= 1) )  ){
        stop("alpha is not in the interval (0, 1). FastCMH aborted. Please ensure alpha is numeric and in interval(0, 1).\n")
    }

    #check Lmax is an integer
    tol <- 1e-5
    #this code checks float part - no easy function in R to check an integer 
    LMaxRemainder <- min(abs(c(Lmax%%1, Lmax%%1-1)))
    if (is.numeric(Lmax)){
        if ((LMaxRemainder > tol) || (LMaxRemainder < 0) ){
            stop("Lmax is not a positive integer. FastCMH aborted.\n")
        }
    } else { 
        stop("Lmax is not a positive integer. FastCMH aborted.\n")
    }



    #check showProcessing:
    #TODO
    #should there be error messages if the booleans are messed up?
    #seems that there is no base is.boolean function; easy way to force non-true
    #to be FALSE
    if (is.null(showProcessing)){
        showProcessing <- FALSE
    } 
    if(isTRUE(showProcessing)){
        showProcessing <- TRUE
    } else {
        showProcessing <- FALSE
    }

    #check saveAllPvals:
    if (is.null(saveAllPvals)){
        saveAllPvals <- FALSE
    } 
    if(isTRUE(saveAllPvals)){
        saveAllPvals <- TRUE
    } else {
        saveAllPvals <- FALSE
    }

    #check doFDR:
    if (is.null(doFDR)){
        doFDR <- FALSE
    } 
    if(isTRUE(doFDR)){
        doFDR <- TRUE
        #if doFDR is true, need saveAllPvals to be true...just a quick fix for now
        saveAllPvals <- TRUE
    } else {
        doFDR <- FALSE
    }

    #check useDependenceFDR:
    if (is.null(useDependenceFDR)){
        useDependenceFDR <- FALSE
    } 
    if(isTRUE(useDependenceFDR)){
        useDependenceFDR <- TRUE
    } else {
        useDependenceFDR <- FALSE
    }

    #filtered intervals filename

    #check saveToFile:
    if (is.null(saveToFile)){
        saveToFile <- FALSE
    } 
    if(isTRUE(saveToFile)){
        useDependenceFDR <- TRUE
    } else {
        useDependenceFDR <- FALSE
    }

    #check saveFilename:
    if (is.null(saveFilename)){
        saveFilename <- "fastcmhresults.RData" 
    } 

    #check saveFolder:
    if (is.null(saveFolder)){
        saveFolder <- "./" 
#        saveFolder <- normalizePath(saveFolder)
        saveFolder <- fixSeparator(saveFolder)
        saveFolder <- checkIsFolder(saveFolder)
    } 


    #now return df
    fastcmhresults <- main_fastcmh_cpp(data, label, cov, alpha, Lmax, showProcessing, saveAllPvals, doFDR, useDependenceFDR)

    #if save to file:
    if (saveToFile){
        saveFolder <- checkIsFolder(saveFolder)
        saveFolder <- checkFolderExists(saveFolder)
        if (saveFolder != "./"){
            saveFilename <- paste0(saveFolder, saveFilename)
        }
        save(fastcmhresults, file=saveFilename)
        cat("saving FastCMH results to: ", saveFilename, "\n", sep="")
    } else {
        #return results:
        return(fastcmhresults)
    }

}




#just wait for a key before continuing
readkey <- function()
{
    cat ("Press [enter] to continue")
    line <- readline()
}

#' Demo of FastCMH
#'
#' This function runs a demo for FastCMH, by first creating a sample data set and then running FastCMH on this data set.
#' @param folder The folder in which the data for the demo will be saved. Default is the current directory, "./". The demo data will created in folder/data and the results will be saved in folder/results as an RData file
#' @section Details:
#' This function will first create a sample data set in the folder/data, and will then run FastCMH on this data set, before saving the results in folder/results. The method runs in several steps, with each step showing the R code that can be used to do the step, then running that R code, and then waiting for the user to press enter before moving onto the next step.
#' @examples
#' DemoFastCMH("../fastcmhdemodata") 
#' DemoFastCMH() #default folder is ./fastcmhdemo
#' @export
DemoFastCMH <- function(folder="fastcmhdemo/"){

    #the amount of delay, when a sleep delay was used instead of readkey()
    sleepDelay <- 2

    #check saveFolder:
    if (is.null(folder)){
        folder <- "./fastcmhdemo/" 
    } 

    #normalise the path
#    folder <- normalizePath(folder)
    folder <- fixSeparator(folder)

    #check it is a folder
    folder <- checkIsFolder(folder)

    cat("\n")
    cat("====================================================\n")
    cat("FastCMH Demo:\n")
    cat("====================================================\n\n")

    #creating folder
    cat("Creating folders for demo: ", folder, "\n", sep="")
    folder <- checkFolderExists(folder)
    #create data folder and results folder
    datafolder <- paste0(folder, "data")
    datafolder <- checkFolderExists(datafolder)
    resultsfolder <- paste0(folder, "results")
    resultsfolder <- checkFolderExists(resultsfolder)
    #add a delay
    cat("\n")
#    Sys.sleep(sleepDelay)
    readkey();


    cat("----------------------------------------------------\n\n")
    cat("Step 1: Creating data files for demo \n", sep="")
    cat("**(Run command)**\n")
    cat("> makeSampleData(\"", datafolder, "\", showOutput=TRUE)\n\n", sep="")
    readkey();

    makeSampleData(datafolder, showOutput=TRUE)
    cat("\n\n")
    cat("Sample data created in ", datafolder, "\n\n", sep="")
    readkey();


    cat("\n\n----------------------------------------------------\n")
    cat("\n")
    cat("Step 2: Now running FastCMH, and saving results to file\n")
    cat("**(Run command)**\n")
    cat("> FastCMH(folder=\"", datafolder, ", saveToFile=TRUE, saveFolder=\"", resultsfolder, "\")\n\n", sep="")
    readkey();


    #Do FastCMH
    FastCMH(folder=datafolder, saveToFile=TRUE, saveFolder=resultsfolder)


    resultsfile <- paste0(resultsfolder, "fastcmhresults.RData")
    cat("\n\n----------------------------------------------------\n")
    cat("\n")
    cat("Step 3: Inspect the results\n")
    cat("**(Run commands)**\n")
    cat("> load(\"", resultsfile, "\")\n", sep="")
    cat("> str(fastcmhresults)\n\n", sep="")
    fastcmhresults <- FastCMH(folder=datafolder, saveToFile=FALSE)
    print(str(fastcmhresults))
    readkey();


    cat("\n\n----------------------------------------------------\n")
    cat("Step 4: Or, access the results directly\n")
    cat("\n\n")
    cat("The data frame of significant intervals (filtered) can be found using:\n")
    cat("**(Run command)**\n")
    cat("> fastcmhresults$sig\n", sep="")
    cat("\n")
    print(fastcmhresults$sig)


    cat("\n")
    cat("In other words, ", sep="") 
    if (nrow(fastcmhresults$sig)==1){
        cat("the only ", sep="")
    } else { 
        cat("the first", sep="")
    }
    cat("significant interval is [", as.numeric(fastcmhresults$sig[1]), ", ", as.numeric(fastcmhresults$sig[2]), "],\n", sep="")
    cat("and the confounded interval [200, 203] is not found/is ignored.\n")


    cat("\n\n")
    cat("====================================================\n")
    cat("End of Demo\n")
    cat("====================================================\n\n")
}





#make sure it is a directory (last character is "/")
checkIsFolder <- function(str){
    returnString <- str
    mysep <- .Platform$file.sep
    if( substr(str, nchar(str), nchar(str)) == mysep ){
        #do nothing
    } else{
        returnString <- paste0(returnString, mysep)
    }
    return(returnString)
}


#make sure it is a directory (last character is "/")
#also, if folder does not exist, create it.
checkFolderExists <- function(str){
    returnString <- str
    mysep <- .Platform$file.sep
    if( substr(str, nchar(str), nchar(str)) == mysep ){
        #do nothing
    } else{
        returnString <- paste0(returnString, mysep)
    }

    #now check if folder exists; if not, create it
    if (dir.exists(returnString)==FALSE){
        cat("folder ", returnString, " does not exist. Creating folder now...\n", sep="")
        dir.create(returnString, showWarnings=FALSE, recursive=TRUE)
    }
    return(returnString)
}




#get a data.frame of default names, with folder prefix
loadDefaultFileNames <- function(folder){
    #make sure last character is "/", otherwise add it
    folder <- checkIsFolder(folder)

    #expand path in case of tilde ~
    folder <- path.expand(folder)

    xfilename <- paste0(folder, "data.txt")
    yfilename <- paste0(folder, "label.txt")
    cfilename <- paste0(folder, "cov.txt")
    Lmax <- 10
    alpha <- 0.05
    basefilename <- paste0(folder, "result")
    pvalfilename <- paste0(folder, "pval.txt")
    filteredfilename <- paste0(folder, "filtered.txt")
    dlist <-list(xfilename=xfilename, yfilename=yfilename, cfilename=cfilename, Lmax=Lmax, alpha=alpha, basefilename=basefilename, pvalfilename=pvalfilename, filteredfilename=filteredfilename)
    return(dlist)
}



#the main function
main_fastcmh_cpp <- 
function(xfilename, yfilename, cfilename, alpha, Lmax, showProcessing, saveAllPvals, doFDR, useDependenceFDR){
    returnList <- main_fastcmh2(xfilename, yfilename, cfilename, alpha, Lmax, showProcessing, saveAllPvals, doFDR, useDependenceFDR)

    return(returnList)
}




