library("jsonlite")
source("benchUtils.r")
source("microBenchmarks.r")
source("normalizedMatrix.r")

# Env settings
set.seed(0)

# Parameter setting
cliArgs = commandArgs(trailingOnly=TRUE)
datasetJSON <- cliArgs[2]
actionParam <- cliArgs[4]
outputDir <- cliArgs[6]
mode <- cliArgs[8]
TR <- cliArgs[10]
FR <- cliArgs[12]

print(datasetJSON)
print(actionParam);

datasetMeta <- fromJSON(txt=datasetJSON)


# Kinds of tasks
algorithmTasks <- c("logisticRegression", "kMeansClustering")
microbenchTasks <- c(
  "scalarAddition", "scalarMultiplication",
  "leftMatrixMultiplication", "rightMatrixMultiplication", 
  "rowWiseSum", "columnWiseSum", "elementWiseSum")

# Task selection
tasks = c()
if(actionParam == "all") {
    tasks = c(microbenchTasks, algorithmTasks)
} else if(actionParam == "micro") {
    tasks = microbenchTasks
} else if(actionParam == "algorithm") {
    tasks = algorithmTasks
} else if(actionParam == "test") {
    tasks = c("leftMatrixMultiplication",
        "rightMatrixMultiplication", "rowWiseSum",
        "columnWiseSum", "elementWiseSum")
} else {
    tasks = c(actionParam)
}

# Benchmarking loop
isDataSynthetic = (datasetMeta$name == "synthesized")
extended = (datasetMeta$outputMeta == "synth_extended")
TRs = c(as.numeric(TR));
FRs = c(as.numeric(FR));


# Create output directory, if it does not exist
workingDir <- getwd()
outputPath <- file.path(workingDir, outputDir)
dir.create(outputPath, outputDir)

print(TRs);
print(FRs)

#TRs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
#FRs <- c(1,2,3,4,5)

for(TR in TRs) {
  for(FR in FRs) {
    gc()

    # Get data
    if(isDataSynthetic){
      if(extended){
        matrices <- genDatasetExtended(mode, datasetMeta$nR, datasetMeta$dS, TR, FR)
      }
      else{
        matrices <- genDataset(mode, datasetMeta$nR, datasetMeta$dS, TR, FR)
      }
    }
    else{
      matrices <- loadDataset(mode, datasetMeta) 
    }

    # Get tasks as anonymous functions
    target <- matrices$target
    matMatrix <- matrices$matMatrix
    taskPairs <- getTasks(datasetMeta, matrices$nMat, matrices$dMat, matMatrix, target, tasks)

    # benchmark and record results
    for(pair in taskPairs) {

      print(paste("Execution mode:", mode, "Experiment:", pair$name, "TR-FR:", "(",TR, ",", FR, ")"))

      # Testing for runtime performance
      if(actionParam != "test") {
        print("Testing performance")
        action <- function() { pair$runner(matrices$data) }
        fname <- file.path(outputPath, sprintf("%s_%s_TR=%i_FR=%i_%s.txt", pair$name, datasetMeta$outputMeta, TR, FR, mode))
        benchmarkIt(action, 5, fname)  
      } 
      else {  # Testing for correctness
        print("Testing correctness")
        action1 <- function() { pair$runner(matrices$data)}
        action2 <- function() { pair$runner(matrices$matMatrix)}
        print("Testing non-transposed")
        result <- checkEquivalence(action2, action1)
        if(pair$name != "leftMatrixMultiplication" && pair$name != "rightMatrixMultiplication"){
          print("Testing transposed")
          action1T <- function() { pair$runner(t(matrices$data))}
          action2T <- function() { pair$runner(t(matrices$matMatrix))}
          result <- checkEquivalence(action1T, action2T)
        }
      }
    }
  }
}

