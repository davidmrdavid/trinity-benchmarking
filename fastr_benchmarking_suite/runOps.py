import sys
import subprocess

CMD = "mx --dynamicimports /compiler,fastr,/tools --cp-sfx ../../mxbuild/dists/jdk1.8/morpheusdsl.jar --jdk jvmci R  --J @'-Xmx220G -da -dsa -agentlib:hprof=cpu=samples' --R.PrintErrorStacktracesToFile=true  --polyglot -f benchmarkRunner.r --args -fpath %s -task %s -outputDir r_raw -mode %s -TR %d -FR %d"

# You need to remove hprof stuff above ^


D = ["movie_metadata.json"]#, 
#D = ["yelp_metadata.json", "walmart_metadata.json",
#     "books_metadata.json", "lastfm_metadata.json", "expedia_metadata.json",
#     "flights_metadata.json"] 
#M = ["linearRegression","logisticRegression"]#, "kMeansClustering"]
#D = ["synthesized.json"]
M = ["crossProduct"]
for d in D:

  d = "./benchparams/" + d
  for m in M:
    for TR in range(3,4):#15,16):
      for FR in range(5,6):
        #print("RUNNING: ", "TR=",TR, "FR=",FR)
        #COMMAND = CMD % ("rowWiseSum", "materialized", TR, FR)
        #pipe = subprocess.Popen(COMMAND, shell=True)
        #pipe.wait()

        #print("RUNNING: ", "TR=",TR, "FR=",FR)
        #COMMAND = CMD % ("linearRegression", "trinity", TR, FR)
        #pipe = subprocess.Popen(COMMAND, shell=True)
        #pipe.wait()

        print("RUNNING: ", "TR=",TR, "FR=",FR)
        COMMAND = CMD % (d, m, "trinity", TR, FR)
        pipe = subprocess.Popen(COMMAND, shell=True)
        pipe.wait()

        print("RUNNING: ", "TR=",TR, "FR=",FR)
        COMMAND = CMD % (d, m, "morpheusR", TR, FR)
        pipe = subprocess.Popen(COMMAND, shell=True)
        pipe.wait()

        print("RUNNING: ", "TR=",TR, "FR=",FR)
        COMMAND = CMD % (d, m, "materialized", TR, FR)
        pipe = subprocess.Popen(COMMAND, shell=True)
        pipe.wait()
