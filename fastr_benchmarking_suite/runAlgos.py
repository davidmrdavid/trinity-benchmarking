import sys
import subprocess

#CMD = "mx --dynamicimports /compiler,fastr --cp-sfx ../../mxbuild/dists/jdk1.8/morpheusdsl.jar --J @'-Xmx220G' --jdk jvmci R --polyglot -f benchmarkRunner.r --args -fpath ./benchparams/synthesized.json -task %s -outputDir r_raw2 -mode %s -TR %d -FR %d"

# THSI IS THE CURRENT ONE!
DATASETS = [
"movie_metadata.json",
"books_metadata.json",
"lastfm_metadata.json",
"walmart_metadata.json",
"expedia_metadata.json",
"movie_metadata.json",
"yelp_metadata.json",
"flights_metadata.json"
]
CMD = "mx --dynamicimports fastr,/compiler --cp-sfx ../../mxbuild/dists/jdk1.8/morpheusdsl.jar --J @'-Xmx220G' --jdk jvmci R --polyglot -f benchmarkRunner.r --args -fpath %s -task %s -outputDir paper_r_algos4 -mode %s -TR %s -FR %s"
for dataset in DATASETS:
    for task in ["GNMFClustering"]:
#         COMMAND = CMD % ("./benchparams/"+dataset, task, "trinity", "1", "1")#(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
#         print(COMMAND)
#         pipe = subprocess.Popen(COMMAND, shell=True)
#         pipe.wait()
        
#         COMMAND = CMD % ("./benchparams/"+dataset, task, "morpheusR", "1", "1")#(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
#         print(COMMAND)
#         pipe = subprocess.Popen(COMMAND, shell=True)
#         pipe.wait()
        
        COMMAND = CMD % ("./benchparams/"+dataset, task, "materialized", "1", "1")#(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        print(COMMAND)
        pipe = subprocess.Popen(COMMAND, shell=True)
        pipe.wait()
