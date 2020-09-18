import subprocess
import sys

TR = 10
FR = [1,2,3,4,5]


CMD = "mx --dy /compiler,js/graal-nodejs --cp-sfx ../../mxbuild/dists/jdk1.8/morpheusdsl.jar  --jdk jvmci node --polyglot morpheus.js ./benchparams/synthesized_logRegJS.json logisticRegression 5000 paper_js_algos materialized T %s %s" 

for fr in FR:
    print("RUNNING: ", "TR=",TR, "FR=",fr)
    COMMAND = CMD % (TR, fr)
    print(COMMAND)
    pipe = subprocess.Popen(COMMAND, shell=True)
    pipe.wait()
