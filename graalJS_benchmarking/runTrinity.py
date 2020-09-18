import subprocess
import sys

CMD = "mx --dy /compiler,js/graal-nodejs --cp-sfx ../../mxbuild/dists/jdk1.8/morpheusdsl.jar  --jdk jvmci node --polyglot morpheus.js ./benchparams/synthesized_logRegJS.json logisticRegression 5000 paper_js_algos22 trinity T %s %s" 

print("RUNNING: ", "TR=",sys.argv[1], "FR=",sys.argv[2])
COMMAND = CMD % (sys.argv[1], sys.argv[2])
print(COMMAND)
pipe = subprocess.Popen(COMMAND, shell=True)
pipe.wait()
