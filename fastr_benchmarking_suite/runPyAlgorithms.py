import subprocess
import sys

CMD = "mx --dynamicimports /compiler,graalpython,fastr --cp-sfx ../../mxbuild/dists/jdk1.8/morpheusdsl.jar:../../../../fastr/mxbuild/dists/jdk1.8/fastr.jar --J @-Xmx220G --jdk jvmci python --polyglot ../graalpython_benchmarking_suite/benchmarkRunner.py --fpath ./benchparams/synthesized2.json --task %s --numWarmups 10000 --mode %s --monolang False --outputDir paper_py_algos --TR %s --FR %s"


print("RUNNING: ")
COMMAND = CMD % ("linearRegression", "materialized", "1", "1")
print(COMMAND)
pipe = subprocess.Popen(COMMAND, shell=True)
pipe.wait()
