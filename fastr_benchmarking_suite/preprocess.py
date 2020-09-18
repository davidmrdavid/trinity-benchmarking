from os import listdir
from json import load, dump
from os.path import isfile, join

import numpy as np
from scipy.io import mmread
from scipy.sparse import hstack, coo_matrix


outputDir = "./datasetsJSON/"
paramsDir = "./benchparams"
paramFiles = [join(paramsDir, f) for f in listdir(paramsDir) if isfile(join(paramsDir, f))] 

json_template = {"mathjs": "SparseMatrix", "values": [], "index": [], "ptr": [], "size": []}
json_template2 = {"mathjs": "DenseMatrix", "data": [], "size": []}

for paramFile in paramFiles:
  if(".swp" in paramFile):
      continue

  params = None
  with open(paramFile, "r") as f:
      params = load(f)
      for key in params.keys():
          params[key] = str(params[key])

  if not("SDir" in params.keys()):
     continue

  s = None
  if(params["SDir"] != ""):
      if params["outputMeta"] == "Walmart":


          print("A")
          s = np.matrix(np.genfromtxt(params["SDir"], skip_header=True, dtype=int)).T
          json_tosave = json_template2
          json_tosave["data"] = s.T.tolist()[0]
          json_tosave["size"] = list(s.shape) 
          fname = outputDir  + params["datasetDir"] + "/"+ params["SDir"].split("/")[-1]
          with open(fname, 'w') as outfile:
            dump(json_tosave, outfile)

      else:
          s = mmread(params['SDir'])# if params['SDir'] != "" else None

          print("B")
          json_tosave = json_template
          json_tosave["values"] = s.tocsr().data.tolist()
          json_tosave["index"] = s.tocsr().indices.tolist()
          json_tosave["ptr"] = s.tocsr().indptr.tolist()
          json_tosave["size"] = list(s.shape)
          fname = outputDir  + params["datasetDir"] + "/"+ params["SDir"].split("/")[-1]
          with open(fname, 'w') as outfile:
            dump(json_tosave, outfile)

  print("!!!")
  join_set1 = np.genfromtxt(params['FK1dir'], skip_header=True, dtype=int) -1 
  print("JSL", len(join_set1))
  row = range(0, len(join_set1))
  col = join_set1
  nrow = len(row)
  ncol = max(join_set1) +1
  K1 = coo_matrix(([1]*nrow, (row, col)), shape=(nrow, ncol))

  join_set2 = np.genfromtxt(params['FK2dir'], skip_header=True, dtype=int) -1
  row = range(0, len(join_set2))
  col = join_set2
  nrow = len(row)
  ncol = max(join_set2) + 1
  K2 = coo_matrix(([1]*nrow, (row, col)), shape=(nrow, ncol))
  
  r1 = mmread(params["R1S"])
  r2 = mmread(params['R2S'],)
  Y = np.matrix(np.genfromtxt(params['Ydir'], skip_header=True, dtype=int)).T
  json_tosave = json_template2
  json_tosave["data"] = Y.T.tolist()[0]
  json_tosave["size"] = list(Y.shape) 
  fname = outputDir  + params["datasetDir"] + "/"+ params["Ydir"].split("/")[-1]
  with open(fname, 'w') as outfile:
    dump(json_tosave, outfile)

  if params["outputMeta"] == "Flights":
    join_set3 = np.genfromtxt(params['FK3dir'], skip_header=True, dtype=int) -1
    row = range(0, len(join_set3))
    col = join_set3
    nrow = len(row)
    ncol = max(join_set3) + 1
    K3 = coo_matrix(([1]*nrow, (row, col)), shape=(nrow, ncol))

    r3 = mmread(params['R3S'],)

    print("C")
    json_tosave = json_template
    json_tosave["values"] = r3.tocsr().data.tolist()
    json_tosave["index"] = r3.tocsr().indices.tolist()
    json_tosave["ptr"] = r3.tocsr().indptr.tolist()
    json_tosave["size"] = list(r3.shape)
    fname = outputDir  + params["datasetDir"] + "/"+ params["R3S"].split("/")[-1]
    with open(fname, 'w') as outfile:
      dump(json_tosave, outfile)


    print("D")
    json_tosave = json_template
    json_tosave["values"] = K3.tocsr().data.tolist()
    json_tosave["index"] = K3.tocsr().indices.tolist()
    json_tosave["ptr"] = K3.tocsr().indptr.tolist()
    json_tosave["size"] = list(K3.shape)
    fname = outputDir  + params["datasetDir"] + "/"+ params["FK3dir"].split("/")[-1]
    with open(fname, 'w') as outfile:
      dump(json_tosave, outfile)

  json_tosave = json_template
  json_tosave["values"] = r1.tocsr().data.tolist()
  json_tosave["index"] = r1.tocsr().indices.tolist()
  json_tosave["ptr"] = r1.tocsr().indptr.tolist()
  json_tosave["size"] = list(r1.shape)
  fname = outputDir + params["datasetDir"] + "/" + params["R1S"].split("/")[-1]
  print(fname)
  with open(fname, 'w') as outfile:
    dump(json_tosave, outfile)
  print(json_tosave["size"], "R1", r1.shape, json_tosave["size"])

  json_tosave = json_template
  json_tosave["values"] = r2.tocsr().data.tolist()
  json_tosave["index"] = r2.tocsr().indices.tolist()
  json_tosave["ptr"] = r2.tocsr().indptr.tolist()
  json_tosave["size"] = list(r2.shape)
  fname = outputDir + params["datasetDir"] + "/"+ params["R2S"].split("/")[-1]
  with open(fname, 'w') as outfile:
    dump(json_tosave, outfile)

  json_tosave = json_template
  json_tosave["values"] = K1.tocsr().data.tolist()
  json_tosave["index"] = K1.tocsr().indices.tolist()
  json_tosave["ptr"] = K1.tocsr().indptr.tolist()
  json_tosave["size"] = list(K1.shape)
  fname = outputDir + params["datasetDir"] + "/" + params["FK1dir"].split("/")[-1]
  print(fname)
  with open(fname, 'w') as outfile:
    dump(json_tosave, outfile)
  print(json_tosave["size"], "K1", K1.shape, json_tosave["size"])

  json_tosave = json_template
  json_tosave["values"] = K2.tocsr().data.tolist()
  json_tosave["index"] = K2.tocsr().indices.tolist()
  json_tosave["ptr"] = K2.tocsr().indptr.tolist()
  json_tosave["size"] = list(K2.shape)
  fname = outputDir + params["datasetDir"] + "/" + params["FK2dir"].split("/")[-1]
  with open(fname, 'w') as outfile:
    dump(json_tosave, outfile)

  print("ROUND: ", params["outputMeta"], params["datasetDir"])


