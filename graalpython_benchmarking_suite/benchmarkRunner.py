from morpheus import morpheus
import numpy as np
import logReg
import kmeans
import statistics
import argparse
import json
import argparse
import timeit
import statistics
from kmeans import NormalizedKMeans
from logReg import NormalizedLogisticRegression
from linReg import NormalizedLinearRegression
from gnmf import GaussianNMF
import polyglot
from time import time
import csv
import os
import gc

def gen_matrices_poly(num_rows_R, num_cols_S, tup_ratio, feat_ratio, mode):
    # Computing matrix dims
    num_rows_S = num_rows_R * tup_ratio
    num_cols_R = num_cols_S * feat_ratio
    num_rows_K = num_rows_S
    num_cols_K = num_rows_R

    #num of elements per matrix
    area_S = num_rows_S * num_cols_S
    area_R = num_rows_R * num_cols_R
    area_K = num_rows_K * num_cols_K

    foo = polyglot.eval(language="R", string="source('benchUtils.r'); function(x,y,z,w){ genMatricesForPy(x,y,z,w); }")
    res = foo(num_rows_R, num_cols_S, tup_ratio, feat_ratio)

    Sarg, Ksarg, Rsarg, matMatrixForeign, target, avatarArg = res
    print(avatarArg)
    print(dir(avatarArg))
    print(avatarArg.scalarAddition(1,2))
    #raise NotImplementedError
    nm = morpheus.NormalizedMatrix(S=Sarg, Ks=Ksarg, Rs=Rsarg, foreign_backend=True, avatar=avatarArg) 
    print("#################################################")
    print(nm.isMorpheus)
    print("#################################################")
    print(avatarArg.scalarAddition(1,2))
    mat = morpheus.NormalizedMatrix(mat=matMatrixForeign, avatar=avatarArg)
    print("---------")
    target = morpheus.NormalizedMatrix(mat=target, avatar=avatarArg)
    print(">>>>>>>>>")
    data = nm if mode == "trinity" else mat
    print("#################################################")
    print(data.isMorpheus)
    print("#################################################")
    n_mat, d_mat = (num_rows_S, num_cols_R + num_cols_S)
    matrices = {
        "data": data,
        "matMatrix": mat,
        "target": target,
        "nMat": n_mat,
        "dMat": d_mat,
        "avatar" : avatarArg
    }
    return matrices

def benchmark_it(action, fname):
    print("BIT")
    times = []
    with open(fname, 'a') as f:
        for i in range(25):
            gc.collect()
            timeStart = int(round(time() * 1000)) 
            action()
            timeEnd = int(round(time() * 1000))
            timeTotal = timeEnd - timeStart 
            times.append(timeTotal)
            f.write(str(timeTotal) + "\n")
            f.flush()
            os.fsync(f.fileno())
            print(i, "__", timeTotal)
    return times

def do_logistic_regression(x, max_iter, winit, gamma, target):
    m1 = NormalizedLogisticRegression() 
    return m1.fit(x,target, winit)

def get_tasks(params, n_mat, d_mat, T, target, tasks, is_monolang):

    lmm_num_rows = d_mat
    lmm_num_cols = 2
    rmm_num_rows = 2
    rmm_num_cols = n_mat
    lmm_arg = lmm_num_rows * lmm_num_cols
    rmm_arg = rmm_num_rows * rmm_num_cols
    
    log_reg_max_iter = 20
    log_reg_gamma = 0.000001

    n_S = n_mat
    center_num = 10
    end = d_mat
    
    if is_monolang:
        #log_reg_winit = np.matrix(np.random.randn(d_mat, 1))
        #gnmf_winit = np.matrix(np.random.randn(n_mat, 5))
        #gnmf_h_init = np.matrix(np.random.rand(5, d_mat))
        #k_center = (T[:center_num, : end]).T 
        raise NotImplementedError

    else:
        print("ELSE")

        foo = polyglot.eval(language="R", string="source('benchUtils.r'); function(x,y){ genRanMatrixForPy(x,y); }")
        log_reg_winit = foo(d_mat, 1)

        #gnmf_winit = morpheus.Faux(urg(n_mat, 5))
        #gnmf_h_init = morpheus.Faux(urg(5, d_mat))
        #k_center = morpheus.Faux(T.obj.splice(0, 10-1, 0, d_mat-1).transpose())
    
    """
    all_tasks = [
        ("scalarAddition", do_scalar_addition),
        ("scalarMultiplication", do_scalar_multiplication),
        ("leftMatrixMultiplication", lambda x: do_left_matrix_multiplication(x, lmm_arg)),
        ("rightMatrixMultiplication", lambda x: do_right_matrix_multiplication(x, rmm_arg)),
        ("rowWiseSum", do_row_wise_sum),
        ("columnWiseSum", do_column_wise_sum),
        ("elementWiseSum", do_element_wise_sum),
        ("logisticRegression", 
           lambda x: do_logistic_regression(x, 
               log_reg_max_iter, log_reg_winit, log_reg_gamma, target)),
        ("linearRegression",
           lambda x: do_linear_regression(x,
               log_reg_max_iter, log_reg_winit, log_reg_gamma, target)),
        ("gaussianNMF", 
            lambda x: do_gnmf(x,
                gnmf_winit, gnmf_h_init)),
        ("kMeansClustering", 
            lambda x: do_kmeans_clustering(x, 
               log_reg_max_iter, center_num, k_center, n_S))
    ]
    """
    print("ALL TASKS")
    all_tasks = [
        ("logisticRegression", 
           lambda x: do_logistic_regression(x, 
               log_reg_max_iter, log_reg_winit, log_reg_gamma, target))
    ]

    chosen_tasks = []
    for name, func in all_tasks:
        if name in tasks:
             chosen_tasks.append((name, func))

    print(len(list(chosen_tasks)))
    return chosen_tasks

def main():

    # parse params
    parser = argparse.ArgumentParser(description='Add some integers.')
    parser.add_argument('--fpath', metavar='fpath', type=str, help='benchparams')
    parser.add_argument('--task', metavar='task', type=str, help='task to run')
    parser.add_argument('--numWarmups', metavar='numWarmups', type=int, help='number of warmups')
    parser.add_argument('--outputDir', metavar='outputDir', type=str, help='what is the output directory')
    parser.add_argument('--mode', metavar='mode', type=str, help='mode of execution')
    parser.add_argument('--monolang', metavar='monolang', type=bool, help='monolang')
    parser.add_argument('--TR', metavar='TR', type=int, help='monolang')
    parser.add_argument('--FR', metavar='FR', type=int, help='monolang')
    args = parser.parse_args()

    action_param = args.task
    dataset_meta = None
    with open(args.fpath) as f:
        dataset_meta = json.load(f)

    algorithm_tasks = ["logisticRegression", "kMeansClustering"]
    microbench_tasks = [
        "scalarAddition", "scalarMultiplication",
        "leftMatrixMultiplication", "rightMatrixMultiplication",
        "rowWiseSum", "columnWiseSum", "elementWiseSum"
    ]
    mode = args.mode
    is_monolang = False

    if not(os.path.exists(args.outputDir)):
        os.makedirs(args.outputDir)
    
    # Task selection
    tasks = []
    if action_param == "all":
        tasks = microbench_tasks + algorithm_tasks
    elif action_param == "micro":
        tasks = microbench_tasks
    elif action_param == "algorithm":
        tasks = algorithm_tasks
    else:
        tasks = [action_param]
    print(tasks)
    
    # Benchmarking loop
    TRs = [1]
    FRs = [1]

    is_data_synthetic = dataset_meta["name"] == "synthesized"
    if is_data_synthetic:
        TRs = dataset_meta["TRs"]
        FRs = dataset_meta["FRs"]

    TRs = [args.TR]
    FRs = [args.FR]
    for TR in TRs:
        for FR in FRs:
            if(is_data_synthetic and is_monolang):
                # matrices = gen_matrices(mode, dataset_meta["nR"], dataset_meta["dS"], TR, FR)
                raise NotImplementException
            elif is_data_synthetic:
                print("MAT")
                matrices = gen_matrices_poly(dataset_meta["nR"], dataset_meta["dS"], TR, FR, mode)
                print("MATS")
            else:
               raise NotImplementedError

            task_pairs = get_tasks(args, matrices["nMat"], matrices["dMat"], matrices["matMatrix"], matrices["target"], tasks, is_monolang)
            print(task_pairs)
            start_time = time()
            print("__________")
            for name, action in task_pairs:
                print("TASK")
                fname = args.outputDir + "/"  + name + "_" + dataset_meta["outputMeta"] + "_" + "TR=" + str(TR) +  "_FR=" + str(FR) + "_"  + mode + ".txt"
                times = benchmark_it(lambda : action(matrices["data"]), fname)
main()


def gen_matrices(mode, num_rows_R, num_cols_S, tup_ratio, feat_ratio):

    # Computing matrix dims
    num_rows_S = num_rows_R * tup_ratio
    num_cols_R = num_cols_S * feat_ratio
    num_rows_K = num_rows_S
    num_cols_K = num_rows_R

    #num of elements per matrix
    area_S = num_rows_S * num_cols_S
    area_R = num_rows_R * num_cols_R
    area_K = num_rows_K * num_cols_K

    #TODO: need sparse K!
    S = np.matrix(np.random.rand(num_rows_S, num_cols_S))
    K = np.matrix(np.random.rand(num_rows_K, num_cols_K))
    R = np.matrix(np.random.rand(num_rows_R, num_cols_R))


    # generate materialized matrix
    KR = np.matmul(K, R)
    materialized_matrix = np.hstack((S, KR))

    data = materialized_matrix
    if mode == "trinity":
        # generate materialized matrix
        data = morpheus.NormalizedMatrix(S, [K], [R])

    # get Y, package them, return
    n_mat, d_mat = materialized_matrix.shape
    Y = np.matrix(np.random.rand(n_mat, 1) > 0.5) 

    matrices = {
        "data": data,
        "matMatrix": materialized_matrix,
        "target": target,
        "nMat": n_mat,
        "dMat": d_mat
    }
     
    return matrices

def do_scalar_addition(x):
    return x + 42

def do_scalar_multiplication(x):
    return x * 42

def do_left_matrix_multiplication(x, lmm_arg):
    return lmm_arg * x

def do_right_matrix_multiplication(x, rmm_arg):
    return x * rmm_arg

def do_row_wise_sum(x):
    return np.sum(x, axis=0)

def do_column_wise_sum(x):
    return np.sum(x, axis=1)

def do_element_wise_sum(x):
    return np.sum(x)


def do_linear_regression(x, max_iter, winit, gamma, target):
    m1 = NormalizedLinearRegression()
    return m1.fit(x, target, winit)

def do_gnmf(x, w_init, h_init):
    m1 = GaussianNMF()
    return m1.fit(x, w_init, h_init)

def do_kmeans_clustering(x, max_iter, center_number, k_center, n_S):
    m1 = NormalizedKMeans()
    return m1.fit(x, k_center, n_S)

