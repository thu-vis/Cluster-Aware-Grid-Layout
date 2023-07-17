import numpy as np
import os
import random
from scipy.spatial.distance import cdist
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import gridlayoutOpt
import time


class gridOptimizer(object):
    def __init__(self):
        super().__init__()

    def solveKM(self, cost_matrix):  # solve bipartite graph assignment with KM algorithm
        # print('start KM')
        # row_asses, _, _ = lapjv.lapjv(cost_matrix)
        row_asses = np.array(gridlayoutOpt.solveKM(cost_matrix))
        N = row_asses.shape[0]
        col_asses = np.zeros(shape=N, dtype='int')
        for i in range(N):
            col_asses[row_asses[i]] = i
        # print('end KM')
        return row_asses, col_asses

    def solveJV(self, cost_matrix):  # solve bipartite graph assignment with JV algorithm
        N = cost_matrix.shape[0]
        row_asses = np.array(gridlayoutOpt.solveLap(cost_matrix, True, max(50, int(0.15 * N))))
        # row_asses = np.array(gridlayoutOpt.solveLap(cost_matrix, True, 100))
        # row_asses = np.array(gridlayoutOpt.solveLap(cost_matrix, True, 50))
        col_asses = np.zeros(shape=N, dtype='int')
        for i in range(N):
            col_asses[row_asses[i]] = i
        # print('end KM')
        return row_asses, col_asses

    def check_cost_type(self, ori_embedded, row_asses, labels, type):  # calculate the cost of grid layout
        # tmp_row_asses = np.array(
        #     gridlayoutOpt.optimizeSwap(ori_embedded, row_asses, labels, change, type, 1, 0,
        #                                0, 10, False, 0))
        N = row_asses.shape[0]
        tmp_row_asses = np.array(
            gridlayoutOpt.checkCostForAll(ori_embedded, row_asses, labels, type, 1, 0))

        new_cost = np.array([tmp_row_asses[N], tmp_row_asses[N + 1], tmp_row_asses[N + 2], tmp_row_asses[N + 3]])
        # print(new_cost)
        return new_cost

    # optimize gridlayout (ours)
    # ori_embedded: grids position of the initial layout
    # row_asses: current layout
    # labels: cluster labels
    # useGlobal: if execute global step
    # useLocal: if execute local step
    # convex_type: type of optimized convexity measure in local step
    # maxit: max number of rounds of bipartite graph(accelerated swap) in local step
    # maxit2: max number of rounds of swap in local step
    # only_compact: if generate a most compactness layout
    # swap_cnt: max times of swap in local step
    # swap_op_order: if execute local step first and then global step
    # choose_k: randomly select from the top k options when swap in local step
    def grid_op(self, ori_embedded, row_asses, labels, useGlobal=True, useLocal=True, convex_type="Triple", maxit=10,
                maxit2=5, only_compact=False, swap_cnt=2147483647, swap_op_order=False, choose_k=1):

        start = time.time()
        N = row_asses.shape[0]

        def solve_op(ori_embedded, row_asses, type, alpha, beta, alter=False, alter_best=None, maxit=5, maxit2=2,
                     swap_cnt=2147483647):
            if alter_best is None:
                alter_best = [1, 1, 1, 1]

            # iterative optimization
            new_row_asses = row_asses.copy()
            new_row_asses2 = row_asses.copy()
            ans_row_asses = row_asses.copy()
            new_cost = np.array([2147483647, 2147483647, 2147483647])
            best_cost = new_cost.copy()
            print("````````````````````")
            print('alpha', alpha, beta)
            print("````````````````````")
            change = np.ones(shape=N, dtype='bool')

            if maxit > 0:  # iteration times for bipartite graph assignment optimization
                tmp_start = time.time()
                # tmp_row_asses = np.array(gridlayoutOpt.optimizeBA(ori_row_asses, row_asses, labels, change, type, alpha, beta, maxit))
                tmp_row_asses = np.array(
                    gridlayoutOpt.optimizeBA(ori_embedded, row_asses, labels, change, type, alpha, beta, alter,
                                             alter_best, maxit))
                for i in range(N):
                    new_row_asses[i] = round(tmp_row_asses[i])
                new_cost = np.array([tmp_row_asses[N], tmp_row_asses[N + 1], tmp_row_asses[N + 2]])

                ans_row_asses = new_row_asses.copy()
                best_cost = new_cost.copy()
                # print("cost1", best_cost)

            # if maxit2 >= 0:
            if maxit2 > 0:  # iteration times for swap optimization
                change = np.ones(shape=N, dtype='bool')
                seed = 10

                tmp_row_asses = np.array(
                    gridlayoutOpt.optimizeSwap(ori_embedded, new_row_asses, labels, change, type, alpha, beta,
                                               maxit2, choose_k, seed, True, swap_cnt))
                for i in range(N):
                    new_row_asses2[i] = round(tmp_row_asses[i])
                new_cost = np.array([tmp_row_asses[N], tmp_row_asses[N + 1], tmp_row_asses[N + 2]])
                best_cost = new_cost.copy()
                ans_row_asses = new_row_asses2.copy()

            if (type != "Global") and (maxit2 == 0):  # assignment optimization in cluster
                change = np.ones(shape=N, dtype='bool')
                new_row_asses2 = np.array(
                    gridlayoutOpt.optimizeInnerCluster(ori_embedded, new_row_asses, labels, change))
                new_cost = np.array([-1, -1, -1])
                best_cost = new_cost.copy()
                ans_row_asses = new_row_asses2.copy()

            return ans_row_asses, best_cost

        compact_it = 3  # iteration times for the most compact layout
        global_it = 5  # max iteration times for the global optimization
        alter = True  # allow to adjust the weight parameter of proximity and compactness

        ori_row_asses = row_asses.copy()

        ans = row_asses
        new_cost = np.array([2147483647, 2147483647, 2147483647])

        if not swap_op_order and useGlobal:  # Global Step

            alter_best = [0, 0, 1, 1]

            if alter:
                # row_asses1, new_cost1 = solve_op(ori_embedded, row_asses, "Global", 0, 0, False, None, 0, 0)
                row_asses1 = ori_row_asses.copy()
                new_cost1 = self.check_cost_type(ori_embedded, row_asses1, labels, "Global")
                # self.show_grid(row_asses1, labels, square_len, "1.png", False, False)
                print("new_cost1", new_cost1)
                print("time1", time.time() - start)
                row_asses2, new_cost2 = solve_op(ori_embedded, ans, "Global", 0, 1, False, None, compact_it, 0)
                # self.show_grid(row_asses2, labels, square_len, "2.png", False, False)
                print("new_cost2", new_cost2)
                print("time2", time.time() - start)

                alter_best[0] = min(new_cost1[0], new_cost2[0])
                alter_best[2] = max(new_cost1[0], new_cost2[0]) - min(new_cost1[0], new_cost2[0])
                alter_best[1] = min(new_cost1[1], new_cost2[1])
                alter_best[3] = max(new_cost1[1], new_cost2[1]) - min(new_cost1[1], new_cost2[1])

                print(alter_best)

                if only_compact:
                    ans, new_cost = row_asses2, new_cost2
                elif (alter_best[3] == 0) or (alter_best[2] == 0):
                    ans, new_cost = row_asses1, new_cost1
                else:
                    # ans, new_cost = solve_op(ori_embedded, ans, "Global", 0, 0.5, True, alter_best, global_it, 0)
                    ans, new_cost = solve_op(ori_embedded, ans, "Global", 0, 0.5, True, alter_best, global_it, 0)
                    # self.show_grid(ans, labels, square_len, "4.png", False, False)
                print("new_cost4", new_cost)
                print("time4", time.time() - start)
            else:
                # ans, new_cost = solve_op(ori_embedded, ans, "Global", 0, 0.5, True, alter_best, global_it, 0)
                ans, new_cost = solve_op(ori_embedded, ans, "Global", 0, 0.5, True, alter_best, 1, 0)
                # self.show_grid(ans, labels, square_len, "4.png", False, False)

        end2 = time.time()
        t2 = end2 - start

        if useLocal:  # Local Step

            ans, new_cost = solve_op(ori_embedded, ans, convex_type, 1, 0, False, None, maxit, maxit2, swap_cnt)

            print("new_cost5", new_cost)
            print("time5", time.time() - start)

        if swap_op_order and useGlobal:  # Global Step after local Step

            alter_best = [0, 0, 1, 1]

            if alter:
                # row_asses1, new_cost1 = solve_op(ori_embedded, row_asses, "Global", 0, 0, False, None, 0, 0)
                row_asses1 = ori_row_asses.copy()
                new_cost1 = self.check_cost_type(ori_embedded, row_asses1, labels, "Global")
                # self.show_grid(row_asses1, labels, square_len, "1.png", False, False)
                print("new_cost1", new_cost1)
                print("time1", time.time() - start)
                row_asses2, new_cost2 = solve_op(ori_embedded, ans, "Global", 0, 1, False, None, compact_it, 0)
                # self.show_grid(row_asses2, labels, square_len, "2.png", False, False)
                print("new_cost2", new_cost2)
                print("time2", time.time() - start)

                alter_best[0] = min(new_cost1[0], new_cost2[0])
                alter_best[2] = max(new_cost1[0], new_cost2[0]) - min(new_cost1[0], new_cost2[0])
                alter_best[1] = min(new_cost1[1], new_cost2[1])
                alter_best[3] = max(new_cost1[1], new_cost2[1]) - min(new_cost1[1], new_cost2[1])

                print(alter_best)

                if only_compact:
                    ans, new_cost = row_asses2, new_cost2
                elif (alter_best[3] == 0) or (alter_best[2] == 0):
                    ans, new_cost = row_asses1, new_cost1
                else:
                    # ans, new_cost = solve_op(ori_embedded, ans, "Global", 0, 0.5, True, alter_best, global_it, 0)
                    ans, new_cost = solve_op(ori_embedded, ans, "Global", 0, 0.5, True, alter_best, global_it, 0)
                    # self.show_grid(ans, labels, square_len, "4.png", False, False)
                print("new_cost4", new_cost)
                print("time4", time.time() - start)
            else:
                # ans, new_cost = solve_op(ori_embedded, ans, "Global", 0, 0.5, True, alter_best, global_it, 0)
                ans, new_cost = solve_op(ori_embedded, ans, "Global", 0, 0.5, True, alter_best, 1, 0)
                # self.show_grid(ans, labels, square_len, "4.png", False, False)

        end = time.time()
        t1 = end - start

        # return ans, t1, t2, new_cost, new_cost2[3]
        return ans, t1, t2

    # baseline + ours
    # X_embedded: t-sne position of elements
    # labels: ground truth labels
    # type : type of optimized convexity measure in local step
    # maxit: max number of rounds of bipartite graph(accelerated swap) in local step
    # maxit2: max number of rounds of swap in local step
    # use_global: if execute global step
    # use_local: if execute local step
    # only_compact: if generate a most compactness layout
    # swap_cnt: max times of swap in local step
    # pred_labels: predict labels
    # swap_op_order: if execute local step first and then global step
    # choose_k: randomly select from the top k options when swap in local step
    def grid(self, X_embedded: np.ndarray, labels: np.ndarray = None, type='Triple', maxit=10, maxit2=5, use_global=True,
             use_local=True, only_compact=False, swap_cnt=2147483647, pred_labels=None, swap_op_order=False,
             choose_k=1):
        if pred_labels is None:
            pred_labels = labels.copy()

        # Initialize
        start = time.time()  # Record start time
        ans = None
        best = 2147483647
        X_embedded -= X_embedded.min(axis=0)  # Normalize data
        X_embedded /= X_embedded.max(axis=0)
        num = X_embedded.shape[0]  # Number of data points
        square_len = math.ceil(np.sqrt(num))  # Grid length
        print('num', num)
        N = square_len * square_len  # Grid size
        maxLabel = 0
        for id in range(num):
            label = labels[id]
            maxLabel = max(maxLabel, label + 1)

        def getDist(x1, y1, x2, y2):
            return np.sqrt(np.square(x1 - x2) + np.square(y1 - y2))  # Calculate distance

        # Baseline layout
        grids = np.dstack(np.meshgrid(np.linspace(0, 1 - 1.0 / square_len, square_len),
                                      np.linspace(0, 1 - 1.0 / square_len, square_len))) \
            .reshape(-1, 2)  # Create grids

        tmp = grids[:, 0].copy()
        grids[:, 0] = grids[:, 1]
        grids[:, 1] = tmp

        # print(grids)

        original_cost_matrix = cdist(grids, X_embedded, "euclidean")  # Compute initial cost matrix
        # knn process
        dummy_points = np.ones((N - original_cost_matrix.shape[1], 2)) * 0.5  # Create dummy data points
        # dummy at [0.5, 0.5]
        dummy_vertices = (1 - cdist(grids, dummy_points, "euclidean")) * 100  # Calculate distances between dummy data points and grids
        cost_matrix = np.concatenate((original_cost_matrix, dummy_vertices), axis=1)  # Merge cost matrix

        cost_matrix = np.power(cost_matrix, 2)  # Square the cost matrix

        # row_asses, col_asses = self.solveKM(cost_matrix)
        row_asses, col_asses = self.solveJV(cost_matrix)  # Solve the assignment problem of the cost matrix
        # col_asses = col_asses[:num]
        # self.show_grid(row_asses, labels, square_len, 'new.png')

        old_X_embedded = X_embedded
        ori_row_asses = row_asses.copy()
        ori_embedded = grids[col_asses]
        ori_labels = labels.copy()

        t0 = time.time() - start

        # Cluster labels
        print("start cluster")
        # labels = np.array(gridlayoutOpt.getClusters2(ori_row_asses, labels, pred_labels))
        labels = pred_labels.copy()
        maxLabel = labels.max() + 1
        print("end cluster")


        # datas = np.load("T-base.npz")
        # labels = datas['labels']
        # row_asses = datas['row_asses']
        # ori_row_asses = row_asses

        # data = np.load('T-global.npz')
        # labels = data['labels']

        # Start optimization
        print("start optimize")
        print("--------------------------------------------------")
        start = time.time()

        ans, t1, t2 = self.grid_op(ori_embedded, row_asses, labels,
                                                 use_global, use_local, type,
                                                 maxit, maxit2, only_compact,
                                                 swap_cnt, swap_op_order,
                                                 choose_k=choose_k)

        end = time.time()
        print("end optimize")
        print("--------------------------------------------------")
        print('time:', end - start)

        if (maxit + maxit2) == 0:
            t1 = t2

        return ans, t1, t0, labels, ori_row_asses

    def translateAdjust(self, ori_embedded, row_asses, labels, useGlobal=True, useLocal=True, convex_type="Triple", maxit=5,
                        maxit2=2, change_list=np.array([]), translate=np.array([0, 0])):
        N = row_asses.shape[0]
        square_len = round(np.sqrt(N))
        num = labels.shape[0]

        # X_embedded = np.zeros((N, 2))
        # for x in range(square_len):
        #     bias = x*square_len
        #     for y in range(square_len):
        #         gid = bias+y
        #         X_embedded[ori_row_asses[gid]][0] = x*1.0/square_len
        #         X_embedded[ori_row_asses[gid]][1] = y*1.0/square_len
        X_embedded = ori_embedded.copy()

        X_embedded[row_asses[change_list]] += translate

        # 生成初始layout
        grids = np.dstack(np.meshgrid(np.linspace(0, 1 - 1.0 / square_len, square_len),
                                      np.linspace(0, 1 - 1.0 / square_len, square_len))) \
            .reshape(-1, 2)

        tmp = grids[:, 0].copy()
        grids[:, 0] = grids[:, 1]
        grids[:, 1] = tmp

        cost_matrix = cdist(grids, X_embedded, "euclidean")

        cost_matrix = np.power(cost_matrix, 2)

        # row_asses, col_asses, info = fastlapjv(cost_matrix, k_value=50 if len(cost_matrix)>50 else len(cost_matrix))
        new_ori_row_asses, col_asses = self.solveJV(cost_matrix)

        # self.show_grid(new_ori_row_asses, labels, square_len, 'new2.png')

        new_row_asses, _, _ = self.grid_op(X_embedded, new_ori_row_asses, labels, useGlobal, useLocal,
                                                 convex_type, maxit, maxit2)
        # print("done")
        change = np.ones(shape=N, dtype='bool')
        new_row_asses = np.array(gridlayoutOpt.optimizeInnerCluster(ori_embedded, new_row_asses, labels, change))
        return new_row_asses

    # draw grids
    def show_grid(self, row_asses, grid_labels, square_len, path='new.png', showNum=True, just_save=False):
        def highlight_cell(x, y, ax=None, **kwargs):
            rect = plt.Rectangle((x - .5, y - .5), 1, 1, fill=False, **kwargs)
            ax = ax or plt.gca()
            ax.add_patch(rect)
            return rect

        from colors import MyColorMap
        cm = MyColorMap()

        data = []
        num = math.ceil(square_len)
        for i in range(num - 1, -1, -1):
            row = []
            for j in range(num):
                # if grid_labels[row_asses[num * i + j]]==-1:
                #     row.append((1, 1, 1, 1))
                # else:
                #     row.append(plt.cm.tab20(grid_labels[row_asses[num * i + j]]))
                row.append(cm.color(grid_labels[row_asses[num * i + j]]))
            data.append(row)
        plt.cla()
        plt.imshow(data)
        for i in range(num - 1, -1, -1):
            for j in range(num):
                highlight_cell(i, j, color="white", linewidth=1)
        if showNum:
            for i in range(num):
                for j in range(num):
                    text = plt.text(j, num - i - 1, row_asses[num * i + j], fontsize=7, ha="center", va="center",
                                    color="w")
                    # text = plt.text(j, num - i - 1, grid_labels[row_asses[num * i + j]], fontsize=7, ha="center", va="center",
                    #                 color="w")
        plt.axis('off')
        plt.savefig(path)
        if not just_save:
            plt.show()