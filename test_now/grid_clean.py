import random

dict = {}
dict2 = {}
dict0 = {}
save_file = True

def body(dataset, grid_width, flag, tt, rid):
    tmp_dict = {}
    full_flag = 0

    import numpy as np
    np.set_printoptions(suppress=True)
    np.set_printoptions(precision=4)
    import os
    tsne_result = "../tsne_result"

    # feature_path = os.path.join(processed_data, "{}_full.npz".format(dataset))
    tsne_path = os.path.join(tsne_result, dataset, "position_pred.npz")
    # features = np.load(feature_path, allow_pickle=True)['features']
    tsne_file = np.load(tsne_path, allow_pickle=True)
    embeddings = tsne_file['positions']
    labels = tsne_file['gtlabels']
    plabels = tsne_file['plabels']
    label_names = tsne_file['label_names']

    import collections

    print(collections.Counter(labels).keys())
    print(collections.Counter(plabels).keys())

    labels_num = len(collections.Counter(labels).keys())

    if dataset == "Indian food":
        labels_num = 11
    if dataset == "Weather":
        labels_num = 4
    if dataset == "Isolet":
        labels_num = 8
    if dataset == "StanfordDog-10":
        labels_num = 7

    name_flag = flag

    if flag == 'all':
        flag = labels_num

    print(label_names)
    print(labels.shape[0])
    # if labels.shape[0] < grid_width*grid_width:
    #     return

    size = labels.shape[0]
    # grid_width = 45
    # flag = 7
    path = "./grid_result"
    if not os.path.exists(path):
        os.mkdir(path)

    path = os.path.join(path, dataset)
    if not os.path.exists(path):
        os.mkdir(path)

    path = os.path.join(path, str(grid_width))
    if not os.path.exists(path):
        os.mkdir(path)

    path = os.path.join(path, str(name_flag) + "-" + str(grid_width) + "-" + str(tt))
    if save_file and not os.path.exists(path):
        os.mkdir(path)

    if save_file:
        for file_name in os.listdir(path):
            os.remove(os.path.join(path, file_name))

    samples_N = grid_width * grid_width

    resample = True
    if rid > 0:
        resample = False

    if resample:
        if flag == 'all':
            if (grid_width-1)*(grid_width-1) >= size:
                return
            if samples_N > size:
                samples_N = size
            samples = np.random.choice(size, samples_N, replace=False)
            np.save("samples_{}0.npy".format(dataset), samples)
        else:
            # if flag >= labels_num:
                # return
            if flag > labels_num:
                flag = random.randint(3, labels_num)

            import collections
            print(collections.Counter(labels))

            cnt = 0
            while cnt < 20:
                chosen_labels = np.random.choice(labels_num, flag, replace=False)
                # chosen_labels = np.array([1, 4, 7, 9])
                print(chosen_labels)

                idx = np.zeros(size, dtype='bool')
                for lb in chosen_labels:
                    idx = (idx | (labels == lb))

                idx2 = np.zeros(size, dtype='bool')
                for lb in chosen_labels:
                    idx2 = (idx2 | (plabels == lb))
                idx = idx & idx2

                print(idx)
                idx = (np.arange(size))[idx]
                chosen_size = idx.shape[0]
                print('chosen_size', chosen_size)

                if (grid_width - 1) * (grid_width - 1) < chosen_size:
                    break

                if (grid_width - 1) * (grid_width - 1) < chosen_size*6:
                    m_idx = idx.copy()
                    while m_idx.shape[0] <= (grid_width - 1) * (grid_width - 1):
                        m_idx = np.concatenate((m_idx, idx))
                    idx = m_idx.copy()
                    chosen_size = idx.shape[0]
                    break

                cnt += 1

            if (grid_width-1)*(grid_width-1) >= chosen_size:
                return

            if samples_N > chosen_size:
                samples_N = chosen_size

            samples = np.random.choice(chosen_size, samples_N, replace=False)
            samples = idx[samples]
            np.save("samples_{}0.npy".format(dataset), samples)

    samples = np.load("samples_{}0.npy".format(dataset))
    # s_features = features[samples]
    s_embeddings = embeddings[samples]
    print(s_embeddings.shape)
    s_labels = labels[samples]
    print(s_labels.shape)
    # print(samples.shape[0], s_features.shape)

    import math
    if rid > 0:
        r = rid*math.pi/16
        for i in range(s_embeddings.shape[0]):
            tmp_x = s_embeddings[i][0]
            tmp_y = s_embeddings[i][1]
            s_embeddings[i][0] = tmp_x*math.cos(r) - tmp_y*math.sin(r)
            s_embeddings[i][1] = tmp_x*math.sin(r) + tmp_y*math.cos(r)

    import numpy as np

    new_labels = s_labels.copy()

    from gridOptimizer_clean import gridOptimizer
    import os

    # type = "T"
    # # showT = ""
    # showT = "-NoneText"

    for type in ["O", "Triple", "PerimeterRatio"]:
        use_type = type
        # for showT in ["", "-NoneText"]:
        showT = "-NoneText"
        for choose_k in [0, 1]:
            if (type != "O") and (choose_k == 0):
                continue
            if (type == "O") and (choose_k != 0):
                continue

            # for op_type in ["base", "compact", "global", "full"]:
            for op_type in ["base", "global", "local", "full", "full2"]:
            # for op_type in ["base", "global", "full"]:
            # for op_type in ["base", "global", "full"]:
                if (type != "O") and ((op_type == "base") or (op_type == 'global')):
                    continue
                if (type == "O") and ((op_type != "base") and (op_type != 'global')):
                    continue
                if ((type != "PerimeterRatio") and (type != "Triple")) and ((op_type == "local") or (op_type == "full2")):
                    continue

                swap_op_order = False
                if op_type == "full2":
                    swap_op_order = True

                m1 = 0
                m2 = 3
                if type == "CutRatio":
                    m1 = 0
                    m2 = 3
                if type == "PerimeterRatio":
                    m1 = 0
                    m2 = 8
                if type == "Triple":
                    m1 = 5
                    m2 = 0
                if type == "AreaRatio":
                    m1 = 0
                    m2 = 3
                if type == "O":
                    m1 = 0
                    m2 = 0

                use_local = True
                if (op_type == "base") or (op_type == "global") or (op_type == "compact"):
                    m1 = 0
                    m2 = 0
                    use_local = False

                use_global = True
                if (op_type == "base") or (op_type == "local"):
                    use_global = False

                only_compact = False
                if op_type == "compact":
                    only_compact = True

                # file_path = os.path.join(path, type + "-" + op_type + showT + ".png")
                # save_path = os.path.join(path, type + "-" + op_type + ".npz")
                file_path = os.path.join(path, type + "-" + op_type + showT + "-" + str(choose_k) + ".png")
                save_path = os.path.join(path, type + "-" + op_type + "-" + str(choose_k) + ".npz")

                Optimizer = gridOptimizer()
                # print("check done", BASolver.checkConvex(np.array(row_asses_c), np.array(s_labels)))
                # row_asses_m, heat = BASolver.grid3(s_embeddings, s_labels, 'E')
                final_it, cv_time, g_time, g_time2 = 0, 0, 0, 0
                row_asses_m, t1, t2, new_labels, new_cost, cc, ori_row = Optimizer.grid(s_embeddings, s_labels, use_type, m1, m2,
                                                                           use_global, use_local, only_compact, swap_op_order=swap_op_order, swap_cnt=2147483647, pred_labels=plabels[samples], choose_k=choose_k)

                show_labels = new_labels
                show_labels = np.array(show_labels)

                if show_labels.max() >= 15:
                    labels_dict = {}
                    dict_num = 0
                    for i in range(show_labels.shape[0]):
                        if show_labels[i] not in labels_dict:
                            labels_dict.update({show_labels[i]: dict_num})
                            dict_num += 1
                        show_labels[i] = labels_dict[show_labels[i]]

                tmp = np.full(row_asses_m.shape[0] - show_labels.shape[0], dtype='int', fill_value=-1)
                show_labels = np.concatenate((show_labels, tmp), axis=0)
                print(row_asses_m.shape[0], show_labels.shape[0])

                print(t1, t2)
                showText = True
                if showT == "-NoneText":
                    showText = False
                    s_samples = samples
                    if dataset=="OoDAnimals":
                        s_samples = tsne_file['true_id'][samples]
                    if dataset=="OoDAnimals3":
                        s_samples = tsne_file['true_id'][samples]
                    np.savez(save_path, row_asses=row_asses_m, labels=new_labels, samples=s_samples, ori_row=ori_row)
                print("new_cost", new_cost)
                name = "\'" + dataset + "\'-" + str(grid_width) + "-" + str(name_flag) + "-" + type + "-" + op_type + "-" + str(choose_k)
                print(name, tt)

                # new_cost = np.append(new_cost, [t1 + t2, t2], None)
                if op_type == "base":
                    new_cost = np.append(new_cost, [t2], None)
                else:
                    # new_cost = np.append(new_cost, [t1+t2], None)
                    new_cost = np.append(new_cost, [t1], None)

                cflag = 0
                new_cost[0] = np.exp(-new_cost[0] / grid_width / grid_width)
                new_cost[1] = np.exp(-new_cost[1] / grid_width / grid_width)
                new_cost[2] = 1 - new_cost[2] / grid_width / grid_width
                new_cost[3] = 1 - new_cost[3] / grid_width / grid_width
                new_cost[4] = 1 - new_cost[4] / grid_width / grid_width
                new_cost[5] = 1 - new_cost[5] / grid_width / grid_width
                new_cost[6] = 1 - new_cost[6] / grid_width / grid_width
                new_cost[7] = 1 - new_cost[7] / grid_width / grid_width
                new_cost[8] = 1 - new_cost[8] / grid_width / grid_width

                if (op_type != "base")and(op_type != "global")and(op_type != "full2")and(cc > 0):
                    cflag = 1
                    full_flag = 1

                print(cc)
                tmp_dict.update({name: new_cost.copy()})

                if cflag == 0:
                    if name not in dict:
                        dict.update({name: new_cost.copy()})
                        dict2.update({name: 1})
                    else:
                        dict[name] += new_cost
                        dict2[name] += 1

                # name0 = name+"-"+str(tt)
                # dict0.update({name0: new_cost.copy()})

                print(show_labels.max())
                if show_labels.max() < 30:
                    Optimizer.show_grid(row_asses_m, show_labels, grid_width, file_path, showText, just_save=True)
                # Optimizer.show_grid(row_asses_m, show_labels, grid_width, "E-full-NoneText.svg", showText)
                # Optimizer.show_grid(row_asses_m, show_labels, grid_width, "test"+name+".png", showText)
                # Optimizer.show_grid(row_asses_m, s_labels, grid_width)


# for dataset in ["MNIST", "STL-10", "CIFAR10", "USPS",  "Weather", "Clothes", "FashionMNIST", "Animals", "Indian food", "Wifi"]:
# for dataset in ["Animals", "Indian food", "Wifi"]:
# for dataset in ["Clothes"]:
# for dataset in ["StanfordDog-10"]:
# for dataset in ["OoDAnimals"]:
# for dataset in ["MNIST", "CIFAR10", "USPS", "Animals", "Weather", "Wifi", "Isolet", "Indian food", "Texture", "StanfordDog-10", "OoDAnimals"]:
# for dataset in ["MNIST", "CIFAR10", "USPS", "Animals", "Wifi"]:
for dataset in ["Indian food"]:
    # for grid_width in [20, 30, 40]:
    # for grid_width in [20, 30, 40]:
    for grid_width in [30]:
        for flag in ['all']:
        # for flag in [7]:
            max_tt = 1
            # if grid_width == 20:
            #     max_tt = 50
            for tt in range(max_tt):
                body(dataset, grid_width, flag, tt+80, tt % 8)
    import numpy as np
    import pickle
    f = open("grid_result/"+dataset+"/data.pkl", 'wb+')
    pickle.dump({"dict": dict, "dict2": dict2}, f, 0)
    f.close()
    # np.savez("grid_result/"+dataset+"/data.npz", dict=dict, dict2=dict2)


for key in dict:
    t = -1
    dict[key] /= dict2[key]
    print(key, "----", "%.3lf"%dict[key][0], "&", "%.3lf"%dict[key][1], "&", "%.3lf"%dict[key][2], "&", "%.3lf"%dict[key][3], "&", "%.3lf"%dict[key][4], "&", "%.3lf"%dict[key][5], "&", "%.3lf"%dict[key][6], "&", "%.3lf"%dict[key][7], "&", "%.3lf"%dict[key][8], "&", "%.3lf"%dict[key][t])

import pickle
f = open("all_data3.pkl", 'wb+')
pickle.dump({"dict": dict}, f, 0)
f.close()