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

    tsne_path = os.path.join(tsne_result, dataset, "position_pred.npz")
    tsne_file = np.load(tsne_path, allow_pickle=True)
    embeddings = tsne_file['positions']
    labels = tsne_file['gtlabels']
    plabels = tsne_file['plabels']
    label_names = tsne_file['label_names']

    import collections

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
    s_embeddings = embeddings[samples]
    print(s_embeddings.shape)
    s_labels = labels[samples]
    print(s_labels.shape)

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

    for type in ["O", "Triple", "PerimeterRatio"]:
        use_type = type
        # showT = "-NoneText"
        for showT in ["-NoneText"]:

            # for op_type in ["base", "global", "local", "full"]:
            for op_type in ["base", "global", "full"]:
                if (type != "O") and ((op_type == "base") or (op_type == 'global')):
                    continue
                if (type == "O") and ((op_type != "base") and (op_type != 'global')):
                    continue
                if ((type != "PerimeterRatio") and (type != "Triple")) and (op_type == "local"):
                    continue

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

                # file_path = os.path.join(path, type + "-" + op_type + showT + ".png")
                # save_path = os.path.join(path, type + "-" + op_type + ".npz")
                file_path = os.path.join(path, type + "-" + op_type + showT + ".png")
                save_path = os.path.join(path, type + "-" + op_type + ".npz")

                Optimizer = gridOptimizer()

                # generate baseline and ours
                row_asses_m, t1, t2, new_labels, ori_row = Optimizer.grid(s_embeddings, s_labels, use_type, m1, m2,
                                                                          use_global, use_local, pred_labels=plabels[samples])

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
                    np.savez(save_path, row_asses=row_asses_m, labels=new_labels, samples=s_samples, ori_row=ori_row)
                name = "\'" + dataset + "\'-" + str(grid_width) + "-" + str(name_flag) + "-" + type + "-" + op_type
                print(name, tt)

                print(show_labels.max())
                if show_labels.max() < 30:
                    Optimizer.show_grid(row_asses_m, show_labels, grid_width, file_path, showText, just_save=True)


# for dataset in ["MNIST", "CIFAR10", "USPS", "Animals", "Weather", "Wifi", "Isolet", "Indian food", "Texture", "StanfordDog-10", "OoDAnimals"]:
for dataset in ["MNIST"]:
    # for grid_width in [20, 30, 40]:
    for grid_width in [30]:
        for flag in ['all']:
        # for flag in [7]:
            max_tt = 1
            for tt in range(max_tt):
                body(dataset, grid_width, flag, tt+80, tt % 8)