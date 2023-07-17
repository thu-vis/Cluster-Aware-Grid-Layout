Cluster-Aware Grid Layout Algorithm
===========================

https://user-images.githubusercontent.com/129478716/229704745-f38ecd8d-d53b-4cf7-a061-42c2d1005834.mp4

=======================================

Codes for grid layout algorithm described in our paper ["Cluster-Aware Grid Layout"](https://xxxx) (TVCG 2023).

Compile
----------
```
pip install pybind11
cd "c module_now"
python setup.py build_ext --inplace
```
Tested on Windows, and Linux with python 3.8.

Usage Example
-----
Run this example to compare our method and baseline method. Runing this example requires the installation of [scipy](https://www.scipy.org/), [numpy](https://www.numpy.org/) and [matplotlib](https://github.com/matplotlib/matplotlib)

Node: Download data in [Google Drive](https://drive.google.com/file/d/11Jas4sRB1uAhCDsYQQBppyc-55SFyl50/view?usp=share_link), and unzip in the root directory.
      Also need to move the compiled c++ package to the same directory as the py file.
```
python test_now/demo.py
```

## Acknowledgement
This code is also used the published code of [fast-lapjv](https://github.com/thu-vis/fast-lapjv), and it is our pleasure to acknowledge their contributions.

## Citation
If you use this code for your research, please consider citing:
```
@article{zhou2023clusteraware,
    title = {Cluster-Aware Grid Layout},
    author = {Zhou, Yuxing and Yang, Weikai and Chen, Jiashu and Chen, Changjian and Shen, Zhiyang and Luo, Xiaonan and Yu, Lingyun and Liu, Shixia},
    journal = {IEEE Transactions on Visualization and Computer Graphics (accepted)},
    year = {2023}
}
```

## Contact
If you have any problem about our code, feel free to contact
- zhouyx23@mails.tsinghua.edu.cn

or describe your problem in Issues.
