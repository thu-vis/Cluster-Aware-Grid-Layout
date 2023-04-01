Cluster-Aware gridlayout Algorithm
==================================================================

Codes for gridlayout algorithm described in our paper ["Cluster-Aware Grid Layout"](https://xxxx) (xxxx).

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

Node: Download data in http://xxxx, and unzip in the root directory.
      Also need to move the compiled c++ package to the same directory as the py file.
```
python test_now/grid_clean.py
```

## Acknowledgement
This code is also used the published code of [fast-lapjv](https://github.com/thu-vis/fast-lapjv), and it is our pleasure to acknowledge their contributions.

## Contact
If you have any problem about our code, feel free to contact
- yx-zhou19@mails.tsinghua.edu.cn

or describe your problem in Issues.