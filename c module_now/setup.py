from setuptools import setup, Extension
import pybind11, numpy
# python3 -m pybind11 --includes
functions_module = Extension(  
    name ='gridlayoutOpt',  
    sources = ['mainop_alter_ac.cpp'],
    include_dirs = [pybind11.get_include()],
    extra_compile_args=["-fopenmp", "-std=c++11"],
    extra_link_args=["-fopenmp", "-std=c++11"],
    language='c++',
)  
  
setup(ext_modules = [functions_module])
