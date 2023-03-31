from setuptools import setup, Extension
import pybind11, numpy

# python3 -m pybind11 --includes
functions_module = Extension(
    name='gridlayoutOpt',
    sources=['mainop_alter_ac.cpp'],
    include_dirs=[pybind11.get_include()],
    extra_compile_args=["/openmp", "-openmp:llvm"],
    extra_link_args=["/openmp"],
    language='c++',
)

setup(ext_modules=[functions_module])
