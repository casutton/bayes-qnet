from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

incdirs=["/opt/local/include/", "/home/eecs/casutton/include", "/exports/home/csutton/usr/include"]
libdirs=["/opt/local/lib/", "/home/eecs/casutton/lib", "/exports/home/csutton/usr/lib"]

incdirs = [ l for l in incdirs if os.path.exists(l) ]
libdirs = [ l for l in libdirs if os.path.exists(l) ]

cflags = [ ]

setup(
  name = "dist",
  ext_modules=[ 
    Extension("cninq", ["cninq.pyx"], extra_compile_args=cflags),
    Extension("misc", ["misc.pyx"], extra_compile_args=cflags),
    Extension("pwfun", ["pwfun.pyx"], extra_compile_args=cflags),
    Extension("sampling", ["sampling.pyx", "randomkit.c", "cdist.c"], extra_compile_args=cflags),
    Extension("arrivals", ["arrivals.pyx"], extra_compile_args=cflags),
    Extension("queues", ["queues.pyx"], extra_compile_args=cflags),
    Extension("qnet", ["qnet.pyx"], extra_compile_args=cflags),
    Extension("distributions", ["distributions.pyx"], extra_compile_args=cflags),
    Extension("pyglpk", ["pyglpk.pyx"], extra_compile_args=cflags, include_dirs=incdirs, libdirs=libdirs, libraries=["glpk"]),
    ],
  cmdclass = {'build_ext': build_ext}
)
