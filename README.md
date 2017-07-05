# bayes-qnet

Code for Bayesian inference for queueing networks with incomplete data.

The statistical methods are described in this paper:

  C. Sutton and M. I. Jordan. 
  Bayesian inference in queueing networks. 
  Annals of Applied Statistics, 5(1):254–282, 2011.
  https://arxiv.org/abs/1001.3355

## System Requirements

This code requires Python 2 and scipy. Also requires
yaml.py

Significant portions
of the time critical portions of the code 
are written in Cython. The C translations
of the Cython modules are included in the git repo.

## Installation

You will need to create the C shared libraries
to use the Python code. If you have Cython properly
installed, this can be done using

```
  cd src/
  make
```

## Getting Started

Of course, there is no real documentation,
but I have created a Jupyter notebook that walks
through an example usage. This is in the `examples/`
directory. 

I've saved the notebook both in native Jupyter format,
so you can run the code one your machine, and as a static
HTML page, in case you do not have Jupyter installed.

Hopefully this is enough to get started.

## Test Suite

(unfortunately this part is from memory)

To test if you have compiled correctly, you can run the unit
tests. These are all in the `src/` directory and their
names begin with `test_`. Why I decided not to put them
in separate directory, who knows.

The most important one is test_qnet.py. So just

```
  cd src/
  python2 test_qnet.py
```

should work. A dozen of the tests will fail because I just
removed some unused (but moderately testing) functionality
to make the package a lot easier to install.


