.. rfbzero documentation master file, created by
   sphinx-quickstart on Sat Oct 14 13:30:52 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

rfbzero.py
=============

:code:`rfbzero.py` is a Python package for zero dimensional simulation of
electrochemical cycling in redox flow batteries (RFBs).

This package contains modules to describe initial flow cell setup, chemical
and electrochemical properties of the redox-active electrolytes being cycled,
cell cycling protocol selection, and optional inputs for capacity degradation
mechanisms and active species crossover.


If you have a feature request or find a bug, please
`file an issue <https://github.com/ericfell/rfbzero/issues>`_
or contribute code improvements and
`submit a pull request <https://help.github.com/articles/creating-a-pull-request-from-a-fork/>`_!

Installation
------------

The easiest way to install :code:`rfbzero.py` is
from PyPI https ????
using pip:

.. code-block:: bash

   pip install ??

See :doc:`./getting-started` for instructions
on getting started from scratch.


Dependencies
~~~~~~~~~~~~

rfbzero.py requires:

-   Python (>=3.10)
-   SciPy


Example Jupyter notebooks are provided in the examples/ directory.

Examples and Documentation
---------------------------

The documentation can be found at
`rfbzero.readthedocs.io <https://ericfell-rfbzero.readthedocs.io/en/latest/>`_.



.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting-started
   examples
   flowcell
   experiment
   degradation
   crossover
   faq

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
