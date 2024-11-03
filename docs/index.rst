Diced |Stars|
=============

.. |Stars| image:: https://img.shields.io/github/stars/althonos/diced.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/diced/stargazers

*A Rust re-implementation of the* `MinCED <https://github.com/ctSkennerton/minced>`_ *algorithm to Detect Instances of* `CRISPRs <https://en.wikipedia.org/wiki/CRISPR>`_ *in Environmental Data.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/diced/python.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/diced/actions

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/diced?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/diced/

.. |PyPI| image:: https://img.shields.io/pypi/v/diced.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/diced

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/diced?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/diced

.. |AUR| image:: https://img.shields.io/aur/version/python-diced?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-diced

.. |Wheel| image:: https://img.shields.io/pypi/wheel/diced?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/diced/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/diced.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/diced/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/diced.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/diced/#files

.. |License| image:: https://img.shields.io/badge/license-GPLv3+-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/gpl-3.0/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/diced/

.. |Mirror| image:: https://img.shields.io/badge/mirror-LUMC-001158?style=flat-square&maxAge=2678400
   :target: https://git.lumc.nl/mflarralde/diced/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/diced.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/diced/issues

.. |Docs| image:: https://img.shields.io/readthedocs/diced?style=flat-square&maxAge=3600
   :target: http://diced.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/diced/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/pypi/dm/diced?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/diced


.. currentmodule:: diced


Overview
--------

MinCED is a method developed by `Connor T. Skennerton <https://github.com/ctSkennerton>`_
to identify `Clustered Regularly Interspaced Short Palindromic Repeats (CRISPRs) <https://en.wikipedia.org/wiki/CRISPR>`_
in isolate and metagenomic-assembled genomes. It was derived from the CRISPR 
Recognition Tool. It uses a fast scanning algorithm to identify candidate 
repeats, combined with an extension step to find maximally spanning regions 
of the genome that feature a CRISPR repeat.

Diced is a Rust reimplementation of the MinCED method, using the original
Java code as a reference. It produces exactly the same results as MinCED,
corrects some bugs, and is much faster. The Diced implementation is available 
as a Rust library for convenience.

.. grid:: 1 2 3 3
   :gutter: 1

   .. grid-item-card:: :fas:`battery-full` Batteries-included

      Diced is a Python package, so you can add it as a dependency to your 
      project, and stop worrying about the `minced` binary invoking the Java 
      Virtual Machine.

   .. grid-item-card:: :fas:`screwdriver-wrench` Flexible I/O

      Directly pass sequence data as Python `str` objects, and retrieve the 
      results using an iterator.

   .. grid-item-card:: :fas:`microchip` Fast

      The Java code uses a handwritten implementation of the `Boyer-Moore algorithm <https://en.wikipedia.org/wiki/Boyer%E2%80%93Moore_string-search_algorithm>`_, 
      while the Rust implementation uses the ``str::find`` method of the standard 
      library, which uses the `Two-way algorithm <https://en.wikipedia.org/wiki/Two-way_string-matching_algorithm>`_. 

   .. grid-item-card:: :fas:`memory` Memory-efficient

      The Rust code powering Diced is *zero-copy*: it will work without 
      copying the sequence data from the Python memory space.

   .. grid-item-card:: :fas:`check` Consistent results 

      Get the same results as MinCED ``v0.4.2+10f0a26e``.

   .. grid-item-card:: :fas:`dolly` Pre-built packages

      Get the pre-built wheels from PyPI for a fast installation on a 
      variety from machines, or compile from source using ``maturin``.
      


Setup
-----

Run ``pip install diced`` in a shell to download the latest release 
from PyPi, or have a look at the :doc:`Installation page <guide/install>` to find 
other ways to install ``diced``.


Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. include:: related.rst


License
-------

This library is provided under the `GNU General Public License v3.0 or later <https://choosealicense.com/licenses/gpl-3.0/>`_.
The code for this implementation was derived from the 
`MinCED source code <https://github.com/ctSkennerton/minced>`_, which is 
available under the GPLv3 as well. See the 
:doc:`Copyright Notice <guide/copyright>` section for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by the original MinCED authors. It was was developed by* 
`Martin Larralde <https://github.com/althonos/>`_ *during his PhD project at the* 
`Leiden University Medical Center <https://www.lumc.nl/en/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.

