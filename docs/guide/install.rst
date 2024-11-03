Installation
============

.. caution::

    Wheels are provided for Linux, MacOS and Windows x86-64, but other machines 
    will have to build the wheel from the source distribution. Building ``diced``
    involves compiling the ``diced`` crate, which requires a Rust compiler to 
    be available on the local machine.


PyPi
^^^^

``diced`` is hosted on GitHub, but the easiest way to install it is to download
the latest release from its `PyPi repository <https://pypi.python.org/pypi/diced>`_.
It will install all dependencies then install ``diced`` either from a wheel if
one is available, or from source after compiling the Rust code :

.. code:: console

    $ pip install --user diced


.. Conda
.. ^^^^^

.. ``diced`` is also available as a `recipe <https://anaconda.org/bioconda/diced>`_
.. in the `bioconda <https://bioconda.github.io/>`_ channel. To install, simply
.. use the ``conda`` installer:

.. .. code:: console

..      $ conda install -c bioconda diced


.. Arch User Repository
.. ^^^^^^^^^^^^^^^^^^^^

.. A package recipe for Arch Linux can be found in the Arch User Repository
.. under the name `python-diced <https://aur.archlinux.org/packages/python-diced>`_.
.. It will always match the latest release from PyPI.

.. Steps to install on ArchLinux depend on your `AUR helper <https://wiki.archlinux.org/title/AUR_helpers>`_
.. (``yaourt``, ``aura``, ``yay``, etc.). For ``aura``, you'll need to run:

.. .. code:: console

..     $ aura -A python-diced


.. BioArchLinux
.. ^^^^^^^^^^^^

.. The `BioArchLinux <https://bioarchlinux.org>`_ project provides pre-compiled packages
.. based on the AUR recipe. Add the BioArchLinux package repository to ``/etc/pacman.conf``:

.. .. code:: ini

.. ..     [bioarchlinux]
..     Server = https://repo.bioarchlinux.org/$arch

.. Then install the latest version of the package and its dependencies with ``pacman``:

.. .. code:: console

..     $ pacman -Sy
..     $ pacman -S python-diced


Piwheels
^^^^^^^^

Diced works on Raspberry Pi computers, and pre-built wheels are compiled 
for `armv7l` on `piwheels <https://www.piwheels.org/project/diced/>`_.
Run the following command to install these instead of compiling from source:

.. code:: console

   $ pip3 install diced --extra-index-url https://www.piwheels.org/simple

Check the `piwheels documentation <https://www.piwheels.org/faq.html>`_ for 
more information.


GitHub + ``pip``
^^^^^^^^^^^^^^^^

If, for any reason, you prefer to download the library from GitHub, you can clone
the repository and install the repository by running (with the admin rights):

.. code:: console

    $ pip install -U git+https://github.com/althonos/diced

.. caution::

    Keep in mind this will install always try to install the latest commit,
    which may not even build, so consider using a versioned release instead.


GitHub + ``build``
^^^^^^^^^^^^^^^^^^

If you do not want to use ``pip``, you can still clone the repository and
use ``build`` and ``installer`` manually:

.. code:: console

    $ git clone --recursive https://github.com/althonos/diced
    $ cd diced
    $ python -m build .
    # python -m installer dist/*.whl

.. Danger::

    Installing packages without ``pip`` is strongly discouraged, as they can
    only be uninstalled manually, and may damage your system.
