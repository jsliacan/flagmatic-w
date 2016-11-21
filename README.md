
A version of [Flagmatic-dev](https://github.com/jsliacan/flagmatic-dev) that works with Sage-7.4.


flagmatic
=============

To install, navigate to the directory with `setup.py` file:

    $ cd path/to/flagmatic/pkg

Then run the following code to install CSDP solver:

    $ sage -i csdp

or alternatively, download the CSDP package from [https://projects.coin-or.org/Csdp/](https://projects.coin-or.org/Csdp/) and copy `csdp` to your PATH.

Finally, install Flagmatic with the code below.

    $ sage -python setup.py install

Have fun.
