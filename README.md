
A version of [Flagmati-dev](https://github.com/jsliacan/flagmatic-dev) that is being made Sage-friendly.


Flagmatic
=============

A package for the Sage.

To install for the first time, type (this will remove other versions of flagmatic and install this one):

    $ cd pkg
    $ cp install-dev-flagmatic /usr/local/bin/
    $ sudo chmod +x /usr/local/bin/install-dev-flagmatic
    $ cd; install-dev-flagmatic

To install again (e.g. after having installed another version of flagmatic in the meantime), type:

    $ install-dev-flagmatic

Then, from the Sage prompt, type:

    sage: from flagmatic.all import *
