
A version of Flagmatic with extra functionality such as assumptions and stability checks. This version was updated to run on SageMath 9.X and several of the core computations have been parallelized using Python's multiprocessing.


flagmatic
=============

To install, navigate to the directory with `setup.py` file:

    $ cd path/to/flagmatic/pkg

Then run the following code to install CSDP solver:

    $ sage -i csdp

or alternatively, download the CSDP package from [https://projects.coin-or.org/Csdp/](https://projects.coin-or.org/Csdp/) and copy `csdp` to your PATH.

Finally, install Flagmatic with the command below.

    $ sage -python setup.py install

**POSSIBLE ISSUES:**

  * Check whether the installation finished without errors. If you're missing permissions, re-run the last command as a superuser. 
  * Make sure you refresh/restart your Terminal after the install for changes to take effect. 
  * If you are getting an error about a missing `gcc` command you are probably on a Mac and  you will need to install Command Line Tools by typing `xcode-select --install` in Terminal.

### Running example scripts ###

Make sure that the script begins with `from flagmatic.all import *` or that you have imported everything yourself. Assuming the script is called `test.sage`, navigate to the directory with the script, start sage, and type

    $ runfile test.sage

Alternatively, you can type in the script into the Sage interpreter one line at a time. 

### Verifying stability ###

**Sage is needed.** Tested with Sage7.4 on Mac OSX Sierra and Ubuntu 16.04.

Assuming your certificates and inspect_certificate.py are in the same directory and that you have writing permissions there, from that directory type:

    $ sage -python inspect_certificate.py <cert1> --stability <bound> <tau> <B> [<cert2> [<cert3>]]

(otherwise use full paths for filenames).
The arguments are:
  
  * `<cert1>` is the filename of the flag algebras certificate for the problem you want to verify stability for, e.g. `c5.js`
  * `<bound>` is the lower bound for this problem (should match the FA upper bound), e.g. `24/625`
  * `<tau>` is the type you chose, e.g. `"3:12"`
  * `<B>` is the B-graph you chose, e.g. `"5:1223344551"`
  * `<cert2>` is the filename of the FA certificate with, additionally, tau forbidden as induced subgraph, e.g. `cert_robust_stab.js`; not needed, if tau = `1:`.
  * `<cert3>` is the filename of the FA certificate with, additionally, B forbidden as induced subgraph, e.g. `cert_perfect_stab.js`; not needed if Q_tau has corank 1.

**EXAMPLE:**
Prove that if triangles are forbidden, C<sub>5</sub> density is at most 24/625. Prove that this problem is perfectly-C<sub>5</sub>-stable. Given are certificates `c5.js` and `c5_stab.js`. 

    $ sage -python inspect_certificate.py c5.js --stability 24/625 "3:12" "5:1223344551" c5_stab.js

Check the output, it indicates if conditions of the Theorem 7.1 in *Strong forms of stability from flag algebra calculations* are met.

**WARNING:**
When verifying stability for Clebsch Graph, make sure you use the following string: "g:12131415162728292a373b3c3d484b4e4f595c5e5g6a6d6f6g7e7f7g8c8d8g9b9d9fabacaebgcfde". It's hardcoded for speed.
