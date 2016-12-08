
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

If you are getting an error about a missing `gcc` command you are probably on a Mac and  you will need to install Command Line Tools by typing `xcode-select --install` in Terminal.

### Verifying stability ###

* make sure you are in a folder where your certificates are and that you have writing permission in that folder
* type

  $ sage -python inspect_certificate.py \<certificateFA\> --stability \<bound\> \<tau\> \<B\> \<certificateTau\>

* where 
  
  * `<certificateFA>` is the filename of the flag algebras certificate for the problem you want to verify stability for, e.g. `cert_flag_alg.js`
  * `<bound>` is the lower bound for this problem (should match the FA upper bound), e.g. `24/625`
  * `<tau>` is the type you chose, e.g. `"3:12"`
  * `<B>` is the B-graph you chose, e.g. `"5:1223344551"`
  * `<certificateTau>` is the filename of the FA certificate with, additionally, tau forbidden as induced subgraph, e.g. `cert_robust_stab.js`


Have fun.
