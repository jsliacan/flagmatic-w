
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

*Sage is needed.* Tested with Sage-7.4 on Mac OSX Sierra and Ubuntu 16.04

From a folder where your certificates are and where you have writing permissions, type:

    $ sage -python inspect_certificate.py <cert1> --stability <bound> <tau> <B> [[<cert2>] <cert3>]

The arguments are:
  
  * `<cert1>` is the filename of the flag algebras certificate for the problem you want to verify stability for, e.g. `cert_flag_alg.js`
  * `<bound>` is the lower bound for this problem (should match the FA upper bound), e.g. `24/625`
  * `<tau>` is the type you chose, e.g. `"3:12"`
  * `<B>` is the B-graph you chose, e.g. `"5:1223344551"`
  * `<cert2>` is the filename of the FA certificate with, additionally, tau forbidden as induced subgraph, e.g. `cert_robust_stab.js`; not needed, if tau = `1:`.
  * `<cert3>` is the filename of the FA certificate with, additionally, B forbidden as induced subgraph, e.g. `cert_perfect_stab.js`; not needed if Q_tau has corank 1.

EXAMPLE:
Prove that the Turán C<sub>5</sub>-density of a triangle is 24/625. Prove that this problem is perfectly-C<sub>5</sub>-stable. Given are certificates `cert.js` and `cert_tau.js`. 

    $ sage -python inspect_certificate.py cert.js --stability 24/625 "3:12" "5:1223344551" cert_tau.js

Check the output, it indicates if conditions of the Theorem 7.1 are met.