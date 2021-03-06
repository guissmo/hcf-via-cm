# hcf-via-cm

This is a toy implementation used to generate examples for the preprint of `Computing the Hilbert Class Fields of Quartic CM Fields Using Complex Multiplication`, as referenced in [this page](https://math.guissmo.com/pubs.php).

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

## Installation

1. Clone the repository.
2. Edit `sage` script to point to the relevant directories.
   - `DOT_SAGE` is where the `init.sage` and `hist.sage` are expected to be found.
   - `LOGFILEDIR` is where you would like SAGE to save its log files.
   - `RECIPDIR` is where Marco Streng's [RECIP](https://github.com/mstreng/recip) installation is found.
   - `FGAGDIR` is where your __finitely generated abelian group__ PARI/GP scripts are found.
   - `SHIMURADIR` is where your [`shimura.gp`](https://gitlab.inria.fr/cmh/cmh/-/blob/master/scripts/shimura.gp) and `shimura.macros.gp` are found.
     - Note that `shimura.gp` is a script in Andreas Enge and Emmanuel Thomé's [cmh](https://gitlab.inria.fr/cmh/). Its license is GPLv3 or later.
     - The `shimura.gp` in this repository is from an old version of the file. The rest of the script has not yet been tested on newer versions, but is expected to work.
   - Replace `~/.magma/2.24-3/` by where your `magma` installation lies.
3. Make `sage` executable.

## Usage

1. Run `.sage`.
2. Load `CMdata.sage` via `load("CMdata.sage")`

## Example 

This is an example usage to reproduce [Example 4.11](https://math.guissmo.com/pubs.php) as found in the preprint _Computing the Hilbert Class Fields of Quartic CM Fields Using Complex Multiplication_.

1. `A = ShGrpData([809, 53, 500],m=2)`
2. `pol = A.find_defining_polynomial_and_verify(rosenhain_invariants(2)[1],prec=400,verbose=True,check_conductor=True,autoretry=5)`
3. `pol[0]` contains a polynomial that defines the extension $H_{K^r}(1) / K^r$.

__Caveat__: When using the polynomial to define an extension of (the SAGE object representing the) reflex field, SAGE attempts a `polredbest`, which is slow. I personally did not wait for SAGE to finish and instead used PARI/GP directly to verify. This is done by doing:

4. `paripol = poly_s2p(pol[0], A.prelt['s2p'])`
5. `gp.rnfconductor(A.Kr_gp, [paripol, 10^6])`

Verifying that the result is `1`, we find that the conductor of the extension has no prime factors below $10^6$. If one is patient enough and wishes for PARI to find the conductor itself, replace `[paripol, 10^6]` with `paripol`.

## To-Do

* Documentation
* Handle cases when neither $\star_1$ nor $\star_2$ is true.
