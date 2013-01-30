.. file:suttonchenpotential class

.. _suttonchenpotential class:



============================
SuttonChenPotential class
============================

This class constructs the Sutton-Chen potential [#]_ as a subclass of :class:`~pysic.interactions.compound.CompoundPotential`.

This potential contains 2-body and many-body terms defined as

.. math::

  U & = \varepsilon \left[ \sum_{i,j} \left( \frac{a}{r_{ij}} \right)^n - c \sum_i \sqrt{\rho_i} \right] \\
  \rho_i & = \sum_j \left( \frac{a}{r_{ij}} \right)^m

where :math:`r_{ij}` is the interatomic distance, and :math:`\varepsilon, a, c, m, n` are parameters.

Load the potential with::

 >>> from pysic.interactions.suttonchen import SuttonChenPotential
 >>> pot = SuttonChenPotential(symbols=...,
 ...                           tags=...,
 ...                           indices=...,
 ...                           parameters=...,
 ...                           cutoff=...,
 ...                           cutoff_margin=...)

The same cutoff is used for both the 2-body and many-body terms.

If the potential is given several target pairs, e.g., ``symbols = [['A', 'B'], ['A', 'C']]``,
all possible combinations are targeted by default. For the above list, the 2-body interaction would be
evaluated for A-B and A-C pairs (not B-C), and the many-body interaction would be evaluated
for elements A with A-B and A-C pairs, for element B with B-A terms, and for element C with C-A
terms. If the many-body term should be only evaluated for certain elements, the method
:meth:`~pysic.interactions.suttonchen.SuttonChenPotential.set_density_symbols` can be used
for specifying the specific target symbols. For instance::

 >>> pot.set_symbols( [['A', 'B'], ['A', 'C']] )
 >>> pot.set_density_symbols( ['A'] )

will have the 2-body terms evaluated for A-B and A-C pairs while the many-body term is only evaluated
for neighborhoods of element A, taking into account the pairs A-B and A-C.


.. [#] Sutton, A. P., and Chen, J., 1990, Philos. Mag. Lett., 61, 139.



Full documentation of the SuttonChenPotential class
----------------------------------------------------

.. currentmodule:: pysic.interactions.suttonchen
.. autoclass:: SuttonChenPotential
   :members:
   :undoc-members:

