.. file:combpotential class

.. _combpotential class:



============================
CombPotential class
============================

This class constructs the Comb potential for silicon and silica [#]_ as a subclass of :class:`~pysic.interactions.compound.CompoundPotential`.

Load the potential with::

 >>> from pysic.interactions.comb import CombPotential
 >>> import pysic     
 >>> calc = pysic.Pysic()
 >>> pot = CombPotential()     
 >>> pot.set_calculator(calc, True)

The potential is only defined for Si and SiO and parameters have been hard coded, so symbols need not be set. The potential also automatically includes screened :class:`~pysic.interactions.coulomb.CoulombSummation`, which is why the calculator should be given to it explicitly with :meth:`~pysic.interactions.CombPotential.set_calculator`.


Comb is a fairly complicated potential with many components such as local attractive and repulsive terms, bond angle terms and electrostatics. The class allows the user to switch terms on and off in order to analyze their contribution.


.. [#] T.-R. Shan, D. Bryce, J. Hawkins, A. Asthagiri, S. Phillpot, and S. Sinnott, Phys Rev B 82, 235302 (2010).



Full documentation of the CombPotential class
----------------------------------------------------

.. currentmodule:: pysic.interactions.comb
.. autoclass:: CombPotential
   :members:
   :undoc-members:

