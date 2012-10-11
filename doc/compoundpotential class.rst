.. file:compoundpotential class

.. _compoundpotential class:



.. file:compoundpotential class - description

.. _compoundpotential class - description:



============================
CompoundPotential class
============================

This class defines an interaction as a collection of :class:`~pysic.interactions.local.Potential`\ s,
:class:`~pysic.interactions.local.ProductPotential`\ s, :class:`~pysic.interactions.bondorder.Coordinator`\ s, and
:class:`~pysic.interactions.bondorder.BondOrderParameters`. It is a wrapper to ease the handling of complicated potentials.



.. file:using compoundpotential

.. _using compoundpotential:



Using CompoundPotential
---------------------------

The :class:`~pysic.interactions.compound.CompoundPotential` class does not define a potential itself, but it defines an abstract superclass which can be used for defining particular potentials as subclasses. One such subclass is the Sutton-Chen potential, which will be used as an example here.

Sutton-Chen potential combines a simple power-law potential with a similar density-like term

.. math::

  U & = \varepsilon \left[ \sum_{i,j} \left( \frac{a}{r_{ij}} \right)^n - c \sum_i \sqrt{\rho_i} \right] \\
  \rho_i & = \sum_j \left( \frac{a}{r_{ij}} \right)^m

In Pysic, this can be defined by a structure like this::

 >>> power_pot = pysic.Potential('power', ...)
 >>> rho = pysic.BondOrderParameters('power_bond', ...)
 >>> root = pysic.BondOrderParameters('sqrt_scale', ...)
 >>> sqrt_rho = pysic.Coordinator([root, rho])
 >>> rho_pot = pysic.Potential('constant', coordinator=sqrt_rho)
 >>> calc = pysic.Pysic()
 >>> calc.set_potentials([power_pot, rho_pot])

but this is already somewhat cumbersome. It is much nicer to wrap the definitions in a container::

 >>> sutton_chen_pot = SuttonChenPotential(...)
 >>> calc = pysic.Pysic()
 >>> calc.set_potentials(sutton_chen_pot)

This is the purpose of :class:`~pysic.interactions.compound.CompoundPotential` and its
subclasses. 

The central classes such as :class:`pysic.interactions.local.Potential` are imported
automatically with the :mod:`pysic` module, which is why you can refer to them with ``pysic.Potential``.
Compound potentials are not, so you need to manually import them, for instance like this::

 >>> from pysic.interactions.suttonchen import SuttonChenPotential
 >>> sutton_chen_pot = SuttonChenPotential(...)

The potentials provided with pysic are all found in `pysic.interactions.x` where `x` is
the name of the submodule containing the potential.



.. file:defining compoundpotentials

.. _defining compoundpotentials:



Defining new wrappers
------------------------

:class:`~pysic.interactions.compound.CompoundPotential` inherits the interface of 
:class:`~pysic.interactions.local.Potential` and it is used in a similar fashion.
A subclass of :class:`~pysic.interactions.compound.CompoundPotential` then also
inherits the interface, and so it is easy to create new potentials this way without 
having to care about the interface. Defining a subclass in Python is done like this::

 >>> class SuttonChenPotential(CompoundPotential):
 ...     def __init__(self, ...):
 ...         super(SuttonChenPotential, self).__init__(...)
 ...

The brackets in the class definition, ``(CompoundPotential)``, tell that the class
derives from another class and inherits its functionality. In the constructor, ``__init__``,
the superclass constructor is called with ``super(SuttonChenPotential, self).__init__(...)``,
but it is also possible to extend the constructor with additional commands.

Compound potentials store their components in a list called ``self.pieces``, and this is
done in the method :meth:`~pysic.interactions.compound.CompoundPotential.define_elements`.
The components are passed to the calculator via the 
:meth:`~pysic.interactions.compound.CompoundPotential.build` method. This is automatically
called a compound potential is given to the :meth:`pysic.calculator.Pysic.add_potential` method,
so that the interface is similar to that of simple potentials.
For the Sutton-Chen example, the definitions could look something like this::

 ...     def define_elements(self):
 ...         power_pot = pysic.Potential('power', ...)
 ...         rho = pysic.BondOrderParameters('power_bond', ...)
 ...         root = pysic.BondOrderParameters('sqrt_scale', ...)
 ...         sqrt_rho = pysic.Coordinator([root, rho])
 ...         rho_pot = pysic.Potential('constant', coordinator=sqrt_rho)
 ...         self.pieces = [power_pot, rho_pot]
 ...         

In principle, this is the only method that a new potential has to override, since the rest is done
automatically. Of course, it is possible
to override other methods as well or write new ones. The potentials can be made very general,
accepting free parameters and target lists, but it is also possible to hard code the parameters.
The latter tends to be easier since handling arbitrary lists of parameters becomes easily a code-intensive 
task.




.. file:compoundpotential class - autogenerated

.. _compoundpotential class - autogenerated:




List of methods
---------------

    
Below is a list of methods in :class:`~pysic.interactions.compound.CompoundPotential`, grouped according to
the type of functionality. Note that most methods are inherited from the :class:`~pysic.interactions.local.Potential` class.



Interaction handling
____________________

- :meth:`~pysic.interactions.compound.CompoundPotential.build`
- :meth:`pysic.interactions.local.Potential.get_cutoff`
- :meth:`pysic.interactions.local.Potential.get_cutoff_margin`
- :meth:`~pysic.interactions.compound.CompoundPotential.get_elements`
- :meth:`~pysic.interactions.compound.CompoundPotential.get_number_of_parameters`
- :meth:`pysic.interactions.local.Potential.get_parameter_names`
- :meth:`pysic.interactions.local.Potential.get_parameter_value`
- :meth:`pysic.interactions.local.Potential.get_parameter_values`
- :meth:`~pysic.interactions.compound.CompoundPotential.set_potential_type` (meant for internal use)
- :meth:`pysic.interactions.local.Potential.get_soft_cutoff`
- :meth:`~pysic.interactions.compound.CompoundPotential.remove`
- :meth:`pysic.interactions.local.Potential.set_cutoff`
- :meth:`pysic.interactions.local.Potential.set_cutoff_margin`
- :meth:`pysic.interactions.local.Potential.set_parameter_value`
- :meth:`pysic.interactions.local.Potential.set_parameter_values`
- :meth:`pysic.interactions.local.Potential.set_soft_cutoff`

Coordinator handling
_____________________

- :meth:`~pysic.interactions.local.Potential.get_coordinator`
- :meth:`~pysic.interactions.local.Potential.set_coordinator`

Target handling
__________________

- :meth:`~pysic.interactions.local.Potential.accepts_target_list`
- :meth:`~pysic.interactions.local.Potential.add_indices`
- :meth:`~pysic.interactions.local.Potential.add_symbols`
- :meth:`~pysic.interactions.local.Potential.add_tags`
- :meth:`~pysic.interactions.local.Potential.get_different_indices`
- :meth:`~pysic.interactions.local.Potential.get_different_symbols`
- :meth:`~pysic.interactions.local.Potential.get_different_tags`
- :meth:`~pysic.interactions.local.Potential.get_indices`
- :meth:`~pysic.interactions.local.Potential.get_number_of_targets`
- :meth:`~pysic.interactions.local.Potential.get_symbols`
- :meth:`~pysic.interactions.local.Potential.get_tags`
- :meth:`~pysic.interactions.local.Potential.set_indices`
- :meth:`~pysic.interactions.local.Potential.set_symbols`
- :meth:`~pysic.interactions.local.Potential.set_tags`


Description
___________

- :meth:`~pysic.interactions.compound.CompoundPotential.describe`
- :meth:`~pysic.interactions.local.Potential.is_multiplier` (meant for internal use)
- :meth:`~pysic.interactions.local.Potential.get_potentials` (meant for internal use)







Full documentation of the CompoundPotential class
-------------------------------------------------

.. currentmodule:: pysic.interactions.compound
.. autoclass:: CompoundPotential
   :members:
   :undoc-members:

