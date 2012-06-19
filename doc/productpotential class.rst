.. file:productpotential class

.. _productpotential class:



.. file:productpotential class - description

.. _productpotential class - description:



============================
ProductPotential class
============================

This class defines an interaction as a product of :class:`~pysic.interactions.local.Potential`\ s.



.. file:multiplying potentials

.. _multiplying potentials:



Adding and multiplying potentials
-----------------------------------

The ordinary :class:`~pysic.interactions.local.Potential` defines a local interaction for :math:`n` bodies (:math:`i,j,\ldots`), 
so that the total potential energy is

.. math::

  V = \sum_{(i,j,\ldots)} v_{i,j,\ldots}.

Defining several such potentials (:math:`p`) for the same set of atoms simply adds them together

.. math::

  V = \sum_{(i,j,\ldots)} \sum_p  v^p_{i,j,\ldots}.

In code, this could be done for instance this way::

  pot1 = pysic.Potential(...)
  pot2 = pysic.Potential(...)
  calc = pysic.Pysic()
  calc.set_potentials( [pot1, pot2] )

However, you may also wish to *multiply* the potentials to obtain

.. math::

  V = \sum_{(i,j,\ldots)} \prod_p v^p_{i,j,\ldots}.

This can be done by wrapping the potentials to be multiplied in a :class:`~pysic.interactions.local.ProductPotential`::

  pot1 = pysic.Potential(...)
  pot2 = pysic.Potential(...)
  prod = pysic.ProductPotential( [pot1, pot2] )
  calc = pysic.Pysic()
  calc.set_potentials( prod )

An immediate benefit of this functionality is the possibility to construct complicated potentials from simple building blocks
without having to hardcode all the different variants.
  

.. file:parameterization

.. _parameterization:



Parameterization
-----------------------------------

The :class:`~pysic.interactions.local.ProductPotential` class has mostly a similar interface as the elemental :class:`~pysic.interactions.local.Potential` class.
However, the object does not store parameters itself. It merely wraps a group of potentials and these
potentials are used for storing the parameters. More precisesly, the *first* 
:class:`~pysic.interactions.local.Potential` object in the list contained by the 
:class:`~pysic.interactions.local.ProductPotential` defines all the general parameters of the potential.

Since the product potential is defined as

.. math::

  V = \sum_{(i,j,\ldots)} \prod_p v^p_{i,j,\ldots},

the potentials :math:`v^p_{i,j,\ldots}` being multiplied need to be defined for the same group of atoms -
both the number of atoms and their types must match. Also the cutoffs should be equal, since if one :math:`v^p = 0`,
the whole product is. Because of this, the targets and cutoffs defined by the *leading* potential 
are used for the entire product::

  # Here, pot1 is the leading potential
  prod = pysic.ProductPotential( [pot1, pot2, pot3] )
  
In fact, since the cutoff and targets of the potentials forming the product besides 
the leading one are ignored, they need not be defined at all. 
Also possible :class:`~pysic.interactions.bondorder.Coordinator` objects defining bond order factors need to
be attached to the leading potential in order to have an effect.

As an example, the following code::

  pot1 = pysic.Potential(..., symbols=[['H','O']])
  pot2 = pysic.Potential(..., symbols=[['H','H']])
  prod1 = pysic.ProductPotential( [pot1, pot2] )
  prod2 = pysic.ProductPotential( [pot2, pot1] )

defines two products, ``prod1`` and ``prod2``, which are both a product of the potentials ``pot1`` and ``pot2``. 
However, the former, ``prod1``, targets H-O pairs (``pot1`` is the leading potential) while the latter,
``prod2``, affects H-H pairs (``pot2`` is the leading potential).


Furthermore, the :class:`~pysic.interactions.local.ProductPotential` defines methods for handling the general parameters
similarly to the :class:`~pysic.interactions.local.Potential` class, and these methods directly affect the leading potential.
That is, methods of the type::

 prod = pysic.ProductPotential( [pot1, pot2] )
 # method call on prod
 prod.set_cutoff(8.0)
 # is equivalent to a method call on the leading potential in prod
 prod.get_potentials()[0].set_cutoff(8.0)
 
are equivalent.

To avoid changing the original potentials, the objects are copied when passed to the product. 
However, this also means that you are unable to modify the product by changing the component potentials after having wrapped them::

 prod = pysic.ProductPotential( [pot1, pot2] )
 # method call on prod
 prod.set_cutoff(8.0)
 # is not the same as a method call on the original pot1
 pot1.set_cutoff(8.0)

For instance::

  calc = pysic.Pysic()
  pot1 = pysic.Potential(..., symbols=[['H','O']])   # Define a potential for H-O
  calc.add_potential( pot1 )                         # Add pot1 to the calculator
  pot2 = pysic.Potential(...)                        # Define another potential
  prod = pysic.ProductPotential( [pot1, pot2] )      # Multiply pot1 and pot2 to make prod
  prod.set_symbols( [['H','H']] )                    # Set prod to target H-H
  calc.add_potential( prod )                         # Add the prod to the calculator

creates two potentials, ``pot1`` and ``prod`` affecting H-O and H-H pairs.
Although ``prod`` is constructed from ``pot1`` and ``pot2``, is actually contains only
copies of these potentials.
If ``prod`` contained the original potential objects ``prod.set_symbols( [['H','H']] )`` 
would affect ``pot1``, the leading potential, and change ``symbols``
in ``pot1`` as well having both potentials end up targeting H-H pairs.
This kind of behaviour could lead to one unintentional changes the properties of other potentials

On the other hand, this example::

  calc = pysic.Pysic()
  pot1 = pysic.Potential(..., symbols=[['H','O']])   # Define a potential for H-O
  calc.add_potential( pot1 )                         # Add pot1 to the calculator
  pot2 = pysic.Potential(...)                        # Define another potential
  prod = pysic.ProductPotential( [pot1, pot2] )      # Multiply pot1 and pot2 to make prod
  pot1.set_symbols( [['H','H']] )                    # Set pot1 to target H-H
  calc.add_potential( prod )                         # Add the prod to the calculator

Leads to ``pot1`` and ``prod`` affecting H-H and H-O pairs, respectively. This is because 
``pot1.set_symbols( [['H','H']] )`` only changes the original potential, not the copy stored in ``prod``.


The physical parameters of the potential components are naturally defined separately for each potential. 
Therefore the :class:`~pysic.interactions.local.ProductPotential` has no methods for accessing these parameters of its constituents.
Instead, one has to directly access the :class:`~pysic.interactions.local.Potential` objects themselves::

  pot1 = pysic.Potential(...)                            # Define a potential
  pot2 = pysic.Potential(...)                            # Define another potential
  prod = pysic.ProductPotential( [pot1, pot2] )          # Multiply pot1 and pot2 to make prod
  prod.get_potentials()[1].set_parameter_value('a', 1.0) # set parameter a to 1.0
  pot2.set_parameter_value('a' 2.0)                      # set parameter a to 2.0
  
Above, the second potential in ``prod`` gets the parameter value ``a = 1.0``. 
The last command changes the parameter value in ``pot2``, but the change does not propagate
to the copy stored in ``prod``.


.. file:productpotential class - autogenerated

.. _productpotential class - autogenerated:




List of methods
---------------

    
Below is a list of methods in :class:`~pysic.interactions.local.ProductPotential`, grouped according to
the type of functionality.


Interaction handling
____________________

- :meth:`~pysic.interactions.local.ProductPotential.get_cutoff`
- :meth:`~pysic.interactions.local.ProductPotential.get_cutoff_margin`
- :meth:`~pysic.interactions.local.ProductPotential.get_soft_cutoff`
- :meth:`~pysic.interactions.local.ProductPotential.set_cutoff`
- :meth:`~pysic.interactions.local.ProductPotential.set_cutoff_margin`
- :meth:`~pysic.interactions.local.ProductPotential.set_soft_cutoff`

Coordinator handling
_____________________

- :meth:`~pysic.interactions.local.ProductPotential.get_coordinator`
- :meth:`~pysic.interactions.local.ProductPotential.set_coordinator`

Target handling
__________________

- :meth:`~pysic.interactions.local.ProductPotential.accepts_target_list`
- :meth:`~pysic.interactions.local.ProductPotential.add_indices`
- :meth:`~pysic.interactions.local.ProductPotential.add_symbols`
- :meth:`~pysic.interactions.local.ProductPotential.add_tags`
- :meth:`~pysic.interactions.local.ProductPotential.get_different_indices`
- :meth:`~pysic.interactions.local.ProductPotential.get_different_symbols`
- :meth:`~pysic.interactions.local.ProductPotential.get_different_tags`
- :meth:`~pysic.interactions.local.ProductPotential.get_indices`
- :meth:`~pysic.interactions.local.ProductPotential.get_number_of_targets`
- :meth:`~pysic.interactions.local.ProductPotential.get_symbols`
- :meth:`~pysic.interactions.local.ProductPotential.get_tags`
- :meth:`~pysic.interactions.local.ProductPotential.set_indices`
- :meth:`~pysic.interactions.local.ProductPotential.set_symbols`
- :meth:`~pysic.interactions.local.ProductPotential.set_tags`


Description
___________

- :meth:`~pysic.interactions.local.ProductPotential.is_multiplier` (meant for internal use)
- :meth:`~pysic.interactions.local.ProductPotential.get_potentials` (meant for internal use)


Full documentation of the ProductPotential class
-------------------------------------------------

.. currentmodule:: pysic.interactions.local
.. autoclass:: ProductPotential
   :members:
   :undoc-members:

