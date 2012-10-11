.. file:bondorderparameters class

.. _bondorderparameters class:



.. file:bondorderparameters class - description

.. _bondorderparameters class - description:





============================
BondOrderParameters class
============================

This class defines a set of parameters for a bond order factor, to be used
in conjunction with the :class:`~pysic.interactions.bondorder.Coordinator` class. 

Similarly to the potentials, the available types of bond order factors 
are always inquired from the
Fortran core to ensure that any changes made to the core are
automatically reflected in the Python interface.

The same utility functions in :mod:`~pysic` for inquiring
keywords and other data needed for creating the potentials also work for
fetching information on bond order factors, if applicable.
The functions check automatically if the inquired name matches a potential
or a bond order factor and gather the correct type of information based on this.

For example:

- Inquire the names of available bond order factors: :meth:`~pysic.list_valid_bond_order_factors`
- Inquire the names of parameters for a bond order factors: :meth:`~pysic.names_of_parameters`

.. file:bond order cutoffs

.. _bond order cutoffs:





Bond order cutoffs
---------------------------

Atomic coordination is an example of a simple bond order factor.
It is calculated by checking all atom pairs
and counting which ones are "close" to each others. Close naturally means
closer than some predefined cutoff distance.
However, in order to make the coordination a continuous and 
differentiable function, a continuous cutoff has to be applied.
This is done similarly to the smooth cutoffs used in :class:`~pysic.interactions.local.Potential`
by defining a proximity function which is 1 for small separations and
0 for large distances.


.. math::

   f(r) = \begin{cases} 1, & r < r_\mathrm{soft} \\ \frac{1}{2}\left(1+\cos \pi\frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}}\right), & r_\mathrm{soft} < r < r_\mathrm{hard} \\ 0, & r > r_\mathrm{hard} \end{cases}.

.. only:: html

 .. plot::
   
   import matplotlib.pyplot as plt
   from math import cos, pi
   x = []; y = []
   start = 0.0
   end = 2.0
   steps = 100
   margin = 0.4
   dx = (end-start)/steps
   for i in range(steps+1):
       xval = start + i*dx
       x.append( xval )
       if xval < 1.0 - margin/2:
           y.append( 1.0 )
       elif xval > 1.0 + margin/2:
           y.append( 0.0 )
       else:
           y.append( 0.5*(1.0 + cos(pi*(xval - 1.0 + margin/2)/margin) ) )
   plt.plot(x,y)
   plt.plot(0)
   plt.xlim(0.0,2.0)
   plt.ylim(-0.1,1.1)
   plt.title(r'Proximity function: $r_\mathrm{soft} = 0.8$, $r_\mathrm{hard} = 1.2$')
   plt.show()

Since bond order factors such as 
atomic coordination need not decay as a function of distance, one must always
define a margin for continuous cutoff in bond order factors.

.. file:defining parameters for bond order factors

.. _defining parameters for bond order factors:





Defining parameters
-------------------

A :class:`~pysic.interactions.bondorder.BondOrderParameters` instance defines the type of the bond order factor,
the cutoffs, and parameters for one set of elements. The parameters are formally split
according to the number of atoms they act on. So, an n-body factor can have parameters
which are applied for 1-body, 2-body, etc. terms. Bond order factors are applied and 
parameterized by atomic element types (chemical symbols). An n-body factor must always have
one or several sets of n symbols to designate the atoms it affects. So, 2- and 3-body factors
could accept for instance the following lists of symbols, respectively::

      >>> two_body_targets = [['H', 'H'], ['H', 'O'], ['O', 'O']]
      >>> three_body_targets = [['Si', 'O', 'H']]

As a rule of thumb, if an n-body bond order factor incorporates parameters for m bodies 
(:math:`m \le n`), then these parameters are targeted at the first m symbols of the target list.
For instance, if a 3-body factor has the targets of the above example (three_body_targets)
and it contains 1- and 2-body parameters, then the 1-body parameters are targeted at Si atoms
and the 2-body parameters at Si-O bonds.

As an example, the Tersoff bond order factor


.. math::

    b_i = \left[ 1 + \left( \beta_i \sum_{j \ne i} \sum_{k \ne i,j} \xi_{ijk} g_{ijk}  \right)^{\eta_i} \right]^{-\frac{1}{2 \eta_i}}

.. math::

    \xi_{ijk} = f(r_{ij}) f(r_{ik}) \exp\left[a_{ij}^{\mu_i} (r_{ij} - r_{ik})^{\mu_i} \right]

.. math::

    g_{ijk} = 1 + \frac{c_{ij}^2}{d_{ij}^2} - \frac{c_{ij}^2}{d_{ij}^2 + (h_{ij} - \cos \theta_{ijk})^2}

is a three-body factor (it includes terms depending on atom triplets :matH:`(i, j, k)`) and therefore requires a set of three
elements as its target. It incorporates three single body and four two body parameters.
Such a bond order factor could be created with the following command::

     >>> bonds = pysic.BondOrderParameters('tersoff', cutoff=3.2, cutoff_margin=0.4,
     ...       	 		           symbols=[['Si', 'Si', 'Si']],
     ...				   parameters=[[beta, eta, mu],
     ...				               [a, c, d, h],
     ...					       [ ]])

or alternatively in pieces by a series of commands::

   >>> bonds = pysic.BondOrderParameters('tersoff')
   >>> bonds.set_cutoff(3.2)
   >>> bonds.set_cutoff_margin(0.4)
   >>> bonds.set_symbols([['Si', 'Si', 'Si']])
   >>> bonds.set_parameter_value('beta', beta)
   >>> bonds.set_parameter_value('eta', eta)
   >>> bonds.set_parameter_value('mu', mu)
   >>> bonds.set_parameter_value('a', a)
   >>> bonds.set_parameter_value('c', c)
   >>> bonds.set_parameter_value('d', d)
   >>> bonds.set_parameter_value('h', h)

To be used in calculations, this is then passed on to a :class:`~pysic.interactions.bondorder.Coordinator`, 
:class:`~pysic.interactions.local.Potential`, and :class:`~pysic.calculator.Pysic` with::

    >>> crd = pysic.Coordinator( bonds )
    >>> pot = pysic.Potential( ... , coordinator=crd )
    >>> cal = pysic.Pysic( potentials=pot )

The above example creates a bond order factor which is applied to all Si triplets (symbols=[['Si','Si','Si']]). The command also assigns 1-body parameters :math:`\beta`, :math:`\eta`, and :math:`\mu`, and 2-body parameters 
:math:`a`, :math:`c`, :math:`d`, and :math:`h`. If there are other elements in the system besides silicon, they will be completely ignored: The bond order factors are calculated as if the other elements do not exist. If one wishes to include, say, Si-O bonds in the bond order factor calculation, the list of symbols needs to be expanded by::

   >>> bonds.add_symbols([['Si', 'Si', 'O'],
   ...			  ['Si', 'O', 'Si'],
   ...			  ['Si', 'O', 'O']])
   >>> bonds.get_symbols()
   [['Si', 'Si', 'Si'], ['Si', 'Si', 'O'], ['Si', 'O', 'Si'], ['Si', 'O', 'O']]

The format of the symbol list is as follows. In each triplet, the first symbol determines the element on which the factor is calculated. Since above the first symbol of each triplet is 'Si', the factor
will only be applied on Si atoms. The other symbols define the other elements in the triplets which are taken in to account. The second and third symbols are not, however, symmetric. As the bond order factor is defined using 2-body parameters (:math:`a_ij` etc.), the first two symbols determine the elements of those two atoms (atoms :math:`i` and :math:`j`). The third symbol determines the element of the third atom (atom k) of the triplet. I.e., in the example above, Si-O bond parameters are included with::

   >>> [['Si', 'O', 'Si'], ['Si', 'O', 'O']]

where the first works for triplets O-Si-Si and the second for O-Si-O.
However, one should note especially that a triplet A-B-C is only taken in to account if both bonds (A-B) and (B-C) have parameters associated with them. Therefore, ['Si', 'O', 'O'] is enough to fully define O-Si-O bond triplets, but to fully describe Si-Si-O, one also has to define the Si-Si bond with O as the third partner, which is given by::

   >>> [['Si', 'Si', 'O']]

Instead of giving a list of symbols to a single :class:`~pysic.interactions.bondorder.BondOrderParameters`, one can define many instances with different symbols and different parameters, and feed a list of these to a :class:`~pysic.interactions.bondorder.Coordinator` object.::

     >>> bond_siosi = pysic.BondOrderParameters('tersoff', cutoff=3.2, cutoff_margin=0.4,
     ...       	      			        symbols=[['Si', 'O', 'Si']],
     ...			                parameters=[[beta_si, eta_si, mu_si],
     ...				                    [a_sio, c_sio, d_sio, h_sio],
     ...					            [ ]])
     >>> bond_sisio = pysic.BondOrderParameters('tersoff', cutoff=3.5, cutoff_margin=0.5,
     ...       	                		symbols=[['Si', 'Si', 'O']],
     ...			                parameters=[[beta_si, eta_si, mu_si],
     ...				                    [a_sisi, c_sisi, d_sisi, h_sisi],
     ...					            []])
     >>> bond_list = [bond_siosi,bond_sisio]
     >>> crd = pysic.Coordinator( bond_list )

The above example would assign the parameter values

.. math::
      
   \beta_\mathrm{Si} = \mathtt{beta\_si} \\
   \eta_\mathrm{Si} = \mathtt{eta\_si} \\
   \mu_\mathrm{Si} = \mathtt{mu\_si} \\
   a_\mathrm{Si-O} = \mathtt{a\_sio} \\
   c_\mathrm{Si-O} = \mathtt{c\_sio} \\
   d_\mathrm{Si-O} = \mathtt{d\_sio} \\
   h_\mathrm{Si-O} = \mathtt{h\_sio} \\
   a_\mathrm{Si-Si} = \mathtt{a\_sisi} \\
   c_\mathrm{Si-Si} = \mathtt{c\_sisi} \\
   d_\mathrm{Si-Si} = \mathtt{d\_sisi} \\
   h_\mathrm{Si-Si} = \mathtt{h\_sisi}

This gives the user the possibility to precisely control the parameters, including cutoffs, for different elements.

Note that the beta, eta, and mu parameters are the same for both :class:`~pysic.interactions.bondorder.BondOrderParameters` objects defined in the above example. They could be different in principle, but when the factors are calculated, the 1-body parameters are taken from the first object in the list of bonds (bond_list) for which the first element is of the correct type. Because of this, the 1-body parameters in bond_sisio are in fact ignored. This feature can be exploited for mixing different types of bond order factors, as explained below.

For three different elements, say C, O, and H, the possible triplets are::

  >>> [ ['H', 'H', 'H'],  # H-H bond in an H-H-H triplet
  ...   ['H', 'H', 'C'],  # H-H bond in an H-H-C triplet
  ...   ['H', 'H', 'O'],  # H-H bond in an H-H-O triplet
  ...   ['H', 'O', 'H'],  # H-O bond in an H-H-O triplet
  ...   ['H', 'O', 'C'],  # H-O bond in an O-H-C triplet
  ...   ['H', 'O', 'O'],  # H-O bond in an O-H-O triplet
  ...   ['H', 'C', 'H'],  # etc.
  ...   ['H', 'C', 'C'],
  ...   ['H', 'C', 'O'],
  ...   ['O', 'H', 'H'],
  ...   ['O', 'H', 'C'],
  ...   ['O', 'H', 'O'],
  ...   ['O', 'O', 'H'],
  ...   ['O', 'O', 'C'],
  ...   ['O', 'O', 'O'],
  ...   ['O', 'C', 'H'],
  ...   ['O', 'C', 'C'],
  ...   ['O', 'C', 'O'],
  ...   ['C', 'H', 'H'],
  ...   ['C', 'H', 'C'],
  ...   ['C', 'H', 'O'],
  ...   ['C', 'O', 'H'],
  ...   ['C', 'O', 'C'],
  ...   ['C', 'O', 'O'],
  ...   ['C', 'C', 'H'],
  ...   ['C', 'C', 'C'],
  ...   ['C', 'C', 'O'] ]

In principle, one can attach a different set of parameters to each of these.
Often though the parameters are mostly the same,
and writing these kinds of lists for all possible combinations is cumbersome. 
To help in generating the tables, the utility method 
:meth:`~pysic.utility.convenience.expand_symbols_table` can be used.
For instance, the full list of triplets above can be created with::

   >>> pysic.utility.convenience.expand_symbols_table([['C', 'O', 'H'], 
   ...                                     ['C', 'O', 'H'], 
   ...                                     ['C', 'O', 'H']])

.. file:list of bond order factors

.. _list of bond order factors:




List of currently available bond order factors
----------------------------------------------

Below is a list of bond order factors currently implemented.

- :ref:`coordination scaling function`
- :ref:`square root scaling function`
- :ref:`tabulated scaling function`
- :ref:`coordination bond order factor`
- :ref:`power decay bond order factor`
- :ref:`tabulated bond order factor`
- :ref:`tersoff bond order factor`


.. file:coordination scaling function

.. _coordination scaling function:





Coordination scaling function
_____________________________

1-body bond order factor defined as

.. math::

   b_i(\Sigma_i) = \varepsilon_i \frac{\Delta \Sigma_i}{1+\exp(\gamma_i \Delta \Sigma_i)}\\
   \Delta \Sigma_i = C_i (\Sigma_i - N_i).

where :math:`\Sigma_i` is the bond order sum.

In other words, this factor only overrides the scaling function of another
bond order factor when mixed. Especially, it is zero if not paired with other bond order factors.

Keywords::

    >>> names_of_parameters('c_scale')
    [['epsilon', 'N', 'C', 'gamma']]

.. file:square root scaling function

.. _square root scaling function:





Square root scaling function
_____________________________

1-body bond order factor defined as

.. math::

   b_i(\Sigma_i) = \varepsilon_i \sqrt{ \Sigma_i }

where :math:`\Sigma_i` is the bond order sum and :math:`\varepsilon` is a scaling constant.

Keywords::

    >>> names_of_parameters('sqrt_scale')
    [['epsilon']]

.. file:tabulated scaling function

.. _tabulated scaling function:





Tabulated scaling function
_____________________________

1-body bond order factor of the type

.. math::

   b_i(\Sigma_i) = s_i(\Sigma_i),

where :math:`s_i(\Sigma)` is a tabulated function. The tabulation works similarly to the :ref:`tabulated potential`.

Keywords::

    >>> names_of_parameters('table_scale')
    [[id, range, scale]]


.. file:coordination bond order factor

.. _coordination bond order factor:





Coordination bond order factor
______________________________

2-body bond order factor defined as

.. math::

   b_i = \sum_{j \ne i} f(r_{ij}).

The coordination of an atom is simply the sum of the proximity functions.
This is a parameterless (besides cutoffs) 2-body bond order factor.


Keywords::

    >>> names_of_parameters('neighbors')
    [[], []]

.. file:power decay bond order factor

.. _power decay bond order factor:





Power decay bond order factor
_____________________________

2-body bond order factor defined as

.. math::

   b_i = \sum_{j \ne i} f(r_{ij})  \left(\frac{a_{ij}}{r_{ij}}\right)^{n_{ij}}.

This is a density-like bond factor, where the contributions of atomic pairs decay with 
interatomic distance according to a power law. In form, it is similar to the :ref:`power decay potential` potential.

Keywords::

    >>> names_of_parameters('power_bond')
    [[], [a, n]]

.. file:tabulated bond order factor

.. _tabulated bond order factor:





Tabulated bond order factor
_____________________________

2-body bond order factor of the type

.. math::

   b_i = \sum_{j \ne i} f(r_{ij}) t(r_{ij}),

where :math:`t(r)` is a tabulated function. The tabulation works similarly to the :ref:`tabulated potential`.

Similarly, to tabulate bond scaling, use the :ref:`tabulated scaling function`.

Keywords::

    >>> names_of_parameters('table_bond')
    [[], [id, range, scale]]


.. file:tersoff bond order factor

.. _tersoff bond order factor:




Tersoff bond order factor
_________________________

3-body bond order factor defined as

.. math::

    b_i = \left[ 1 + \left( \beta_i \sum_{j \ne i} \sum_{k \ne i,j} \xi_{ijk} g_{ijk}  \right)^{\eta_i} \right]^{-\frac{1}{2 \eta_i}}

.. math::

    \xi_{ijk} = f(r_{ij}) f(r_{ik}) \exp\left[a_{ij}^{\mu_i} (r_{ij} - r_{ik})^{\mu_i} \right]

.. math::

    g_{ijk} = 1 + \frac{c_{ij}^2}{d_{ij}^2} - \frac{c_{ij}^2}{d_{ij}^2 + (h_{ij} - \cos \theta_{ijk})^2}

where r and theta are distances and angles between the atoms.
This rather complicated bond factor takes also into account the directionality of bonds in its angle dependency.

Keywords::

    >>> names_of_parameters('tersoff')
    [['beta', 'eta', 'mu'], ['a', 'c', 'd', 'h'], []]

.. file:bondorderparameters class - autogenerated

.. _bondorderparameters class - autogenerated:





List of methods
---------------

Below is a list of methods in :class:`~pysic.interactions.bondorder.BondOrderParameters`, grouped according to
the type of functionality.


Parameter handling
__________________

- :meth:`~pysic.interactions.bondorder.BondOrderParameters.accepts_parameters`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_bond_order_type`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_cutoff`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_cutoff_margin`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_number_of_parameters`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_parameter_names`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_parameter_value`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_parameter_values`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_parameters_as_list`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_soft_cutoff`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.set_cutoff`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.set_cutoff_margin`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.set_parameter_value`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.set_parameter_values`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.set_parameters`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.set_soft_cutoff`

Target handling
_______________

- :meth:`~pysic.interactions.bondorder.BondOrderParameters.accepts_target_list`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.add_symbols`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_different_symbols`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_number_of_targets`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.get_symbols`
- :meth:`~pysic.interactions.bondorder.BondOrderParameters.set_symbols`


Full documentation of the BondOrderParameters class
---------------------------------------------------

.. currentmodule:: pysic.interactions.bondorder
.. autoclass:: BondOrderParameters
   :members:
   :undoc-members:

