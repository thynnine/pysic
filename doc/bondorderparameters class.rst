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

.. file:atomic and pairwise factors

.. _atomic and pairwise factors:



Atomic and pairwise factors
---------------------------------------

Two types of factors can be defined: atomic or pairwise (per bond) factors. 
Let us first give formal definitions for these types and then discuss the differences in their use and behavior.

Atomic factors
__________________

Atomic factors are of the form

.. math::

  b_i = s_i( \sum_{j,\ldots} c_{ij\ldots}),

where the local factors :math:`c_{ij\ldots}` may have 2, 3, or more indices depending on how many bodies
affect the factor. When this kind of a factor is applied on an :math:`n`-body potential, :math:`v_{ij\ldots}`, 
an :math:`n`-body factor is created as the average of the atomic factors

.. math::

 b_{ij\ldots} = \frac{1}{n}(b_i + b_j + \ldots).

The resulting potential is then

.. math::

 U = \sum_{i,j,\ldots} b_{ij\ldots} v_{ij\ldots},

where the summation goes over all :math:`n`-chains.

Pairwise factors
___________________

Pairwise factors, on the other hand are only defined for pair potentials :math:`v_{ij}`. These factors scale
the interaction by

.. math::

 U = \sum_{i,j} b_{ij} v_{ij},

where the summation goes over all pairs :math:`(i,j)`. Note that this form requires :math:`b_{ij}` to be
symmetric with respect to :math:`i` and :math:`j`.
It would also be possible to define the factor through

.. math::

 U = \frac{1}{2} \sum_{i \ne j} \tilde{b}_{ij} v_{ij},

where the sum goes over all indices --- and thus over all pairs twice. Using the notation above, this is equivalent to

.. math::

 U = \sum_{i,j} \frac{1}{2}(\tilde{b}_{ij}+\tilde{b}_{ji}) v_{ij}.

Clearly, :math:`b_{ij} = \frac{1}{2}(\tilde{b}_{ij}+\tilde{b}_{ji})`.
The factor :math:`\tilde{b}_{ij}` defined in this manner need not be symmetric, since the summation automatically
leads to the symmetric form. Because of this, to avoid the need to force symmetry on :math:`b_{ij}`, 
Pysic calculates the pairwise factors using :math:`b_{ij} = \frac{1}{2}(\tilde{b}_{ij}+\tilde{b}_{ji})`. 
Therefore it is only necessary to implement the non-symmetric factors :math:`\tilde{b}_{ij}`.

It is important to note that any scaling functions are applied on the factors :math:`\tilde{b}_{ij}`,
not on :math:`b_{ij}`. The difference can be seen with a simple example. Let our factor be :math:`\tilde{b}_{ij} = \sqrt{\sum_{k} c_{ijk}}`. Then, :math:`b_{ij} = \frac{1}{2} (\sqrt{\sum_{k} c_{ijk}} + \sqrt{\sum_{k} c_{jik}}) \ne \sqrt{\frac{1}{2}(\sum_{k} c_{ijk} + c_{jik})}`. Also note that :ref:`mixing bond order types` is not possible between atomic and pairwise factors.

Only use a pairwise factor with a pair potential!
If a pairwise factor is applied on an :math:`n`-body potential where :math:`n \ne 2`, it will automatically be zero. However, applying such a factor will result in the potential being multiplied by the factor, and so the potential becomes zero also.

Pairwise and atomic factors can be used on one and the same pair potential, in which case they are simply added together:

.. math::
  
  b_{ij} = \frac{1}{2}(b_i + b_j + \tilde{b}_{ij} + \tilde{b}_{ji}).

The usefulness of this is probably limited --- this behavior is adopted to avoid conflicts.

.. file:defining parameters for bond order factors

.. _defining parameters for bond order factors:





Defining parameters
-------------------

A :class:`~pysic.interactions.bondorder.BondOrderParameters` instance defines the type of the bond order factor,
the cutoffs, and parameters for one set of elements. The parameters are formally split to scaling function parameters
and local summation parameters, where for a factor of type

.. math::

  b_i = s_i(\sum_j c_{ij})

:math:`s_i` is the scaling function and :math:`\sum_j c_{ij}` is the local sum (similarly for many-body factors). 

This division is mostly cosmetic, though, and the parameters could just as well be defined
as a single list.

Bond order factors are applied and 
parameterized by atomic element types (chemical symbols). An :math:`n`-body factor must always have
one or several sets of :math:`n` symbols to designate the atoms it affects. So, 2- and 3-body factors
could accept for instance the following lists of symbols, respectively::

      >>> two_body_targets = [['H', 'H'], ['H', 'O'], ['O', 'O']]
      >>> three_body_targets = [['Si', 'O', 'H']]

Atomic factors
__________________

As an example, the :ref:`power decay bond order factor`

.. math::

   b_i = \sum_{j \ne i} f(r_{ij})  \left(\frac{a}{r_{ij}}\right)^{n}
 
is a two-body factor and therefore requires two elements as its target. It includes no scaling and two local summation factors.

Such a bond order factor could be created with the following command::

     >>> bonds = pysic.BondOrderParameters('power_bond', cutoff=3.2, cutoff_margin=0.4,
     ...                                   symbols=[['Si', 'Si']],
     ...                                   parameters=[[],[a, n]])

or alternatively in pieces by a series of commands::

   >>> bonds = pysic.BondOrderParameters('power_bond')
   >>> bonds.set_cutoff(3.2)
   >>> bonds.set_cutoff_margin(0.4)
   >>> bonds.set_symbols([['Si', 'Si']])
   >>> bonds.set_parameter_value('a', a)
   >>> bonds.set_parameter_value('n', n)

To be used in calculations, this is then passed on to a :class:`~pysic.interactions.bondorder.Coordinator`, 
:class:`~pysic.interactions.local.Potential`, and :class:`~pysic.calculator.Pysic` with::

    >>> crd = pysic.Coordinator( bonds )
    >>> pot = pysic.Potential( ... , coordinator=crd )
    >>> cal = pysic.Pysic( potentials=pot )

The above factor will only consider Si-Si pairs in the local summation :math:`b_i = \sum_j c_{ij}`, i.e., the atoms :math:`i` and :math:`j` must both be silicons. 
If also, say Si-O pairs should be taken into account,
the list needs to be expanded:

    >>> bonds.add_symbols([['Si','O']])
    >>> bonds.get_symbols()
    [['Si', 'Si'], ['Si', 'O']]

This will apply the factor on Si atoms so that the summation :math:`b_i = \sum_j c_{ij}` goes over Si-Si and Si-O pairs, i.e., atom :math:`i` is Si but atom :math:`j` may be either Si or O.

If you want to define a bond factor also for oxygens, this can be done separately with for instance::

     >>> bonds2 = pysic.BondOrderParameters('power_bond', cutoff=3.2, cutoff_margin=0.4,
     ...                                    symbols=[['O', 'O'], ['O', 'Si']],
     ...                                    parameters=[[],[a, n]])
     >>> crd2 = pysic.Coordinator( bonds2 )
     >>> pot2 = pysic.Potential( ... , coordinator=crd2 )
     >>> cal.add_potential( pot2 )

Here one needs to be careful. If you apply a bond order factor on a potential, say a Si-O pair potential, the factor is applied on both Si and O atoms. However, if no parameters are applied for O, the factor for it is zero. That is, in the case of a Si-O pair (Si is :math:`i` and O is :math:`j`) the bond factor :math:`b_{ij} = \frac{1}{2}(b_i + b_j)=\frac{1}{2}b_i`, which may be not the intended result.

If you want to have different local parameters for the different pairs O-O and O-Si, you must define two 
:class:`~pysic.interactions.bondorder.BondOrderParameters` objects and wrap them in a :class:`~pysic.interactions.bondorder.Coordinator`::

   >>> bonds_oo  = pysic.BondOrderParameters('power_bond', cutoff=3.2, cutoff_margin=0.4,
   ...                                       symbols=[['O', 'O']],
   ...                                       parameters=[[],[a_oo, n_oo]])
   >>> bonds_osi = pysic.BondOrderParameters('power_bond', cutoff=3.2, cutoff_margin=0.4,
   ...                                       symbols=[['O', 'Si']],
   ...                                       parameters=[[],[a_osi, n_osi]])
   >>> crd3 = pysic.Coordinator( [bonds_oo, bonds_osi] )
   >>> pot3 = pysic.Potential( ... , coordinator=crd3 )
   >>> cal.set_potentials( [pot, pot3] )

That is, local summation is done using all the given parameters summed together.
If you apply a scaling factor on a bond order factor, however, it is applied only once. The scaling is determined by the first
:class:`~pysic.interactions.bondorder.BondOrderParameters` object in the :class:`~pysic.interactions.bondorder.Coordinator`
which defines a scaling function and whose first target is the atom for which the factor is being calculated. This can be used
for :ref:`mixing bond order types`.

As a rule of thumb, remember that one :class:`~pysic.interactions.local.Potential` can contain only one :class:`~pysic.interactions.bondorder.Coordinator`,
but a :class:`~pysic.interactions.bondorder.Coordinator` can contain many :class:`~pysic.interactions.bondorder.BondOrderParameters`. 
So if your bond order factor requires several sets of parameters due to the different element pairs, it is safest to define each set
of parameters using its own :class:`~pysic.interactions.bondorder.BondOrderParameters` object and wrap the parameters involved in one
local summation in a :class:`~pysic.interactions.bondorder.Coordinator`.

Pairwise factors
___________________

As an example, the :ref:`tersoff bond order factor`


.. math::

    \tilde{b}_{ij} = \left[ 1 + \left( \beta \sum_{k \ne i,j} \xi_{ijk} g_{ijk}  \right)^{\eta} \right]^{-\frac{1}{2 \eta}}

.. math::

    \xi_{ijk} = f(r_{ik}) \exp\left[a^{\mu} (r_{ij} - r_{ik})^{\mu} \right]

.. math::

    g_{ijk} = 1 + \frac{c^2}{d^2} - \frac{c^2}{d^2 + (h - \cos \theta_{ijk})^2}

is a three-body factor (it includes terms depending on atom triplets :math:`(i, j, k)`) and therefore requires a set of three
elements as its target. It incorporates two scaling and five local sum parameters.
Such a bond order factor could be created with the following command::

     >>> bonds = pysic.BondOrderParameters('tersoff', cutoff=3.2, cutoff_margin=0.4,
     ...       	 		           symbols=[['Si', 'Si', 'Si']],
     ...				   parameters=[[beta, eta],
     ...				               [a, c, d, h, mu]])

or alternatively in pieces by a series of commands::

   >>> bonds = pysic.BondOrderParameters('tersoff')
   >>> bonds.set_cutoff(3.2)
   >>> bonds.set_cutoff_margin(0.4)
   >>> bonds.set_symbols([['Si', 'Si', 'Si']])
   >>> bonds.set_parameter_value('beta', beta)
   >>> bonds.set_parameter_value('eta', eta)
   >>> bonds.set_parameter_value('a', a)
   >>> bonds.set_parameter_value('c', c)
   >>> bonds.set_parameter_value('d', d)
   >>> bonds.set_parameter_value('h', h)
   >>> bonds.set_parameter_value('mu', mu)


The above example creates a bond order factor which is applied to all Si triplets (symbols=[['Si','Si','Si']]). The command also assigns scaling parameters :math:`\beta`, :math:`\eta`, and :math:`\mu`, and local summation parameters 
:math:`a`, :math:`c`, :math:`d`, and :math:`h`. If there are other elements in the system besides silicon, they will be completely ignored: The bond order factors are calculated as if the other elements do not exist. If one wishes to include, say, Si-O bonds in the bond order factor calculation, the list of symbols needs to be expanded by::

   >>> bonds.add_symbols([['Si', 'Si', 'O'],
   ...			  ['Si', 'O', 'Si'],
   ...			  ['O', 'Si', 'Si'],
   ...			  ['Si', 'O', 'O'],
   ...			  ['O', 'Si', 'O']])
   >>> bonds.get_symbols()
   [['Si', 'Si', 'Si'], ['Si', 'Si', 'O'], ['Si', 'O', 'Si'], ['O', 'Si', 'Si'], ['Si', 'O', 'O'], ['O', 'Si', 'O']]

The format of the symbol list is as follows. In each triplet, the first two symbols determine the bond on which the factor is calculated (atoms :math:`i` and :math:`j`). 
(For atomic factors, the first symbol determines the element on which the factor is applied.)
The third symbol defines the other atoms in the triplets which are taken in to account (atom :math:`k`). 
That is, in the example above, Si-(Si-O) bond parameters are included with::

   >>> ['Si', 'O', 'Si']

O-(Si-O) with::

   >>> ['Si', 'O', 'O']

Si-(O-Si) with::

   >>> ['O', 'Si', 'Si']

and O-(O-Si) with::

   >>> ['O', 'Si', 'O']

The definition is complicated like this to enable the tuning of parameters of all the various bond combinations separately.

Instead of giving a list of symbols to a single :class:`~pysic.interactions.bondorder.BondOrderParameters`, one can define many instances with different symbols and different parameters, and feed a list of these to a :class:`~pysic.interactions.bondorder.Coordinator` object.::

     >>> bond_sioo = pysic.BondOrderParameters('tersoff', cutoff=3.2, cutoff_margin=0.4,
     ...       	      			        symbols=[['Si', 'O', 'O']],
     ...			                parameters=[[beta_si, eta_si],
     ...				                    [a_sio, c_sio, d_sio, h_sio, mu_si]])
     >>> bond_sisio = pysic.BondOrderParameters('tersoff', cutoff=3.5, cutoff_margin=0.5,
     ...       	                		symbols=[['Si', 'Si', 'O']],
     ...			                parameters=[[beta_si, eta_si],
     ...				                    [a_sisi, c_sisi, d_sisi, h_sisi, mu_si]])
     >>> bond_list = [bond_sioo,bond_sisio]
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

Note that the beta, eta, and mu parameters are the same for both :class:`~pysic.interactions.bondorder.BondOrderParameters` objects defined in the above example. They could be different in principle, but when the factors are calculated, the scaling parameters are taken from the first object in the list of bonds (`bond_list`) for which the first element is of the correct type. Because of this, the scaling parameters in `bond_sisio` are in fact ignored. This feature can be exploited for :ref:`mixing bond order types`.

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
   ...                                                 ['C', 'O', 'H'], 
   ...                                                 ['C', 'O', 'H']])

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
    [['epsilon', 'N', 'C', 'gamma'], []]

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
    [['epsilon'], []]

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
    [[id, range, scale], []]


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

Pairwise 3-body bond order factor defined as

.. math::

    \tilde{b}_{ij} = \left[ 1 + \left( \beta_{ij} \sum_{k \ne i,j} \xi_{ijk} g_{ijk}  \right)^{\eta_{ij}} \right]^{-\frac{1}{2 \eta_{ij}}}

.. math::

    \xi_{ijk} = f(r_{ik}) \exp\left[a_{ijk}^{\mu_{ijk}} (r_{ij} - r_{ik})^{\mu_{ijk}} \right]

.. math::

    g_{ijk} = 1 + \frac{c_{ijk}^2}{d_{ijk}^2} - \frac{c_{ijk}^2}{d_{ijk}^2 + (h_{ijk} - \cos \theta_{ijk})^2}

where r and theta are distances and angles between the atoms.
This rather complicated bond factor takes also into account the directionality of bonds in its angle dependency.

Keywords::

    >>> names_of_parameters('tersoff')
    [['beta', 'eta'], ['a', 'c', 'd', 'h', 'mu']]

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

