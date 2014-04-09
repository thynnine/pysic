.. file:coulombsummation class

.. _coulombsummation class:



.. file:coulombsummation class - description

.. _coulombsummation class - description:



.. module:: pysic.interactions.coulomb

======================
CoulombSummation class
======================

If a periodic system contains charges interacting via the :math:`\frac{1}{r}` Coulomb potential, direct summation of the interactions
    
    .. math::
       :label: direct_sum
    
       E = \sum_{(i,j)} \frac{1}{4\pi\epsilon_0}\frac{q_i q_j}{r_{ij}},
    
where the sum is over pairs of charges :math:`q_i, q_j` (charges of the entire system, not just the simulation cell) and the distance between the charges is :math:`r_{ij} = |\mathbf{r}_j - \mathbf{r}_i|`, does not work in general because the sum :eq:`direct_sum` converges very slowly (it actually converges only conditionally). Therefore truncating the sum may lead to severe errors. More advanced techniques must be used in order to accurately evaluate such sums.

This class represents the algorithms used for evaluating the :math:`1/r` sums. It wraps the summation parameters and activates the summation of Coulomb interactions. If an instance of :class:`~pysic.interactions.coulomb.CoulombSummation` is given to the :class:`~pysic.calculator.Pysic` calculator, Coulomb interactions between all charged atoms are automatically included in the calculations, regardless of possible :class:`~pysic.interactions.local.Potential` potentials the calculator may also contain. Otherwise the charges do not directly interact. This is due to two reasons: First, the direct Coulomb interaction is usually always required and it is convenient that it is easily enabled. Second, the specific potentials described by :class:`~pysic.interactions.local.Potential` are evaluated by direct summation and so the Coulomb summation is separate also on algorithm level in the core.


.. file:coulomb scaling

.. _coulomb scaling:



Charge scaling
--------------

Sometimes, you may want to scale the effective charges before calculating the Coulomb sum. Especially, you may want to exclude some atoms from the long range summation. This can be done by giving the :class:`~pysic.interactions.coulomb.CoulombSummation` a list of scaling values, one per atom. The actual charges of the atoms are then multiplied by the given scaling values before the Coulomb potential is calculated. If a scaling value is 0, the corresponding atom is always treated as if it had no charge. Note though, that  scaling with unequal scaling constants may lead to the cell being effectively charged.

Using the Python `map <http://docs.python.org/tutorial/datastructures.html#functional-programming-tools>`_ function is a convenient way to generate such atom-by-atom lists. For instance, if you want to generate a list by element::

  >>> atoms = ase.Atoms('H2O')  
  >>> def charge_by_elem(elem):
  ...     if elem == 'H':
  ...         return 0.1
  ...     elif elem == 'O':
  ...         return -0.2
  ...     else:
  ...         return 0.0
  ... 
  >>> system.set_initial_charges(map(charge_by_elem, system.get_chemical_symbols()))
  >>> charges
  [0.1, 0.1, -0.2]



.. file:coulomb parameters

.. _coulomb parameters:



Automatic parameters
------------------------

As the summation algorithms are parameter dependent, one should always check numeric convergence before real simulations.
As a first guess, the utility function :func:`pysic.interactions.coulomb.estimate_ewald_parameters` can be used for estimating the parameters of the Ewald method.


.. file:list of summation modes

.. _list of summation modes:




List of currently available summation algorithms
------------------------------------------------

Below is a list of summation algorithms currently implemented.

.. file:Ewald summation

.. _Ewald summation:





Ewald summation
_______________

The standard technique for overcoming the problem of summing long ranged periodic potentials is the so called Ewald summation method. The idea is to split the long ranged and singular Coulomb potential to a short ranged singular and long ranged smooth parts, and calculate the long ranged part in reciprocal space via Fourier transformations. This is possible (for a smooth potential) since the system is periodic and the same supercell repeats infinitely in all directions. In practice the calculation can be done by adding (and subtracting) Gaussian charge densities over the point charges to screen the potential in real space. That is, the original charge density :math:`\rho(\mathbf{r}) = \sum_i q_i \delta(\mathbf{r} - \mathbf{r}_i)` is split by
    
    .. math::
      :nowrap:
    
      \begin{eqnarray}
      \rho(\mathbf{r}) & = & \rho_s(\mathbf{r}) + \rho_l(\mathbf{r}) \\
      \rho_s(\mathbf{r}) & = & \sum_i \left[ q_i \delta(\mathbf{r} - \mathbf{r}_i) - q_i G_\sigma(\mathbf{r} - \mathbf{r}_i) \right] \\
      \rho_l(\mathbf{r}) & = & \sum_i q_i G_\sigma(\mathbf{r} - \mathbf{r}_i) \\
      G_\sigma(\mathbf{r}) & = & \frac{1}{(2 \pi \sigma^2)^{3/2}} \exp\left( -\frac{|\mathbf{r}|^2}{2 \sigma^2} \right)
      \end{eqnarray}
    
Here :math:`\rho_l` generates a long range interaction since at large distances the Gaussian densities :math:`G_\sigma` appear the same as point charges (:math:`\lim_{\sigma/r \to 0} G_\sigma(\mathbf{r}) = \delta(\mathbf{r})`). Since the charge density is smooth, so will be the potential it creates. The density :math:`\rho_s` exhibits short ranged interactions for the same reason: At distances longer than the width of the Gaussians the point charges are screened by the Gaussians which exactly cancel them (:math:`\lim_{\sigma/r \to 0} \delta(\mathbf{r}) - G_\sigma(\mathbf{r}) = 0`).
    
The short ranged interactions are directly calculated in real space
    
    .. math::
       :nowrap:
    
       \begin{eqnarray}
       E_s & = & \frac{1}{4 \pi \varepsilon_0} \int \frac{\rho_s(\mathbf{r}) \rho_s(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} \mathrm{d}^3 r \mathrm{d}^3 r' \\
           & = & \frac{1}{4 \pi \varepsilon_0} \sum_{(i,j)} \frac{q_i q_j}{r_{ij}} \mathrm{erfc} \left( \frac{r_{ij}}{\sigma \sqrt{2}} \right) - \frac{1}{4 \pi \varepsilon_0} \frac{1}{\sqrt{2 \pi} \sigma} \sum_i^N q_i^2.
       \end{eqnarray}
    
The complementary error function :math:`\mathrm{erfc}(r) = 1 - \mathrm{erf}(r) = 1 - \frac{2}{\sqrt{\pi}} \int_0^r e^{-t^2/2} \mathrm{d}t` makes the sum converge rapidly as :math:`\frac{r_{ij}}{\sigma} \to \infty`. The latter sum is the self energy of each point charge in the potential of the particular Gaussian that screens the charge, and the sum runs over all charges in the supercell spanning the periodic system. (The self energy is cancelled by the long range part of the energy.)
    
The long ranged interaction 

.. math::
	
	E_l = \frac{1}{4 \pi \varepsilon_0} \int \frac{\rho_l(\mathbf{r}) \rho_l(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} \mathrm{d}^3 r \mathrm{d}^3 r' 

can be calculated in reciprocal space by Fourier transformation. The result is
    
    .. math::
       :nowrap:
    
       \begin{eqnarray}
       E_l & = & \frac{1}{2 V \varepsilon_0} \sum_{\mathbf{k} \ne 0} \frac{e^{-\sigma^2 k^2 / 2}}{k^2} |S(\mathbf{k})|^2 \\
       S(\mathbf{k}) & = & \sum_i^N q_i e^{\mathrm{i} \mathbf{k} \cdot \mathbf{r}_i}
       \end{eqnarray}
    
The sum in :math:`E_l` runs over the reciprocal lattice :math:`\mathbf{k} = k_1 \mathbf{b}_1 + k_2 \mathbf{b}_2 + k_3 \mathbf{b}_3` where :math:`\mathbf{b}_i` are the vectors spanning the reciprocal cell (:math:`[\mathbf{b}_1 \mathbf{b}_2 \mathbf{b}_3] = ([\mathbf{v}_1 \mathbf{v}_2 \mathbf{v}_3]^{-1})^T` where :math:`\mathbf{v}_i` are the real space cell vectors). Likewise the sum in the structure factor :math:`S(\mathbf{k})` runs over all charges in the supercell.
    
The total energy is then the sum of the short and long range energies
    
    .. math::
    
       E = E_s + E_l.
    
    
If the system carries a net charge, the total Coulomb potential of the infinite periodic system is infinite. Excess charge can be neutralized by imposing a uniform background charge of opposite sign, which results in the correction term

.. math::

   E_c = - \frac{\sigma^2}{4 V \varepsilon_0} \left| \sum_i q_i \right|^2.

This correction is applied automatically.

Forces are obtained as the gradient of the total energy. For atom :math:`\alpha`, the force is

.. math::

  \mathbf{F}_\alpha = - \nabla_\alpha E = - \nabla_\alpha E_s - \nabla_\alpha E_l.

(There is no contribution from :math:`E_c`.)
The short ranged interactions are easily calculated in real space

.. math::

  - \nabla_\alpha E_s = \frac{q_\alpha}{4 \pi \varepsilon_0} \sum_{j} q_j \left[ \mathrm{erfc} \left( \frac{r_{\alpha j}}{\sigma \sqrt{2}} \right) \frac{1}{r_{\alpha j}^2}  + \frac{1}{\sigma} \sqrt{\frac{2}{\pi}} \exp\left( -\frac{r_{\alpha j}^2}{2 \sigma^2} \right) \frac{1}{r_{\alpha j}} \right] \hat{r}_{\alpha j},

where :math:`\hat{r}_{\alpha j} = \mathbf{r}_{\alpha j} / r_{\alpha j}` is the unit vector pointing from atom :math:`\alpha` to :math:`j`.
The long range forces are obtained by differentiating the structure factor

.. math::
  :nowrap:

  \begin{eqnarray}
  - \nabla_\alpha E_l & = &- \frac{1}{2 V \varepsilon_0} \sum_{\mathbf{k} \ne 0} \frac{e^{-\sigma^2 k^2 / 2}}{k^2} 2 \mathrm{Re} [ S^*(\mathbf{k}) \nabla_\alpha S(\mathbf{k}) ] \\
  \nabla_\alpha S(\mathbf{k}) & = & q_\alpha \mathbf{k} (- \sin \mathbf{k} \cdot \mathbf{r}_\alpha + i \cos \mathbf{k} \cdot \mathbf{r}_\alpha ).
  \end{eqnarray}






.. file:coulombsummation class - autogenerated

.. _coulombsummation class - autogenerated:




List of methods
---------------
  
Below is a list of methods in :class:`~pysic.interactions.coulomb.CoulombSummation`, grouped according to
the type of functionality.

Initialization
______________

- :meth:`~pysic.interactions.coulomb.CoulombSummation.initialize_parameters` (meant for internal use)
- :meth:`~pysic.interactions.coulomb.CoulombSummation.get_summation`
- :meth:`~pysic.interactions.coulomb.CoulombSummation.set_summation`
- :data:`~pysic.interactions.coulomb.CoulombSummation.summation_modes`
- :data:`~pysic.interactions.coulomb.CoulombSummation.summation_parameter_descriptions`
- :data:`~pysic.interactions.coulomb.CoulombSummation.summation_parameters`

Parameter handling
__________________

- :meth:`~pysic.interactions.coulomb.CoulombSummation.get_parameters`
- :meth:`~pysic.interactions.coulomb.CoulombSummation.set_parameter_value`
- :meth:`~pysic.interactions.coulomb.CoulombSummation.set_parameter_values`
- :meth:`~pysic.interactions.coulomb.CoulombSummation.set_parameters`


Miscellaneous
______________

- :meth:`~pysic.interactions.coulomb.CoulombSummation.get_realspace_cutoff`
- :meth:`~pysic.interactions.coulomb.CoulombSummation.get_scaling_factors`
- :meth:`~pysic.interactions.coulomb.CoulombSummation.set_scaling_factors`




Full documentation of the CoulombSummation class
-------------------------------------------------

.. currentmodule:: pysic.interactions.coulomb
.. autoclass:: CoulombSummation
   :members:
   :undoc-members:

.. autofunction:: estimate_ewald_parameters


