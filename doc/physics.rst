.. file:physics

.. _physics:



.. file:physics forewords

.. _physics forewords:



Physics of Pysic
----------------

In this section we briefly describe the physical motivation and main ideas behind Pysic. The details of actual algorithms, functions, and implementations included in Pysic are discussed later in sections :ref:`run forewords` and :ref:`syntax forewords`.

.. file:target problems

.. _target problems:



Applications and goals
-----------------------

Phenomena such as Si covalent bonding and charge redistribution at defects and interfaces make semiconductors a very challenging group of materials to describe computationally. Usually quantum mechanical methods such as density functional theory are used, but these approaches are limited by their high computational cost. More sophisticated empirical potentials are being developed for these materials, however, combining fairly accurate precision with a reasonable cost. Notably, the ReaxFF [1]_ and COMB [2]_ [3]_ variable charge potentials have been recently introduced and shown to reproduce the structural properties of for instance Si, SiO\ :sub:`2` and HfO\ :sub:`2`. Still these methods are efficient enough to allow the simulation of semiconductor systems at size and time scales relevant for actual device performance (hundreds of thousands of atoms). 



.. file:interactions

.. _interactions:



Atomistic interactions
-----------------------

To accurately model condensed matter systems, one often needs to know how the electrons behave, as electrons determine the chemical properies of atoms. Even if the primary interest is the atomic structure, not the electronic one, it may be necessary to include electrons in the simulations in order to calculate the atomistic interactions. The electronic structure must be treated on a quantum mechanical level making these kinds of calculations quite heavy, and so one would often wish to bypass the electronic level of detail and only work with the atomistic structure and a classical description. If we assume that the atomic and electronic degrees of freedom are separate (this is the Born-Oppenheimer approximation, and it is often justified), then in principle the electronic ground state :math:`| \psi \rangle` is determined by the atomic coordinates :math:`\{ \mathbf{R}_i \}` and nuclear charges :math:`\{ Q_i \}` through the Schrödinger equation

.. math::
  :nowrap:

  \begin{eqnarray}
  H(\{ \mathbf{R}_i \}, \{ Q_i \}) | \psi \rangle  & = & E | \psi \rangle \\
  \Rightarrow | \psi \rangle & = & \psi(\{ \mathbf{R}_i \}, \{ Q_i \}),
  \end{eqnarray}

where :math:`H = H(\{ \mathbf{R}_i \}, \{ Q_i \})` is the Hamiltonian.
As the total energy of the system, :math:`E`, is determined by the electronic state, it too is a function of the atomic configuration

.. math::

 E = \langle \psi | H | \psi \rangle = E(\{ \mathbf{R}_i \}, \{ Q_i \}).

So, in principle, the energy of the system and the atomic forces :math:`F_\alpha = -\nabla_{\mathbf{R}_\alpha} E` can be determined from the atomic configuration. Construction of the energy function :math:`E(\{ \mathbf{R}_i \}, \{ Q_i \})` is very challenging though. It is a very complicated function of the coordinates of *all* the atoms in the system, and in practice one needs to further assume that the system can be split in local substructures for which the energy function can be parameterized, and that the total energy is obtained as a sum of the local contributions. For :math:`n`-body local interactions, this could be written

.. math::
 :label: local_energy

 E \approx \sum_{(i_1,i_2,\ldots,i_n)} E_\mathrm{local}(\mathbf{R}_{i_1}, \ldots, \mathbf{R}_{i_n}, Q_{i_1}, \ldots, Q_{i_n}),

where :math:`(i_1,i_2,\ldots,i_n)` are all the sets of :math:`n` atoms in the system.



.. file:charge equilibration

.. _charge equilibration:



Dynamic charges
----------------

In practice, calculations get exponentially more complex when the number of bodies included in the local :math:`n`-body energy :math:`E_\mathrm{local}(\mathbf{R}_{i_1}, \ldots, \mathbf{R}_{i_n}, Q_{i_1}, \ldots, Q_{i_n})` in :eq:`local_energy` increases. Therefore usually only few particles close to each other are included in the local energy function, and all configurations with distant atoms are given a zero energy. This approach does not work, however, if the system exhibits long-ranged collective electronic reconstruction such as local charging or polarization, for instance, due to the presence of defects or interfaces. Changes in charge distribution can drastically alter the local energy function :math:`E_\mathrm{local}`. Increasing the number of bodies :math:`n` in the local interaction function does not solve this problem in general, since the charge redistribution may be a system wide phenomenon that no local function can properly capture.

One way to treat the redistribution of charge is to make the local energy also a function of the total atomic charge :math:`q_i = Q_i - \eta_i e`, where :math:`\eta_i` is the (possibly fractional) number of electrons associated with the atom

.. math::

  E_\mathrm{local} = E_\mathrm{local}(\mathbf{R}_{i_1}, \ldots, \mathbf{R}_{i_n}, Q_{i_1}, \ldots, Q_{i_n}, q_{i_1}, \ldots, q_{i_n}).

Although electrons do not specifically belong to any atom making :math:`\eta_i` ambiguous, this approach allows the inclusion of long ranged charge distribution in the model with reasonable computational cost. 

In addition, having the local charge as a parameter of the energy function allows for the optimization of the energy with respect to the local charges, allowing one to search for an equilibrium charge distribution. This is analogous to finding equilibrium structures by optimizing the energy with respect to the atomic coordinates :math:`\mathbf{R}_i`.  



.. file:hybrid calculations

.. _hybrid calculations:



Hybrid QM/MM calculations in Pysic
------------------------------------

Pysic provides a framework for creating and running hybrid QM/MM simulations. Within this framework it is possible to calculate the potential energy and forces in an atomistic configuration which has multiple subsystems and interactions between them. Any external calculator with an ASE interface can be assigned to a subsystem or the classical potentials provided by Pysic can be used. The QM/MM implementation in Pysic uses the mechanical embedding scheme with hydrogen link atoms. It is possible to enable any Pysic-supported interaction potentials between the subsystems, with special functions provided for the easy use of Coulomb and COMB interactions.

.. file:pysic approach

.. _pysic approach:



Pysic approach
----------------

The immediate aim in the development of Pysic is to implement atomistic potentials with variable local atomic charges and apply them in the study of semiconductor interfaces, e.g., silicon-hafnia. In the long term, Pysic will include a full library of different atomistic potentials.

Pysic implements a very general framework for calculating energies and forces due to arbitrary local atomistic pair or many body potentials. It is straightforward to implement new types of interactions in the code, and mixing different potentials during the simulations is simple. Furthermore, one can easily evaluate the contribution of different interactions on the total energy and forces by switching on and off specific interactions. So called bond order, or Tersoff, potentials [4]_ are also supported, and the user is free to scale any potential with a bond-dependent factor. In addition, in a system with local charges, long ranged Coulomb interactions need to be evaluated. Such :math:`1/r`-potentials are calculated with the standard Ewald summation algorithm. Implementation of other algorithms such as Particle mesh Ewald [5]_ or Wolf summation [6]_ is also planned.

In addition, it is planned that various advanced analysis tools are included with the Pysic package. These would include tools for tasks such as potential parametrization or structural analysis using techniques like evolutionary algorithms, machine learning, or Bayesian mehods.

.. file:references

.. _references:



.. [1] \ A. van Duin, S. Dasgupta, F. Lorant, and W. Goddard, J Phys Chem A 105, 9396 (2001). 
.. [2] \ T.-R. Shan, B. D. Devine, T. W. Kemper, S. B. Sinnott, and S. R. Phillpot, Phys Rev B 81, 125328 (2010).   
.. [3] \ T.-R. Shan, D. Bryce, J. Hawkins, A. Asthagiri, S. Phillpot, and S. Sinnott, Phys Rev B 82, 235302 (2010).  

.. [4] \ J. Tersoff, Phys Rev B 37, 6991 (1988).   
.. [5] \ T. Darden, D. York, and L. Pedersen, Journal of Chemical Physics 98, 10089 (1993). 
.. [6] \ D. Wolf, P. Keblinski, S. Phillpot, and J. Eggebrecht, Journal of Chemical Physics 110, 8254 (1999).
