.. file:chargerelaxation class

.. _chargerelaxation class:



.. file:chargerelaxation class - description

.. _chargerelaxation class - description:



.. module:: pysic.charges
.. module:: pysic.charges.relaxation

=======================
ChargeRelaxation class
=======================

This class controls equilibration of atomic charges in the system.

It is possible for the user to define the charges of atoms in ASE.
If a system exhibits charge transfer, polarization, charged defects etc., one may not know the charges beforehand or the charges may change dynamically during simulation. To handle such systems, it is possible to let the charges in the system develop dynamically.

Since charge dynamics are usually much faster than dynamics of the ions, it is usually reasonable to allow the charges to equilibrate between ionic steps. This does not conserve energy exactly, however, since the charge equilibration drives the system charge distribution towards a lower energy. The energy change in charge redistribution is lost unless it is fed back to the system.

.. file:connecting to Pysic

.. _connecting to Pysic:




Connecting the structure, calculator and relaxation algorithm
-------------------------------------------------------------


Special care must be taken when setting up links between the atomic structure (`ASE Atoms`_), the calculator (:class:`~pysic.calculator.Pysic`), and the charge relaxation algorithm (:class:`~pysic.charges.relaxation.ChargeRelaxation`). While some of the objects must know the others, in some cases the behavior of the simulator changes depending on whether or not they have access to the other objects.

The atoms and the calculator are linked as required in the `ASE calculator interface`_: One can link the two by either the :meth:`~pysic.calculator.Pysic.set_atoms` method of :class:`~pysic.calculator.Pysic`, or the `set_calculator`_ method of `ASE Atoms`_. In either case, the atomic structure is given a link to the calculator, and a **copy** of the structure is stored in the calculator. This must be done in order to do any calculations on the system.

Also the relaxation algorithm has to know the :class:`~pysic.calculator.Pysic` calculator, since the relaxation is done according to the :class:`~pysic.interactions.local.Potential` interactions stored in the calculator. The algorithm can be made to know the calculator via the :meth:`~pysic.ChargeRelaxation.set_calculator` method of :class:`~pysic.charges.relaxation.ChargeRelaxation`. By default, this does not make the calculator know the relaxation algorithm, however. Only if the optional argument ``reciprocal=True`` is given the backwards link is also made. :class:`~pysic.calculator.Pysic` can be made to know the relaxation algorithm also by calling the method :meth:`~pysic.calculator.Pysic.set_charge_relaxation`. Unlike the opposite case, by making the link from the calculator, the backwards link from the relaxation algorithm is always made automatically. In fact, even though linking an algorithm to a calculator does not automatically link the calculator to the algorithm, if a different calculator was linked to the algorithm, the link is automatically removed.

This slightly complicated behavior is summarized as follows: The algorithm should always have a link to a calculator, but a calculator need not have a link to an algorithm. If a calculator does link to an algorithm, the algorithm must link back to the same calculator. Clearly one does not always want to perform charge relaxation on the system and so it makes sense that the calculator need not have a link to a charge relaxation algorithm. If such a link does exist, then the relaxation is *automatically* invoked prior to each energy and force evaluation. This is necessary in simulations such as molecular dynamics (MD). A charge relaxation can be linked to a calculator in order to do charge equilibration, but if one does not wish to trigger the charge relaxation automatically, then it is enough to just not let the calculator know the relaxation algorithm.

The atomic structure cannot be given a link to the relaxation algorithm since the charge relaxation is not part of the ASE API and so the atoms object does not know how to interact with it. In essence, from the point of view of the structure, the charge relaxation is fully contained in the calculator.

The charge relaxation algorithm always acts on the structure contained in the calculator. The atomic charges of this structure are automatically updated during the relaxation. Since the calculator only stores a copy of the original structure, the original is not updated. This may be desired if, for instance, one wishes to revert back to the original charges. However, during strucural dynamics simulations such as MD, it is necessary that the relaxed charges are saved between structural steps. This is a problem, since structural dynamics are handled by ASE, and ASE invokes the calculation of forces with the original `ASE Atoms`_ object. Therefore, if the relaxed charges are not saved, charge relaxation is always started from the original charges, which may be very inefficient. In order to have also the original structure updated automatically, the charge relaxation can be made to know the original structure with :meth:`~pysic.ChargeRelaxation.set_atoms`. Note that the structure given to the algorithm is not used in the actual relaxation; the algorithm always works on the structure in the calculator, which may be different. The given structure is merely updated according to the calculation results.


.. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
.. _set_calculator: https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.atoms.Atoms.set_calculator
.. _ASE calculator interface: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#calculator-interface

.. file:piping chargerelaxations

.. _piping chargerelaxations:



Piping ChargeRelaxation algorithms
---------------------------------------

Sometimes the most efficient way to optimize a system is to combine several algorithms or several parameterizations of a single algorithm in a sequence. First, one can run a robust but less efficient algorithm to find a rough solution and follow with a more efficient method to pinpoint the optimum at high precision. For a plain relaxation, this can be done manually. However, if the optimization needs to happen automatically, e.g., between molecular dynamics steps, one needs to combine several ChargeRelaxation objects together. This can be done through the method :meth:`~pysic.charges.relaxation.ChargeRelaxation.add_relaxation_pipe`.

Adding a piped ChargeRelaxation object simply calls the :meth:`~pysic.charges.relaxation.ChargeRelaxation.charge_relaxation` methods of both objects one after another when relaxation is invoked. Several objects can be piped by recursive piping. For instance, the following::

  rel1 = ChargeRelaxation(...)
  rel2 = ChargeRelaxation(...)
  rel3 = ChargeRelaxation(...)

  rel1.add_relaxation_pipe(rel2)
  rel2.add_relaxation_pipe(rel3)

creates a pipe ``rel1``-``rel2``-``rel3``, where ``rel1`` is executed first and ``rel3`` last. Creating looped pipes where the same relaxation appears twice is not allowed.

By default, when a pipe is created, the calculator and structure objects of the first relaxation are copied also to all the piped ChargeRelaxation objects (including recursively defined ones).

.. file:observers

.. _observers:



Observers
----------------

Similarly to `ASE molecular dynamics <https://wiki.fysik.dtu.dk/ase/tutorials/md/md.html>`_, also the ChargeRelaxation allows one to attach callback functions to the dynamics. That is, one can have a specific functions be automatically invoked at given intervals. This can be used for instance for printing statistics or saving data during the optimization run. For instance, the following example saves and plots the charge of the atom with index 0::

  system = Atoms(...)
  calc = Pysic(...)
  rel = ChargeRelaxation(calculator = calc, system = system)

  charge0 = []
  def save_charge():
      charge0.append(calc.get_atoms().get_charges()[0])

  rel.add_observer(save_charge)
  rel.charge_relaxation()

  import matplotlib.pyplot as plt
  plt.plot(np.array(range(len(charge0))), thecharges)
  plt.show()

.. file:list of relaxation modes

.. _list of relaxation modes:




List of currently available relaxation methods
----------------------------------------------

Below is a list of the charge relaxation methods currently implemented.


.. file:damped dynamics

.. _damped dynamics:




Damped dynamics
_______________

Assigning an inertia, :math:`M_q`, on the atomic charges, :math:`q_i`, we can describe the system with the Lagrangian

.. math::

    L = \sum_i \frac{1}{2} m_i \dot{\mathbf{r}}_i^2 + \sum_i \frac{1}{2} M_q \dot{q}_i^2 - U(\{q\},\{\mathbf{r}\}) - \nu \sum_i q_i,

where :math:`m_i, \mathbf{r}_i` are the mass and position of atom :math:`i`, respectively. The last term is a Lagrange multiplier corresponding to the constraint of fixed total charge, i.e., :math:`\sum_i q_i = Q_\mathrm{tot}` being constant. The total potential energy :math:`U` is a function of all charges and positions.

The equations of motion for this system are

.. math::
    :nowrap:

    \begin{eqnarray}
    m_i \ddot{\mathbf{r}}_i & = & -\nabla_i U \\
    M_q \ddot{q}_i & = & -\frac{\partial U}{\partial q_i} - \nu.
    \end{eqnarray}

In the charge equation, the Lagrange multiplier can be shown to equal the average electronegativity of the system, :math:`\nu = \bar{\chi}`, and the derivative is the effective electronegativity of atom :math:`i`, :math:`-\frac{\partial U}{\partial q_i} = \chi_i`. Thus, the effective force driving the change in atomic charges is the electronegativity difference from the avarage

.. math::

    M_q \ddot{q}_i =  -\frac{\partial U}{\partial q_i} - \nu = \chi_i - \bar{\chi} = \Delta \chi_i.

In the damped dynamic equilibration, the charges are developed dynamically according to the equation of motion with an added damping (friction) term :math:`- \eta \dot{q}_i`

.. math::
    :label: charge_equation

    M_q \ddot{q}_i = \Delta \chi_i - \eta \dot{q}_i.

This leads to the charges being driven towards a state where the driving force vanishes :math:`\Delta \chi_i = 0`, i.e., the electronegativities are equal.

During simulation such as molecular dynamics or geometry optimization, charge equilibration is done by running the damped charge dynamics :eq:`charge_equation` before each force or energy evaluation.

Keywords::

    >>> names_of_parameters('dynamic')
    ['n_steps', 'timestep', 'inertia', 'friction', 'tolerance']


.. file:potentiostat

.. _potentiostat:




Potentiostat
_______________

Similarly to :ref:`damped dynamics`, a potentiostat drives the atomic charges :math:`q_i`
with the dynamics

.. math::
    :nowrap:

    \begin{eqnarray}
    m_i \ddot{\mathbf{r}}_i & = & -\nabla_i U \\
    M_q \ddot{q}_i & = & -\frac{\partial U}{\partial q_i} - \Phi.
    \end{eqnarray}

Charge equilibration tries to optimize the charge distribution (with fixed total charge) so that the total energy is minimized. The potentiostat, however, connects the system to an electrode (a charge resevoir) at constant potential :math:`\Phi`. This means that the total charge of the system is allowed to change so that the electronegativity of each atom equals the external potential :math:`\chi_i = -\frac{\partial U}{\partial q_i} = \Phi`.


Keywords::

    >>> names_of_parameters('potentiostat')
    ['n_steps', 'timestep', 'inertia', 'potential']

.. file:optimization

.. _optimization:



Optimization
_______________

This charge optimization method uses the 
`<SciPy fmin_slsqp http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_slsqp.html>`_ function
for constrained optimization using `Sequential Least Squares Programming (SLSQP) <http://en.wikipedia.org/wiki/Sequential_quadratic_programming>`_. 

Instead of dynamic simulation, the algorithm searches for energy minimum :math:`U(\mathbf{q})` with the constraint of constant charge :math:`\sum_i q_i = Q`.

As the external function does not support callback functions, observers cannot be attached to this algorthm.


Keywords::

    >>> names_of_parameters('optimize')
    ['n_steps', 'tolerance']

.. file:chargerelaxation class - autogenerated

.. _chargerelaxation class - autogenerated:



List of methods
---------------
  
Below is a list of methods in :class:`~pysic.charges.relaxation.ChargeRelaxation`, grouped according to
the type of functionality.
  
Initialization
______________

- :meth:`~pysic.charges.relaxation.ChargeRelaxation.initialize_parameters`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.get_relaxation`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_relaxation`
- :data:`~pysic.charges.relaxation.ChargeRelaxation.relaxation_modes`
- :data:`~pysic.charges.relaxation.ChargeRelaxation.relaxation_parameter_descriptions`
- :data:`~pysic.charges.relaxation.ChargeRelaxation.relaxation_parameters`

Atoms handling
______________

- :meth:`~pysic.charges.relaxation.ChargeRelaxation.get_atoms`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_atoms`

Calculator handling
___________________

- :meth:`~pysic.charges.relaxation.ChargeRelaxation.get_calculator`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_calculator`


Parameter handling
__________________

- :meth:`~pysic.charges.relaxation.ChargeRelaxation.get_parameters`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_parameter_value`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_parameter_values`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_parameters`


Pipe handling
______________
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.clear_relaxation_pipe`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.get_relaxation_pipe`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_relaxation_pipe`


Observer handling
__________________
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.add_observer`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.call_observers`
- :meth:`~pysic.charges.relaxation.ChargeRelaxation.clear_observers`



Charge relaxation
_________________

- :meth:`~pysic.charges.relaxation.ChargeRelaxation.charge_relaxation`


Full documentation of the ChargeRelaxation class
------------------------------------------------

.. currentmodule:: pysic.charges.relaxation
.. autoclass:: ChargeRelaxation
   :members:
   :undoc-members:

