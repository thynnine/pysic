.. file:Mersenne_f90


.. _mt95:
        
=============================================
mt95 (Mersenne.f90)
=============================================


The module mt95 contains a random number generator.
Currently 'Mersenne twister' is used, but it could in principle be replaced with any generator. This is an external module not written as part of the Pysic project.

The generator is initialized with the routine :func:`genrand_init` and random real numbers  in [0,1] and [0,1) are extracted with :func:`genrand_real1` and :func:`genrand_real2`, respectively. 


Routines of the mt95 module
---------------------------

.. function:: genrand_init(seed)

.. function:: genrand_real1(random)

.. function:: genrand_real2(random)

