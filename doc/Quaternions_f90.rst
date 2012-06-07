
.. _quaternions:
        
=====================================================
quaternions (Quaternions.f90)
=====================================================



A module for basic quarternion operations
and 3D spatial rotations using quaternion representation
Quaternions are a 4-component analogue to complex numbers

.. math::

 \mathbf{q} = w + x \mathbf{i} + y \mathbf{j} + z \mathbf{k} = [w,x,y,z]

where the three imaginary components obey :math:`\mathbf{ijk} = \mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = -1`.
This leads to a structure similar to that of complex numbers, except that the quaternion product defined
according to the above rules is non-commutative :math:`\mathbf{q}_1\mathbf{q}_2 \ne \mathbf{q}_2\mathbf{q}_1`.

It turns out that unit quaternions :math:`||\mathbf{q}|| = \sqrt{w^2+x^2+y^2+z^2} = 1` represent the space
of 3D rotations so that the rotation by angle :math:`\alpha` around unit axis :math:`\mathbf{u} = [x,y,z]`
is represented by the quaternion

.. math::

  \mathbf{q} = [\cos \frac{\alpha}{2}, x \sin \frac{\alpha}{2}, y \sin \frac{\alpha}{2}, z \sin \frac{\alpha}{2}]

and joining rotations is represented as a a quaternion product
(rotating first by :math:`\mathbf{q}_1`, then by :math:`\mathbf{q}_2` yields the combined rotation
of :math:`\mathbf{q}_{12} = \mathbf{q}_2 \mathbf{q}_1`).

.. only:: html


    List of global variables in quaternions
    ---------------------------------------
    - :data:`norm_tolerance`

    List of custom types in quaternions
    -----------------------------------
    - :data:`qtrn`

    List of subroutines in quaternions
    ----------------------------------
        
    - :func:`norm_quaternion`

    List of functions in quaternions
    --------------------------------
        
    - :func:`cross`
    - :func:`dot`
    - :func:`q2angle`
    - :func:`q2axis`
    - :func:`q2matrix`
    - :func:`qconj`
    - :func:`qdiv`
    - :func:`qinv`
    - :func:`qminus`
    - :func:`qnorm`
    - :func:`qplus`
    - :func:`qprod`
    - :func:`qtimes`
    - :func:`qtimesB`
    - :func:`rot2q`
    - :func:`rotate_a`
    - :func:`rotate_au`
    - :func:`rotate_q`
    - :func:`vec2q`
    - :func:`vnorm`


Full documentation of global variables in quaternions
-----------------------------------------------------
        
        
  .. data:: norm_tolerance

    double precision    *scalar*  *parameter*  

    *initial value* = 1.0d-8
    
    the threshold value for the norm for treating vectors as zero vectors
    

Full documentation of custom types in quaternions
-------------------------------------------------
        
        
  .. data:: qtrn

    The quarternion type. It only contains four
    real components, but the main advantage for
    defining it as a custom type is the possibility
    to write routines and operators for quaternion
    algebra.

    Contained data:

    y: double precision    *scalar*
        an "imaginary" component of the quaternion
    x: double precision    *scalar*
        an "imaginary" component of the quaternion
    z: double precision    *scalar*
        an "imaginary" component of the quaternion
    w: double precision    *scalar*
        the "real" component of the quaternion

Full documentation of subroutines in quaternions
------------------------------------------------
        
        
            
  .. function:: norm_quaternion(qq)

    norms the given quaternion

    Parameters:

    qq: TYPE(qtrn)  *intent()*    *scalar*  
        quaternion to be normed to unity

Full documentation of functions in quaternions
------------------------------------------------
        
        
            
  .. function:: cross(v, u)

    Normal cross product of vectors :math:`\mathbf{v} \times \mathbf{u}` (Note: for 3-vectors only!)

    Parameters:

    v: double precision  *intent()*    *size(3)*  
        vector
    u: double precision  *intent()*    *size(3)*  
        vector
            
  .. function:: dot(v, u)

    Normal dot product of vectors :math:`\mathbf{v}\cdot\mathbf{u}` (Note: for 3-vectors only!)

    Parameters:

    v: double precision  *intent()*    *size(3)*  
        vector
    u: double precision  *intent()*    *size(3)*  
        vector
            
  .. function:: q2angle(q)

    Returns the angle of rotation described by the
    UNIT quarternion :math:`\mathbf{q}`. Note that the unity of :math:`\mathbf{q}`
    is not checked (as it would be time consuming to
    calculate the norm all the time if we know
    the quaternions used have unit length).

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion representation of rotation
            
  .. function:: q2axis(q)

    Returns the axis of rotation described by the
    UNIT quarternion :math:`\mathbf{q}`. Note that the unity of :math:`\mathbf{q}`
    is not checked (as it would be time consuming to
    calculate the norm all the time if we know
    the quaternions used have unit length).

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion representation of rotation
            
  .. function:: q2matrix(q)

    Returns the rotation matrix described by the
    UNIT quarternion :math:`\mathbf{q}`. Note that the unity of :math:`\mathbf{q}`
    is not checked (as it would be time consuming to
    calculate the norm all the time if we know
    the quaternions used have unit length).

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion representation of rotation
            
  .. function:: qconj(q)

    Returns the quarternion conjugate of :math:`\mathbf{q}`: :math:`\mathbf{q}^* = w+x\mathbf{i}+y\mathbf{j}+z\mathbf{k} \to w-x\mathbf{i}-y\mathbf{j}-z\mathbf{k}`

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
            
  .. function:: qdiv(q, r)

    Returns the quarternion :math:`\mathbf{q}` divided by scalar :math:`r`
    component-wise

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
    r: double precision  *intent()*    *scalar*  
        a real scalar
            
  .. function:: qinv(q)

    Returns the quarternion inverse of :math:`\mathbf{q}`: :math:`\mathbf{q}^*/||\mathbf{q}||`

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
            
  .. function:: qminus(q, r)

    Returns the quarternion :math:`\mathbf{q}` subtracted by scalar :math:`r`
    component-wise

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
    r: double precision  *intent()*    *scalar*  
        a real scalar
            
  .. function:: qnorm(q)

    Returns the quarternion norm of :math:`\mathbf{q}`: :math:`||\mathbf{q}|| = \sqrt{w^2+x^2+y^2+z^2}`

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
            
  .. function:: qplus(q, r)

    Returns the quarternion :math:`\mathbf{q}` added by scalar :math:`r`
    component-wise

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
    r: double precision  *intent()*    *scalar*  
        a real scalar
            
  .. function:: qprod(q1, q2)

    Returns the quarternion product :math:`\mathbf{q}_1\mathbf{q}_2`
    Note that the product is non-commutative:  :math:`\mathbf{q}_1\mathbf{q}_2 \ne \mathbf{q}_2\mathbf{q}_1`

    Parameters:

    q1: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
    q2: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
            
  .. function:: qtimes(r, q)

    Returns the quarternion :math:`\mathbf{q}` multiplied by scalar :math:`r`
    component-wise

    Parameters:

    r: double precision  *intent()*    *scalar*  
        a real scalar
    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
            
  .. function:: qtimesB(q, r)

    Returns the quarternion :math:`\mathbf{q}` multiplied by scalar :math:`r`
    component-wise

    Parameters:

    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion
    r: double precision  *intent()*    *scalar*  
        a real scalar
            
  .. function:: rot2q(a, u)

    Returns the quarternion representing a rotation
    around axis :math:`\mathbf{u}` by angle :math:`\alpha`

    Parameters:

    a: double precision  *intent()*    *scalar*  
        angle in radians
    u: double precision  *intent()*    *size(3)*  
        3D vector, defining an axis of rotation
            
  .. function:: rotate_a(vec, da)

    Returns the vector rotated according to
    the vector :math:`\mathbf{d}`. The axis of rotation is given by
    the direction of :math:`\mathbf{d}` and the angle by :math:`||\mathbf{d}||`.

    Parameters:

    vec: double precision  *intent()*    *size(3)*  
        vector to be rotated
    da: double precision  *intent()*    *size(3)*  
        rotation vector (e.g., angular velocity x time :math:`\mathbf{\omega} t`)
            
  .. function:: rotate_au(vec, a, u)

    Returns the vector rotated according to
    the axis :math:`\mathbf{u}` and angle :math:`\alpha`.

    Parameters:

    vec: double precision  *intent()*    *size(3)*  
        vector to be rotated
    a: double precision  *intent()*    *scalar*  
        angle of rotation
    u: double precision  *intent()*    *size(3)*  
        axis of rotation
            
  .. function:: rotate_q(vec, q)

    Returns the 3D vector rotated according to
    the UNIT quarternion :math:`\mathbf{q}`. Note that the unity of :math:`\mathbf{q}`
    is not checked (as it would be time consuming to
    calculate the norm all the time if we know
    the quaternions used have unit length).

    Parameters:

    vec: double precision  *intent()*    *size(3)*  
        vector to be rotated
    q: TYPE(qtrn)  *intent()*    *scalar*  
        a quaternion representation of rotation
            
  .. function:: vec2q(v)

    Returns the quarternion representing a rotation
    around axis :math:`\mathbf{v}` by angle :math:`||\mathbf{v}||`. If :math:`\mathbf{v} = 0`,
    the quaternion :math:`\mathbf{q} = [1 0 0 0]` will be returned.

    Parameters:

    v: double precision  *intent()*    *size(3)*  
        3D vector, defining both the angle and axis of rotation
            
  .. function:: vnorm(v)

    Norm of a vector, :math:`||\mathbf{v}||`

    Parameters:

    v: double precision  *intent()*    *size(3)*  
        vector