!
! A module for basic quarternion operations
! and 3D spatial rotations using quaternion representation
!
MODULE quaternions
  IMPLICIT NONE
  double precision, parameter :: norm_tolerance = 1.0d-8
  
  ! The quarternion type. It only contains four
  ! real components, but the main advantage for
  ! defining it as a custom type is the possibility
  ! to write routines and operators for quaternion
  ! algebra.
  ! *w the "real" component of the quaternion
  ! *x an "imaginary" component of the quaternion
  ! *y an "imaginary" component of the quaternion
  ! *z an "imaginary" component of the quaternion
  TYPE qtrn
     double precision :: w
     double precision :: x
     double precision :: y
     double precision :: z
  END TYPE qtrn

! Quarternion product operator
INTERFACE operator(.qp.)
   MODULE PROCEDURE qprod
END INTERFACE

! Quarternion product operator
INTERFACE operator(*)
   MODULE PROCEDURE qprod
END INTERFACE

! Quarternion norm operator
INTERFACE operator(.norm.)
   MODULE PROCEDURE qnorm
END INTERFACE

! Quarternion conjugate operator
INTERFACE operator(.c.)
   MODULE PROCEDURE qconj
END INTERFACE

! Quarternion conjugate operator
INTERFACE operator(.conj.)
   MODULE PROCEDURE qconj
END INTERFACE

! Quarternion inverse operator
INTERFACE operator(.inv.)
   MODULE PROCEDURE qinv
END INTERFACE

! Quarternion to angle operator
INTERFACE operator(.angle.)
   MODULE PROCEDURE q2angle
END INTERFACE

! Quarternion to axis operator
INTERFACE operator(.axis.)
   MODULE PROCEDURE q2axis
END INTERFACE

! Quarternion to matrix operator
INTERFACE operator(.mat.)
   MODULE PROCEDURE q2matrix
END INTERFACE

! Quarternion + scalar operator
INTERFACE operator(+)
   MODULE PROCEDURE qplus
END INTERFACE

! Quarternion - scalar operator
INTERFACE operator(-)
   MODULE PROCEDURE qminus
END INTERFACE

! Scalar * quarternion operator
INTERFACE operator(*)
   MODULE PROCEDURE qtimes
END INTERFACE

! Quarternion * scalar operator
INTERFACE operator(*)
   MODULE PROCEDURE qtimesB
END INTERFACE

! Quarternion / scalar operator
INTERFACE operator(/)
   MODULE PROCEDURE qdiv
END INTERFACE

! Vector dot product
INTERFACE operator(.o.)
   MODULE PROCEDURE dot
END INTERFACE

! Vector cross product
INTERFACE operator(.x.)
   MODULE PROCEDURE cross
END INTERFACE

! Vector norm
INTERFACE operator(.norm.)
   MODULE PROCEDURE vnorm
END INTERFACE

! Overloading the "rotate" subroutine
INTERFACE rotate
   MODULE PROCEDURE rotate_q, rotate_au, rotate_a
END INTERFACE

CONTAINS

  ! Returns the quarternion q added by scalar r
  ! component-wise
  ! *q a quaternion
  ! *r a real scalar
  ! *qn q+r
  FUNCTION qplus(q,r) &
       RESULT(qn)
    IMPLICIT NONE
    double precision, INTENT(IN) :: r
    TYPE(qtrn), INTENT(IN) :: q
    TYPE(qtrn) :: qn

    qn%w=q%w+r
    qn%x=q%x+r
    qn%y=q%y+r
    qn%z=q%z+r

  END FUNCTION qplus

  ! Returns the quarternion q subtracted by scalar r
  ! component-wise
  ! *q a quaternion
  ! *r a real scalar
  ! *qn q-r
  FUNCTION qminus(q,r) &
       RESULT(qn)
    IMPLICIT NONE
    double precision, INTENT(IN) :: r
    TYPE(qtrn), INTENT(IN) :: q
    TYPE(qtrn) :: qn

    qn%w=q%w-r
    qn%x=q%x-r
    qn%y=q%y-r
    qn%z=q%z-r

  END FUNCTION qminus

  ! Returns the quarternion q multiplied by scalar r
  ! component-wise
  ! *q a quaternion
  ! *r a real scalar
  ! *qn r*q
  FUNCTION qtimes(r,q) &
       RESULT(qn)
    IMPLICIT NONE
    double precision, INTENT(IN) :: r
    TYPE(qtrn), INTENT(IN) :: q
    TYPE(qtrn) :: qn

    qn%w=q%w*r
    qn%x=q%x*r
    qn%y=q%y*r
    qn%z=q%z*r

  END FUNCTION qtimes

  ! Returns the quarternion q multiplied by scalar r
  ! component-wise
  ! *q a quaternion
  ! *r a real scalar
  ! *qn q*r
  FUNCTION qtimesB(q,r) &
       RESULT(qn)
    IMPLICIT NONE
    double precision, INTENT(IN) :: r
    TYPE(qtrn), INTENT(IN) :: q
    TYPE(qtrn) :: qn

    qn%w=q%w*r
    qn%x=q%x*r
    qn%y=q%y*r
    qn%z=q%z*r

  END FUNCTION qtimesB


  ! Returns the quarternion q divided by scalar r
  ! component-wise
  ! *q a quaternion
  ! *r a real scalar
  ! *qn q/r
  FUNCTION qdiv(q,r) &
       RESULT(qn)
    IMPLICIT NONE
    double precision, INTENT(IN) :: r
    TYPE(qtrn), INTENT(IN) :: q
    TYPE(qtrn) :: qn

    qn%w=q%w/r
    qn%x=q%x/r
    qn%y=q%y/r
    qn%z=q%z/r

  END FUNCTION qdiv

  ! Returns the quarternion product q1*q2
  ! Note that the product is non-commutative: q1*q2 /= q2*q1
  ! *q1 a quaternion
  ! *q2 a quaternion
  ! *qn q1*q2
  FUNCTION qprod(q1,q2) &
       RESULT(qn)
    IMPLICIT NONE
    TYPE(qtrn), INTENT(IN) :: q1, q2
    TYPE(qtrn) :: qn
    
    qn%w = q1%w*q2%w-q1%x*q2%x-q1%y*q2%y-q1%z*q2%z
    qn%x = q1%w*q2%x+q1%x*q2%w+q1%y*q2%z-q1%z*q2%y
    qn%y = q1%w*q2%y-q1%x*q2%z+q1%y*q2%w+q1%z*q2%x
    qn%z = q1%w*q2%z+q1%x*q2%y-q1%y*q2%x+q1%z*q2%w

  END FUNCTION qprod

  ! Returns the quarternion conjugate of q: w+xi+yj+zk -> w-xi-yj-zk
  ! *q a quaternion
  ! *cq conjugate of q
  FUNCTION qconj(q) &
       RESULT(cq)
    IMPLICIT NONE
    TYPE(qtrn), INTENT(IN) :: q
    TYPE(qtrn) :: cq

    cq%w = q%w
    cq%x = -q%x
    cq%y = -q%y
    cq%z = -q%z

  END FUNCTION qconj

  ! Returns the quarternion norm of q: sqrt(w**2+x**2+y**2+z**2)
  ! *q a quaternion
  ! *norm the norm of q
  FUNCTION qnorm(q) &
       RESULT(norm)
    IMPLICIT NONE
    TYPE(qtrn), INTENT(IN):: q
    double precision :: norm

    norm = sqrt(q%w*q%w+q%x*q%x+q%y*q%y+q%z*q%z)

  END FUNCTION qnorm

  ! Returns the quarternion inverse of q: conj(q)/norm(q)
  ! *q a quaternion
  ! *iq inverse of q  
  FUNCTION qinv(q) &
       RESULT(iq)
    IMPLICIT NONE
    TYPE(qtrn), INTENT(IN) :: q
    TYPE(qtrn) :: iq

    iq = qconj(q)/qnorm(q)

  END FUNCTION qinv

  ! Returns the quarternion representing a rotation
  ! around axis u by angle a
  ! *a angle in radians
  ! *u 3D vector, defining an axis of rotation
  ! *q representing the rotation
  FUNCTION rot2q(a,u) &
       RESULT(q)
    IMPLICIT NONE
    double precision, INTENT(IN) :: a, u(3)
    TYPE(qtrn) :: q
    double precision :: normer

    q%w = cos(a/2)
    normer = sin(a/2)/sqrt(u(1)**2+u(2)**2+u(3)**2)
    q%x = u(1)*normer
    q%y = u(2)*normer
    q%z = u(3)*normer

  END FUNCTION rot2q

  ! Returns the quarternion representing a rotation
  ! around axis v by angle |v|. If v is 0,
  ! the quaternion [1 0 0 0] will be returned.
  ! *v 3D vector, defining both the angle and axis of rotation
  ! *q representing the rotation
  FUNCTION vec2q(v) &
       RESULT(q)
    IMPLICIT NONE
    double precision, INTENT(IN) :: v(3)
    TYPE(qtrn) :: q
    double precision :: normer, a

    !a = sqrt(v(1)**2+v(2)**2+v(3)**2)
    a = .norm.v
    IF(ABS(a) < norm_tolerance)THEN ! zero rotation
       q%w = 1.0
       q%x = 0.0
       q%y = 0.0
       q%z = 0.0
       RETURN
    END IF
    q%w = cos(a/2)
    normer = sin(a/2)/a
    q%x = v(1)*normer
    q%y = v(2)*normer
    q%z = v(3)*normer    

  END FUNCTION vec2q

  ! Returns the angle of rotation described by the
  ! UNIT quarternion q. Note that the unity of q
  ! is not checked (as it would be time consuming to
  ! calculate the norm all the time if we know
  ! the quaternions used have unit length).
  ! *q a quaternion representation of rotation
  ! *a angle of rotation
  FUNCTION q2angle(q) &
       RESULT(a)
    IMPLICIT NONE
    double precision :: a
    TYPE(qtrn), INTENT(IN) :: q

    a = 2.d0*acos(q%w)

  END FUNCTION q2angle

  ! Returns the axis of rotation described by the
  ! UNIT quarternion q. Note that the unity of q
  ! is not checked (as it would be time consuming to
  ! calculate the norm all the time if we know
  ! the quaternions used have unit length).
  ! *q a quaternion representation of rotation
  ! *u axis of rotation
  FUNCTION q2axis(q) &
       RESULT(u)
    IMPLICIT NONE
    double precision :: u(3), normer
    TYPE(qtrn), INTENT(IN) :: q

    normer = 1.d0/sin(q2angle(q))
    u(1) = q%x*normer
    u(2) = q%y*normer
    u(3) = q%z*normer

  END FUNCTION q2axis

  ! Returns the rotation matrix described by the
  ! UNIT quarternion q. Note that the unity of q
  ! is not checked (as it would be time consuming to
  ! calculate the norm all the time if we know
  ! the quaternions used have unit length).
  ! *q a quaternion representation of rotation
  ! *mat rotation matrix
  FUNCTION q2matrix(q) &
       RESULT(mat)
    IMPLICIT NONE
    double precision :: mat(3,3), w,x,y,z
    TYPE(qtrn), INTENT(IN) :: q
  
    w=q%w
    x=q%x
    y=q%y
    z=q%z
    mat(1,1) = 1.d0-2.d0*(y**2+z**2)
    mat(1,2) = 2.d0*(x*y-w*z)
    mat(1,3) = 2.d0*(x*z+w*y)
    mat(2,1) = 2.d0*(x*y+w*z)
    mat(2,2) = 1.d0-2.d0*(x**2+z**2)
    mat(2,3) = 2.d0*(y*z-w*x)
    mat(3,1) = 2.d0*(x*z-w*y)
    mat(3,2) = 2.d0*(y*z+w*x)
    mat(3,3) = 1.d0-2.d0*(x**2+y**2)

  END FUNCTION q2matrix

  ! Returns the 3D vector vec rotated according to
  ! the UNIT quarternion q. Note that the unity of q
  ! is not checked (as it would be time consuming to
  ! calculate the norm all the time if we know
  ! the quaternions used have unit length).
  ! *q a quaternion representation of rotation
  ! *vec vector to be rotated
  ! *v the rotated vector
  FUNCTION rotate_q(vec,q) &
       RESULT(v)
    IMPLICIT NONE
    double precision, INTENT(IN) :: vec(3)
    TYPE(qtrn), INTENT(IN) :: q
    double precision :: v(3)
    TYPE(qtrn) :: temp, newq
    
    temp%w=0.d0
    temp%x=vec(1)
    temp%y=vec(2)
    temp%z=vec(3)

    newq = q .qp. temp .qp. (.c.q)
    v(1) = newq%x
    v(2) = newq%y
    v(3) = newq%z

  END FUNCTION rotate_q

  ! Returns the vector vec rotated according to
  ! the axis u and angle a.
  ! *u axis of rotation
  ! *a angle of rotation
  ! *vec vector to be rotated
  ! *v the rotated vector
  FUNCTION rotate_au(vec,a,u) &
       RESULT(v)
    IMPLICIT NONE
    double precision, INTENT(IN) :: vec(3),a,u(3)
    double precision :: v(3)
    
    v = rotate_q(vec,rot2q(a,u))

  END FUNCTION rotate_au

  ! Returns the vector vec rotated according to
  ! the vector da. The axis of rotation is given by
  ! the direction of da and the angle by |da|.
  ! *da rotation vector (angular velocity x time)
  ! *vec vector to be rotated
  ! *v the rotated vector  
  FUNCTION rotate_a(vec,da) &
       RESULT(v)
    IMPLICIT NONE
    double precision, INTENT(IN) :: vec(3),da(3)
    double precision :: v(3)
    
    v = rotate_q(vec,rot2q( &
         sqrt(da(1)**2+da(2)**2+da(3)**2),da) )

  END FUNCTION rotate_a

  ! Normal dot product of vectors (DOT_PRODUCT for 3-vectors only)
  ! *v vector
  ! *u vector
  ! *product v . u
  FUNCTION dot(v,u) &
       RESULT(product)
    IMPLICIT NONE
    double precision, INTENT(IN) :: v(3), u(3)
    double precision :: product

    product = v(1)*u(1)+v(2)*u(2)+v(3)*u(3)

  END FUNCTION dot

  ! Normal cross product of vectors
  ! *v vector
  ! *u vector
  ! *product v x u
  FUNCTION cross(v,u) &
       RESULT(product)
    IMPLICIT NONE
    double precision, INTENT(IN) :: v(3), u(3)
    double precision :: product(3)

    product(1) = v(2)*u(3)-v(3)*u(2)
    product(2) = v(3)*u(1)-v(1)*u(3)
    product(3) = v(1)*u(2)-v(2)*u(1)

  END FUNCTION cross

  ! Norm of a vector, |v|
  ! *v vector
  ! *normi norm of v
  FUNCTION vnorm(v) &
       RESULT(normi)
    IMPLICIT NONE
    double precision, INTENT(IN) :: v(3)
    double precision :: normi

    normi = sqrt(v.o.v)

  END FUNCTION vnorm

  ! norms the given quaternion
  ! *qq quaternion to be normed to unity
  SUBROUTINE norm_quaternion(qq)
    IMPLICIT NONE
    TYPE(qtrn) :: qq

    IF(ABS(.norm.qq) > norm_tolerance)THEN
       qq = qq/.norm.qq
    ELSE
       qq%w = 1.d0
       qq%x = 0.d0
       qq%y = 0.d0
       qq%z = 0.d0
    END IF
    
    RETURN
  END SUBROUTINE norm_quaternion

END MODULE quaternions
