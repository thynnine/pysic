! A C-program for MT19937, with initialization improved 2002/1/26.
! Coded by Takuji Nishimura and Makoto Matsumoto.                 

! Code converted to Fortran 95 by José Rui Faustino de Sousa
! Date: 2002-02-01

! Enhanced version by José Rui Faustino de Sousa
! Date: 2003-04-30

! Interface:
!
! Kinds:
!   genrand_intg
!     Integer kind used must be at least 32 bits.
!   genrand_real
!     Real kind used 
!
! Types:
!   genrand_state
!     Internal representation of the RNG state.
!   genrand_srepr
!     Public representation of the RNG state.
!     Should be used to save the RNG state.
!
! Procedures:
!   assignment(=)
!     Converts from type genrand_state to genrand_srepr and vice versa.
!   genrand_init
!     Internal RNG state initialization subroutine accepts either 
!     an genrand_intg integer or a vector as seed or a new state
!     using "put=" returns the present state using "get=". 
!     If it is called with "get=" before being seeded with "put=" 
!     returns a state initialized with a default seed.
!
!   genrand_int32
!  Subroutine returns an array or scalar whose elements are random
!  integer on the [0,0xffffffff] interval.
!
!   genrand_int31
!  Subroutine returns an array or scalar whose elements are random
!  integer on the [0,0x7fffffff] interval.
!
!   genrand_real1
!     Subroutine returns an array or scalar whose elements are random real on the
!     [0,1] interval.
!   genrand_real2
!     Subroutine returns an array or scalar whose elements are random real on the
!     [0,1[ interval.
!   genrand_real3
!     Subroutine returns an array or scalar whose elements are random real on the
!     ]0,1[ interval.
!   genrand_res53
!     Subroutine returns an array or scalar whose elements are random real on the
!     [0,1[ interval with 53-bit resolution.

! Before using, initialize the state by using genrand_init( put=seed )  

! This library is free software.                                  
! This library is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of  
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            

! Copyright (C) 1997, 2002 Makoto Matsumoto and Takuji Nishimura. 
! Any feedback is very welcome.                                   
! http://www.math.keio.ac.jp/matumoto/emt.html                    
! email: matumoto@math.keio.ac.jp                                 
MODULE mt95
  
  IMPLICIT NONE

  PUBLIC  :: genrand_init, ASSIGNMENT(=)
  PUBLIC  :: genrand_int32, genrand_int31, genrand_real1
  PUBLIC  :: genrand_real2, genrand_real3, genrand_res53
  PRIVATE :: uiadd, uisub, uimlt, uidiv, uimod
  PRIVATE :: init_by_type, init_by_scalar, init_by_array, next_state
  PRIVATE :: genrand_encode, genrand_decode, genrand_load_state, genrand_dump_state
  PRIVATE :: genrand_int32_0d, genrand_int32_1d, genrand_int32_2d, genrand_int32_3d
  PRIVATE :: genrand_int32_4d, genrand_int32_5d, genrand_int32_6d, genrand_int32_7d
  PRIVATE :: genrand_int31_0d, genrand_int31_1d, genrand_int31_2d, genrand_int31_3d
  PRIVATE :: genrand_int31_4d, genrand_int31_5d, genrand_int31_6d, genrand_int31_7d
  PRIVATE :: genrand_real1_0d, genrand_real1_1d, genrand_real1_2d, genrand_real1_3d
  PRIVATE :: genrand_real1_4d, genrand_real1_5d, genrand_real1_6d, genrand_real1_7d
  PRIVATE :: genrand_real2_0d, genrand_real2_1d, genrand_real2_2d, genrand_real2_3d
  PRIVATE :: genrand_real2_4d, genrand_real2_5d, genrand_real2_6d, genrand_real2_7d
  PRIVATE :: genrand_real3_0d, genrand_real3_1d, genrand_real3_2d, genrand_real3_3d
  PRIVATE :: genrand_real3_4d, genrand_real3_5d, genrand_real3_6d, genrand_real3_7d
  PRIVATE :: genrand_res53_0d, genrand_res53_1d, genrand_res53_2d, genrand_res53_3d
  PRIVATE :: genrand_res53_4d, genrand_res53_5d, genrand_res53_6d, genrand_res53_7d

  INTRINSIC :: selected_int_kind, selected_real_kind

  INTEGER, PUBLIC, PARAMETER  :: genrand_intg = SELECTED_INT_KIND( 9 )
  INTEGER, PUBLIC, PARAMETER  :: genrand_real = SELECTED_REAL_KIND( 15 )

  INTEGER, PRIVATE, PARAMETER :: wi = genrand_intg
  INTEGER, PRIVATE, PARAMETER :: wr = genrand_real

  ! Period parameters   
  INTEGER(kind=wi), PRIVATE, PARAMETER :: n = 624_wi
  INTEGER(kind=wi), PRIVATE, PARAMETER :: m = 397_wi

  INTEGER(kind=wi), PRIVATE, PARAMETER :: default_seed = 5489_wi

  INTEGER(kind=wi), PRIVATE, PARAMETER :: fbs = 32_wi
  INTEGER(kind=wi), PRIVATE, PARAMETER :: hbs = fbs / 2_wi
  INTEGER(kind=wi), PRIVATE, PARAMETER :: qbs = hbs / 2_wi
  INTEGER(kind=wi), PRIVATE, PARAMETER :: tbs = 3_wi * qbs

  REAL(kind=wr), PRIVATE, PARAMETER :: p231       = 2147483648.0_wr
  REAL(kind=wr), PRIVATE, PARAMETER :: p232       = 4294967296.0_wr
  REAL(kind=wr), PRIVATE, PARAMETER :: p232_1     = p232 - 1.0_wr
  REAL(kind=wr), PRIVATE, PARAMETER :: pi232      = 1.0_wr / p232
  REAL(kind=wr), PRIVATE, PARAMETER :: pi232_1    = 1.0_wr / p232_1
  REAL(kind=wr), PRIVATE, PARAMETER :: pi227      = 1.0_wr / 134217728.0_wr
  REAL(kind=wr), PRIVATE, PARAMETER :: pi253      = 1.0_wr / 9007199254740992.0_wr
  REAL(kind=wr), PRIVATE, PARAMETER :: p231d232_1 = p231 / p232_1
  REAL(kind=wr), PRIVATE, PARAMETER :: p231_5d232 = ( p231 + 0.5_wr ) / p232

  CHARACTER(len=*), PRIVATE, PARAMETER  :: alph = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
  CHARACTER(len=*), PRIVATE, PARAMETER  :: sepr = "&"
  INTEGER(kind=wi), PRIVATE, PARAMETER  :: alps = 62_wi
  INTEGER(kind=wi), PRIVATE, PARAMETER  :: clen = ( n + 1_wi ) * 7_wi !n * ( ceiling( fbs * log( 2.0 ) / log( alps ) ) + 1 )

  TYPE, PUBLIC :: genrand_state
    PRIVATE
    LOGICAL(kind=wi)                :: ini = .FALSE._wi
    INTEGER(kind=wi)                :: cnt = n+1_wi
    INTEGER(kind=wi), DIMENSION(n)  :: val = 0_wi
  END TYPE genrand_state

  TYPE, PUBLIC :: genrand_srepr
    CHARACTER(len=clen) :: repr
  END TYPE genrand_srepr

  TYPE(genrand_state), PRIVATE, SAVE  :: state

! aph@iki.fi:
  LOGICAL, PARAMETER :: debug=.TRUE.
  INTERFACE random_number
    MODULE PROCEDURE genrand_real1_0d
    MODULE PROCEDURE genrand_real1_1d
    MODULE PROCEDURE genrand_real1_2d
    MODULE PROCEDURE genrand_real1_3d
    MODULE PROCEDURE genrand_real1_4d
    MODULE PROCEDURE genrand_real1_5d
    MODULE PROCEDURE genrand_real1_6d
    MODULE PROCEDURE genrand_real1_7d
  END INTERFACE random_number
! aph until here

  INTERFACE ASSIGNMENT( = )
    MODULE PROCEDURE genrand_load_state
    MODULE PROCEDURE genrand_dump_state
  END INTERFACE ASSIGNMENT( = )

  INTERFACE genrand_init
    MODULE PROCEDURE init_by_type
    MODULE PROCEDURE init_by_scalar
    MODULE PROCEDURE init_by_array
  END INTERFACE genrand_init

  INTERFACE genrand_int32
    MODULE PROCEDURE genrand_int32_0d
    MODULE PROCEDURE genrand_int32_1d
    MODULE PROCEDURE genrand_int32_2d
    MODULE PROCEDURE genrand_int32_3d
    MODULE PROCEDURE genrand_int32_4d
    MODULE PROCEDURE genrand_int32_5d
    MODULE PROCEDURE genrand_int32_6d
    MODULE PROCEDURE genrand_int32_7d
  END INTERFACE genrand_int32

  INTERFACE genrand_int31
    MODULE PROCEDURE genrand_int31_0d
    MODULE PROCEDURE genrand_int31_1d
    MODULE PROCEDURE genrand_int31_2d
    MODULE PROCEDURE genrand_int31_3d
    MODULE PROCEDURE genrand_int31_4d
    MODULE PROCEDURE genrand_int31_5d
    MODULE PROCEDURE genrand_int31_6d
    MODULE PROCEDURE genrand_int31_7d
  END INTERFACE genrand_int31

  INTERFACE genrand_real1
    MODULE PROCEDURE genrand_real1_0d
    MODULE PROCEDURE genrand_real1_1d
    MODULE PROCEDURE genrand_real1_2d
    MODULE PROCEDURE genrand_real1_3d
    MODULE PROCEDURE genrand_real1_4d
    MODULE PROCEDURE genrand_real1_5d
    MODULE PROCEDURE genrand_real1_6d
    MODULE PROCEDURE genrand_real1_7d
  END INTERFACE genrand_real1

  INTERFACE genrand_real2
    MODULE PROCEDURE genrand_real2_0d
    MODULE PROCEDURE genrand_real2_1d
    MODULE PROCEDURE genrand_real2_2d
    MODULE PROCEDURE genrand_real2_3d
    MODULE PROCEDURE genrand_real2_4d
    MODULE PROCEDURE genrand_real2_5d
    MODULE PROCEDURE genrand_real2_6d
    MODULE PROCEDURE genrand_real2_7d
  END INTERFACE genrand_real2

  INTERFACE genrand_real3
    MODULE PROCEDURE genrand_real3_0d
    MODULE PROCEDURE genrand_real3_1d
    MODULE PROCEDURE genrand_real3_2d
    MODULE PROCEDURE genrand_real3_3d
    MODULE PROCEDURE genrand_real3_4d
    MODULE PROCEDURE genrand_real3_5d
    MODULE PROCEDURE genrand_real3_6d
    MODULE PROCEDURE genrand_real3_7d
  END INTERFACE genrand_real3

  INTERFACE genrand_res53
    MODULE PROCEDURE genrand_res53_0d
    MODULE PROCEDURE genrand_res53_1d
    MODULE PROCEDURE genrand_res53_2d
    MODULE PROCEDURE genrand_res53_3d
    MODULE PROCEDURE genrand_res53_4d
    MODULE PROCEDURE genrand_res53_5d
    MODULE PROCEDURE genrand_res53_6d
    MODULE PROCEDURE genrand_res53_7d
  END INTERFACE genrand_res53

  CONTAINS

  elemental FUNCTION uiadd( a, b ) RESULT( c )

    INTRINSIC :: ibits, ior, ishft

    INTEGER( kind = wi ), INTENT( in )  :: a, b

    INTEGER( kind = wi )  :: c

    INTEGER( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = IBITS( a, 0, hbs )
    a2 = IBITS( a, hbs, hbs )
    b1 = IBITS( b, 0, hbs )
    b2 = IBITS( b, hbs, hbs )
    s1 = a1 + b1
    s2 = a2 + b2 + IBITS( s1, hbs, hbs )
    c  = IOR( ISHFT( s2, hbs ), IBITS( s1, 0, hbs ) )
    RETURN

  END FUNCTION uiadd
  
  elemental FUNCTION uisub( a, b ) RESULT( c )

    INTRINSIC :: ibits, ior, ishft

    INTEGER( kind = wi ), INTENT( in )  :: a, b

    INTEGER( kind = wi )  :: c

    INTEGER( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = IBITS( a, 0, hbs )
    a2 = IBITS( a, hbs, hbs )
    b1 = IBITS( b, 0, hbs )
    b2 = IBITS( b, hbs, hbs )
    s1 = a1 - b1
    s2 = a2 - b2 + IBITS( s1, hbs, hbs )
    c  = IOR( ISHFT( s2, hbs ), IBITS( s1, 0, hbs ) )
    RETURN

  END FUNCTION uisub
  
  elemental FUNCTION uimlt( a, b ) RESULT( c )

    INTRINSIC :: ibits, ior, ishft

    INTEGER(kind=wi), INTENT(in)  :: a, b

    INTEGER(kind=wi)  :: c

    INTEGER(kind=wi)  :: a0, a1, a2, a3
    INTEGER(kind=wi)  :: b0, b1, b2, b3
    INTEGER(kind=wi)  :: p0, p1, p2, p3

    a0 = IBITS( a, 0, qbs )
    a1 = IBITS( a, qbs, qbs )
    a2 = IBITS( a, hbs, qbs )
    a3 = IBITS( a, tbs, qbs )
    b0 = IBITS( b, 0, qbs )
    b1 = IBITS( b, qbs, qbs )
    b2 = IBITS( b, hbs, qbs )
    b3 = IBITS( b, tbs, qbs )
    p0 = a0 * b0
    p1 = a1 * b0 + a0 * b1 + IBITS( p0, qbs, tbs )
    p2 = a2 * b0 + a1 * b1 + a0 * b2 + IBITS( p1, qbs, tbs )
    p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + IBITS( p2, qbs, tbs )
    c  = IOR( ISHFT( p1, qbs ), IBITS( p0, 0, qbs ) )
    c  = IOR( ISHFT( p2, hbs ), IBITS( c, 0, hbs ) )
    c  = IOR( ISHFT( p3, tbs ), IBITS( c, 0, tbs ) )
    RETURN

  END FUNCTION uimlt

  elemental FUNCTION uidiv( a, b ) RESULT( c )
    
    INTRINSIC :: btest, ishft

    INTEGER(kind=wi), INTENT(in)  :: a, b

    INTEGER(kind=wi)  :: c

    INTEGER(kind=wi)  :: dl, rl

    IF ( BTEST( a, fbs-1 ) ) THEN
      IF ( BTEST( b, fbs-1 ) ) THEN
        IF ( a < b ) THEN
          c = 0
        ELSE
          c = 1
        END IF
      ELSE
        dl = ISHFT( ISHFT( a, -1 ) / b, 1 )
        rl = uisub( a, uimlt( b, dl ) )
        IF ( rl < b ) THEN
          c = dl
        ELSE
          c = uiadd( dl, 1 )
        END IF
      END IF
    ELSE
      IF ( BTEST( b, fbs-1 ) ) THEN
        c = 0
      ELSE
        c = a / b
      END IF
    END IF
    RETURN

  END FUNCTION uidiv

  elemental FUNCTION uimod( a, b ) RESULT( c )
    
    INTRINSIC :: modulo, btest, ishft

    INTEGER(kind=wi), INTENT(in)  :: a, b

    INTEGER(kind=wi)  :: c

    INTEGER(kind=wi)  :: dl, rl

    IF ( BTEST( a, fbs-1 ) ) THEN
      IF ( BTEST( b, fbs-1 ) ) THEN
        IF ( a < b ) THEN
          c = a
        ELSE
          c = uisub( a, b )
        END IF
      ELSE
        dl = ISHFT( ISHFT( a, -1 ) / b, 1 )
        rl = uisub( a, uimlt( b, dl ) )
        IF ( rl < b ) THEN
          c = rl
        ELSE
          c = uisub( rl, b )
        END IF
      END IF
    ELSE
      IF ( BTEST( b, fbs-1 ) ) THEN
        c = a
      ELSE
        c = MODULO( a, b )
      END IF
    END IF
    RETURN

  END FUNCTION uimod

  SUBROUTINE init_by_type( put, get )

    INTRINSIC :: present

    TYPE(genrand_state), OPTIONAL, INTENT(in ) :: put
    TYPE(genrand_state), OPTIONAL, INTENT(out) :: get

    IF ( PRESENT( put ) ) THEN
      IF ( put%ini ) state = put
    ELSE IF ( PRESENT( get ) ) THEN
      IF ( .NOT. state%ini ) CALL init_by_scalar( default_seed )
      get = state
    ELSE
      CALL init_by_scalar( default_seed )
    END IF
    RETURN

  END SUBROUTINE init_by_type

  ! initializes mt[N] with a seed
  SUBROUTINE init_by_scalar( put )

    INTRINSIC :: ishft, ieor, ibits

    INTEGER(kind=wi), PARAMETER :: mult_a = 1812433253_wi !z'6C078965'

    INTEGER(kind=wi), INTENT(in)  :: put

    INTEGER(kind=wi)  :: i

    state%ini = .TRUE._wi
    state%val(1) = IBITS( put, 0, fbs )
    DO i = 2, n, 1
      state%val(i) = IEOR( state%val(i-1), ISHFT( state%val(i-1), -30 ) )
      state%val(i) = uimlt( state%val(i), mult_a )
      state%val(i) = uiadd( state%val(i), i-1_wi )
      ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. 
      ! In the previous versions, MSBs of the seed affect   
      ! only MSBs of the array mt[].                        
      ! 2002/01/09 modified by Makoto Matsumoto             
      state%val(i) = IBITS( state%val(i), 0, fbs )
      ! for >32 bit machines 
    END DO
    state%cnt = n + 1_wi
!    WRITE(*,*) ''
!    WRITE(*,*) 'Initialized by a scalar ', put
!    WRITE(*,*) ''
    RETURN

  END SUBROUTINE init_by_scalar

  ! initialize by an array with array-length 
  ! init_key is the array for initializing keys 
  ! key_length is its length
  SUBROUTINE init_by_array( put )

    INTRINSIC :: size, max, ishft, ieor, ibits
    
    INTEGER(kind=wi), PARAMETER :: seed_d =    19650218_wi !z'12BD6AA'
    INTEGER(kind=wi), PARAMETER :: mult_a =     1664525_wi !z'19660D'
    INTEGER(kind=wi), PARAMETER :: mult_b =  1566083941_wi !z'5D588B65'
    INTEGER(kind=wi), PARAMETER :: msb1_d = ISHFT( 1_wi, fbs-1 ) !z'80000000'

    INTEGER(kind=wi), DIMENSION(:), INTENT(in)  :: put

    INTEGER(kind=wi)  :: i, j, k, tp, key_length

    CALL init_by_scalar( seed_d )
    key_length = SIZE( put, dim=1 )
    i = 2_wi
    j = 1_wi
    DO k = MAX( n, key_length ), 1, -1
      tp = IEOR( state%val(i-1), ISHFT( state%val(i-1), -30 ) )
      tp = uimlt( tp, mult_a )
      state%val(i) = IEOR( state%val(i), tp )
      state%val(i) = uiadd( state%val(i), uiadd( put(j), j-1_wi ) ) ! non linear 
      state%val(i) = IBITS( state%val(i), 0, fbs ) ! for WORDSIZE > 32 machines 
      i = i + 1_wi
      j = j + 1_wi
      IF ( i > n ) THEN
        state%val(1) = state%val(n)
        i = 2_wi
      END IF
      IF ( j > key_length) j = 1_wi
    END DO
    DO k = n-1, 1, -1
      tp = IEOR( state%val(i-1), ISHFT( state%val(i-1), -30 ) )
      tp = uimlt( tp, mult_b )
      state%val(i) = IEOR( state%val(i), tp )
      state%val(i) = uisub( state%val(i), i-1_wi ) ! non linear 
      state%val(i) = IBITS( state%val(i), 0, fbs ) ! for WORDSIZE > 32 machines 
      i = i + 1_wi
      IF ( i > n ) THEN
        state%val(1) = state%val(n)
        i = 2_wi
      END IF
    END DO
    state%val(1) = msb1_d ! MSB is 1; assuring non-zero initial array
    RETURN

  END SUBROUTINE init_by_array

  SUBROUTINE next_state( )

    INTRINSIC :: ishft, ieor, btest, ibits, mvbits

    INTEGER(kind=wi), PARAMETER :: matrix_a = -1727483681_wi !z'9908b0df'

    INTEGER(kind=wi)  :: i, mld

    IF ( .NOT. state%ini ) CALL init_by_scalar( default_seed )
    DO i = 1, n-m, 1
      mld = IBITS( state%val(i+1), 0, 31 )
      CALL MVBITS( state%val(i), 31, 1, mld, 31 )
      state%val(i) = IEOR( state%val(i+m), ISHFT( mld, -1 ) )
      IF ( BTEST( state%val(i+1), 0 ) ) state%val(i) = IEOR( state%val(i), matrix_a )
    END DO
    DO i = n-m+1, n-1, 1
      mld = IBITS( state%val(i+1), 0, 31 )
      CALL MVBITS( state%val(i), 31, 1, mld, 31 )
      state%val(i) = IEOR( state%val(i+m-n), ISHFT( mld, -1 ) )
      IF ( BTEST( state%val(i+1), 0 ) ) state%val(i) = IEOR( state%val(i), matrix_a )
    END DO
    mld = IBITS( state%val(1), 0, 31 )
    CALL MVBITS( state%val(n), 31, 1, mld, 31 )
    state%val(n) = IEOR( state%val(m), ISHFT( mld, -1 ) )
    IF ( BTEST( state%val(1), 0 ) ) state%val(n) = IEOR( state%val(n), matrix_a )
    state%cnt = 1_wi
    RETURN

  END SUBROUTINE next_state

  elemental SUBROUTINE genrand_encode( chr, val )
    
    INTRINSIC :: len

    CHARACTER(len=*), INTENT(out) :: chr
    INTEGER(kind=wi), INTENT(in ) :: val
    
    INTEGER(kind=wi)  :: i, m, d

    d = val
    chr = ""
    DO i = 1, LEN( chr ), 1
      m = uimod( d, alps ) + 1
      chr(i:i) = alph(m:m)
      d = uidiv( d, alps )
      IF ( d == 0 ) EXIT
    END DO
    RETURN

  END SUBROUTINE genrand_encode

  elemental SUBROUTINE genrand_decode( val, chr )
    
    INTRINSIC :: len, len_trim, trim, adjustl, scan

    INTEGER(kind=wi), INTENT(out) :: val
    CHARACTER(len=*), INTENT(in ) :: chr
    
    INTEGER(kind=wi)        :: i, e, p
    CHARACTER(len=LEN(chr)) :: c

    e = 1
    c = TRIM( ADJUSTL( chr ) )
    val = 0
    DO i = 1, LEN_TRIM( c ), 1
      p = SCAN( alph, c(i:i) ) - 1
      IF( p >= 0 ) THEN
        val = uiadd( val, uimlt( p, e ) )
        e = uimlt( e, alps )
      END IF
    END DO
    RETURN

  END SUBROUTINE genrand_decode

  elemental SUBROUTINE genrand_load_state( stt, rpr )

    INTRINSIC :: scan

    TYPE(genrand_state), INTENT(out)  :: stt
    TYPE(genrand_srepr), INTENT(in )  :: rpr

    INTEGER(kind=wi)    :: i, j
    CHARACTER(len=clen) :: c

    i = 1
    c = rpr%repr
    DO
      j = SCAN( c, sepr )
      IF ( j /= 0 ) THEN
        CALL genrand_decode( stt%val(i), c(:j-1) )
        i = i + 1
        c = c(j+1:)
      ELSE
        EXIT
      END IF
    END DO
    CALL genrand_decode( stt%cnt, c )
    stt%ini = .TRUE._wi
    RETURN

  END SUBROUTINE genrand_load_state

  elemental SUBROUTINE genrand_dump_state( rpr, stt )

    INTRINSIC :: len_trim

    TYPE(genrand_srepr), INTENT(out) :: rpr
    TYPE(genrand_state), INTENT(in ) :: stt

    INTEGER(kind=wi)  :: i, j

    j = 1
    rpr%repr = ""
    DO i = 1, n, 1
      CALL genrand_encode( rpr%repr(j:), stt%val(i) )
      j = LEN_TRIM( rpr%repr ) + 1
      rpr%repr(j:j) = sepr
      j = j + 1
    END DO
    CALL genrand_encode( rpr%repr(j:), stt%cnt )
    RETURN

  END SUBROUTINE genrand_dump_state

  ! generates a random number on [0,0xffffffff]-interval
  SUBROUTINE genrand_int32_0d( y )

    INTRINSIC :: ieor, iand, ishft

    INTEGER(kind=wi), PARAMETER :: temper_a = -1658038656_wi !z'9D2C5680'
    INTEGER(kind=wi), PARAMETER :: temper_b =  -272236544_wi !z'EFC60000'

    INTEGER(kind=wi), INTENT(out) :: y
    
    IF ( state%cnt > n ) CALL next_state( )
    y = state%val(state%cnt)
    state%cnt = state%cnt + 1_wi
    ! Tempering 
    y = IEOR( y, ISHFT( y, -11 ) )
    y = IEOR( y, IAND( ISHFT( y,  7 ), temper_a ) )
    y = IEOR( y, IAND( ISHFT( y, 15 ), temper_b ) )
    y = IEOR( y, ISHFT( y, -18 ) )
    RETURN

  END SUBROUTINE genrand_int32_0d

  SUBROUTINE genrand_int32_1d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:), INTENT(out) :: y

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 1 ), 1
      CALL genrand_int32_0d( y(i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int32_1d

  SUBROUTINE genrand_int32_2d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:), INTENT(out) :: y
    
    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 2 ), 1
      CALL genrand_int32_1d( y(:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int32_2d

  SUBROUTINE genrand_int32_3d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:), INTENT(out) :: y
    
    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 3 ), 1
      CALL genrand_int32_2d( y(:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int32_3d

  SUBROUTINE genrand_int32_4d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:,:), INTENT(out) :: y
    
    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 4 ), 1
      CALL genrand_int32_3d( y(:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int32_4d

  SUBROUTINE genrand_int32_5d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:,:,:), INTENT(out) :: y
    
    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 5 ), 1
      CALL genrand_int32_4d( y(:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int32_5d

  SUBROUTINE genrand_int32_6d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:,:,:,:), INTENT(out) :: y
    
    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 6 ), 1
      CALL genrand_int32_5d( y(:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int32_6d

  SUBROUTINE genrand_int32_7d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:,:,:,:,:), INTENT(out) :: y
    
    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 7 ), 1
      CALL genrand_int32_6d( y(:,:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int32_7d

  ! generates a random number on [0,0x7fffffff]-interval
  SUBROUTINE genrand_int31_0d( y )

    INTRINSIC :: ishft

    INTEGER(kind=wi), INTENT(out) :: y

    CALL genrand_int32_0d( y )
    y = ISHFT( y, -1 )
    RETURN

  END SUBROUTINE genrand_int31_0d

  SUBROUTINE genrand_int31_1d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:), INTENT(out) :: y

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 1 ), 1
      CALL genrand_int31_0d( y(i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int31_1d

  SUBROUTINE genrand_int31_2d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:), INTENT(out) :: y

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 2 ), 1
      CALL genrand_int31_1d( y(:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int31_2d

  SUBROUTINE genrand_int31_3d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:), INTENT(out) :: y

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 3 ), 1
      CALL genrand_int31_2d( y(:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int31_3d

  SUBROUTINE genrand_int31_4d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:,:), INTENT(out) :: y

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 4 ), 1
      CALL genrand_int31_3d( y(:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int31_4d

  SUBROUTINE genrand_int31_5d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:,:,:), INTENT(out) :: y

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 5 ), 1
      CALL genrand_int31_4d( y(:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int31_5d

  SUBROUTINE genrand_int31_6d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:,:,:,:), INTENT(out) :: y

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 6 ), 1
      CALL genrand_int31_5d( y(:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int31_6d

  SUBROUTINE genrand_int31_7d( y )

    INTRINSIC :: size

    INTEGER(kind=wi), DIMENSION(:,:,:,:,:,:,:), INTENT(out) :: y

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( y, 7 ), 1
      CALL genrand_int31_6d( y(:,:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_int31_7d

  ! generates a random number on [0,1]-real-interval
  SUBROUTINE genrand_real1_0d( r )

    INTRINSIC :: REAL

    REAL(kind=wr), INTENT(out)  :: r

    INTEGER(kind=wi)  :: a

    CALL genrand_int32_0d( a )
    r = REAL( a, kind=wr ) * pi232_1 + p231d232_1
    ! divided by 2^32-1

    RETURN

  END SUBROUTINE genrand_real1_0d

  SUBROUTINE genrand_real1_1d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 1 ), 1
      CALL genrand_real1_0d( r(i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real1_1d

  SUBROUTINE genrand_real1_2d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 2 ), 1
      CALL genrand_real1_1d( r(:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real1_2d

  SUBROUTINE genrand_real1_3d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 3 ), 1
      CALL genrand_real1_2d( r(:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real1_3d

  SUBROUTINE genrand_real1_4d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 4 ), 1
      CALL genrand_real1_3d( r(:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real1_4d

  SUBROUTINE genrand_real1_5d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 5 ), 1
      CALL genrand_real1_4d( r(:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real1_5d

  SUBROUTINE genrand_real1_6d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 6 ), 1
      CALL genrand_real1_5d( r(:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real1_6d

  SUBROUTINE genrand_real1_7d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 7 ), 1
      CALL genrand_real1_6d( r(:,:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real1_7d

  ! generates a random number on [0,1)-real-interval
  SUBROUTINE genrand_real2_0d( r )

    INTRINSIC :: REAL

    REAL(kind=wr), INTENT(out)  :: r

    INTEGER(kind=wi)  :: a

    CALL genrand_int32_0d( a )
    r = REAL( a, kind=wr ) * pi232 + 0.5_wr
    ! divided by 2^32

    RETURN

  END SUBROUTINE genrand_real2_0d

  SUBROUTINE genrand_real2_1d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 1 ), 1
      CALL genrand_real2_0d( r(i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real2_1d

  SUBROUTINE genrand_real2_2d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 2 ), 1
      CALL genrand_real2_1d( r(:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real2_2d

  SUBROUTINE genrand_real2_3d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 3 ), 1
      CALL genrand_real2_2d( r(:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real2_3d

  SUBROUTINE genrand_real2_4d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 4 ), 1
      CALL genrand_real2_3d( r(:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real2_4d

  SUBROUTINE genrand_real2_5d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 5 ), 1
      CALL genrand_real2_4d( r(:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real2_5d

  SUBROUTINE genrand_real2_6d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 6 ), 1
      CALL genrand_real2_5d( r(:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real2_6d

  SUBROUTINE genrand_real2_7d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 7 ), 1
      CALL genrand_real2_6d( r(:,:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real2_7d

  ! generates a random number on (0,1)-real-interval
  SUBROUTINE genrand_real3_0d( r )

    INTRINSIC :: REAL

    REAL(kind=wr), INTENT(out)  :: r

    INTEGER(kind=wi)  :: a

    CALL genrand_int32_0d( a )
    r = REAL( a, kind=wr ) * pi232 + p231_5d232
    ! divided by 2^32 
    RETURN

  END SUBROUTINE genrand_real3_0d

  SUBROUTINE genrand_real3_1d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 1 ), 1
      CALL genrand_real3_0d( r(i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real3_1d

  SUBROUTINE genrand_real3_2d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 2 ), 1
      CALL genrand_real3_1d( r(:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real3_2d

  SUBROUTINE genrand_real3_3d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 3 ), 1
      CALL genrand_real3_2d( r(:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real3_3d

  SUBROUTINE genrand_real3_4d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 4 ), 1
      CALL genrand_real3_3d( r(:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real3_4d

  SUBROUTINE genrand_real3_5d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 5 ), 1
      CALL genrand_real3_4d( r(:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real3_5d

  SUBROUTINE genrand_real3_6d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 6 ), 1
      CALL genrand_real3_5d( r(:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real3_6d

  SUBROUTINE genrand_real3_7d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 7 ), 1
      CALL genrand_real3_6d( r(:,:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_real3_7d

  ! generates a random number on [0,1) with 53-bit resolution
  SUBROUTINE genrand_res53_0d( r )

    INTRINSIC :: ishft, REAL

    REAL(kind=wr), INTENT(out)  :: r

    INTEGER(kind=wi)  :: a, b

    CALL genrand_int32_0d( a )
    CALL genrand_int32_0d( b )
    a = ISHFT( a, -5 )
    b = ISHFT( b, -6 )
    r = REAL( a, kind=wr ) * pi227 + REAL( b, kind=wr ) * pi253
    RETURN

  END SUBROUTINE genrand_res53_0d

  SUBROUTINE genrand_res53_1d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 1 ), 1
      CALL genrand_res53_0d( r(i) )
    END DO
    RETURN

  END SUBROUTINE genrand_res53_1d

  SUBROUTINE genrand_res53_2d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 2 ), 1
      CALL genrand_res53_1d( r(:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_res53_2d

  SUBROUTINE genrand_res53_3d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 3 ), 1
      CALL genrand_res53_2d( r(:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_res53_3d

  SUBROUTINE genrand_res53_4d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 4 ), 1
      CALL genrand_res53_3d( r(:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_res53_4d

  SUBROUTINE genrand_res53_5d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 5 ), 1
      CALL genrand_res53_4d( r(:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_res53_5d

  SUBROUTINE genrand_res53_6d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 6 ), 1
      CALL genrand_res53_5d( r(:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_res53_6d

  SUBROUTINE genrand_res53_7d( r )

    INTRINSIC :: size

    REAL(kind=wr), DIMENSION(:,:,:,:,:,:,:), INTENT(out)  :: r

    INTEGER(kind=wi)  :: i

    DO i = 1, SIZE( r, 7 ), 1
      CALL genrand_res53_6d( r(:,:,:,:,:,:,i) )
    END DO
    RETURN

  END SUBROUTINE genrand_res53_7d
  ! These real versions are due to Isaku Wada, 2002/01/09 added 
  ! Altered by José Sousa genrand_real[1-3] will not return exactely
  ! the same values but should have the same properties and are faster

END MODULE mt95

