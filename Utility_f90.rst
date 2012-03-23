
.. _utility:
        
=============================================
utility (Utility.f90)
=============================================



A module containing utility functionas and constants.

The module is a collection of standalone helper tools,
not physically relevant functionality.


.. only:: html


    List of global variables in utility
    -----------------------------------
    - :data:`pi`

    List of subroutines in utility
    ------------------------------
        
    - :func:`int2str`
    - :func:`pad_string`
    - :func:`str2int`


Full documentation of global variables in utility
-------------------------------------------------
        
        
  .. data:: pi

    double precision    *scalar*  *parameter*  

    *initial value* = 4.d0 * datan(1.d0)
    
    the constant :math:`\pi` calculated as :math:`\pi = 4 \arctan 1`
    

Full documentation of subroutines in utility
--------------------------------------------
        
        
            
  .. function:: int2str(length, ints, string_out)

    Transforms an integer array to string through a codec::
    
       1 - a
       2 - b
       ...
       101 - A
       102 - B
       ...
       -1 - 1
       -2 - 2
       ...
    
    Unrecognized numbers are treated as white spaces.
    
    The function is used for communicating string arrays between
    Python and Fortran oer f2py.
    

    Parameters:

    length: integer  *intent(in)*    *scalar*  
        string length
    ints: integer  *intent(in)*    *size(length)*  
        the integers
    **string_out**: character(len=length)  **intent(out)**    *scalar*  
        the string
            
  .. function:: pad_string(str_in, str_length, str_out)

    Adds spaces after the given string to
    create a string of a certain length.
    If the given string is longer than the
    specified length, it is truncated.
    
    This is used to ensure strings are of a
    certain length, since character arrays in
    Fortran may be forced to a certain length.
    

    Parameters:

    str_in: character(len=*)  *intent(in)*    *scalar*  
        the original string
    str_length: integer  *intent(in)*    *scalar*  
        the required string length
    **str_out**: character(len=str_length)  **intent(out)**    *scalar*  
        the padded string
            
  .. function:: str2int(length, string, ints)

    Transforms a string to an integer array through a codec::
    
       1 - a
       2 - b
       ...
       101 - A
       102 - B
       ...
       -1 - 1
       -2 - 2
       ...
    
    Unrecognized characters are mapped to 0.
    
    The function is used for communicating string arrays between
    Python and Fortran oer f2py.
    

    Parameters:

    length: integer  *intent(in)*    *scalar*  
        string length
    string: character(len=length)  *intent(in)*    *scalar*  
        the string
    **ints**: integer  **intent(out)**    *size(length)*  
        the integers