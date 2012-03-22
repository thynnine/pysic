
.. _utility:
        
=============================================
utility (Utility.f90)
=============================================



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
    
    
    

Full documentation of subroutines in utility
--------------------------------------------
        
        
            
  .. function:: int2str(length, ints, string_out)


    Parameters:

    length: integer  *intent(in)*    *scalar*  
        
    ints: integer  *intent(in)*    *size(length)*  
        
    **string_out**: character(len=length)  **intent(out)**    *scalar*  
        
            
  .. function:: pad_string(str_in, str_length, str_out)


    Parameters:

    str_in: character(len=*)  *intent(in)*    *scalar*  
        
    str_length: integer  *intent(in)*    *scalar*  
        
    **str_out**: character(len=str_length)  **intent(out)**    *scalar*  
        
            
  .. function:: str2int(length, string, ints)


    Parameters:

    length: integer  *intent(in)*    *scalar*  
        
    string: character(len=length)  *intent(in)*    *scalar*  
        
    **ints**: integer  **intent(out)**    *size(length)*  
        