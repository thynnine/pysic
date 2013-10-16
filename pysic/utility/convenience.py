#! /usr/bin/env python

def expand_symbols_string(symbol_string):
    """Expands a string of chemical symbols to list.
    
    The function parses a string of chemical symbols and turns
    it to a list such as those expected by 
    :class:`~pysic.interactions.local.Potential`.
    
    Examples::
    
      >>> print expand_symbols_string("HH")
      [['H', 'H']]
      >>> print expand_symbols_string("SiSi,SiO,SiH")
      [['Si', 'Si'], ['Si', 'O'], ['Si', 'H']]
    
    Parameters:
    
    symbol_string: string
        the string to be expanded
    """

    parts = symbol_string.split(',')
    symbol_list = []
    for part in parts:
        newsymbol = ""
        newset = []
        for letter in part:
            if letter.isupper():
                if newsymbol != "":
                    newset.append(newsymbol)
                newsymbol = letter
            else:
                newsymbol += letter
        newset.append(newsymbol)
        symbol_list.append(newset)

    return symbol_list


def expand_symbols_table(symbol_list,type=None):
    """Creates a table of symbols object.
            
            The syntax for defining the targets is precise
            but somewhat cumbersome due to the large number of permutations one gets
            when the number of bodies increases. Oftentimes one does not need such
            fine control over all the parameters since many of them have the same
            numerical values. Therefore it is convenient to be able to define
            the targets in a more compact way.
            
            This method generates the detailed target tables from compact syntax.
            By default, the method takes a list of list and multiplies each list
            with the others (note the call for a static method)::
            
             >>> pysic.utility.convenience.expand_symbols_table([  'Si',
             ...                                                  ['O', 'C'],
             ...                                                  ['H', 'O'] ])
             [['Si', 'O', 'H'],
              ['Si', 'C', 'H'],
              ['Si', 'O', 'O'],
              ['Si', 'C', 'O']]
            
            Other custom types of formatting can be defined with the type parameter.
            
            For type 'triplet', the target list is created for triplets A-B-C from an input list of the form::
            
             ['A', 'B', 'C']
            
            Remember that in the symbol table accepted by the BondOrderParameters, one needs to define
            the B-A and B-C bonds separately and so B appears as the first symbol in the output and the other
            two appear as second and third (both cases)::
            
             [['B', 'A', 'C'],
              ['B', 'C', 'A']]
            
            However, for an A-B-A triplet, the A-B bond should only be defined once to prevent double counting.
            Like the default function, also here several triplets can be defined at once::
            
             >>> pysic.BondOrderParameters.expand_symbols_table([ ['H', 'O'],
             ...                                                   'Si',
             ...                                                  ['O', 'C'] ],
             ...                                                type='triplet')
             [['Si', 'H', 'O'],
              ['Si', 'O', 'H'],
              ['Si', 'H', 'C'],
              ['Si', 'C', 'H'],
              ['Si', 'O', 'O'],
              ['Si', 'O', 'C'],
              ['Si', 'C', 'O']]
            
            
            Parameters:
            
            symbol_list: list of strings
            list to be expanded to a table
            type: string
            specifies a custom way of generating the table
            """
    if not isinstance(symbol_list,list):
        return [symbol_list]
    n_slots = len(symbol_list)
    n_subslots = []
    for sub in symbol_list:
        if isinstance(sub,list):
            n_subslots.append(len(sub))
        else:
            n_subslots.append(1)
        
    table = []
    if(type == None):
        slot_indices = n_slots*[0]
        while (slot_indices[n_slots-1] < n_subslots[n_slots-1]):
            row = []
            for i in range(n_slots):
                if isinstance(symbol_list[i],list):
                    symb = symbol_list[i][slot_indices[i]]
                else:
                    symb = symbol_list[i]
                row.append(symb)
            table.append(row)
            slot_indices[0] += 1
            for i in range(n_slots):
                if slot_indices[i] >= n_subslots[i]:
                    if(i < n_slots-1):
                        slot_indices[i] = 0
                        slot_indices[i+1] += 1
                    else:
                        pass
    elif(type == 'triplet'):
        if n_slots != 3:
            return None
        for i in range(n_subslots[1]):
            for j in range(n_subslots[0]):
                for k in range(n_subslots[2]):
                    
                    if isinstance(symbol_list[1],list):
                        symb1 = symbol_list[1][i]
                    else:
                        symb1 = symbol_list[1]
                    if isinstance(symbol_list[0],list):
                        symb0 = symbol_list[0][j]
                    else:
                        symb0 = symbol_list[0]
                    if isinstance(symbol_list[2],list):
                        symb2 = symbol_list[2][k]
                    else:
                        symb2 = symbol_list[2]
                        
                    table.append([symb1,symb0,symb2])
                    if symb0 != symb2:
                        table.append([symb1,symb2,symb0])
        
        
    return table

