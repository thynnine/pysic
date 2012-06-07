#! /usr/bin/env python


codec_s2i = {'1' : -1,
             '2' : -2,
             '3' : -3,
             '4' : -4,
             '5' : -5,
             '6' : -6,
             '7' : -7,
             '8' : -8,
             '9' : -9,
             '0' : -10,
             'a' : 1,
         'b' : 2, 
         'c' : 3, 
         'd' : 4, 
         'e' : 5, 
         'f' : 6, 
         'g' : 7, 
         'h' : 8, 
         'i' : 9, 
         'j' : 10,
         'k' : 11,
         'l' : 12,
         'm' : 13,
         'n' : 14,
         'o' : 15,
         'p' : 16,
         'q' : 17,
         'r' : 18,
         's' : 19,
         't' : 20,
         'u' : 21,
         'v' : 22,
         'w' : 23,
         'x' : 24,
         'y' : 25,
         'z' : 26,
         'A' : 101,
         'B' : 102,
         'C' : 103,
         'D' : 104,
         'E' : 105,
         'F' : 106,
         'G' : 107,
         'H' : 108,
         'I' : 109,
         'J' : 110,
         'K' : 111,
         'L' : 112,
         'M' : 113,
         'N' : 114,
         'O' : 115,
         'P' : 116,
         'Q' : 117,
         'R' : 118,
         'S' : 119,
         'T' : 120,
         'U' : 121,
         'V' : 122,
         'W' : 123,
         'X' : 124,
         'Y' : 125,
         'Z' : 126,
         '_' : 201 }

codec_i2s = {-1 : '1',
             -2 : '2',
             -3 : '3',
             -4 : '4',
             -5 : '5',
             -6 : '6',
             -7 : '7',
             -8 : '8',
             -9 : '9',
             -10 : '0',
             1 : 'a',
         2 : 'b',
         3 : 'c',
         4 : 'd',
         5 : 'e',
         6 : 'f',
         7 : 'g',
         8 : 'h',
         9 : 'i',
         10 : 'j',
         11 : 'k',
         12 : 'l',
         13 : 'm',
         14 : 'n',
         15 : 'o',
         16 : 'p',
         17 : 'q',
         18 : 'r',
         19 : 's',
         20 : 't',
         21 : 'u',
         22 : 'v',
         23 : 'w',
         24 : 'x',
         25 : 'y',
         26 : 'z',
         101 : 'A',
         102 : 'B',
         103 : 'C',
         104 : 'D',
         105 : 'E',
         106 : 'F',
         107 : 'G',
         108 : 'H',
         109 : 'I',
         110 : 'J',
         111 : 'K',
         112 : 'L',
         113 : 'M',
         114 : 'N',
         115 : 'O',
         116 : 'P',
         117 : 'Q',
         118 : 'R',
         119 : 'S',
         120 : 'T',
         121 : 'U',
         122 : 'V',
         123 : 'W',
         124 : 'X',
         125 : 'Y',
         126 : 'Z',
         201 : '_' }




def char2int(char_in):
    """Codes a single character to an integer.
    """
    output = codec_s2i.get(char_in)
    if output == None:
        return 0
    else:
        return output

def int2char(int_in):
    """Decodes an integer to a single character.
    """
    output = codec_i2s.get(int_in)
    if output == None:
        return ' '
    else:
        return output

def str2ints(string_in,target_length=0):
    """Codes a string to a list of integers.

    Turns a string to a list of integers for f2py interfacing.
    If required, the length of the list can be specified and trailing spaces
    will be added to the end.
    """
    ints_out = []
    for char in string_in:
        ints_out.append(char2int(char))
    while len(ints_out) < target_length:
        ints_out.append(char2int(' '))
    return ints_out

def ints2str(ints_in):
    """Decodes a list of integers to a string.
    """
    string_out = ""
    for number in ints_in:
        string_out += int2char(number)
    return string_out


