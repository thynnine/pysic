module utility
  implicit none
  double precision, parameter :: pi = 4.d0 * datan(1.d0)

contains
  subroutine int2str(length,ints,string_out)
    implicit none
    integer, intent(in) :: length
    integer, intent(in) :: ints(length)
    character(len=length), intent(out) :: string_out
    character(len=length+1) :: string
    integer :: i

    string = ''

    do i = 1, length
       select case(ints(i))
       case(-1)
          string = string(1:i)//'1'
       case(-2)
          string = string(1:i)//'2'
       case(-3)
          string = string(1:i)//'3'
       case(-4)
          string = string(1:i)//'4'
       case(-5)
          string = string(1:i)//'5'
       case(-6)
          string = string(1:i)//'6'
       case(-7)
          string = string(1:i)//'7'
       case(-8)
          string = string(1:i)//'8'
       case(-9)
          string = string(1:i)//'9'
       case(-10)
          string = string(1:i)//'0'

       case(1)
          string = string(1:i)//'a'
       case(2)
          string = string(1:i)//'b'
       case(3)
          string = string(1:i)//'c'
       case(4)
          string = string(1:i)//'d'
       case(5)
          string = string(1:i)//'e'
       case(6)
          string = string(1:i)//'f'
       case(7)
          string = string(1:i)//'g'
       case(8)
          string = string(1:i)//'h'
       case(9)
          string = string(1:i)//'i'
       case(10)
          string = string(1:i)//'j'
       case(11)
          string = string(1:i)//'k'
       case(12)
          string = string(1:i)//'l'
       case(13)
          string = string(1:i)//'m'
       case(14)
          string = string(1:i)//'n'
       case(15)
          string = string(1:i)//'o'
       case(16)
          string = string(1:i)//'p'
       case(17)
          string = string(1:i)//'q'
       case(18)
          string = string(1:i)//'r'
       case(19)
          string = string(1:i)//'s'
       case(20)
          string = string(1:i)//'t'
       case(21)
          string = string(1:i)//'u'
       case(22)
          string = string(1:i)//'v'
       case(23)
          string = string(1:i)//'w'
       case(24)
          string = string(1:i)//'x'
       case(25)
          string = string(1:i)//'y'
       case(26)
          string = string(1:i)//'z'

       case(101)
          string = string(1:i)//'A'
       case(102)
          string = string(1:i)//'B'
       case(103)
          string = string(1:i)//'C'
       case(104)
          string = string(1:i)//'D'
       case(105)
          string = string(1:i)//'E'
       case(106)
          string = string(1:i)//'F'
       case(107)
          string = string(1:i)//'G'
       case(108)
          string = string(1:i)//'H'
       case(109)
          string = string(1:i)//'I'
       case(110)
          string = string(1:i)//'J'
       case(111)
          string = string(1:i)//'K'
       case(112)
          string = string(1:i)//'L'
       case(113)
          string = string(1:i)//'M'
       case(114)
          string = string(1:i)//'N'
       case(115)
          string = string(1:i)//'O'
       case(116)
          string = string(1:i)//'P'
       case(117)
          string = string(1:i)//'Q'
       case(118)
          string = string(1:i)//'R'
       case(119)
          string = string(1:i)//'S'
       case(120)
          string = string(1:i)//'T'
       case(121)
          string = string(1:i)//'U'
       case(122)
          string = string(1:i)//'V'
       case(123)
          string = string(1:i)//'W'
       case(124)
          string = string(1:i)//'X'
       case(125)
          string = string(1:i)//'Y'
       case(126)
          string = string(1:i)//'Z'
       case(201)
          string = string(1:i)//'_'
       case default
          string = string(1:i)//' '
       end select
    end do

    string_out = string(2:length+1)

  end subroutine int2str

  subroutine str2int(length,string,ints)
    implicit none
    integer, intent(in) :: length
    integer, intent(out) :: ints(length)
    character(len=length), intent(in) :: string
    integer :: i

    do i = 1, length
       select case(string(i:i))
       case('1')
          ints(i) = -1
       case('2')
          ints(i) = -2
       case('3')
          ints(i) = -3
       case('4')
          ints(i) = -4
       case('5')
          ints(i) = -5
       case('6')
          ints(i) = -6
       case('7')
          ints(i) = -7
       case('8')
          ints(i) = -8
       case('9')
          ints(i) = -9
       case('0')
          ints(i) = -10

       case('a')
          ints(i) = 1
       case('b')
          ints(i) = 2
       case('c')
          ints(i) = 3
       case('d')
          ints(i) = 4
       case('e')
          ints(i) = 5
       case('f')
          ints(i) = 6
       case('g')
          ints(i) = 7
       case('h')
          ints(i) = 8
       case('i')
          ints(i) = 9
       case('j')
          ints(i) = 10
       case('k')
          ints(i) = 11
       case('l')
          ints(i) = 12
       case('m')
          ints(i) = 13
       case('n')
          ints(i) = 14
       case('o')
          ints(i) = 15
       case('p')
          ints(i) = 16
       case('q')
          ints(i) = 17
       case('r')
          ints(i) = 18
       case('s')
          ints(i) = 19
       case('t')
          ints(i) = 20
       case('u')
          ints(i) = 21
       case('v')
          ints(i) = 22
       case('w')
          ints(i) = 23
       case('x')
          ints(i) = 24
       case('y')
          ints(i) = 25
       case('z')
          ints(i) = 26

       case('A')
          ints(i) = 101
       case('B')
          ints(i) = 102
       case('C')
          ints(i) = 103
       case('D')
          ints(i) = 104
       case('E')
          ints(i) = 105
       case('F')
          ints(i) = 106
       case('G')
          ints(i) = 107
       case('H')
          ints(i) = 108
       case('I')
          ints(i) = 109
       case('J')
          ints(i) = 110
       case('K')
          ints(i) = 111
       case('L')
          ints(i) = 112
       case('M')
          ints(i) = 113
       case('N')
          ints(i) = 114
       case('O')
          ints(i) = 115
       case('P')
          ints(i) = 116
       case('Q')
          ints(i) = 117
       case('R')
          ints(i) = 118
       case('S')
          ints(i) = 119
       case('T')
          ints(i) = 120
       case('U')
          ints(i) = 121
       case('V')
          ints(i) = 122
       case('W')
          ints(i) = 123
       case('X')
          ints(i) = 124
       case('Y')
          ints(i) = 125
       case('Z')
          ints(i) = 126

       case('_')
          ints(i) = 201
       case default
          ints(i) = 0
       end select
    end do

  end subroutine 

  subroutine pad_string(str_in,str_length,str_out)
    implicit none
    character(len=*), intent(in) :: str_in
    integer, intent(in) :: str_length
    character(len=str_length), intent(out) :: str_out
    !integer :: i

    if(len(str_in) > str_length)then
       str_out = str_in(1:str_length)
    !else if(len(str_in) < str_length)then
    !   str_out
    else
       str_out = str_in
    end if

  end subroutine pad_string


end module utility
