#! /usr/bin/env python
import inspect
from ase.parallel import rank


class Warning:
    """A warning raised due to a potentially dangerous action.

    The warning is not an exception, so it doesn't by default interrupt execution.
    It will, however, display a message or perform an action depending on the
    warning settings.

    The warning levels are:

        1: a condition that leads to unwanted behaviour (e.g., creating redundant potential in core)
        2: a condition that is likely unwanted, but possibly a hack
        3: a condition that is likely unwanted, but is a featured hack (e.g., bond order mixing)
        4: a condition that is harmless but may lead to errors later (e.g., defining a potential without targets - often you'll specify the targets later)
        5: a condition that may call for attention (notes to user)

    Parameters:

    message: string
        information describing the cause for the warning
    level: int
        severity of the warning (1-5, 1 being most severe)
    """

    warning_level = 5  # 0-6, 0: no warnings, 5: all warnings, 6: also interrupt
    headers = ['', 'WARNING', 'Warning', 'NOTE', 'Note', 'note']

    def __init__(self, message, level):
        self.message = message
        self.level = level

    def display(self):
        """Displays the message related to this warning, but only if the level of
        the warning is smaller (more severe) than the global warning_level.
        """
        if Warning.warning_level >= self.level:
            warn = self.message + "\nfrom function: " + str(inspect.stack()[2][3])
            if rank == 0:
                print style_message("PYSIC " + self.headers[self.level], warn)

        if Warning.warning_level == 6:
            if self.level <= 3:
                raise WarningInterruptException(self.message)


class WarningInterruptException(Exception):
    """An error raised automatically for any warning when
    strict warnings are in use (warning_level = 6).
    """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return """Interrupted execution due to a warning.

        Strict warnings are on. To control the warnings level, use
        pysic.utility.error.set_warning_level(level), where 'level'
        is an integer between 0 (no warnings) and 6 (treat warnings as errors).
        """


def set_warning_level(level):
    """Set the warning level.

    Parameters:

    level: int
        The level of warnings displayed. Should be an integer between 0 (no warnings) and 6 (warnings are treated as errors).
    """
    if level < 0 or level > 6:
        print "warning level must be between 0 (no warnings) and 6 (warnings are errors)"
    else:
        Warning.warning_level = level


def style_message(header, message, width=80, x_pad=2, y_pad=1):
    """Styles a message to be printed into console.

    The message can be a list of string, where each list item corresponds to a
    row.  If a single string is provided, it is converted to a list. Rows are
    cut into pieces so that they will fit into the defined width.
    """

    if isinstance(message, str):
        message = [message]

    # The rows are split if newlines are present
    broken_message = []
    for line in message:
        broken_message += line.split('\n')

    # The rows are wrapped to fit the message width
    wrapped_message = []
    for line in broken_message:
        if line != "":
            wrapped_message += [line[x:x+(width-2-2*x_pad)] for x in range(0, len(line), (width-2-2*x_pad))]
        else:
            wrapped_message.append("")

    str_message = "\n"
    str_message += "|" + str((width-2)*"=")+"|\n"
    str_message += "|"
    space = width-len(header)-2
    pre_space = space/2
    post_space = space-pre_space
    str_message += str((pre_space)*" ")
    str_message += header
    str_message += str((post_space)*" ")
    str_message += "|\n"
    str_message += "|" + str((width-2)*"=")+"|\n"
    for row in range(y_pad):
        str_message += "|"+str((width-2)*" ")+"|\n"

    for line in wrapped_message:
        line = "|" + str((x_pad)*" ") + line + str((x_pad)*" ")
        length = len(line)
        space = width - length
        if space > 0:
            line += str((space-1)*" ")
        line += "|\n"
        str_message += line

    for row in range(y_pad):
        str_message += "|"+str((width-2)*" ")+"|\n"

    str_message += "|" + str((width-2)*"=")+"|\n"

    return str_message


def warn(message, level):
    """Raise and display a warning.

    Parameters:

    message: string
        information describing the cause for the warning
    level: int
        severity of the warning (1-5, 1 being most severe)
    """

    wrn = Warning(message, level)
    wrn.display()
    return wrn


def error(message):
    """Raises an error and displays the reason. Independent of the warning level.

    Parameters:

    message: string
        information describing the cause for the warning
    """
    if rank == 0:
        print style_message("PYSIC ERROR", message)
    raise Exception()


class InvalidPotentialError(Exception):
    """An error raised when an invalid potential is about to be created or used.

    Parameters:

    message: string
        information describing why the error occurred
    potential: :class:`~pysic.interactions.local.Potential`
        the errorneous potential
    """
    def __init__(self, message='', potential=None):
        self.message = message
        self.potential = potential

    def __str__(self):
        if(self.potential is None):
            return self.message
        else:
            return self.messsage + " \n the Potential: " + str(self.potential)


class InvalidCoordinatorError(Exception):
    """An error raised when an invalid coordination calculator is about to be created or used.

    Parameters:

    message: string
        information describing why the error occurred
    coordinator: :class:`~pysic.interactions.bondorder.Coordinator`
        the errorneous coordinator
    """
    def __init__(self, message='', coordinator=None):
        self.message = message
        self.coordinator = coordinator

    def __str__(self):
        if(self.coordinator is None):
            return self.message
        else:
            return self.messsage + " \n  the Coordinator: " + str(self.coordinator)


class InvalidParametersError(Exception):
    """An error raised when an invalid set of parameters is about to be created.

    Parameters:

    message: string
        information describing why the error occurred
    params: :class:`~pysic.interactions.bondorder.BondOrderParameters`
        the errorneous parameters
    """
    def __init__(self, message='', params=None):
        self.message = message
        self.params = params

    def __str__(self):
        if(self.params is None):
            return self.message
        else:
            return self.messsage + " \n  the Parameters: " + str(self.params)


class InvalidSummationError(Exception):
    """An error raised when an invalid coulomb summation is about to be created.

    Parameters:

    message: string
        information describing why the error occurred
    params: :class:`~pysic.interactions.coulomb.CoulombSummation`
        the errorneous summation
    """
    def __init__(self, message='', summer=None):
        self.message = message
        self.summer = summer

    def __str__(self):
        if(self.summer is None):
            return self.message
        else:
            return self.messsage + " \n  the Summation: " + str(self.summer)


class InvalidRelaxationError(Exception):
    """An error raised when an invalid charge relaxation is about to be created.

        Parameters:

        message: string
            information describing why the error occurred
        params: :class:`~pysic.charges.relaxation.ChargeRelaxation`
            the errorneous parameters
        """
    def __init__(self, message='', relaxation=None):
        self.message = message
        self.relaxation = relaxation

    def __str__(self):
        if(self.params is None):
            return self.message
        else:
            return self.messsage + " \n  the Relaxation: " + str(self.relaxation)


class MissingAtomsError(Exception):
    """An error raised when the core is being updated with per atom information before updating the atoms.

    Parameters:

    message: string
        information describing why the error occurred
    """
    def __init__(self, message=''):
        self.message = message

    def __str__(self):
        return self.message


class MissingNeighborsError(Exception):
    """An error raised when a calculation is initiated without initializing the neighbor lists.

    In principle :class:`~pysic.calculator.Pysic` should always take care of handling the neighbors automatically.
    This error is an indication that there is loophole in the built-in preparations.

    Parameters:

    message: string
        information describing why the error occurred
    """
    def __init__(self, message=''):
        self.message = message

    def __str__(self):
        return self.message


class LockedCoreError(Exception):
    """An error raised when a :class:`~pysic.calculator.Pysic` tries to access the core which is locked
    by another calculator.

    Parameters:

    message: string
        information describing why the error occurred
    """
    def __init__(self, message=''):
        self.message = message

    def __str__(self):
        return self.message
