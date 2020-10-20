"""
This module defines a Particle class used to define the properties
of an element used in Synchrotrons. The particle can be composed of any
number of protons, neutrons and electrons.

@copyright: CERN
@date: 08/03/2018
@author: David Amorim
"""

from scipy.constants import (speed_of_light, elementary_charge,
                             proton_mass, neutron_mass, electron_mass)


class Particle(object):
    """ Define a particle

    It is defined by its name and its components

    The object contains the fundamental constants e and c

    A method allows to return the rest energy in MeV

    Attributes
    ----------
    e, c: floats
        elementary charge and speed of light, from scipy.constants
    species: str
        particle name
    A: int or float
        total number of nucleons, can be decimal for standard
        atomic weight
    Z: int
        number of protons
    C: int
        electrical charge
    mass: float
        total particle mass
    E_rest: float
        particle rest energy in eV

    Example
    -------
    import particle_parameters
    16_8_O_2 = particle_parameters.Particle('O', 16, 8, 2)
    239_94_Pu_2- = particle_parameters.Particle('Pu', 239, 94, -2)
    """
    def __init__(self, species, A, Z, C):
        self.e = elementary_charge
        self.c = speed_of_light
        self.species = species
        self.A = A
        self.Z = Z
        self.C = C

    @property
    def mass(self):
        """
        Compute a particle mass from its components:
            - number of nucleons A
            - number of protons Z
            - number of electron C
        """
        particle_mass = ((self.A-self.Z)*neutron_mass
                         + self.Z*(proton_mass + electron_mass)
                         - self.C*electron_mass)
        return particle_mass

    @property
    def charge(self):
        return self.C * self.e

    @property
    def name(self):
        """ Create the particle full name from its properties and ist name

        A       C
         element
        Z
        """
        return '{:.2f}_{:d}_{:s}_{:d}'.format(self.A, self.Z,
                                              self.species, self.C)

    @property
    def E_rest(self):
        return self.mass * self.c ** 2 / self.e


class Proton(Particle):
    """ Define a proton particle

    Example
    -------
    import particle_parameters
    proton = particle_parameters.Proton()
    """
    def __init__(self):
        Particle.__init__(self, 'proton', 1, 1, 1)


class Electron(Particle):
    """ Define an electron particle

    Example
    -------
    import particle_parameters
    electron = particle_parameters.Electron()
    """
    def __init__(self):
        Particle.__init__(self, 'electron', 0, 0, -1)


class Pb54(Particle):
    """ Define an Pb54 particle

    The number of nucleons is the standard atomic weight [1]

    References
    ----------
    [1] "Atomic weights of the elements 2013
    (IUPAC Technical Report)", J.Meija et al.,
    Pure and Applied Chemistry. 88 (3): 265-91

    Example
    -------
    import particle_parameters
    Pb54 = particle_parameters.Pb54()
    """
    def __init__(self):
        Particle.__init__(self, 'Pb54', 207.2, 82, 54)

