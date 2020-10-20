"""
This module defines a Synchrotron class which includes the parameters
useful for stability simulations. Some basic parameters such as the beam
energy or the RF voltage are loaded from a YAML configuration file
and used to compute all the needed parameters such as the synchrotron tune
or the revolution frequency.

Some of these parameters can be modified. The values derived from this will
then be recomputed accordingly.

@copyright: CERN
@date: 08/03/2018
@author: David Amorim
"""

from __future__ import division

import warnings
import yaml
import numpy as np

import particle_parameters


class Synchrotron(object):
    """
    Define the parameters of a synchrotron

    The parameters can be loaded from a YAML file. If no file is provided,
    a default file is used.

    Methods are defined to compute various parameters from the ones in the file

    Parameters
    ----------
    parameter_file: str, default './beams/Basic_synchrotron.yaml'
        Path to the YAML file containing the beam and ring parameters
        A default synchrotron based on the PS is provided

    Attributes
    ----------
    circumference, radius: floats
        Accelerator circumference and radius. One of the value is taken from
        the YAML file and the other one is computed from it. Can be
        overwritten by the user.
    particle: object of class Particle from particle_parameters
        Object storing the particle properties and physical constants
    E_kinetic, gamma: floats
        Beam kinetic energy and Lorentz factor. The kinetic energy is given
        in MeV. One of the value is taken from the YAML file and the other one
        is computed from it. They can be overwritten by the user.
    f0, omega0: float
        Beam revolution frequency and angular revolution frequency derived
        from the beam energy and the machine circumference.
    Qs: float
        Synchrotron tune. Use the value given in the YAML file if
        it is present. Otherwise it is computed from RF voltage, harmonic
        number, beam kinetic parameters, synchronous phase and slippage factor.
        It can be overwritten by the user.
    omegas: float
        Synchrotron angular frequency. Computed from Qs and beam angluar
        revolution frequency.
    alphap: float
        Momentum compaction factor of the machine. It is taken from the
        YAML file.
    Qx, Qy: floats
        Transverse betatron tunes taken from the YAML file. They can be
        overwritten by the user.
    Qxfrac, Qyfrac: floats
        Fractionnal part of the transverse betatron tunes computed from
        the transverse tunes

    Example
    -------
    >>> import synchrotron
    >>> SPS = synchrotron.Synchrotron(
                 parameter_file=('./machine_configuration/
                                 'SPS/SPS_Q26_injection_MOSES.yaml'))
    >>> print SPS.Qs
    0.004179
    >>> SPS.Qs = 0.0043
    """
    def __init__(self,
                 parameter_file=('./machine_configuration/'
                                 'Basic_synchrotron.yaml')):

        # Warn the user if the Basic synchrotron is loaded
        if parameter_file is './machine_configuration/Basic_synchrotron.yaml':
            warnings.warn('Using the example synchrotron', UserWarning)

        self._parameters = self._read_parameters_file(parameter_file)

        self.particle = getattr(particle_parameters,
                                self._parameters['Beam Parameters']
                                                ['particle_name'])()

        self._compute_ring_geometry()
        self._compute_kinetic_parameters()
        self._set_alphap()
        self._set_synchrotron_tune()
        self._set_bunch_length()

    @property
    def circumference(self):
        return self._circumference

    @circumference.setter
    def circumference(self, value):
        self._circumference = value
        self._radius = self.circumference / (2*np.pi)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        self._radius = value
        self._circumference = 2*np.pi*self.radius

    @property
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self, value):
        """ Compute the beam Lorentz factor and derived parameters
        such as the relativistic beta, the kinetic and total energy.
        """
        self._gamma = value
        self._beta = np.sqrt(1. - 1./(self.gamma**2))
        self._p0 = (self._beta * self.gamma
                    * self.particle.mass * self.particle.c)
        self._E_kinetic = (self.gamma - 1)*self.particle.E_rest
        self._E_total = self.E_kinetic + self.particle.E_rest
        try:
            self._eta = self.alphap - 1/self.gamma**2
        except AttributeError:
            print('Momentum compaction factor not attributed yet, '
                  'pass slippage factor computation')
            pass

    @property
    def beta(self):
        return self._beta

    @property
    def E_kinetic(self):
        return self._E_kinetic

    @E_kinetic.setter
    def E_kinetic(self, value):
        E_rest = self.particle.E_rest
        self.gamma = (value + E_rest) / E_rest

    @property
    def p0(self):
        return self._p0

    @p0.setter
    def p0(self, value):
        self.gamma = np.sqrt(1 + (value / (self.particle.mass
                                           * self.particle.c)) ** 2)

    @property
    def alphap(self):
        return self._alphap

    @alphap.setter
    def alphap(self, value):
        self._alphap = value
        self._eta = self.alphap - 1/self._gamma**2

    @property
    def Qs(self):
        """ Compute the synchrotron tune from the beam parameters,
        the RF voltage, harmonic number and synchronous phase.
        If the parameter has been overwitten, the computation is ignored
        """
        if hasattr(self, '_Qs'):
            return self._Qs
        else:
            beta_rel = self._beta
            c = self.particle.c
            RF_voltage = self._parameters['Beam Parameters']['RF_voltage']
            harmonic_number = self._parameters['Beam Parameters']['harmonic']
            phi_s = self._parameters['Beam Parameters']['synchrotron_phase']

            Q_s = np.sqrt(self.particle.e
                          * RF_voltage
                          * np.abs(self._eta)
                          * harmonic_number
                          * np.cos(phi_s)
                          / (2*np.pi*beta_rel*c*self.p0))
            return Q_s

    @Qs.setter
    def Qs(self, value):
        warnings.warn('Overwriting the synchrotron tune parameters: '
                      'RF voltage will not be updated', UserWarning)
        self._Qs = value

    @property
    def f0(self):
        return self._beta * self.particle.c / self.circumference

    @property
    def omega0(self):
        return self.f0 * 2 * np.pi

    @property
    def omegas(self):
        return self.Qs * self.omega0

    @property
    def Qx(self):
        if hasattr(self, '_Qx'):
            return self._Qx
        else:
            return self._parameters['Beam Parameters']['Qx']

    @Qx.setter
    def Qx(self, value):
        warnings.warn('Overwriting the H transverse tune', UserWarning)
        self._Qx = value

    @property
    def Qy(self):
        if hasattr(self, '_Qy'):
            return self._Qy
        else:
            return self._parameters['Beam Parameters']['Qy']

    @Qy.setter
    def Qy(self, value):
        warnings.warn('Overwriting the V transverse tune', UserWarning)
        self._Qy = value

    @property
    def Qxfrac(self):
        return self.Qx - np.floor(self.Qx)

    @property
    def Qyfrac(self):
        return self.Qy - np.floor(self.Qy)

    @property
    def taub(self):
        if hasattr(self, '_taub'):
            return self._taub
        else:
            return self._parameters['Beam Parameters']['taub']

    @taub.setter
    def taub(self, value):
        self._taub = value

    @property
    def sigmaz(self):
        return self.taub * self._beta * self.particle.c / 4

    @sigmaz.setter
    def sigmaz(self, value):
        self.taub = 4 * value / (self._beta * self.particle.c)

    def _read_parameters_file(self, file_path):
        with open(file_path) as yaml_file:
            machine_parameters = yaml.load(yaml_file)
        return machine_parameters

    def _compute_ring_geometry(self):
        try:
            self.circumference = (self._parameters['Ring Parameters']
                                                  ['circumference'])
        except KeyError:
            print('No circumference parameter, using radius')
            try:
                self.radius = self._parameters['Ring Parameters']['radius']
            except KeyError:
                raise ValueError('No radius or circumference given')

    def _compute_kinetic_parameters(self):
        try:
            self.E_kinetic = self._parameters['Beam Parameters']['E_kinetic']
        except KeyError:
            print('No kinetic energy parameter, using gamma')
            try:
                self.gamma = self._parameters['Beam Parameters']['gamma']
            except KeyError:
                raise ValueError('No kinetic energy or gamma given')

    def _set_alphap(self):
        """Check if momentum compaction factor is present at initialization"""
        try:
            self.alphap = self._parameters['Beam Parameters']['alphap']
        except KeyError:
            raise ValueError('No momentum compaction factor found')

    def _set_synchrotron_tune(self):
        """Set the synchrotron tune at the object initialization."""
        try:
            self.Qs = self._parameters['Beam Parameters']['Qs']
        except KeyError:
            print('Synchrotron tune will be computed from RF voltage, '
                  'harmonic number...')
            pass

    def _set_bunch_length(self):
        """Set the bunch length at the object initialization."""
        try:
            self.taub = self._parameters['Beam Parameters']['taub']
        except KeyError:
            print('No full bunch length (in seconds) given in '
                  'the paramters file, trying with '
                  'RMS bunch length (in meters)')
            try:
                self.sigmaz = self._parameters['Beam Parameters']['sigmaz']
            except KeyError:
                raise ValueError('No full bunch length (in seconds) or '
                                 'RMS bunch length (in meters)'
                                 'given in parameters')
