# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.

import collections
import datetime as dt
import itertools as it
import os
import warnings

import quantities as pq
import numpy as np

import e21.core

#-Loader-----------------------------------------------------------------------
def _to_quantity(x):
    if not isinstance(x, pq.Quantity):
        val, unit = x
        x = pq.Quantity(val, unit)
    return x


class ParsingError(Exception):
    """Raised by the Loader classes when a parsing error occures."""


class Loader(object):
    def __init__(self, **kw):
        if 'm_ref' in kw:
            kw['m_ref'] = _to_quantity(kw['m_ref'])
        if 'U_ref' in kw:
            kw['U_ref'] = _to_quantity(kw['U_ref'])
        self.params = kw

    def __call__(self, path, **kw):
        data, file_params = self.parse(path)
        params = dict(self.params)
        params.update(kw)
        params.update(file_params)
        if 'Tset' in data and 'H' in data:
            # > 0 if target temperature varies
            dT = np.diff(data['Tset']).any()
            # > 0 if target field varies
            dH = np.diff(data['H']).any()
            if dT and not dH:
                return TScan(data, params)
            if not dT and dH:
                return BScan(data, params)
        return Magnetisation(data, params)

    def parse(self, path):
        """Importer for VSM datafiles.

        The VSM files have the following syntax:

            COMMENTS
            DATETIME

            COLNAMES
            COLUNITS

            PRECISSION
            DATA

        `COMMENTS` start with a `;` character. Multiple `COMMENTS` per line are
        possible. The datetime has the following syntax:
        `day month date hh:mm:ss yyyy` e.g. `Wed May 23 16:11:30 2012`.

        """
        def assert_blank(line):
            """Tiny helper function, checking for a blank line."""
            if line.strip():
                warnings.warn('Line is not blank:{0}'.format(line.strip()))

        params = {
            'path': path,
            'filename': os.path.split(path)[1],
        }

        with open(path, 'rb') as f:
            for line in iter(f.readline, ''):
                if line.startswith(';'):
                    # TODO parse comments
                    continue
                else:
                    params['datetime'] = self._parse_datetime(line)
                    break

            assert_blank(f.readline())
            variables = [x.strip() for x in f.readline().split(',')]
            units = [x.strip() for x in f.readline().split(',')]
            assert_blank(f.readline())

            precission = f.readline()  # The precission isn't used.
            M = np.genfromtxt(f)
        if not (len(variables) == len(units) == M.shape[1]):
            raise ParsingError('Column mismatch')

        data = {}
        for varname, col, unit in it.izip(variables, M.T, units):
            data[varname] = pq.Quantity(col, unit)
        return data, params

    def _parse_datetime(self, x):
        try:
            return dt.datetime.strptime(x.strip(), '%a %b %d %H:%M:%S %Y')
        except ValueError:
            raise ParsingError('Invalid Date.')

#-Measurement-classes----------------------------------------------------------
class Magnetisation(e21.core.Measurement):
    @property
    def x(self):
        return self.data['X']

    def y(self):
        return self.data['Y']

    @property
    def m(self):
        """Calculates the magnetic moment of the X signal."""
        m_ref = self.params['m_ref']
        U_ref = self.params['U_ref']
        return self.data['X'] * m_ref / U_ref

    @property
    def field(self):
        return self.data['H']

    @property
    def control_temperature(self):
        return self.data['Tcontrol']

    @property
    def target_temperature(self):
        return self.data['Tset']

    @property
    def temperature(self):
        return self.data['Tsample']

    def temperature_stability(self):
        dT = self.data['Tcontrol'] - self._data['Tset']
        return np.mean(dT), np.std(dT)


class FieldScan(Magnetisation, e21.core.Plottable):
    def plot(self, y='m', x='field', axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        subplot_default = {
            'xlabel': 'B ({0})'.format(x.dimensionality),
            'ylabel': 'm ({0})'.format(y.dimensionality),
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', str(np.mean(self.control_temperature)))
        return super(FieldScan, self).plot(y, x, axes, subplot_default, fig_kw, label, **kw)


class TemperatureScan(Magnetisation, e21.core.Plottable):
    def plot(self, y='m', x='temperature', axes=None, subplot_kw={}, fig_kw={}, **kw):
        """A default implementation for a generic fieldscan plot."""
        subplot_default = {
            'xlabel': 'T ({0})'.format(x.dimensionality),
            'ylabel': 'm ({0})'.format(y.dimensionality),
        }
        subplot_default.update(subplot_kw)
        label = kw.pop('label', '{0}'.format(self.field[0]))
        return super(TemperatureScan, self).plot(y, x, axes, subplot_default, fig_kw, label, **kw)

