# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.

import quantities as pq
import numpy as np

import datetime as dt
import itertools as it
import os

import e21.core
from e21.core import env


def _to_quantity(x):
    if not isinstance(x, pq.Quantity):
        val, unit = x
        x = pq.Quantity(val, unit)
    return x


class ParsingError(Exception):
    pass


class ListView(object):
    def __init__(self, source, indices=None):
        self._source = source
        self._indices = indices or []

    def __len__(self):
        return len(self._indices)

    def __getitem__(self, item):
        return self._source[self._indices[item]]

    def __setitem__(self, item, value):
        self._source[self._indices[item]] = value

    def append(self, value):
        idx = len(self._source)
        self._source.append(value)
        self._indices.append(idx)


class Experiment(e21.core.Experiment):
    def __init__(self, date, description=None, measurements=None):
        super(Experiment, self).__init__(date, description)
        self._measurements = []
        self.field_scans = ListView(self._measurements)
        self.temperature_scans = ListView(self._measurements)

    def __getitem__(self, item):
        return self._measurements[item]

    def append(self, measurement):
        if isinstance(measurement, TScan):
            self.temperature_scans.append(measurement)
        elif isinstance(measurement, BScan):
            self.field_scans.append(measurement)
        else:
            self._measurements.append(measurement)

    def _repr_html_(self):
        template = e21.core.env.get_template('vsm.Experiment.html')
        return template.render(obj=self)


class MagnetisationMeasurement(e21.core.Measurement):
    @property
    def x(self):
        return self._data['X']

    def y(self):
        return self._data['Y']

    @property
    def m(self):
        """Calculates the magnetic moment of the X signal."""
        m_ref = self.params['m_ref']
        U_ref = self.params['U_ref']
        return self._data['X'] * m_ref / U_ref

    @property
    def field(self):
        return self._data['H']

    @property
    def control_temperature(self):
        return self._data['Tcontrol']

    @property
    def target_temperature(self):
        return self._data['Tset']

    def temperature_stability(self):
        dT = self._data['Tcontrol'] - self._data['Tset']
        return np.mean(dT), np.std(dT)


class BScan(MagnetisationMeasurement):
    def __init__(self, datas, params):
        super(BScan, self).__init__(datas, params)


class TScan(MagnetisationMeasurement):
    def __init__(self, datas, params):
        super(TScan, self).__init__(datas, params)


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

        return MagnetisationMeasurement(data, params)

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
                raise ParsingError('Line is not blank:{0}'.format(line.strip()))

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
