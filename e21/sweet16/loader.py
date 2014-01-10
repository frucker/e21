# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.core import Measurement
from e21.utility import merge_dicts


class Loader(object):
    """A sweet16 file loader.

    E.g. loading a susceptibility measurement::

        >>>from e21.sweet16 import Loader
        >>>s16 = Loader(mode='Susceptibility')
        >>>measurement = s16('path/to/susceptibility_measurement.dat')

    """
    def __init__(self, **kw):
        self.kw = kw

    def __call__(self, path, **kw):
        data, params = self.parse(path)
        params = merge_dicts(kw, self.kw, params)
        try:
            return Loader.CREATOR[params['mode']](data, params)
        except KeyError:
            return Measurement(data, params)

    def parse(self, path):
        raise NotImplementedError()
