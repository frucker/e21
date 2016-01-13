# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.core import Measurement
from e21.utility import merge_dicts
import datetime as dt
import os
import matplotlib as mpl
import quantities as pq
import numpy as np


class Loader(object):
    """A ppms file loader.

    E.g. loading a susceptibility measurement::

        >>>from e21.sweet16 import Loader
        >>>s16 = Loader(mode='susceptibility')
        >>>measurement = s16('path/to/susceptibility_measurement.dat')

    """
    CREATOR = {}

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
            """
            A PPMS Data Parser. 
            
            Ignores ppms header. Converts empty data points to nan. Conservses
            all data. 
            Returns data dictionary and params dictionary.

            """

            with open(path) as f:
                
                # Skip every line until Data starts
                for line in f:
                    if line.startswith('[Data]'):
                        break
                
                for line in f:
                    if line.startswith('Comment'):
                        header = line.split(',')
                        data = {k: [] for k in header}
                    else:
                        tokens = line.strip().split(',')
                        row = {key: val for key, val in zip(header, tokens)}
                        for key, val in row.iteritems():
                            data[key].append(val)
            
            for k, v in data.iteritems():
                val = []
                if data[k] is not 'Comment':
                    for i, j in enumerate(data[k]):
                        try:
                            val.append(float(data[k][i]))
                        except ValueError:
                            val.append(float('NaN'))


                    data[k] = val
            
            # Generate Parameter Dict to be able to execute loader            
            parameter = ['info']
            params = {k: {} for k in parameter}
            params['info']['mode'] = 'BSWEEP'
            params['general']={}
            params['info']['filename'] = path.split('/')[-1]
            params['info']['filepath'] = path
                       
            return data, params

    
    
