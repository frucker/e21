# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.core import Measurement
import e21.utility 
import datetime as dt
import os
import matplotlib as mpl
import quantities as pq
import numpy as np

class Loader(object):
    """A mira file loader.

    E.g. loading a mira measurement::

        >>>from e21.mira import Loader
        >>>s16 = Loader(mode='mira')
        >>>measurement = s16('path/to/???.dat')

    """
    CREATOR = {}

    def __init__(self, **kw):
        self.kw = kw

    def __call__(self, path, **kw):
        data, params = self.parse(path)
        params = e21.utility.merge_dicts(kw, self.kw, params)
        try:
            return Loader.CREATOR[params['mode']](data, params)
        except KeyError:
            return Measurement(data, params)

    def parse(self,path):
        """
        Parser for mira measurement files. Header files start with '#' and are parsed in
        blocks. Each block starts with '###'+block_title+''.

        Returns two dicitonaries. Header information is stored in 'params'
        dicitonary, data is stored in 'data' dictionary.

        Unit information is restored from data file and added to dictionary
        via quantities package.

        Input parameters:

            path: 'string'         file path to valid mira file

        Output parameters:

            data: 'dict'           dictionary containing data
            params: 'dict'         dictionary containing header information
                                   as blocks (each block a dictionary)

        """
        with open(path) as f:
            units = []
            variables = []
            params = {}
            linenum = 0
            for line in f:
                linenum = linenum + 1
                if line.startswith('### NICOS data file'):
                    block_title = 'general'
                    params[block_title] = {}
                    params[block_title]['date'] = line.split(' at ')[1]
                elif line.startswith('###'):
                    block_title = self.parse_title(line).strip()
                    params[block_title] = {}
                elif line.startswith('# '):
                    if block_title == 'Scan data':
                        key, value, variables = self.parse_units(line, params)
                        data = {k:[] for k in variables}
                    else:
                        key, value = self.parse_params(line)
                    params[block_title][key] = value
                else:
                    row = self.parse_data(line, variables)
                    for key, val in row.iteritems():
                        data[key].append(val)

        try: 
            if data:
                pass 
        except UnboundLocalError:
            data={}
        return data, params

    def parse_title(self,line):
        '''

        '''
        key = line.strip('###').lstrip()
        if key.startswith('Device positions and sample environment state'):
            key = 'devices'
        elif key.startswith('Devices'):
            key = 'error'
        return key

    def parse_units(self,line, params):
        try:
            if params['Scan data']['devices']:
                return 'units', line.strip('#').split(), params['Scan data']['devices']
        except KeyError:
            return 'devices', line.strip('#').split(), line.strip('#').split()

    def parse_params(self,line):
        '''

        '''
        key = line.split(':')[0]
        key = key.strip('#').strip()
        #customize params dict for easier access
        value = ':'.join(line.split(':')[1:])
        return key, value

    def parse_data(self,line, variables):
        tokens = line.strip().split()
        for i,t in enumerate(tokens):
            try:
                float(t)
            except:
                tokens[i] = np.NaN
        row = {key: float(val) for key, val in zip(variables, tokens)}
        return row
    












