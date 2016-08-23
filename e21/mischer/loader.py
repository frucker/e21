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
    """A sweet16 file loader.

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
        params = e21.utility.merge_dicts(kw, self.kw, params)
        try:
            return Loader.CREATOR[params['mode']](data, params)
        except KeyError:
            return Measurement(data, params)

    def parse(self, path):
        """
        Parser for mischer measurement files. Only works for files younger
        than september 2013. Header files start with '#' and are parsed in
        blocks. Each block starts with '#'+block_title+''.

        Returns two dicitonaries. Header information is stored in 'params'
        dicitonary, data is stored in 'data' dictionary.

        Converts the first two columns of sweet 16 data files (date, time)
        into two columns, 'datetime' in datetime format, and 'time_numpy'
        in float, which can be plotted directly by matplotlib.

        Unit information is restored from data file and added to dictionary
        via quantities package.

        Input parameters:

            path: 'string'         file path to valid sweet 16 file

        Output parameters:

            data: 'dict'           dictionary containing data
            params: 'dict'         dictionary containing header information
                                   as blocks (each block a dictionary)

        """
        with open(path) as f:
            block_title = 'info'
            units = []
            variables = []
            params = {}
            params[block_title] = {}
            linenum = 0
            for line in f:
                linenum = linenum + 1
                if line.startswith('#'):
                    if line[1] is not '#':
                        params, block_title = parse_header(line[1:],
                                                            params,
                                                            block_title)
                # splitting variables and units
                elif line.startswith('datetime'):
                    temp = line.strip().split('\t')
                    for i in temp:
                        try:
                            tok = []
                            tok = i.split('(')
                            variables.append(tok[0].strip('(').strip())
                            units.append(tok[1].strip(')'))
                        except IndexError:
                            #variables.append(tok[0])
                            units.append('NaN')

                    data = {k:[] for k in variables[0:]}
                    data['time_numpy'] = []

                    # TODO: Es tritt ein Fehler auf da Ohm am ende in Klammern
                    # The Sweet 16 uses 'Deg', quantities uses 'deg'.
                    units = ['deg' if u == 'Deg' else u for u in units]
                    units = ['ohm' if u == 'Ohm' else u for u in units]

                elif line.strip():
                    row = self.parse_data(line, variables, path, linenum)
                    for key, val in row.iteritems():
                        data[key].append(val)
            # convert python datetime to plottable float
            # (matplotlib: plot_date function) and store in data dictionary
            data['time_numpy'] = mpl.dates.date2num(data['datetime'])

            # each device gets its unit (except for...)
            for col, unit in zip(variables, units):
                if col not in ['datetime', 'capacity', 'loss']:
                    data[col] = pq.Quantity(data[col], unit)
        return data, params

    
    def parse_data(self, line, variables, path, linenum):
        warning = ''
        tokens = line.strip().split('\t')
        row = {key: float(val.replace(',','.')) for key, val in zip(variables[1:], tokens[1:])}
        # (convert date and time to one) datetime and insert in dict
        row['datetime'] = dt.datetime.strptime(tokens[0],
                                               '%d.%m.%Y %H:%M:%S,%f') #%f up to nanoseconds
        # TODO: From time to time, the LS340 creates misreadings which result
        # in either 100K or 0K readings.
        # Currently we simply ignore these lines.
        # For th future it would be better to check if a value of 100K is
        # actually an outlier an should be removed or a valid temperature.
        try:       
            if row['sample_temp_1'] == 0.:
                return {}
            elif row['sample_temp_1'] == 100.:
                return {}
        
            # check if number of data headers fits to number of data columns
            # (happens e.g. at last line after measurement abort)
            elif (len(tokens) != len(variables)):
                while len(tokens) < len(variables):
                    tokens.append(np.NaN)
            # TODO: Should raise a warning (per file) when inserting NaNs 
        except KeyError:
            pass
        return row

def parse_commands(arguments):
    """ Parses the command line of a mischer measurement file and returns
        a dicitionary of command string.

        Input parameters:

            arguments:     list, list of command string entries,
                           (like one line in mischer command_file.txt)

        """
    command_line = ['mode', 'init_temperature', 'init_temperature_rate',
                    'target_temperature', 'target_temperature_rate',
                    'init_field', 'init_field_rate', 'target_field',
                    'target_field_rate', 'lockin1_sensitivity',
                    'lockin2_sensitivity', 'lockin3_sensitivity', 'delay',
                    'steps', 'target_current', 'current_rate']
    command_line_units = ['', ' K', ' K/min', ' K', ' K/min', ' T',
                          ' T/min', ' T', ' T/min', '', '', '',
                          ' min', '', ' A', ' A/min']
    args = []
    for i in range(len(arguments.split('\t'))):
        args.append(arguments.split('\t')[i] + command_line_units[i])
    modes = {'1': 'BSWEEP', '2': 'TSWEEP', '3': 'CONST', '4': 'BSTEP',
             '6': 'TSTEP', 'Tramp': 'TSWEEP', 'Bramp': 'BSWEEP',
             '7': 'ASWEEP', '8':'FSWEEP','Bstep':'BSTEP', '10':'ASTEP', 'Const':'CONST'}
    # Parse sensitivities
    sens = {'0':'2e-9 V','00': '2e-9 V', '01': '5e-9 V', '02': '1e-8 V',
            '03': '2e-8 V', '04': '5e-8 V', '05': '1e-7 V',
            '06': '2e-7 V', '07': '5e-7 V', '08': '1e-6 V',
            '09': '2e-6 V', '0':'2e-9 V','1': '2e-9 V', '2': '1e-8 V',
            '3': '2e-8 V', '4': '5e-8 V', '5': '1e-7 V',
            '6': '2e-7 V', '7': '5e-7 V', '8': '1e-6 V',
            '9': '2e-6 V', '10': '5e-6 V', '11': '1e-5 V',
            '12': '2e-5 V', '13': '5e-5 V', '14': '1e-4 V',
            '15': '2e-4 V', '16': '5e-4 V', '17': '1e-3 V',
            '18': '2e-3 V', '19': '5e-3 V', '20': '0.01 V',
            '21': '0.02 V', '22': '0.05 V', '23': '0.1 V',
            '24': '0.2 V', '25': '0.5 V', '26': '1 V'}
    args = dict(zip(command_line, args))
    # Convert Lockin Sensitivity Units
    args['lockin1_sensitivity'] = sens[args['lockin1_sensitivity']]
    args['lockin2_sensitivity'] = sens[args['lockin2_sensitivity']]
    args['lockin3_sensitivity'] = sens[args['lockin3_sensitivity']]
    args['mode'] = modes[args['mode']]
    return args

def parse_header(line, params, block_title):
        """
        Parses header of sweet 16 files. Function is called within parse,
        but can be called seperately

        input parameters:

            line:       'string'        line of data file
            params:     'dict'          dictionary, empty dictionary containing
                                        only one key, value pair:
                                        params[block_title] = {}
            block_title:'string'        initial block tile

        output:

            params:      'dict'         containing header information
                                        (key, value pair)of this line, stored
                                        in params[block_title][kw] = argument

            block_title  'string'       new block title if line defines new
                                        block in data file,
                                        i.e. '#'+block_title+'' (no argument)
        """
        # True, wenn Zeile nicht nur aus leerzeichen besteht.
        if line.strip():
            try:
                line = line.split(':')
                kw, arg = line[0].strip(), line[1].strip()                

                # check if block title or information
                if (arg == ''):
                    block_title = line[0].strip()
                    params[block_title] = {}
                # check if time information
                # (containing an additional ':', thus treated seperately)
                if (kw == 'time'):
                    arg = ':'.join((line[1].strip(), line[2].strip()))
                # check if 'command' line. Then argument is dictionary containing
                # infromation of command line of command_file.txt of measurement
                if (kw == 'command'):
                    arg = parse_commands(arg)
                # Corret bad data file output (files befor 08/14)
                if (kw == 'offset.'):
                    kw = 'offset'
                # if not block title, then store information
                if not (arg == ''):
                    params[block_title][kw] = arg
            except IndexError:
                params[line[0].replace(' ','_').strip()] = {}
                block_title = line[0].replace(' ','_').strip()
        return params, block_title



