# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
from e21.core import Measurement
from e21.utility import merge_dicts
import datetime as dt
import os

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
        params = merge_dicts(kw, self.kw, params)
        try:
            return Loader.CREATOR[params['mode']](data, params)
        except KeyError:
            return Measurement(data, params)

    def parse(self,path):  
        """
        Parser for sweet 16 measurement files. Only works for files younger than september 2013.
        Header files start with '#' and are parsed in blocks. Each block starts with '#'+block_title+''.
        
        Returns two dicitonaries. Header information is stored in 'params' dicitonary, data is stored in 'data' dictionary.
        
        Converts the first two columns of sweet 16 data files (date, time) into two columns, 'datetime' in datetime format, and 'time_numpy' 
        in float, which can be plotted directly by matplotlib.
        
        Unit information is restored from data file and added to dictionary via quantities package.
        
        Input parameters:
        
            path: 'string'         file path to valid sweet 16 file
            
        Output parameters:
        
            data: 'dict'           dictionary containing data
            params: 'dict'         dictionary containing header information as blocks (each block a dictionary)
            
        """
        with open(path) as f:
            block_title = 'info'
            units = []
            variables = []
            params = {}
            params[block_title] = {}
            for line in f:
                if line.startswith('#'):
                    params, block_title = parse_header(line, params, block_title)
                elif line.startswith('date'):
                    variables = line.strip().split('\t')
                    data = { k:[] for k in variables[2:]} 
                    data['datetime'] = []
                    data['time_numpy'] = []
                elif line.startswith('(dd'):
                    units = line.strip().split('\t')
                elif line.strip():
                    data = parse_data(line, variables, data)
        return  params,data

    def parse_header(self, line, params, block_title):
        """
        Parses header of sweet 16 files. Function is called within parse, but can be called seperately
        
        input parameters:
        
            line:       'string'        line of data file
            params:     'dict'          dictionary, empty dictionary containing only one key, value pair: params[block_title] = {}
            block_title:'string'        initial block tile
        
        output:
        
            params:      'dict'         containing header information (key, value pair)of this line, 
                                        stored in params[block_title][kw] = argument
            block_title  'string'       new block title if line defines new block in data file, i.e. '#'+block_title+'' (no argument)  
        """
        
        line = line.strip('#')
        if line.strip():                          # True, wenn Zeile nicht nur aus leerzeichen besteht.   
            line = line.split(':')
            kw, arg = line[0].strip(), line[1].strip()
            # check if block title or information
            if (arg == ''):                              
                block_title = line[0].strip()
                params[block_title] = {}
            # check if time information (containing an additional ':', thus treated seperately
            if (kw == 'time'):
                arg = line[1].strip()+':'+line[2].strip()
            # check if 'command' line. Then argument is dictionary containing infromation of command line of command_file.txt of measurement
            if (kw == 'command'):
                arg = parse_commands(arg) 
            # if not block title, then store information
            if not (arg == ''):
                params[block_title][kw] = arg
        return params, block_title
        
    def parse_data(self, line,variables, data):        
        columns = line.strip().split('\t')
        #iterate over values in line and append to corresponding column in data dictionary
        for k,v in zip(variables[2:],columns[2:]):
            data[k].append(float(v))
        # convert date and time row to python datetime
        data['datetime'].append(dt.datetime.strptime(columns[0]+columns[1], '%d.%m.%Y%H:%M:%S')) # convert date and time to one datetime and insert in dict
        # convert python datetime to plottable float (matplotlib: plot_date function) and store in data dictionary
        data['time_numpy'] = matplotlib.dates.date2num(data['datetime'])   
        return data

    def parse_commands(self, arguments):
        """ Parses the command line of a Sweet 16 measurement file and returns a dicitionary of command string.
        
            Input parameters: 
            
                arguments:     list, list of command string entries, (like one line in S16 command_file.txt
            
            """            
        command_line = ['mode', 'init_temperature','init_temperature_rate', 'target_temperature','target_temperature_rate',
                        'init_field', 'init_field_rate', 'target_field', 'target_field_rate', 
                        'lockin1_sensitivity', 'lockin2_sensitivity',
                        'delay', 'needle_valve_const', 'needle_valve_percentage', 'steps', 'target_current', 'target_rate','init_angle', 'target_angle']
        args =  dict(zip(command_line,arguments.split('\t')))
        return args
