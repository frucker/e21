# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
import e21.core
from e21.core import lookup
import numpy as np
from e21.Nicos import Loader
import quantities as pq
from e21.Nicos.Mira.core import Mira
from operator import itemgetter
from IPython.display import display, HTML
from e21.utility import ProgressBar



def _correct_field(field, correction):
    return field - correction * pq.T

# Susceptibility-classes ------------------------------------------------------

    

class miraHe3(Mira):

    @property
    def experiment_type(self):
        return 'Mira He3'

    @property
    def om_offset(self):
        return self.params['Offsets']['om_offset']

    


def check_empty_files(filelist, mode = 'susceptibility'):
    """ Checks for empty files in list of files.
        Empty files are files containing only header lines.
        input: filelist
        output: filelist, list of non-empty files
    """ 

    s16 = Loader(mode=mode)
    fl = []
    for n, files in enumerate(filelist):
        Test = s16(files) 
        try:
            if Test.data['datetime'][0]:
                fl.append(files)
        except IndexError:
            print 'file empty: {}'.format(files)
    return fl

def create(data, params):
    """Takes data and param and creates a FieldScan or TemperatureScan."""
    return miraHe3(data, params)




class Experiment(e21.core.Experiment):

    @property
    def temperature_scans(self):
        return [
            x for x in sorted(
                self._measurements,
                key=lambda x: x.field[0]) if isinstance(
                x,
                TemperatureScan)]

    @property
    def field_scans(self):
        return [
            x for x in sorted(
                self._measurements,
                key=lambda x: x.mean_temperature) if isinstance(
                x,
                FieldScan)]          

    def add_measurements(self, filelist):
        """ Adds measurements to Experiment

            filelist: list of measurement files
        """
        p = ProgressBar(len(filelist))
        s16 = Loader(mode='miraHe3')
        meas_len = len(self._measurements)
        for n, files in enumerate(filelist):
            p.animate(n+1)
            index = n + meas_len
            self._measurements[index] = s16(files)

def measurement_details(Exp, key='devices'):
    '''
    returns a function g(i) which can be used with ipython.widges.interact

    input parameters:

        Exp: e21 Experiment class dicitonary

    output:
    
        g(i): function that displays a HTML table containing comprehensive
              measurement information on Measurement Exp[i]

    ''' 


    def g(i = 0):
        s = '<h2> Measurement Details for Measurement Number {}</h2>'.format(i)
        s += '<h3> General information: </h3>'
        s += '<table class="table table-hover">'
        keys = Exp[i].params[key].items()
        keys_sorted = sorted(keys, key=itemgetter(0))
        for j in range(len(keys_sorted)):
            s += '<tr> <td><strong> {} </strong></td> <td> {} </td></tr>'.format(keys_sorted[j][0],keys_sorted[j][1])  
        s += '</table>'    
        s = display(HTML(s))
        
    return g

def important_info(Exp):
    '''
    returns a function g(i) which can be used with ipython.widges.interact

    input parameters:

        Exp: e21 Experiment class dicitonary

    output:
    
        g(i): function that displays a HTML table containing comprehensive
              measurement information on Measurement Exp[i]

    ''' 


    def g(i = 0):
        s = '<h2> Measurement Details for Measurement Number {}</h2>'.format(i)
        s += '<h3> General information: </h3>'
        s += '<table class="table table-hover">'
        s += '<tr> <td><strong> Sample </strong></td> <td> {} </td></tr>'.format(Exp[i].sample) 
        s += '<tr> <td><strong> B offset </strong></td> <td> {} </td></tr>'.format(Exp[i].params['Offsets']['B_offset']) 
        s += '</table>'
        s = display(HTML(s))
    return g

def data_details(Exp):
    '''
    returns a function g(i) which can be used with ipython.widges.interact

    input parameters:

        Exp: e21 Experiment class dicitonary

    output:
    
        g(i): function that displays a HTML table containing comprehensive
              measurement information on Measurement Exp[i]

    ''' 


    def g(i = 0):
        s = '<h2> Measurement Details for Measurement Number {}</h2>'.format(i)
        s += '<h3> General information: </h3>'
        s += '<table class="table table-hover">'

        for j in range(len(Exp[i].params['Scan data']['devices'])):
            s += '<tr> <td><strong> {} </strong></td> <td> {} </td></tr>'.format(Exp[i].params['Scan data']['devices'][j], Exp[i].params['Scan data']['units'][j])  
        s += '</table>'    
        s = display(HTML(s))
        
    return g

