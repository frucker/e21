import itertools as it
import numpy as np
import quantities as pq
import matplotlib
from IPython.display import display, HTML
# colors taken from bootstrap
COLORS = {
    'black': '#000000',
    'darker gray': '#222222',
    'dark gray': '#333333',  #text color
    'gray': '#555555',
    'light gray': '#999999',
    'lighter gray': '#eeeeee',
    'white': '#ffffff',
    'blue': '#049cdb',
    'dark blue': '#0064cd',
    'light blue': '#3a87ad',
    'lighter blue': '#d9edf7',
    'green': '#46a546',
    'light green': '#468847',
    'lighter green': '#dff0d8',
    'red': '#9d261d',
    'light red': '#b94a48',
    'lighter red': '#f2dede',
    'yellow': '#ffc40d',
    'light yellow': '#fcf8e3',
    'orange': '#f89406',
    'pink': '#c3325f',
    'purple': '#7a43b6',
}

def Felix_colormap(n):
    cdict = {'red': ((0.0, 1.0, 1.0),
                     (0.7, 1.0, 1.0),
                     (1.0, 0.24, 0.24)),
         'green': ((0.0, 0.8, 0.8),
                   (0.7, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
         'blue': ((0.0, 0.0, 0.0),
                  (0.7, 0.0, 0.0), 
                  (1.0, 1.0, 1.0))}
    cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap(n)

def merge_dicts(*dicts):
    """Merges multiple dicts together.
    
    ::
        >>> x = dict(a=1, b=1)
        >>> y = dict(b=2)
        >>> merge_dicts(x, y, {'c': 1})
        {'a': 1, 'b': 2, 'c': 1}
    
    """
    return dict(it.chain(*[x.iteritems() for x in dicts]))



def correct_temp(temps,T):
    """
    correct temperatures (temps) for constant offset (T)
    
    """
    return temps+T*pq.K
    
def calibrate(data,scaling):

    """
    calibrate data by multiplying a scaling factor
    
    """
    return data*scaling

def export_FS(data, numbers, path, option = 'susceptibility',
              scailing = 1, t_off = 0):
    """
    Export Field Sweeps to txt files. 
    Outputfiles are stored in 'path' directory.

    Options: 

        Susceptibility: - One .txt file for each FS, 
                          named after mean temperature: "FSxxxxK"
                          Temperature given in 4 digits:  4.38 K  --> 0483
                                                          12.80 K --> 1280    
                                                          320.1 K --> 3201
                          File header: Temperature Field SUS SUSIM
                        - One .txt file containing all FS, horizontally stacked

    Input parameters:

        data:       Experiment/Susceptibility Class Data Set
        numbers:    list of desired measurements for output
        path:       output folder (must exist)

    optional:        

        options:    - 'susceptibility' (default) for susceptiblity measurements    
        scailing:   - scailing number to calibrate susceptibility
        t_off:      - temperature offset    
    """
    if (option == 'susceptibility'):

        for i in numbers:
            temp = '{:0>4d}'.format(
                   int(correct_temp(data[i].mean_temperature,t_off)*100))[:4]
            temp = temp.replace('.','')
            f = open(path+'FS{}K.txt'.format(temp),'w')
            f.write('Temperature\tField\tSUS\tSUSIM\n')
            for j in range(len(data[i].field)):
                f.write('{}\t{}\t{}\t{}\n'.format(
                         float(correct_temp(data[i].temperature[j],t_off)),
                         float(data[i].field[j]), 
                         float(calibrate(data[i].real[j],scailing)),
                         float(calibrate(data[i].imag[j],scailing))))
            f.close()
        f = open(path+'FS_gesamt.txt','w')
        f.write('Temperature\tField\tSUS\tSUSIM\n')
        temps = []
        fields = []
        sus = []
        susim = []
        for i in numbers:
            T = correct_temp(data[i].temperature,t_off)
            B = data[i].field
            S = calibrate(data[i].real,scailing)
            SI = calibrate(data[i].imag,scailing)
            temps = np.hstack([temps, T])
            fields = np.hstack([fields, B])
            sus = np.hstack([sus, S])
            susim = np.hstack([susim, SI])
        for i in range(len(temps)):
            f.write('{}\t{}\t{}\t{}\n'.format(float(temps[i]),
                    float(fields[i]),float(sus[i]),
                    float(susim[i])))
        f.close()
    else:
        print "{} option not implemented yet. Do nothing".format(option)
        
def measurement_details(Exp):
    def g(i = 0):
        s = '<h2> Measurement Details for Measurement Number {}</h2>'.format(i)
        s += '<h3> General information: </h3>'
        s += '<table class="table table-hover">'
        s += '<tr> <td><strong> Date </strong></td> <td> {} </td></tr>'.format(Exp[i].params['info']['date'])
        s += '<tr> <td><strong> Sample </strong></td> <td> {} </td></tr>'.format(Exp[i].params['info']['sample'])
        s += '<tr> <td><strong> Time </strong></td> <td> {} </td></tr>'.format(Exp[i].params['info']['time'])
        s += '<tr> <td><strong> User </strong></td> <td> {} </td></tr>'.format(Exp[i].params['info']['user'])
        s += '<tr> <td><strong> Mode </strong></td> <td> {} </td></tr>'.format(Exp[i].params['mode'])
        s += '<tr> <td><strong> Amplification </strong></td> <td> {} </td></tr>'.format(Exp[i].params['general']['amplification'])
        s += '<tr> <td><strong> Dropping Resistance </strong></td> <td> {} k$\Omega$ </td></tr>'.format(Exp[i].params['general']['dropping_resistance'])
        s += '<tr> <td><strong> Filename </strong></td> <td> {} </td></tr>'.format(Exp[i].params['general']['filename'])
        if Exp[i].params['contact_distance']:
            s += '<tr> <td><strong> Contact Distance </strong></td> <td> {} </td></tr>'.format(Exp[i].params['contact_distance'])
            s += '<tr> <td><strong> Sample Thickness </strong></td> <td> {} </td></tr>'.format(Exp[i].params['sample_thickness'])
            s += '<tr> <td><strong> Sample Width </strong></td> <td> {} </td></tr></table>'.format(Exp[i].params['sample_width'])
        else:
            s += '</table>'
        
        s += '<h3> Lockin Information: </h3>'
        s += '<table class="table table-hover">' 
        for k,v in Exp[i].params['lock_in_1'].items():
            s += '<tr> <td><strong> {} </strong></td> <td> {} </td>'.format(k, v)
            if 'lock_in_2' in Exp[i].params.keys():
                s +='<td> {} </td>'.format(Exp[i].params['lock_in_2'][k])
            if 'lock_in_3' in Exp[i].params.keys():
                s +='<td> {} </td></tr> '.format(Exp[i].params['lock_in_3'][k])
            else:
                s += '</tr>' 
        s += '</table>'
        s += '<h3> Command File: </h3>'
        s += '<table class="table table-hover">'
        for k,v in Exp[i].params['info']['command'].items():
            s += '<tr> <td><strong> {} </strong></td> <td> {} </td></tr>'.format(k,v)
   
        
        s = display(HTML(s))
    return g
