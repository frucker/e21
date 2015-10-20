import itertools as it
import numpy as np
import quantities as pq
import matplotlib
import matplotlib.pyplot as plt
from IPython.display import display, HTML
import e21.sweet16
from operator import itemgetter
from tempfile import mkstemp
from shutil import move
from os import remove, close
import os
import pickle

def Calc_Sus_PPMS_SI(input_emu, rho_g_cm3, mass_in_mg, field_in_mT):
    input_SI = np.array(input_emu)*1e-3
    rho_SI = rho_g_cm3*1e-3/1e-6
    field_SI = field_in_mT*1e-3/(4*np.pi*1e-7)
    mass_SI = mass_in_mg*1e-6
    output=input_SI*rho_SI/mass_SI/field_SI
    return output

def Calc_Sus_PPMS_SI_volume(input_emu, volume_in_mm3, field_in_mT):
    input_SI = np.array(input_emu)*1e-3
    volume_SI = volume_in_mm3*1e-9
    field_SI = field_in_mT*1e-3/(4*np.pi*1e-7)
    output=input_SI/volume_SI/field_SI
    return output



def Felix_colormap():
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
    return cmap

def STT_colormap():
    cpool = [ '#1b1da3', '#ff0000', '#e52b2b', '#0051ff','#00868b']
    cmap = matplotlib.colors.ListedColormap(cpool[0:6], 'indexed')
    return cmap

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
    if type(scaling) is list:
        if (len(scaling) ==2):
            m = scaling[0]
            b = scaling[1]
            z = 0
        elif (len(scaling) == 3):
            m = scaling[0]
            b = scaling[1]
            z = scaling[2]
        else:
            raise ValueError('Scaling must be list of length 2 not {}'.format(len(scaling))) 
    elif type(scaling) is float or int:
        m = scaling
        b = 0
        z = 0
    else:
        raise ValueError('Scaling must be list of length 2 or float not {}'.format(type(scaling)))  
    return (data+z*pq.V)*m+b*pq.V

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
        s += ('<tr> <td><strong> Date </strong></td>'
             ' <td> {} </td></tr>'.format(Exp[i].params['info']['date']))
        s += ('<tr> <td><strong> Sample </strong></td>'
             '<td> {} </td></tr>'.format(Exp[i].params['info']['sample']))
        s += ('<tr> <td><strong> Time </strong></td>'
             '<td> {} </td></tr>'.format(Exp[i].params['info']['time']))
        s += ('<tr> <td><strong> User </strong></td>'
             '<td> {} </td></tr>'.format(Exp[i].params['info']['user']))
        s += ('<tr> <td><strong> Mode </strong></td>'
             '<td> {} </td></tr>'.format(Exp[i].params['mode']))
        s += ('<tr> <td><strong> Amplification </strong></td> <td> {} </td>'
             '</tr>'.format(Exp[i].params['general']['amplification']))
        s += ('<tr> <td><strong> Dropping Resistance </strong></td> <td>'
             ' {} k$\Omega$ </td></tr>'.format(
                 Exp[i].params['general']['dropping_resistance']))
        s += ('<tr> <td><strong> Filename </strong></td> <td> {} </td>'
             '</tr>'.format(Exp[i].params['general']['filename']))
        try: 
            if Exp[i].params['general']['contact_distance']:
                s += ('<tr> <td><strong> Contact Distance </strong></td>'
                     '<td> {} </td></tr>'.format(Exp[i].params['general']['contact_distance']))
                s += ('<tr> <td><strong> Sample Thickness </strong></td>'
                     '<td> {} </td></tr>'.format(Exp[i].params['general']['sample_thickness']))
                s += ('<tr> <td><strong> Sample Width </strong></td> <td> {} </td>'
                     '</tr>'.format(Exp[i].params['general']['sample_width']))
        except:
            pass
        s += '</table>'     
        s += '<h3> Lockin Information: </h3>'
        s += '<table class="table table-hover">' 
        keys = Exp[i].params['lock_in_1'].items()
        keys_sorted = sorted(keys, key=itemgetter(0))       
        for j in range(len(keys_sorted)):
            s += '<tr> <td><strong> {} </strong></td> <td> {} </td>'.format(keys_sorted[j][0],keys_sorted[j][1])
            if 'lock_in_2' in Exp[i].params.keys():
                s +='<td> {} </td>'.format(Exp[i].params['lock_in_2'][keys_sorted[j][0]])
            if 'lock_in_3' in Exp[i].params.keys():
                s +='<td> {} </td></tr> '.format(Exp[i].params['lock_in_3'][keys_sorted[j][0]])
            else:
                s += '</tr>' 
        s += '</table>'
        s += '<h3> Command File: </h3>'
        s += '<table class="table table-hover">'
        keys = Exp[i].params['info']['command'].items()
        keys_sorted = sorted(keys, key=itemgetter(0))
        for j in range(len(keys_sorted)):
            s += '<tr> <td><strong> {} </strong></td> <td> {} </td></tr>'.format(keys_sorted[j][0],keys_sorted[j][1])
   
        
        s = display(HTML(s))
    return g

def symmetrize(Meas, Meas2):
        """ symmetrize transport data. 

        input:
            Meas: transport measurement class, upsweep
            Meas2:  transport measurement class, downsweep
        
        output:
            list, symmetrized data        

        """

        dat = Meas.res
        rev = Meas2.res
        data = []
        for i in range(len(dat)):
            data.append(float(dat[i]+rev[i]))
        data = data*pq.ohm*pq.m
        return 0.5*data

def antisymmetrize(Meas, Meas2):
        """ symmetrize transport data. 

        input:
            Meas: transport measurement class, upsweep
            Meas2:  transport measurement class, downsweep
        
        output:
            list, symmetrized data        

        """

        dat = Meas.hall
        rev = Meas2.hall
        data = []
        for i in range(len(dat)):
            data.append(float(dat[i]-rev[i]))
        data = data*pq.ohm*pq.m
        return -1*0.5*data

def MakeOverview(Exp, *args):
        """
        Overview over all measurements including every type of measurement
        (so far B/T Sweeps)

        args options append a coulmn with the
        corresponding measurement property:

                files, current, needle valve
        """
        # Check measurement Type
        if (isinstance(Exp, e21.sweet16.susceptibility.Experiment)):
            experiment = 'susceptibility'
            Tscan = e21.sweet16.susceptibility.TemperatureScan
            Bscan = e21.sweet16.susceptibility.FieldScan
        elif (isinstance(Exp, e21.sweet16.transport.Experiment)):
            experiment = 'transport'
            Tscan = e21.sweet16.transport.TemperatureScan
            Bscan = e21.sweet16.transport.FieldScan
        else:
            experiment = 'not defined'
            Tscan = e21.sweet16.transport.TemperatureScan
            Bscan = e21.sweet16.transport.FieldScan

        # Create Table Header
        html_table = '<h1> {} - Measurement Overview </h1>' .format(
            str(Exp[0].params['info']['sample']).decode('utf-8', 'ignore')) 
        html_table += '<strong> Sample: </strong> {} <br>'.format(
            Exp[0].params['info']['sample']).decode('utf-8', 'ignore')
        html_table += ('<strong> Measurement Option:'
                      '</strong> {} <br> <br>'.format(experiment))
        html_table += '<table class="table table-hover">'
        html_table += ('<tr> <th> Number </th> <th> Sweep Type </th> <th>'
                       'Start </th> <th> Stop </th> <th> T / B </th> <th>'
                       'B_init </th> <th> T_init </th><th> Ampl </th><th>'
                       'dropping res </th><th> sweep </th>')
        if 'angle' in args:
            html_table += '<th> Angle </th>'
        if 'current' in args:
            html_table += '<th> Current </th>'
        if 'NV' in args:
            html_table += '<th> Needle Valve </th>'
        if 'filename' in args:
            html_table += '<th> Filename </th>'
        else:
            html_table += '</tr>'
        for num in range(len(Exp)):
            if (isinstance(Exp[num], Tscan)):
                T_init = float(np.round(Exp[num].temperature[0], 6))
                T_final = float(
                    np.round(Exp[num].temperature[-1], 6))
                B_field = float(
                    np.round(np.median(Exp[num].data['B_field']), 3))
                sigma = float(np.round(np.std(Exp[num].data['B_field']), 3))
                html_table += ('<tr> <td>{}</td> <td> Tsweep </td>'
                               '<td> {} K </td>'
                               '<td> {} K </td>'
                               '<td> {}+-{} T </td>'
                               '<td> {} </td>'
                               '<td> {} </td>'
                               '<td> {} </td>'
                               '<td> {} k&Omega; </td>'
                               '<td> {} </td>').format(
                                    num,
                                    T_init,
                                    T_final,
                                    B_field,
                                    sigma,
                                    Exp[num].params['info']['command']['init_field'],
                                    Exp[num].params['info']['command']['init_temperature'],
                                    Exp[num].params['general']['amplification'],
                                    Exp[num].params['general']['dropping_resistance'],
                                    Exp[num].params['info']['command']['target_temperature_rate'])
                if 'reserve' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].params['lock_in_1']['dyn_reserve'])
                if 'angle' in args:
                    html_table += '<td> {} </td>'.format(
                        round(Exp[num].data['angle'][0], 2))
                if 'current' in args:
                    html_table += '<td> {} </td>'.format(Exp[num].mean_current)
                if 'NV' in args:
                    html_table += '<td> {} </td>'.format(Exp[num].NV)             
                if 'filename' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].params['general']['filename'].split('/')[-1])
                html_table += '</tr>'   

            elif (isinstance(Exp[num], Bscan)):
                B_init = float(np.round(Exp[num].data['B_field'][0], 6))
                B_final = float(np.round(Exp[num].data['B_field'][-1], 6))
                Temp = float(
                    np.round(np.median(Exp[num].temperature), 3))
                sigma = float(
                    np.round(np.std(Exp[num].temperature), 3))
                html_table += ('<tr> <td>{}</td> <td> Field </td>'
                               '<td> {} T </td>'
                               '<td> {} T </td>'
                               '<td> {}+-{} K </td>'
                               '<td> {} </td>'
                               '<td> {} </td>'
                               '<td> {} </td>'
                               '<td> {} k&Omega; </td>'
                               '<td> {} </td>').format(
                                    num,
                                    B_init,
                                    B_final,
                                    Temp,
                                    sigma,
                                    Exp[num].params['info']['command']['init_field'],
                                    Exp[num].params['info']['command']['init_temperature'],
                                    Exp[num].params['general']['amplification'],
                                    Exp[num].params['general']['dropping_resistance'],
                                    Exp[num].params['info']['command']['target_field_rate'])
                if 'reserve' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].params['lock_in_1']['dyn_reserve'])
                if 'angle' in args:
                    html_table += '<td> {} </td>'.format(
                        round(
                            Exp[num].data['angle'][0],
                            2))
                if 'current' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].mean_current)
                if 'NV' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].NV)
                if 'filename' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].params['general']['filename'].split('/')[-1])
                html_table += '</tr>'   
            else:
                B_init = float(np.round(Exp[num].data['B_field'][0], 6))
                B_final = float(np.round(Exp[num].data['B_field'][-1], 6))
                Temp = float(
                    np.round(np.median(Exp[num].temperature), 3))
                sigma = float(
                    np.round(np.std(Exp[num].temperature), 3))
                html_table += ('<tr> <td>{}</td> <td> {} </td>'
                               '<td> {} T </td>'
                               '<td> {} T </td>'
                               '<td> {}+-{} K </td>'
                               '<td> {} </td>'
                               '<td> {} </td>'
                               '<td> {} </td>'
                               '<td> {} k&Omega; </td>'
                               '<td> {} </td>').format(
                                    num,
                                    Exp[num].params['info']['command']['mode'],
                                    B_init,
                                    B_final,
                                    Temp,
                                    sigma,
                                    Exp[num].params['info']['command']['init_field'],
                                    Exp[num].params['info']['command']['init_temperature'],
                                    Exp[num].params['general']['amplification'],
                                    Exp[num].params['general']['dropping_resistance'],
                                    Exp[num].params['info']['command']['target_field_rate'])
                if 'reserve' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].params['lock_in_1']['dyn_reserve'])
                if 'angle' in args:
                    html_table += '<td> {} </td>'.format(
                        round(
                            Exp[num].data['angle'][0],
                            2))
                if 'current' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].mean_current)
                if 'NV' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].NV)
                if 'filename' in args:
                    html_table += '<td> {} </td>'.format(
                        Exp[num].params['general']['filename'].split('/')[-1])
                html_table += '</tr>' 
        html_table += '</table>'
        return html_table

def replace(file_path, pattern, subst):
    ''' Replace string in file

    input: 
        file_path:  str, file path
        pattern:    str, pattern to be replaced
        subst:      str, substitution

    '''
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    for line in old_file:
        new_file.write(line.replace(pattern, subst))
    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

def replace_zeros(only_y = False):
    """ Replaces 0.0 with 0 for Chrisitan.
    
    E.g.::

        f, ax = figure()
        replace_zeros()

    :param f: matplotlib figure.

    """
    f = plt.gcf()
    f.canvas.draw()
    ax = plt.gca()
    xlabels = [item.get_text() for item in ax.get_xticklabels()]
    ylabels = [item.get_text() for item in ax.get_yticklabels()]
    for i, j in enumerate(xlabels):
        if not j == '':
            if float(j) == 0:
                xlabels[i] = '0'
    for i, j in enumerate(ylabels):
        if not j == '':
            if float(j) == 0:
                ylabels[i] = '0'
    if not only_y:    
        ax.set_xticklabels(xlabels)
    ax.set_yticklabels(ylabels)

def replace_zeros_legend():
    """ Replaces 0.0 with 0 for Chrisitan in legends
    
    """
    ax = plt.gca()
    lin = ax.lines
    lab = [l.get_label() for l in lin]
    for i, j in enumerate(lab):
        if j == '0.0':
            lab[i] = '0'
    plt.legend(lin, lab)
    return lin, lab

def order_measurements(Exp, meas, param = 'temp'):
    """ Orders a list of measurements after a given parameter ascending
    
    
    E.g.::
    
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'), key=os.path.getmtime))
        meas = [0,1,2,3]
        order_measurements(MnSi_T, meas, param = 'temp')
    
    :param dict Exp: Sweet 16 Experiment class dictionary
    :param list meas: List of measurements
    :param str param: parameter to sort after:
                      * 'temp' = temperature, 
    
    :raises: NotImplementedError if given order parameter not implemented  
       
    """
    
    val = []
    if param == 'temp':
        for j, i in enumerate(meas):
            T = np.round(Exp[i].mean_temperature,2)
            val.append((T,i))
        meas = sorted(val, key=lambda x: x[0])  
        meas = [meas[i][1] for i in range(len(meas))]#
    elif param == 'field':
        for j, i in enumerate(meas):
            B = np.round(Exp[i].mean_field,3)
            val.append((B,i))
        meas = sorted(val, key=lambda x: x[0])  
        meas = [meas[i][1] for i in range(len(meas))]#
    else:
        raise NotImplementedError('Method not implemented yet: {}'.format(param))
        
    return meas

def save(Obj, path, name, binary = False):
    if not os.path.exists(path):
        os.makedirs(path)
    with open(path+name, 'wb') as f:
        if binary:
            pickle.dump(Obj, f, pickle.HIGHEST_PROTOCOL)
        else:
            pickle.dump(Obj, f)

def read(path, name):
    with open(path+name, 'rb') as f:
        return pickle.load(f)
        

