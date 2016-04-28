import e21.utility as e21u
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import quantities as pq
import matplotlib.cm as cm
from pylab import get_cmap
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter

def get_current(Exp, meas, tol = 0.01):
    """ Returns current of a list of measurement. 
      
    Exception raised if measurements have different currents.
    
    E.g.::plot_grid
        
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                key=os.path.getmtime))
        current, current_string = get_current(MnSi_T, meas)
    
    :param dict Exp: sweet16 Experiment class dictionary
    :param list meas: list of numbers of desired measurements
    :param float tol: Tolerance in which all currents have to be equal
    
    :raises: ValueError, Measuremetns don't have the same current within tol
    
    :return (float, str): current of measurements as float and string, (e.g. '0.0 A')
                          returns '0.0 A' if current atrribute not available
    
    """
    # check if measurements given
    if len(meas)==0:
        raise ValueError('No measurements given')
    
    try:
        #get current of first measurement in meas
        current = abs(Exp[meas[0]].mean_current)
        # format into string
        current_string = str(current).replace('.','_').replace(' ','')
        # check if all measurements have the same current within a tolerance
        for i, j in enumerate(meas):
            if (Exp[j].mean_current-current) > 0.01:
                raise ValueError("Current of measurement number {} not " 
                      "the same ({})"
                      "as standard: {}!".format(j,Exp[j].mean_current, current))

        return current, current_string
    except AttributeError:
        return 0, '0A'

def lookup(obj, name):
    """Helper function to provide function and string lookup of attributes.

    E.g.::

        >>>class Test(object):
        ...    pass
        >>>t = Test()
        >>>t.value = 'myvalue'
        >>>lookup(t, 'value')
        myvalue
        >>>lookup(t, lambda x: x.value)
        myvalue
        >>lookup(t, [1, 2, 3])
        [1, 2, 3]

    """
    if isinstance(name, (str,unicode)):
        return getattr(obj, name)
    elif isinstance(name, collections.Callable):
        return name(obj)
    else:
        raise ValueError('No valid data')

def make_grid(Exp, meas = [],data = 'temperature', min_y = 0, max_y = 0.6,
              scaling = [-1,0], N = 400j, **kwargs):
    
    """ Creates an extended grid for 2D plotting of Data.

    The desired data set is interpolated on a regular gird of size (NxN). X and
    Y axis of the grid temperature and magnetic field. The date can be chosen 
    freely.
    
    
    E.g.::
    
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                       key=os.path.getmtime))
        meas = [0,1,2,3]
        g, c = make_grid(MnSi_T, meas)
        Grids[c] = g
    
    :param dict Exp: sweet16 Experiment class
    :param list meas: list of measurements
    :param str,callable data: Desired Measurement data
                              (property of Susc. Measurement),
                              e.g. 'imag' or callable
    :param float min_y: min field value
    :param float max_y: max field value
    :param list scaling: scaling factors m, b (m*x+ b) of data
    :param complex N: Gridsize (NxN)
            
    :return (dict, str): Returns a dictionary containing:

        * `temps`: list of temperature values
        * 'grid':  actual grid
        * 'fields':list of field values

    """
    
    Grd = {}
    if len(meas) == 0:
        raise Exception("No measurements given!")
    
    # get current
    current, current_string = get_current(Exp, meas)
       
    # get temperature boundaries
    if not 'min_x' in kwargs.keys():
        min_x = min([np.median(Exp[i].temperature) for i in meas])
        #min_x = min_x + 0.05*min_x
    else:
        min_x = kwargs['min_x']*pq.K
    
    if not 'max_x' in kwargs.keys():
        max_x = max([np.median(Exp[i].temperature) for i in meas])
        #max_x = max_x - 0.05*max_x
    else:
        max_x = kwargs['max_x']*pq.K
    #generate grid
    grid_x, grid_y = np.mgrid[min_x:max_x:N, min_y:max_y:N]
    # x coordinates (list)
    x = []
    # y cooordinates
    y = []
    # values
    values = []
    for i in meas: 
        x1 = Exp[i].temperature[0::1]
        y1 = Exp[i].field[0::1]
        val = e21u.calibrate(lookup(Exp[i],data), scaling)[0::1]
        x = np.hstack((x,x1))
        y = np.hstack((y,y1))
        values = np.hstack((values, val))
    #scatter(x,y)
    # scale x and y scale for interpolation. Else bad interpolation
    # due do difference in absolute values
    xscale = max_x - min_x
    yscale = max_y - min_y
    scale = np.array([xscale, yscale])
    punkte = np.vstack((x,y)).T
    grid = griddata(punkte/scale, values, (grid_x/xscale, grid_y/yscale),
                    method = 'linear', fill_value = 0)
    
    Grd['temps'] = np.linspace(min_x, max_x, int(N.imag))
    Grd['grid'] = grid
    Grd['fields'] = np.linspace(min_y, max_y, int(N.imag))
    
    return Grd, '{}'.format(current)

def plot_grid(Grids, current, vmin = 0, vmax = 1e-6, nbins = 4, **kwargs): 
       
    """ Makes a 2D Contourplot of a given Grid. The input grid has 
        to be stored in a dicitionary (result of make_grid)
    
    
    E.g.::
    
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                       key=os.path.getmtime))
        meas = [0,1,2,3]
        g, c = make_grid_real(MnSi_T, meas)
        Grids[c] = g
        current = 0.0 # enter your measurement current here
        plot_grid(Grids, current)
    
    :param dict Grids: Dicitonary containing a grid of data as well
                       as x and y axis values (result of make_grid)
    :param float current: current of measurements to plot
    :param float vmin: min susceptiblity value
    :param float vamx: max susceptiblity value
    :param float nbins: number of bins
    
    """
    # get grid information
    grid = Grids['{}'.format(current)]['grid']
    temps = Grids['{}'.format(current)]['temps']
    fields = Grids['{}'.format(current)]['fields']
    # get temperature boundaries
    min_x = min(temps)
    max_x = max(temps)
    min_y = min(fields)
    max_y = max(fields)
        
    
    # Gernerate Plot
    if not 'cmap' in kwargs.keys():
        cmap = cm.jet
    else:
        cmap = kwargs['cmap']  
    plt.imshow(grid.T, origin='lower', extent=[min_x,max_x,min_y,max_y],
           aspect='auto', vmin=vmin, vmax=vmax, cmap = cmap )
    plt.title('Bsweeps @ I = {}\n'.format(current))
    plt.locator_params(nbins=nbins)
    plt.xlabel('$T$ (K)')
    plt.ylabel(r'$\mu_0H\,\rm{(T)}$')
    if not 'clabel' in kwargs.keys():
        clabel = "$\chi'_T$"
    else:
        clabel = kwargs['clabel']
    if not 'cbin' in kwargs.keys():
        cbin = 4
    else:
        cbin = kwargs['cbin']
    if not 'remove_cb' in kwargs.keys():
        plt.colorbar(label = clabel)
    current_string = str(current).replace('.','_').replace(' ','')
    #savefig(plot_path + 'Interpolation_{}.png'.format(current_string))


def plot_data(Exp, meas = [], data = 'real', scaling = -1, N = 400j, nbins = 4,
              legopts = [], **kwargs):
    """Plots a set of data which corresponds
    to :func:`.make_grid` function as a probe.

    Usually only used together with :func:`~.make_grid`
       
    E.g.::
    
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                       key=os.path.getmtime))
        meas = [0,1,2,3]
        plot_data(MnSi_T, meas)
    
    :param dict Exp: Sweet 16 Experiment class dictionary
    :param list meas: List of measurements
    :param float scaling: Scales the susceptibility with the factor scaling
    :param complex N: Gridsize of corresponding grid
    :param int nbins: Number of bins in plot
    :param dict legopts: Valid options are:

                * 'small': Small legend (min and max value) within plot
                * 'short': Smmall legend (min, max) outside of plot
                * 'remove': no legend
    :param dict kwargs: Valid kwargs are:

                * min_x: float, minimal field value
                * max_x: float, maximum field value
              
    """
    
    # get temperature boundaries
    if not 'min_x' in kwargs.keys():
        min_x = min([min(Exp[i].field) for i in meas])
        #min_x = min_x + 0.05*min_x
    else:
        min_x = kwargs['max_x']
    
    if not 'max_x' in kwargs.keys():
        max_x = max([max(Exp[i].field) for i in meas])
        #max_x = max_x - 0.05*max_x
    else:
        max_x = kwargs['max_x']
    
    current, current_string = get_current(Exp, meas)
    
    # Order Measuerments after Temperature
    meas = e21u.order_measurements(Exp, meas, param = 'temp')
    ax = plt.gca()
    for j, i in enumerate(meas):
        c = e21u.Felix_colormap()
        color = c(j/float(len(meas)))
        T = Exp[i].init_temp.strip('K')
        mT = round(Exp[i].mean_temperature,2)
        Exp[i].plot(y = e21u.calibrate(lookup(Exp[i],data), scaling),color=color, 
                    linewidth=3, axes = ax,
                    label = '{} / {}'.format('{:.2f}'.format(mT),
                                            float(T)))
    plt.xlabel('B (T)')
    plt.ylabel('U ($\mu$V)')
    plt.title('Bsweeps @ I = {}\n'.format(current))
    plt.locator_params(nbins=nbins)
    if 'small' in legopts:
        lines = ax.lines[0], ax.lines[-1]
        labels = [l.get_label() for l in lines]
        l = plt.legend(lines, labels, loc = 'lower right')
    elif 'short' in legopts:
        lines = ax.lines[0], ax.lines[-1]
        labels = [l.get_label() for l in lines]
        l = plt.legend(lines, labels, title = 'T (K)\n meas / set',
                   bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    else:
        lines = ax.lines
        labels = [l.get_label() for l in lines]
        l = plt.legend(lines, labels, title = 'T (K)\n meas / set',
                   bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if 'remove' in legopts:
        l.remove()

def plot_interpolated_sweeps(Grids, current, N = 400j):
    """Plots interpolated sweeps
       
    E.g.::
    
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                       key=os.path.getmtime))
        meas = [0,1,2,3]
        g, c = make_grid(MnSi_T, meas)
        Grids[c] = g
        plot_interpolated_sweeps(Grids, c)
    
    :param dict Grids: Dicitonary containing a grid of data as well as x 
                       and y axis values (result of make_grid)
    :param float current: Current for which the grid should be plotted.
    :param complex N: Gridsize (NxN)
    
    """
    grid = Grids['{} A'.format(current)]['grid'] 
    temps = Grids['{} A'.format(current)]['temps']
    fields = Grids['{} A'.format(current)]['fields']
    plt.figure()
    plt.figsize(16,10)
    n = int(N.imag)
    steps = int(n/80)
    nmb = range(0,n,steps)[2:-4]
    for j, i in enumerate(nmb):
        c = e21u.Felix_colormap()
        color = c(j/float(n/steps))
        plt.plot(fields,grid[i], color = color)
        plt.xlabel('B (T)')
        plt.ylabel('U ($\mu$V)')

def get_current_meas(Exp, current):
    """Returns a list of measurement numbers with given current.
       
    E.g.::
    
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                key=os.path.getmtime))
        meas = get_current_meas(MnSi_T, 1.0)
           
    :param dict Exp: Sweet 16 Experiment class
    :param float current: Current
    :return (list): List of measurements
        
    """
    
    meas = []
    for i in range(len(Exp)):
        if Exp[i].mean_current == current :

            meas.append(i)
    return meas

def find_index(Grids,value, parameter = 'fields', current = 0):
    """Returns the index closest to the given field at a given 
    current in a given Grid
       
    E.g.::
    
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                key=os.path.getmtime))
        find_field_index(MnSi_T[0].field, 2)
    
    :param dict Grids: Dictionary, output of :func:`~.make_grid`
    :param list Grids: Arabitrary measurement, eg. Exp[meas].field
    :param float value: Value to look for closest index
    :param str parameter: 'fields' or 'temps', irrelevant if Grids is list
    :return (int): Returns index
    
        
    """
    if isinstance(Grids, dict):
        a = Grids['{} A'.format(current)][parameter]
    elif isinstance(Grids, list):
        a = Grids
    else:
        raise ValueError('No valid input format, must be dict or list')
    ind = min(range(len(a)), key=lambda i: abs(float(a[i])-float(value)))
    return ind

def find_sweep(Grids, value, parameter = 'fields', current = 0):
    """Find a certain B odr T sweep in the Grids dictionary
    
    E.g.::
    
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                       key=os.path.getmtime))
        meas = [0,1,2,3]
        g, c = make_grid(MnSi_T, meas)
        Grids[c] = g
        find_sweep(Grids, 1, parameter = 'fields')
        
    
    :param dict Grids: Dictionary, output of :func:`~.make_grid`
    :param float value: Value at which sweep was performed
    :param str parameter: 'fields' for Tsweep or 'temps' for Bsweep
    :param float current: Current at which to look for
    :return (list): Returns corresponding sweep
    
    """
    
    num = find_index(Grids, value, parameter, current)
    return Grids['{} A'.format(current)]['grid'][num]

def find_res(Grids, Grids_zero, T=28, current=0.0, B_min = 0,
             B_max = 0.55, full_out=False):
    """Finds residuum of two curves. 
    
    Compares a Bsweep at given T and 0 A stored in Grids_zero with a Bsweep
    stored in Grids at T and I in the region between B_min and B_max.
    Returns the residuum of the two best fitting curves and
    the corresponding index. If full_out is True, then returns the
    residua of the difference between every curve at I and
    the curve at I = 0.
    
    E.g.::
    
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                       key=os.path.getmtime))
        meas = get_current_meas(MnSi_T, 0)
        g, c = make_grid(MnSi_T, meas)
        Grids[c] = g
        res, ind = find_res(Grids)
        
    
    :param dict Grids: Dictionary, output of :func:`~.make_grid`
    :param dict Grids_zero: Dictionary, output of :func:`~.make_grid`
    :param float T: Temperature at which to compare Tsweeps
    :param float I: Current of sweeps to compare with zero current sweeps
    :param float B_min: Lower boundary of comparison region
    :param float B_max: Upper boundary of comparison region
    :param bool full_out: False returns the minimal res and the
                          corresponding ind,True returns all residua
    :return (float, int): Returns the residuum and the corresponding index
    
    """
    I0 = find_index(Grids_zero, B_min, parameter = 'fields', current = 0.0)
    I1 = find_index(Grids_zero, B_max, parameter = 'fields', current = 0.0)
    residua = []
    y1 = find_sweep(Grids_zero, T, parameter = 'temps', current = 0.0)[I0:I1]
    for i in range(len(Grids_zero['{} A'.format(current)]['grid'])):
        y2 = Grids_zero['{} A'.format(current)]['grid'][i][I0:I1]
        res = np.sum((y1-y2)**2)
        residua.append((res,i))
    if full_out:
        return min(residua, key=lambda x: x[0]), residua
    return min(residua, key=lambda x: x[0])  

def get_temperature_shift(Grids, Grids_zero, I=0.0, T = 28, B_min = 0.49, B_max = 0.55):
    """ gets the temperature shift of two sweeps at same 
    Temperature but different currents
       
    :param dict Grids: Dictionary, output of :func:`~.make_grid`
    :param float T: Temperature at which to compare Bsweeps
    :param float I: Current of sweep to compare with zero current sweep
    :param float B_min: Lower boundary of comparison region
    :param float B_max: Upper boundary of comparison region
    :return (float): Returns temperature shift 
    
    """
    
    res, ind = find_res(Grids, Grids_zero, T=T, current=I, B_min=B_min, B_max=B_max)
    T0 = Grids['{} A'.format(I)]['temps'][ind]
    return T*pq.K-T0

def get_all_Tshifts(Grids, Grids_zero, I = [], T = 28, B_min = 0.49, B_max = 0.55):
    """ gets the temperature shifts for different currents
       
    :param dict Grids: Dictionary, output of :func:`~.make_grid`
    :param float T: Temperature at which to compare Bsweeps
    :param list I: list of currents
    :param float B_min: Lower boundary of comparison region
    :param float B_max: Upper boundary of comparison region
    :return (list): List of temperature differences

    """
    return [get_temperature_shift(Grids, Grids_zero, I=i, T=T, B_min=B_min, B_max=B_max) for i in I]

def get_valid_currents(Grids):
    """ Returns all valid currents in Grids

    """
    
    currents = Grids.keys()
    currents = [i.strip('A').strip() for i in currents]
    currents = sorted(currents)
    return currents



def chi_I_neu(Data, Data_im, T = 28.0, B = [0.2], offset = 0, y_out = 'real',
              delta = False, **kwargs):


    # Possibility to define own color map
    if not 'cmap' in kwargs.keys():
        cmap = cm.jet
    else:
        cmap = kwargs['cmap']

    if 'fontsize' in kwargs.keys():
        fsize = kwargs['fontsize']
    else:
        fsize = 18.
    
    # get all valid currents of measurement set
    currents = [ float(i.strip('A')) for i in Data[T]['data'].keys()]
    currents = sorted(currents)

    # Iterate over all given field values
    for j, b in enumerate(B):
        
        # find index of magentic field B
        field_index = find_index(list(Data[T]['fields']), b)

        # define color map
        color = cmap((j)/float(len(B)))

        values = [] # data point (chi)
        x_axis = []  # x-coordinate (j)

        # Get chi_0 (@zero current and given T)
        x0 = Data[T]['data']['0.0 A'][field_index]
        x0_im = Data_im[T]['data']['0.0 A'][field_index]
        
        # Iterate over currents
        
        for i, I in enumerate(currents):
            x1 = Data[T]['data']['{} A'.format(I)][field_index]
            x2 = Data_im[T]['data']['{} A'.format(I)][field_index]
            
            points, zlab = calculate_output(x0,x0_im,x1,x2, y_out,delta) 

            # Convert current to current density
            if 'A' in kwargs.keys():
                A = kwargs['A']
                I = float(I)/A
                xlab = r'$\vec{j}$ ($\frac{MA}{m^2}$)'
            else:
                xlab = r'$I$ (A)'
            x_axis.append(I)
            
            # make default title
            values.append(points+(i*offset))        
        plt.title('T = {} K'.format(T))
        plt.plot(x_axis, values, '-o', label = '{}'.format(b), color=color)
    plt.ylabel(zlab)
    plt.xlabel(xlab)

def chi_I(Grids, Grids_zero, T, B=[0.2], B_min=0.3, B_max=0.4, **kwargs):
    """ Produces a plot of Chi vs. I. 
        
        Takes a Dictionary containing various grids of data for different
        currents. Checks which currents are available in grid.

    
    :param dict Grids: Dictionary, output of :func:`~.make_grid`
    :param dict Grids_zero: Dictionary, output of :func:`~.make_grid`
                            usually Grids
    :param float T: Temperature
    :param list B: List of Field values at which to evaluate Chi Vs I
    :param float B_min: Lower boundary of comparison region
    :param float B_max: Upper boundary of comparison region
    :param bool full_out: False returns the minimal res and the
                          corresponding ind,True returns all residua

    """
    
    def get_chi(B, current, ind):
        '''Helper function to get Chi value out of grid

        :param float B: Field Value
        :param float current: Current Values
        :param int ind: Index of desired temperature sweep

        '''
        return Grids['{} A'.format(current)]['grid'][ind][find_index(Grids, B, parameter = 'fields', current = current)]    

    # Possibility to define own color map
    if not 'cmap' in kwargs.keys():
        cmap = cm.jet
    else:
        cmap = kwargs['cmap']
    
    # get all valid currents of measurement set
    currents = get_valid_currents(Grids) 

    # Iterate over all given field values
    for j, b in enumerate(B):

        # define color map
        color = cmap((j)/float(len(B)))

        points = [] # data point (chi)
        x = []  # x-coordinate (j)

        # Get chi_0 (@zero current and given T)
        res, ind = find_res(Grids, Grids_zero,T, 0., B_min, B_max)
        x0 = get_chi(b, 0., ind)
        min_I = min(currents)
        max_I = max(currents)
        
        # Iterate over currents
        for i, I in enumerate(currents):
            res, ind = find_res(Grids,Grids_zero, T, I, B_min, B_max)

            # Deside if output is in percent or in absolute difference Delta
            if not 'output' in kwargs.keys():            
                y = (get_chi(b, I, ind)/x0-1)*100
                ylab = r'$\Delta\chi$ (%)'
            elif kwargs['output'] == 'Delta':
                y = get_chi(b, I, ind)-x0
                ylab = r'$\Delta\chi$ (Si)'
            points.append(y)

            # Convert current to current density
            if 'A' in kwargs.keys():
                A = kwargs['A']
                I = float(I)/A
            x.append(I)

            # make default title
            plt.title('T = {} K'.format(T))
        plt.plot(x, points, '-o', label = '{}'.format(b), color=color)
    plt.ylabel(ylab)
    plt.xlabel('I (A)')

def plot_anpassung_auto(Grids, Grids_zero, T = 28, current = 0.0, B_min = 0, B_max = 0.55):
    """ Muss noch kommentiert werden

    """
    f, ax = plt.subplots(1,1)
    plt.plot(Grids['0.0 A']['fields'], find_sweep(Grids, T, parameter = 'temps',
             current = 0.0),label = '0 A')
    res, ind = find_res(Grids, Grids_zero, T, current=current, B_min=B_min,
                        B_max=B_max)
    plt.plot(Grids['{} A'.format(current)]['fields'], 
             Grids['{} A'.format(current)]['grid'][ind],
             label = '{} A'.format(current))
    plt.xlabel('Magnetic Field (T)')
    plt.ylabel('$\chi$ (a.u.)')

def current_dep(Grids, Grids_zero, T = 28.4, B_min = 0.45, B_max = 0.55, **kwargs):
    """ Muss noch kommentiert werden

    """    
    for i in kwargs.keys():
        if i not in ['currents', 'cmap', 'output', 'col_label']:
            raise KeyError('Invalid Key "{}"'.format(i))
    if 'currents' in kwargs.keys():
        currents = kwargs['currents']
    else:    
        currents = get_valid_currents(Grids)
    if 'col_label' in kwargs.keys():
        col_label = kwargs['col_label']
    else:    
        col_label = 'Re_Chi'
    #f, ax = plt.subplots(1,1)
    if 'cmap' in kwargs.keys():   
        c = kwargs['cmap'] 
    else:    
        c = e21u.Felix_colormap()
    color = c(0./float(len(currents)+1))
    #plot(Grids['0.0 A']['fields'], 
    #     find_temperature(Grids, T, 0.0),label = '0', color = color)
    data = {}
    data['fields'] = Grids['{} A'.format(currents[0])]['fields']
    data['data'] = {}
    for i, I in enumerate(currents):
        color = c((i+1)/float(len(currents)+1))
        res, ind = find_res(Grids, Grids_zero, T, I, B_min, B_max)
        
        
        data['data']['{} A'.format(I)]=Grids['{} A'.format(I)]['grid'][ind]
        if 'output' not in kwargs.keys():
            plt.plot(Grids['{} A'.format(I)]['fields'], 
                     Grids['{} A'.format(I)]['grid'][ind],
                     label = '{}'.format(currents[i]), color = color)
        else:
            if kwargs['output'] == False:
                plt.plot(Grids['{} A'.format(I)]['fields'], 
                         Grids['{} A'.format(I)]['grid'][ind],
                         label = '{}'.format(currents[i]), color = color)
            elif kwargs['output'] is not True:
                raise ValueError('Output needs to be True or False, not {}'.format(kwargs['output']))
    if 'output' in kwargs.keys():
        if kwargs['output'] == True:
            return data
            

def current_dep_neu(Data, T = 28.4, leg = True, **kwargs):
    """ Muss noch kommentiert werden

    """    
    for i in kwargs.keys():
        if i not in ['currents', 'cmap', 'output', 'col_label', 'A', 'color']:
            raise KeyError('Invalid Key "{}"'.format(i))
    if 'currents' in kwargs.keys():
        currents = kwargs['currents']
    else:    
        currents = Data[T]['data'].keys()
    if 'col_label' in kwargs.keys():
        col_label = kwargs['col_label']
    else:    
        col_label = 'Re_Chi'
    #f, ax = plt.subplots(1,1)
    if 'cmap' in kwargs.keys():   
        c = kwargs['cmap'] 
    else:    
        c = e21u.Felix_colormap()
    color = c(0./float(len(currents)+1))
    #plot(Grids['0.0 A']['fields'], 
    #     find_temperature(Grids, T, 0.0),label = '0', color = color)
  
    for i, I in enumerate(sorted(currents)):
        if 'color' in kwargs.keys():
            color = c(kwargs['color'])
        else: 
            color = c((i+1)/float(len(currents)+1))
        if 'A' in kwargs.keys():
            J = round(float(I.strip(' A'))/kwargs['A']*1e-6,2)
            label = '{}'.format(J)
        else:
            label = '{}'.format(I.strip(' A'))
            
        plt.plot(Data[T]['fields'], Data[T]['data'][I], color = color, label = label)
    if leg == True:    
        if 'A' in kwargs.keys():    
            l = plt.legend(title=r'$\vec{j}$ ($\frac{MA}{m^2}$)',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            return l
        else:
            l = plt.legend(title=r'I (A)',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            return l

def mark_critical_current(I = 0.426, pos = .1):
    """Places a dashed line in plot at I with name Ic
    
    :param float I: Critical current
    :param float pos: Horizontal position of Text in percent of total height

    """
    f = plt.gcf()
    ax = plt.gca()
    xmin, xmax, ymin, ymax = ax.axis()
    plt.axvline(x=I, ymin=ymin, ymax=ymax, ls = '--', color = 'grey')
    plt.text(I-(np.abs(xmin-xmax))*0.02, ymin+pos*np.abs(ymin-ymax), '$I_c$')

def calculate_output(x0,x0_im,x1,x2, y_out,delta):
    x0 = np.array(x0)
    x0_im = np.array(x0_im)
    x1 = np.array(x1)
    x2 = np.array(x2)
    if delta == False:
        if y_out == 'real':
            x = x1
            points = (x-x0)/x0*100
            zlab = '\n'+r'$\frac{\Delta \rm{Re}\chi^\perp}{\rm{Re} \chi_0^\perp}$ (%)'
        elif y_out == 'imag':
            x = x2
            points = (x-x0_im)/x0_im*100
            zlab = '\n'+r'$\frac{\Delta  \rm{Im}  \chi^\perp}{\rm{Im} \chi_0^\perp}$ (%)'
        elif y_out == 'abs':
            x = np.sqrt(x1**2 + x2**2)
            x0_abs = np.sqrt(x0**2+x0_im**2)
            points = (x-x0_abs)/x0_abs*100
            zlab = '\n'+r'$\frac{\Delta |\chi^\perp|}{|\chi_0^\perp|}$ (%)'
        elif y_out == 'phi':
            x = np.arctan(x1/x2)
            x0_phi = np.arctan(x0/x0_im)
            points = (x-x0_phi)/x0_phi*100
            zlab = '\n'+r'$\frac{\Delta \varphi}{\varphi_0}$ (%)'
        else:
            raise ValueError('{} invalid for y_out)'.format(y_out))
            
    elif delta == True:
        if y_out =='real':
            x = x1
            points = x-x0
            zlab = '\n'+r'$\Delta \rm{Re} \chi^\perp$'
        elif y_out == 'imag':
            x = x2
            points = x-x0_im
            zlab = '\n'+r'$\Delta \rm{Im} \chi^\perp$'
        elif y_out == 'abs':
            x = np.sqrt(x1**2 + x2**2)
            x0_abs = np.sqrt(x0**2+x0_im**2)
            points = x-x0_abs
            zlab = '\n'+r'$\Delta$ |$\chi^\perp$|'
        elif y_out == 'phi':
            x = np.arctan(x1/x2)
            x0_phi = np.arctan(x0/x0_im)
            points = x-x0_phi
            zlab = '\n'+r'$\Delta \varphi$ (Rad)'
        else:
            raise ValueError('{} invalid for y_out)'.format(y_out))
    else:
        raise ValueError('delta needs to be true or false, not {}'.format(type(delta)))

    return points, zlab

            
            
 


def chi_B(Data, Data_im, T = 28.0, offset = 0, y_out = 'real',
          plot_3d = False, delta = False, **kwargs):
    """ Produces a plot of Chi vs. Magnetic Field 
        
        Takes a Dictionary containing various grids of data for different
        currents. Checks which currents are available in grid.

    
    :param dict Grids: Dictionary, output of :func:`~.make_grid`
    :param dict Grids_imag: Dictionary, output of :func:`~.make_grid`
                            usually Grids
    :param float T: Temperature
    :param list B: List of Field values at which to evaluate Chi Vs I
    :param float B_min: Lower boundary of comparison region
    :param float B_max: Upper boundary of comparison region
    :param bool full_out: False returns the minimal res and the
                          corresponding ind,True returns all residua

    """
    
 
    # Possibility to define own color map
    if not 'cmap' in kwargs.keys():
        cmap = cm.jet
    else:
        cmap = kwargs['cmap']      
    
    # get all valid currents of measurement set
    if 'currents' in kwargs.keys():
        currents = kwargs['currents']
    else:
        currents = [ float(i.strip('A')) for i in Data[T]['data'].keys()]
        currents = sorted(currents)
    points = [] # data points (chi)

    # Get chi_0 (@zero current and given T)
    x0 = Data[T]['data']['0.0 A']
    x0_im = Data_im[T]['data']['0.0 A']
    min_I = min(currents)
    max_I = max(currents)
    
    # For 3D output
    verts = []
    cl = []

    # Iterate over currents 
    for i, I in enumerate(currents):
        # define color map
        color = cmap((i)/float(len(currents)))
        
        #print x
        field = Data[T]['fields']
        #print field
        x1 = Data[T]['data']['{} A'.format(I)]
        x2 = Data_im[T]['data']['{} A'.format(I)]
        xlab = '\n'+r'$\mu_0 H$ (T)'
        ylab = '\n'+r'$I$ (A)' 
        points, zlab = calculate_output(x0,x0_im,x1,x2, y_out,delta) 
        
        for j, val in enumerate(points):
            points[j]=points[j]+(i*offset)

        # Convert current to current density
        if 'A' in kwargs.keys():
            A = kwargs['A']
            j = round(float(I)/A*1e-6,2)
           
        else:
            j = I

        if plot_3d == True:
            points[0], points[-1] = 0, 0
            verts.append(list(zip(field, points)))
            cl.append(color)
        else:
            # Generate 2D Waterfall
            # make default title
            plt.title('T = {} K'.format(T))
            zeroline = [0+i*offset]*len(field)
            plt.plot(field, points, '-', label = '{}'.format(j), color = color)
            #plt.plot(field, zeroline, color = 'g')
            ax = plt.gca()
            ax.fill_between(field, points, zeroline, where=points>=zeroline, facecolor='grey', interpolate=True, alpha = 0.3)
            ax.fill_between(field, points, zeroline, where=points<=zeroline, facecolor='red', interpolate=True, alpha = 0.2)
            
            plt.ylabel(zlab)
            plt.xlabel(xlab)
            #plt.legend(title = r'$\vec{j} (\frac{MA}{m^2})$')
            if 'A' in kwargs.keys():
                plt.legend(title=r'$\vec{j}$ ($\frac{MA}{m^2}$)',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                #e21u.replace_zeros_legend()
                #plt.legend(title=r'$\vec{j}$ ($\frac{MA}{m^2}$)',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    # Generate 3D Watefall         
    if plot_3d == True:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)
        zs = [float(i) for i in currents]
        poly = PolyCollection(verts, facecolor=cl)
        poly.set_alpha(0.7)
        ax.add_collection3d(poly, zs=zs, zdir='y')
        ax.set_xlabel(xlab)
        if 'XL' in kwargs.keys():
            ax.set_xlim3d(kwargs['XL'])
        ax.set_ylabel(ylab)
        if 'YL' in kwargs.keys():
            ax.set_ylim3d(kwargs['YL'])
        ax.set_zlabel('\n'+zlab)
        if 'ZL' in kwargs.keys():
            ax.set_zlim3d(kwargs['ZL'])

def chi_T(Data, Data_im, T = [28.0], current = '0.0', offset = 0, y_out = 'real',
          plot_3d = False, delta = False, **kwargs):
    
    # Possibility to define own color map
    if not 'cmap' in kwargs.keys():
        cmap = cm.jet
    else:
        cmap = kwargs['cmap']      
    
    points = [] # data points (chi)

    # Get chi_0 (@zero current and given T)

    I = current
    # For 3D output
    verts = []
    cl = []

    # Iterate over currents 
    for i, t in enumerate(T):
        # define color map
        x0 = Data[t]['data']['0.0 A']
        x0_im = Data_im[t]['data']['0.0 A']
        color = cmap((i)/float(len(T)))
        
        #print x
        field = Data[t]['fields']
        #print field
        x1 = Data[t]['data']['{} A'.format(I)]
        x2 = Data_im[t]['data']['{} A'.format(I)]
        xlab = '\n'+r'$\mu_0 H$ (T)'
        ylab = '\n'+r'$I$ (A)' 
        points, zlab = calculate_output(x0,x0_im,x1,x2, y_out,delta) 
        
        for j, val in enumerate(points):
            points[j]=points[j]+(i*offset)

        # Convert current to current density
        if 'A' in kwargs.keys():
            A = kwargs['A']
            j = round(float(I)/A*1e-6,2)
           
        else:
            j = I

        if plot_3d == True:
            points[0], points[-1] = 0, 0
            verts.append(list(zip(field, points)))
            cl.append(color)
        else:
            # Generate 2D Waterfall
            # make default title
            plt.title('T = {} K'.format(T))
            zeroline = [0+i*offset]*len(field)
            plt.plot(field, points, '-', label = '{}'.format(t), color = color)
            #plt.plot(field, zeroline, color = 'g')
            ax = plt.gca()
            ax.fill_between(field, points, zeroline, where=points>=zeroline, facecolor='grey', interpolate=True, alpha = 0.3)
            ax.fill_between(field, points, zeroline, where=points<=zeroline, facecolor='red', interpolate=True, alpha = 0.2)
            
            plt.ylabel(zlab)
            plt.xlabel(xlab)
            #plt.legend(title = r'$\vec{j} (\frac{MA}{m^2})$')
            if 'A' in kwargs.keys():
                plt.legend(title=r'T (K)',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                #e21u.replace_zeros_legend()
                #plt.legend(title=r'$\vec{j}$ ($\frac{MA}{m^2}$)',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    # Generate 3D Watefall         
    if plot_3d == True:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)
        zs = [float(i) for i in T]
        poly = PolyCollection(verts, facecolor=cl)
        poly.set_alpha(0.7)
        ax.add_collection3d(poly, zs=zs, zdir='y')
        ax.set_xlabel(xlab)
        if 'XL' in kwargs.keys():
            ax.set_xlim3d(kwargs['XL'])
        ax.set_ylabel(ylab)
        if 'YL' in kwargs.keys():
            ax.set_ylim3d(kwargs['YL'])
        ax.set_zlabel('\n'+zlab)
        if 'ZL' in kwargs.keys():
            ax.set_zlim3d(kwargs['ZL'])

