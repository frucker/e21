import e21.utility as e21u
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import quantities as pq
import matplotlib.cm as cm

def get_current(Exp, meas):
    """ Returns current of a list of measurement. 
      
    Exception raised if measurements have different currents.
    
    E.g.::
        
        path = '/Your_Path/'
        MnSi_T = e21.sweet16.susceptibility.Experiment()
        MnSi_T.add_measurements(sorted(glob.glob(path+'*.dat'),
                                key=os.path.getmtime))
        current, current_string = get_current(MnSi_T, meas)
    
    :param dict Exp: sweet16 Experiment class dictionary
    :param list meas: list of numbers of desired measurements
    
    :raises: ValueError, Measuremetns don't have the same current
    
    :return (float, str): current of measurements as float and string, (e.g. '0.0 A')
    
    """
    if len(meas)==0:
        raise ValueError('No measurements given')
    
    try:
        #get current
        current = abs(Exp[meas[0]].mean_current)
        # format string
        current_string = str(current).replace('.','_').replace(' ','')
        # check consitency
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
              scaling = -1, N = 400j, **kwargs):
    
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
    :param float scaling: scaling factor of data
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
    f, ax = plt.subplots(1,1)    
    plt.imshow(grid.T, origin='lower', extent=[min_x,max_x,min_y,max_y],
           aspect='auto', vmin=vmin, vmax=vmax)
    plt.title('Bsweeps @ I = {}\n'.format(current))
    plt.locator_params(nbins=nbins)
    plt.xlabel('T (K)')
    plt.ylabel('B (T)')
    if not 'clabel' in kwargs.keys():
        clabel = "$\chi'_T$"
    else:
        clabel = kwargs['clabel']
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
    f= plt.figure()
    ax = plt.gca()
    for j, i in enumerate(meas):
        color = e21u.Felix_colormap(j/float(len(meas)))
        T = Exp[i].params['info']['command']['init_temperature'].strip('K')
        mT = round(Exp[i].mean_temperature,2)
        Exp[i].plot(y = lookup(Exp[i],data)*scaling ,color=color, 
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
        color = e21u.Felix_colormap(j/float(n/steps))
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

def chi_I(Grids, Grids_zero, T, B=[0.2], B_min=0.3, B_max=0.4):
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
        return Grids['{} A'.format(current)]['grid'][ind][find_index(Grids, B, parameter = 'fields', current = 0.0)]

    currents = get_valid_currents(Grids)
    f, ax = plt.subplots(1,1)
    
    for j, b in enumerate(B):
        color = cm.jet((j)/float(len(B)))
        #color = e21u.Felix_colormap((j)/float(len(B)))
        points = []
        x = []
        I = 0.
        res, ind = find_res(Grids, Grids_zero,T, I, B_min, B_max)
        x0 = get_chi(b, I, ind)
        min_I = min(currents)
        max_I = max(currents)
        for i, I in enumerate(currents):
            res, ind = find_res(Grids,Grids_zero, T, I, B_min, B_max)
            points.append((get_chi(b, I, ind)/x0-1)*100)
            x.append(I)
            plt.title('T = {} K'.format(T))
        #s = UnivariateSpline(x, points)
        #xs = np.linspace(min_I, max_I, 1000)
        #ys = s(xs)
        plt.plot(x, points, '-o', label = '{}'.format(b), color=color)
        #plot(xs, ys)
    plt.ylabel(r'$\Delta\chi$ (%)')
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

def current_dep(Grids, Grids_zero, T = 28.4, B_min = 0, B_max = 0.55):
    """ Muss noch kommentiert werden

    """    
    currents = get_valid_currents(Grids)
    f, ax = plt.subplots(1,1)
    color = e21u.Felix_colormap(0./float(len(currents)+1))
    #plot(Grids['0.0 A']['fields'], 
    #     find_temperature(Grids, T, 0.0),label = '0', color = color)
    for i, I in enumerate(currents):
        color = e21u.Felix_colormap((i+1)/float(len(currents)+1))
        res, ind = find_res(Grids, Grids_zero, T, I, B_min, B_max)
        plt.plot(Grids['{} A'.format(I)]['fields'], 
                 Grids['{} A'.format(I)]['grid'][ind],
                 label = '{}'.format(currents[i]), color = color)

def mark_critical_current(I = 0.426, pos = .1):
    """Places a dashed line in plot at I with name Ic
    
    :param float I: Critical current
    :param float pos: Horizontal position of Text in percent of total height

    """
    f = plt.gcf()
    ax = plt.gca()
    xmin, xmax, ymin, ymax = ax.axis()
    plt.axvline(x=0.426, ymin=ymin, ymax=ymax, ls = '--', color = 'grey')
    plt.text(I-(np.abs(xmin-xmax))*0.02, ymin+pos*np.abs(ymin-ymax), '$I_c$')


