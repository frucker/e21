# -*- coding: utf-8 -*-
#
# e21, (c) 2013, see AUTHORS. Licensed under the GNU GPL.
'''

STT package includes utility to refine Spin Transfer Torque Resistivity data
take by sweet 16

'''
# Import E21 uitlity
import e21.utility

# Import Quantities 
import quantities as pq

# Add quantity microohms
uohm = pq.unitquantity.UnitQuantity('muohm', pq.ohm*1e-6, u_symbol='muOhm', aliases=['uOhm', 'uohm', 'muohms'])

# Import Python Utility
import scipy
from pylab import *
import numpy

def check_current(Exp = {}, meas = []):

    '''  Get SpinToqrue current of list of measurements and
         check if all the same.

    input:
        Exp: e21 Experiment class. Current Experiment
        meas: list, List of measurements within current Experiment
    
    output:
        current: float. Current of measurements.

    '''    
    # Get current of Measurement
    current = abs(Exp[meas[0]].mean_current)
    
    # Check if all currents are equal
    for i, j in enumerate(meas):
        if (Exp[j].mean_current-current) > 0.01:
            raise Exception("Current of measurement" 
                            "number {} not the same ({}) as standard:"
                            " {}!".format(j,Exp[j].mean_current, current))
            #Todo: Raise warning instead!
    return current

def plot_BS(Exp = {}, meas = [], opt = 'res'):

    ''' Plots rho/hall over B at defined current and different Temperatures.
    
    Reads the Spin-Torque current and raises Exception if not equal for 
    every given measurement.
    Plots symmetrized resistivity or antisymmetrized hall resistivity
    over Magnetic Field for different Temperatures. 
    
    input: 
        Exp: e21.Experiment class, current Experiment
        meas: list, indices of measurements. IMPORTANT: always Up and down sweeps after each other. No single sweeps.
        opt: string, 'res' for resistivity and 'hall' for Hall resisitvity
    
    
    TODO: Check if up-down sweep is correct.
    '''
    # get and check current
    current = check_current(Exp, meas)
    
    # plot data
    
    f = figure()
    ax = gca()
    for j, i in enumerate(meas[::2]):
        color = e21.utility.Felix_colormap()(j/float(len(meas)))
        if opt is 'res':
            sym = e21.utility.symmetrize
        elif opt is 'hall':
            sym = e21.utility.antisymmetrize
        else:
            raise KeyError('Unknown Option: {}'.format(opt))
        Exp[i].plot_res(y = sym(Exp[i], Exp[i+1]).rescale(uohm*pq.cm), color=color, linewidth=3, axes = ax,
                    label = '{} / {} K'.format(np.round(Exp[i].mean_temperature,2),
                                               float(Exp[i].params['info']['command']['init_temperature'].strip('K'))))
        xlabel('B (T)')
        if opt is 'hall':
            ylabel(r'$\rho_{xy}\,(\mu\Omega cm)$', labelpad=15)
        elif opt is 'res':
            ylabel(r'$\rho_{xx}\,(\mu\Omega cm)$', labelpad=15)
        suptitle('Bsweeps @ I = {}'.format(current) , fontsize=20)
        legend(fontsize = 12)


def Interp_Bsweeps_res(Grids = {}, Exp = {}, 
                       meas = [], extent = [7,7.2, -1, 1],
                       N = 400j, opt = 'res', unit = uohm*pq.cm, **kwargs):

    
    ''' Plot and Interpolate B-sweeps of one selected current.
    
    Plot all given B-sweeps or Measruements in meas. Makes a contourplot and interpolates linarily. 
    Plots a given number of equally spaced interpolated B-sweeps.
    
    Furthermore corrects the Temperature for heating effects of STT currents. 
    This is not universally doable, therfore careful with this script!
    Requirements: All curves have to have the length. If not interpolation 
    has to be added for raw data.
    
    input:
        Grids:      dict, Containing grid data and corresponding information. 
                    Keys are currents of measurements
        Exp:        Sweet 16 Experiment class. Current Experiment
        meas:       list, List of Measurements within current Experiment
        extent:     [min_x, max_x, min_y, max_y]: list of float. 
                    window size of interpolation, 
        N:          complex, Grid size (NxN)
        opt:        str, 'res' for resisitvity, 'hall' for Hall-Effect
        unit:       output Unit. Has to be multiples
                    of muohm*cm (which is standard)
        
    output:
    
        Grids:      dict, Dictionary entry at given current 
                    value of interpolated Data
        
    
    '''
    
    # Check input to be an even number and if i and i+1 are 'opposite' sweeps
    if len(meas)%2 != 0:
        raise ValueError('Odd number of measurements. Symmetrize not possible') 
    for i in range(len(meas))[::2]:
        B1 = Exp[meas[i]].params['info']['command']['init_field'].strip('T')
        B2 = Exp[meas[i+1]].params['info']['command']['init_field'].strip('T')
        T1 = Exp[meas[i]].params['info']['command']['init_temperature'].strip('K')
        T2 = Exp[meas[i+1]].params['info']['command']['init_temperature'].strip('K')
        if float(B1) != -1*float(B2): 
            raise ValueError('Sweeps {} and {} are not up and down sweeps'.format(meas[i], meas[i+1]))
        if float(T1) != float(T2):
            raise ValueError('Up and down sweeps {} and {} not at same T'.format(meas[i], meas[i+1]))
        
    # Import scipy utility
    from scipy import interpolate
    from scipy.interpolate import griddata
    
    # Get and check measurement current
    current = check_current(Exp, meas)
    
    # get extent
    
    min_x = extent[0]
    max_x = extent[1]
    min_y = extent[2]
    max_y = extent[3]
    
    # Make color map
    figure()
           
    # make grid
    grid_x, grid_y = np.mgrid[min_x:max_x:N, min_y:max_y:N]
    # x coordinates (list)
    x = []
    # y cooordinates
    y = []
    # values
    values = []
    for i in meas[::2]:    
        #x1 = e21.utility.calibrate(MnSi_13_05T[i].temperature, 2.665e3)
        x1 = Exp[i].temperature
        y1 = Exp[i].field
        if opt is 'hall':
            val = e21.utility.antisymmetrize(Exp[i],Exp[i+1]).rescale(unit)
        elif opt is 'res':
            val = e21.utility.symmetrize(Exp[i],Exp[i+1]).rescale(unit)
        else:
            raise Exception('unknown option: {}'.format(opt))
        x = hstack((x,x1))
        y = hstack((y,y1))
        values = hstack((values, val))
    punkte = vstack((x,y)).T 
    grid = griddata(punkte, values, (grid_x, grid_y), method = 'linear') # make grid
    imshow(grid.T, origin='lower', extent=extent, aspect = 'auto')
    suptitle('Bsweeps @ I = {}'.format(current), fontsize=20)
    xlabel('T (K)')
    ylabel('B (T)')
    cb = colorbar()
    if opt is 'hall':
        cb.set_label(r'$\rho_{xy}\,(\mu\Omega cm)$', labelpad=15)
        ct = '{}_hall'.format(current)
    elif opt is 'res':
        cb.set_label(r'$\rho_{xx}\,(\mu\Omega cm)$', labelpad=15)
        ct = '{}'.format(current)
    Grids[ct] = {}
    Grids[ct]['extent'] = extent
    Grids[ct]['data'] = grid*unit
    Grids[ct]['N'] = int(N.imag)
    
    return Grids

def plot_interpolated_data(Grids, current = 0, smoothing_factor = 11, opt = 'res', ncurves = 40):
    ''' Plot interpolated Bsweeps. 
    
    This function is supposed to be executed after Interp_Bsweeps_res().
    Plots a desired number of interpolated curves. Interpolation is done in
    Interp_Bsweeps_res()

    input:
        Grids:      dict, Interpolation Matrix and specifications. Output
                    of Interp_Bsweeps_res()
        current:    float, measurement current of desired interpolation
        smoothing_factor: smoothing factor of smooth utility
        opt:        str, 'res' for resisitvity, 'hall' for Hall-Effect
        ncurves:    int, number of desired interpolated curves

    '''
        
    
    # Plot interpolated Data
    if opt is 'hall':
        current = '{}_hall'.format(float(current))
    elif opt is 'res':
        current = '{}'.format(float(current))
    else:
        raise KeyError('Unknown option: {}'.format(opt))
        
    # Get interpolation extent
    min_x = Grids[current]['extent'][0]
    max_x = Grids[current]['extent'][1]
    min_y = Grids[current]['extent'][2]
    max_y = Grids[current]['extent'][3]
    
    #start plotting
    figure()
    n = Grids[current]['N']
    steps = int(n/ncurves)
    for j, i in enumerate(range(0,n,steps)):
        color = e21.utility.Felix_colormap()(j/float(n/steps))
        plot(linspace(min_y,max_y,len(smooth(Grids[current]['data'][i],
             window_len = smoothing_factor))), smooth(Grids[current]['data'][i],
             window_len = smoothing_factor), color = color)
        xlabel('B (T)')
        if opt is 'hall':
            ylabel(r'$\rho_{xy}\,(\mu\Omega cm)$', labelpad=15)
        elif opt is 'res':
            ylabel(r'$\rho_{xx}\,(\mu\Omega cm)$', labelpad=15)
        xlim(min_y, max_y)
        rcParams['lines.linewidth'] = 3

def corrected_BS(Grids, path, ref_current = 0, sttcurrent = 579,
                 T_i = 7.01, opt = 'res', plotopt = 1, B_min = -1, 
                 B_max = 1):
    
    ''' Calculate temperatrue drift corrected Bsweep at sttcurrernt.
    
    Takes a table of temperatur differences at specified fields.
    Adds the Delta Ts to T_init and generates a new list of temperatures, 
    which is now corrected for dirfts due to changing magnetoresistance.
    
    Requirements:
    Datafile of temperature differences has to be stored in the path folder
    and has to be named "'enter your current'.txt". File has to contain values
    between B_min and B_max. May be of arbitraty length.

    input:
        Grids:      dict, Interpolation Matrix and specifications. Output
                    of Interp_Bsweeps_res()
        path:       str, path of temperature differences file
        ref_current: float, reference current. Usually 0
        sttcurrent: spin transfer torque current of measurement that should be
                    corrected
        T_i:     float, Temperature to which offsett is added
        opt:        str, 'res' for resisitvity, 'hall' for Hall-Effect
        plotopt:    0 or 1. 0: no plot. 1: results are plotted
        B_min:      float, minimum magnetic field of calibration file
        B_max:      float, maximum magnetic field of calibration file

    output:
        original:   quantity list, original interpolated data
        corrected:  quantity list, corrected data

    '''

    # check if hall or res
    if opt is 'hall':
        current = '{}_hall'.format(float(ref_current))
    elif opt is 'res':
        current = '{}'.format(float(ref_current))
    else:
        raise KeyError('Unknown option: {}'.format(opt))
    
    # get extent and interpolation dimension
    min_x = Grids[current]['extent'][0]
    max_x = Grids[current]['extent'][1]
    min_y = Grids[current]['extent'][2]
    max_y = Grids[current]['extent'][3]
    N = Grids[current]['N']
    
    # open calibration file
    with open(path + '{}.txt'.format(sttcurrent)) as f:
        Temps = []
        for line in f:
            Temps.append(float(line.strip()))

    # Make list of init Temperatures    
    T_init = [T_i]*len(Temps)
    T_real = []

    # Add calibration to init Temperatures
    for i in range(len(T_init)):
        T_real.append(T_init[i] + Temps[i])
    value = []
    
    # set x axis. 
    x = linspace(B_min,B_max,len(Temps))
    
    # loop over all temperatures
    for i, j in enumerate(T_real):

        # get Grids index of Temperature i
        ind = calc_temp_ind(j, min_x, max_x, N)

        # interpolate original data to mach grid
        f = scipy.interpolate.interp1d(linspace(B_min,B_max,N),Grids[current]['data'][ind])
        data = f(x)
        value.append(data[i])

    # get data units
    units = Grids[current]['data'][0].units

    # get original interpolated data, match grid by interp1d and add unit
    init_ind = calc_temp_ind(T_i, min_x, max_x, N)
    f = scipy.interpolate.interp1d(linspace(B_min,B_max,N),Grids[current]['data'][init_ind])
    original = f(x)
        
    # add unit to corrected data    
    corrected = value*units
    
    # plot corrected and original data if plotoption is 1
    if plotopt == 1:
        figure()
        plot(x, data, label = 'original')
        plot(x, value, label = 'new')
        xlabel('B (T)')
        if opt is 'hall':
            ylabel(r'$\rho_{xy}\,(\mu\Omega cm)$', labelpad=15)
        elif opt is 'res':
            ylabel(r'$\rho_{xx}\,(\mu\Omega cm)$', labelpad=15)
        suptitle('At {} mA'.format(sttcurrent), fontsize = 20)
        legend()
    
    return original, corrected



def calc_temp_ind(Temp, min_x = 0, max_x = 1, N = 200):
    ''' Calculate Temperature index of Interpolation Matrix (Grids)
    
    input:
        Temp:           float, Temperature
        min_x, max_x:   Boarders of interpolation
        N:              Dimension of interpolation (NxN)

    output:
        index: int, index of Measurement corresponding to Temp

    '''
    steps = (max_x-min_x)/N
    return int((Temp-min_x)/steps)

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

