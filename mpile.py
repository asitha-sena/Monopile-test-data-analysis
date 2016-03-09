#Copyright Asitha Senanayake (2015)

from numpy import *
from pylab import *


###################################################
###  CONVERT DATA FILE TO READABLE FORMATS      ###
###################################################
    
def calibration_curves(degree=3):
    """
    Produce calibration curves for each sensel in pressure sensor. This function is called by 'convert_raw_psi()' to 
    convert raw readings to psi.
    
    Input:
    -----
    degree - Degree of polynomial fitted into calibration data
    
    Output:
    ------
    A 32 by 32 2-d list populated with poly1d objects. Each poly1d object represents the
    calibration curve for the corresponding sensel. 
    """
    pressure = array([0.0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0])
    
    raw = loadtxt(r'D:\Dropbox\UT\Research\Fall 2014\DAQ\Original Calibration Files from Sensor Products Inc\UT5010-8802\0 - 5psi calibration.txt')
    raw = raw.reshape(18,32,32)
    
    cal_curves = [[] for i in range(0,32)]
    
    for i in range(0,32):
        for j in range(0,32):
            c = polyfit(raw[:,i,j], pressure, degree)   #i - width/circumference, j-height/depth
            cal_curves[i].append(poly1d(c))
            
    return cal_curves


def convert_pressure_sensor_file(filename):
    """Converts pressure sensor data files exported in matrix form by the Tactilus software to a format readable by 
    'loadtxt()' by commenting out the lines with non-float type data. Use 'loadtext()' to obtain a 3D array with the 
    pressure sensor data from this converted text file.

    Input:
    -----
    filename - Path to data file

    Output:
    ------
    Text file readable by loadtxt() with '.N' appended to name to name of input file
    
    Example:
    -------
    convert_pressure_sensor_file(r'D:\Box Sync\Lab Tests\2014_11_25\Pressure data, matrix, before load test.txt')
    data = loadtxt(r'D:\Box Sync\Lab Tests\2014_11_25\Pressure data, matrix, before load test.N.txt').reshape(58,32,32)
    
    data[Frame][width][height]
    """
    
    f_old = open(filename, 'r')
        
    #Open new file with '_N' appended to save the cleaned up data
    temp     = filename.split('.')
    temp[-1] = 'N.txt'        
    filename_new     = '.'.join(temp)
    f_new = open (filename_new, 'w')
    
    
    for line in f_old:
        temp = line.split()
        
        if len(temp) > 0:
            
            #Write line if first item is a float
            try: 
                float(temp[0])
                f_new.write(line)
                
            #Comment out lines which are not floats
            except ValueError:
                temp = '#' + line
                f_new.write(temp)
            
    f_old.close()
    f_new.close()


    
def convert_raw_to_psi(filename1,filename2,degree=3):
    """
    Input:
    -----
    filename1   - Input file with raw data, Text file with pressure matrices from Tactilus
    filename2   - Output file with pressure in psi
    degree      - Degree of polynomial fitted into calibration data (default is 3 i.e cubic polynomial)
    
    Output:
    ------
    File with converted pressure readings in psi and readable by loadtxt() function
    """
    
    #Load calibration curves
    cal_curves = calibration_curves(3)
    
    #Convert input file into format readable by loadtxt() function i.e. comment out strings
    convert_pressure_sensor_file(filename1)
    
    #Load raw data into 3D array from loadtxt() readable text file
    temp     = filename1.split('.')
    temp[-1] = 'N.txt'        
    file1     = '.'.join(temp)
    raw_data = loadtxt(file1)
    
    x,y      = raw_data.shape
    frames   = x/32
    raw_data = raw_data.reshape(frames,32,32)            #p[frame][width][height]
    
    #Initialize 3D array for pressure readings
    pressure = zeros((frames,32,32))
    
    #Populate pressure array
    
    for i in range(0,frames): #frame/time
        for j in range(0,32): #width/circumference
            for k in range(0,32): #column/height/depth
                pressure[i,j,k] = cal_curves[j][k](raw_data[i,j,k])
    
    #Write 3D pressure array into new file
    file_new = open(filename2,'w')
    
    frame=0
    for data_slice in pressure:
        # Description of slice
        file_new.write('# Frame = %d \n' %frame)
        
        # The formatting string indicates that I'm writing out
        # the values in left-justified columns 6 characters in width
        # with 2 decimal places.  
        np.savetxt(file_new, data_slice, fmt='%6.3f')
        
        frame+=1
    
    file_new.close()
    

def convert_labview_datafile(filename):
    """Covert datafile from Yunhan's VI with data in engineering units and voltages to proper format.
    
    Input:
    ------
    filename - Data file from Yunhan's Labview VI
    
    Output:
    -------
    Datafile with following format:
    Date    Time   Load (lb)    Load (V)   LVDT (in)    LVDT (V)   LMT (in)   LMT (V)
    """
    
    import datetime as datatime
    
    f_old = open(filename, 'r')
    
    #Open new file with '_N' appended to save the cleaned up data
    temp     = filename.split('.')
    temp[-1] = 'N.txt'        
    filename_new     = '.'.join(temp)
    f_new = open (filename_new, 'w')
    f_new.write('\tDate\t\tTime\t\tLoad(lb) Load(V) LVDT(in) LVDT(V) LMT(in) LMT(V)\n')

    #Read raw datafile and write extracted data to new file              
    for line in f_old:
        if line.startswith('D')==True:
            temp = line.split()
            temp1 = temp[0:4]       #Initialize temp1
            
            temp1.append(temp[6]) #Load (V)
            temp1.append(temp[4]) #LVDT (in)
            temp1.append(temp[7]) #LVDT (V)
            temp1.append(temp[5]) #LMT (in)
            temp1.append(temp[8]) #LMT (V)
            
            temp1 = '\t'.join(temp1)
            f_new.write(temp1 + '\n')
    
    f_old.close()
    f_new.close()


###################################################
###              READ DATA FILES                ###
###################################################

def read_data_file_1(filename, cal_load=-16.285, cal_disp=0.1313, method='LMT'):
    """Read data from OTRC-SCF and OTRC-MPILE output files with a single displacement measurement 

    Input:
    -----
    'filename' - Location of raw data file 
    'cal_load' - Load cell calibration factor (16.285 by default)
    'cal_disp_1' - LVDT1/LMT calibration factor (-0.1313 by default, change it to 7.78 for LMT)
    'cal_disp_2' - LVDT2 calibration factor (-0.1354 by default)

    Output:
    ------
    data_time  - Array with time data
    data_load  - Array with load data
    data_disp  - Array with displacement data
    zero_load  - Zero error in load
    loc        - Number of entries in time, load, and displacement arrays
    """
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename, 'r')
    
    if method == 'LMT':     column = 4; #cal_disp = -7.806
    elif method == 'LVDT':  column = 4; #cal_disp = 0.1313
    elif method == 'LVDT_2': column= 4; #cal_disp = -0.1354
    elif method == 'Motor': column = 22; #cal_disp = 0.000062
        
    #Initialize arrays in which to store extracted data
    time_stamp = array([datetime.timedelta(hours=i) for i in xrange(100000)])
    data_load  = zeros((100000))
    data_disp  = zeros((100000))
    data_time  = zeros((100000))
    
    loc  = 0               
    for line in file1:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[loc] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell and LVDT/LMT Voltage
            data_load[loc] = float(temp[3])
            data_disp[loc] = float(temp[column])
            loc += 1                           #Counter to keep track of the last data entry into the data arrays
    
            
    #Identify zero errors 
    zero_time = time_stamp[0]
    zero_load = data_load[0]
    zero_disp = data_disp[0]
    
    #Convert Timestamps in to Seconds, Remove zero errors in Load & Disp, Apply Calibration Factors
    for i in range(loc):
        data_time[i] = (time_stamp[i] - time_stamp[0]).total_seconds()
        data_load[i] = (data_load[i] - zero_load)*cal_load #- 1.6
        data_disp[i] = (data_disp[i] - zero_disp)*cal_disp
        
    file1.close()
    return data_time, data_load, data_disp, zero_load, loc


def read_data_file_2(filename, cal_load=-1, cal_disp=1, method='LMT'):
    """Read data from Yunhan's new DAQ VI output files with a single displacement measurement 

    Input:
    -----
    'filename' - Location of raw data file 
    'cal_load' - Load cell calibration factor 
    'cal_disp_1' - LVDT1/LMT calibration factor 
    'cal_disp_2' - LVDT2 calibration factor

    Output:
    ------
    data_time  - Array with time data
    data_load  - Array with load data
    data_disp  - Array with displacement data
    zero_load  - Zero error in load
    loc        - Number of entries in time, load, and displacement arrays

    """
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename, 'r')
    
    if method == 'LMT':     column = 5; #cal_disp = 1
    elif method == 'LVDT':  column = 3; #cal_disp = 1
    elif method == 'LVDT_2': column= 4; #cal_disp = 1
        
    #Initialize arrays in which to store extracted data
    time_stamp = array([datetime.timedelta(hours=i) for i in xrange(100000)])
    data_load  = zeros((100000))
    data_disp  = zeros((100000))
    data_time  = zeros((100000))
    
    loc  = 0               
    for line in file1:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[loc] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell and LVDT/LMT Voltage
            data_load[loc] = float(temp[4])
            data_disp[loc] = float(temp[column])
            loc += 1                           #Counter to keep track of the last data entry into the data arrays
    
            
    #Identify zero errors 
    zero_time = time_stamp[0]
    zero_load = data_load[0]
    zero_disp = data_disp[0]
    
    #Convert Timestamps in to Seconds, Remove zero errors in Load & Disp, Apply Calibration Factors
    for i in range(loc):
        data_time[i] = (time_stamp[i] - time_stamp[0]).total_seconds()
        data_load[i] = (data_load[i] - zero_load)*cal_load #- 1.6
        data_disp[i] = (data_disp[i] - zero_disp)*cal_disp
        
    file1.close()
    return data_time, data_load, data_disp, zero_load, loc

def read_data_file_3(filename, cal_load=-1, cal_disp=1, method='LVDT'):
    """Read data from Yunhan's new DAQ VI output files with both engineering units and raw voltages. 

    Input:
    -----
    'filename' - Location of raw data file 
    'cal_load' - Load cell calibration factor 
    'cal_disp_1' - LVDT1/LMT calibration factor 
    'cal_disp_2' - LVDT2 calibration factor

    Output:
    ------
    data_time  - Array with time data
    data_load  - Array with load data
    data_disp  - Array with displacement data
    zero_load  - Zero error in load
    loc        - Number of entries in time, load, and displacement arrays

    """
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename, 'r')
    
    if   method == 'LVDT':   column = 5;
    elif method == 'LMT':    column = 7;
    elif method == 'LVDT_2': column = 7;
        
    #Initialize arrays in which to store extracted data
    time_stamp = array([datetime.timedelta(hours=i) for i in xrange(100000)])
    data_load  = zeros((100000))
    data_disp  = zeros((100000))
    data_time  = zeros((100000))
    
    loc  = 0               
    for line in file1:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[loc] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell and LVDT/LMT Voltage
            data_load[loc] = float(temp[3])
            data_disp[loc] = float(temp[column])
            loc += 1                           #Counter to keep track of the last data entry into the data arrays
    
            
    #Identify zero errors 
    zero_time = time_stamp[0]
    zero_load = data_load[0]
    zero_disp = data_disp[0]
    
    #Convert Timestamps in to Seconds, Remove zero errors in Load & Disp, Apply Calibration Factors
    for i in range(loc):
        data_time[i] = (time_stamp[i] - time_stamp[0]).total_seconds()
        data_load[i] = (data_load[i] - zero_load)*cal_load #- 1.6
        data_disp[i] = (data_disp[i] - zero_disp)*cal_disp
        
    file1.close()
    return data_time, data_load, data_disp, zero_load, loc


def read_pressure_sensor_file_1(filename, D):
    """Reads data files exported from the Tactilus software with the format of statistics over two regions of 
    interest (ROI) which have been pre-selected. Use 'convert_pressure_sensor_file()' in conjuction with 'loadtext()' 
    in order to extract data from files exported in the form of matrices.

    Input:
    -----
    filename - Path to file with pressure sensor data
    D        - Diameter of pile

    Output:
    ------
    p_net    - Array with load/unit length data
    loc      - Number of entries in p_net array
    """

    #Pressure sensor data
    file1 = open(filename, 'r')
    
    
    P_ROI_1 = zeros(200000) #Front face of pile
    P_ROI_2 = zeros(200000) #Back face of pile
    
    loc=0
    for line in file1:
        temp = line.split()  
        if line.startswith('Average pressure(ROI 1)')==True:
            P_ROI_1[loc] = float(temp[3])
        elif line.startswith('Average pressure(ROI 2)')==True:
            P_ROI_2[loc] = float(temp[3])     
            loc += 1
    
    #plot(P_ROI_1[0:loc], label='ROI 1')
    #plot(P_ROI_2[0:loc], label='ROI 2')
    #grid(True), legend(), xlabel('Number of Frames'), ylabel('Pressure (psi)')
    
    #Remove zero errors
    P_ROI_1 = P_ROI_1 - P_ROI_1[0]
    P_ROI_2 = P_ROI_2 - P_ROI_2[0]
    
    P_net = P_ROI_1 - P_ROI_2
    p_net = -(P_net - P_net[0])*D                    #Remove zero reading and convert pressure to a line load
    
    file1.close()
    return p_net, loc



def read_pressure_sensor_file_2(filename,alpha=0):
    """Calculate resultant lateral load acting on the pile based on pressure sensor data. The input file should be a pressure
    sensor matrix file exported from Tactilus and converted by 'convert_pressure_sensor_file()' to be readable by 'loadtxt()'.
    
    Input:
    -----
    filename - Path and name of converted pressure sensor matrix datafile
    frames   - Number of frames recorded in the file
    alpha    - Initial offset angle. This can be changed to find out the front and the back facing sensels. 
    
    Output:
    ------
    Returns 'p', a 2D array of shape 'frames by height', with the force per unit length (lb/in) acting at each height ranging 
    from 0 to 32 inches.
    """
    
    pressure = loadtxt(filename)
    x,y      = pressure.shape
    frames   = x/32
    pressure = pressure.reshape(frames,32,32)            #p[frame][width][height]
    
    dA      = 32*12.6/(32*32)                #Area of a single sensel
    d_theta = 2*pi/32                        #Angle subtended by a single sensel (radians)
    theta = arange(alpha + d_theta/2, alpha + 2*pi, d_theta)   #Vector of average angles to sensels, measured from axis of symmetry
    
    #Resolving forces on individual sensels and integrating over the pile face
    p = zeros((frames,32))
    
    for t in range(0,frames):
        for h in range(0,32):
            p[t,h] = sum(pressure[t,:,h]*sin(theta)*dA)  #Sum pressure around circumference
    
    return p


###################################################
###   ANALYZE MONOTONIC TESTS                   ###
###################################################

def plot_monotonic_test(filename, method='LVDT', x1=0, x2=350, dP=0.0, name='Model Test',DAQ="OLD", graphs='Load/Disp',
                        cal_load=16.285, cal_disp=0.1313, ls='-'):
    """Plots the results from monotonic lateral load tests (to failure) from raw or cleaned up data. Uses the motor displacement
    instead of the LMT/LVDT data due to shortcomings in the sensors (i.e. LVDT stroke insufficient and LMT sensitivity inadequate.

    Input:
    -----
    filename - Location of data file
    x1, x2   - start and end points for plotting
    dP       - Manually add imbalance load due to counter-weights being attached. Compare counter-weight and initial load cell
                 reading to find dP
    DAQ      - 'OLD', 'NEW', 'NEW_2'
    graphs   - 'Complete' to plot Load-Time, Disp-Time, and Load-Disp graphs
               'Load/Disp' to plot plot only Load-Disp graph (Default)
    cal_load - Load calibration. Depends on the type of sensor and test date. Use 1 if engineering units are been directly read from datafile.
    cal_disp - Load calibration. Depends on the type of sensor and test date. Use 1 if engineering units are been directly read from datafile.
    ls       - Linestyle in plots. Follow PyPlot format from linestyles eg: 'b-', 'r:' 

    Output:
    ------
    Plot of load versus displacement curve.
    """
    
    import datetime as datetime
    import time as time
    
    #Read Data
    if DAQ=="OLD":
        [data_time, data_load, data_disp, zero_load, loc] = read_data_file_1(filename, cal_load=cal_load, cal_disp=cal_disp, method=method)
    elif DAQ=="NEW":
        [data_time, data_load, data_disp, zero_load, loc] = read_data_file_2(filename, cal_load=cal_load, cal_disp=cal_disp, method=method)
    elif DAQ=="NEW_2":
        [data_time, data_load, data_disp, zero_load, loc] = read_data_file_3(filename, cal_load=cal_load, cal_disp=cal_disp, method=method)

        
    print "Number of data points = %d" %loc
    print "Approximate initial load = %3.1flb\n" %(zero_load*16.285 + 1.6)
    
    #Plot Data

    #x1 = 0; x2 = loc;          #5 readings/sec, 10 sec/cycle

    if graphs=='Complete':
        #rcParams['figure.figsize'] = 18, 12  #Change the size of the inline figures/graphs
        #rcParams['lines.linewidth'] = 2
    
        #Load time history
        subplot(3,1,1)
        plot(data_time[x1:x2],data_load[x1:x2], ls)
        xlabel('Time (seconds)'), ylabel('Load (lb)'), grid(True)
        
        #Displacement time history
        subplot(3,1,2)
        plot(data_time[x1:x2],data_disp[x1:x2], ls)
        xlabel('Time (seconds)'), ylabel('Displacement (in)'),grid(True)
        
        #Hysterisis loops
        subplot(3,1,3)
        plot(data_disp[x1:x2], data_load[x1:x2], ls, label=name), #xlim(xmin=0)
        xlabel('Displacement (in)'), ylabel('Load (lb)'), grid(True)

    else:
        #rcParams['figure.figsize'] = 12, 8  #Change the size of the inline figures/graphs
        #rcParams['lines.linewidth'] = 2
        
        plot(data_disp[x1:x2], data_load[x1:x2], ls, label=name), #xlim(xmin=0)
        xlabel('Displacement (in)'), ylabel('Load (lb)'), grid(True)



def plot_pressure_sensor_monotonic_1(filename1, filename2, D, L, e, z_rot, z, x, dx=0, A=1,B=1, DAQ='OLD', method='LVDT'):
    """Plots p-y curves using a load-displacement data file and the corresponding pressure sensor data file in the form of ROI statistics.

    Input:
    -----
    filename1 - file with displacement and load cell data
    filename2 - file with pressure sensor data as ROI statistics exported from the Tactilus software
    D - Pile diameter (inches)
    L - Length of embedment (inches)
    e - Load eccectricity from mudline (inches)
    z_rot - Depth to the center of rotation from mudline (inches)
    z - Depth of the p-y curve to plot
    x - Number of datapoints to plot (depends on the duration of the test)
    dx - Variable to synchronize the start times of the pressure sensor and displacement measurements
    A,B - Assign either +1 or -1 in order to control the sign when plotting
    Note: This function takes into account two ROIs from the pressure sensor data. Always check if these ROIs and the averaging method matches with
    those selected for each test. The ROIs will depend on the orientation of the pressure sensor relative to the direction of movement.

    Output:
    ------
    Plot p-y curves at from specified depth.
    """
    
    #Load and LVDT data
    if DAQ=='OLD':
        [t, load, disp, zero_load, x] = read_data_file_1(filename1, method=method)
    elif DAQ=='NEW':
        [t, load, disp, zero_load, x] = read_data_file_2(filename1, method=method)
    
    #Load pressure sensor data
    [p_net,loc] = read_pressure_sensor_file_1(filename2,D)
    
    #Remove initial offset (or zero error)
    p_net = p_net - p_net[0]
    
    #Distance to p-y curve of interest from the point of load application
    z_py = z + e #inches
    
    #Calculate lateral displacement (y)
    y = disp*(e+z_rot - z_py)/(e+z_rot)             #Convert displacement at load application point to displacement at depth z
    y = y - y[0]                                    #Remove zero reading
    
    #Control +/- for plotting 
    p_net = A*p_net
    y = B*y
        

    rcParams['figure.figsize'] = 12,8
    """
    plot(y[0:x],-p_net[dx:x+dx], label='z = %.1f-in' %z)  
    xlabel('Displacement (in)'), ylabel('p (lb/in)'), grid(True), legend(loc='lower right')
    #ylim((-0.5,4.5))
    
    print '\nTime duration over which load/disp data was measured       = %.1f seconds' %t[x-1]
    print 'Number of load/disp data points                            = %d' %x
    print 'Time duration over which pressure sensor data was measured = %.1f seconds' %(loc/5.0)
    print 'Number of load/disp data points                            = %d\n' %loc
    """
    
    plot(disp[0:x], abs(p_net[0:x]*32), label='Pressure Sensor')



def plot_pressure_sensor_monotonic_2(filename1, filename2, D, L, e, z_rot, DAQ='OLD', method='LVDT',
                                     delta_z=0, depths=[2,12,24], alpha=0, x=1000, A=1):
    """Plots p-y curves using a load-displacement data file and the corresponding pressure sensor data file in matrix form.

    Input:
    -----
    filename1 - file with displacement and load cell data
    filename2 - file with pressure sensor data in matrix form and converted by 'convert_pressure_sensor_file()'
    D         - Pile diameter (inches)
    L         - Length of embedment (inches)
    e         - Height of LVDT measured from mudline (inches)
    z_rot     - Depth to the center of rotation from mudline (inches)
    DAQ       - The output files from the new and old DAQ systems are different therefore they require different 
                functions to read them. Enter 'OLD' or 'NEW' in order to select the correct read function.
    method    - Displacement measurement method. 'LVDT', 'LMT', or 'LVDT2'.
    depths    - Depths of p-y curves to be plotted (in) eg: [2,4,10,30]
    delta_z   - Height of top row of sensels above the mudline (in)
    alpha     - Pile orientation (rad)
    x         - Number of data points to plot. The max data points in the input file wil be chosen if x is too big
    A         - Coefficient to control the moment equilibrium output (A=1 by default)
    
    Output:
    ------
    Figure 1  - p versus depth
    Figure 2  - p-y curves at specified depths
    """
    
    #Load and LVDT data
    if DAQ=='OLD':
        [t, load, disp, zero_load, loc] = read_data_file_1(filename1, method=method)
    elif DAQ=='NEW':
        [t, load, disp, zero_load, loc] = read_data_file_2(filename1, method=method)
    elif DAQ=='NEW_2':
        [t, load, disp, zero_load, loc] = read_data_file_3(filename1, method=method)
        
    #Load pressure sensor data
    p = read_pressure_sensor_file_2(filename2, alpha=alpha)
    
    frames = p.shape[0]
    x = min(x,loc,frames)
    
    #Remove initial offset (or zero error)
    p = p - p[1,:]
    
    z    = arange(0,32)          #Depth of each row of sensors (inches)
    h    = 32 - z - delta_z      #Distance vertically upwards from pile tip to each row of sensors (inches)
    z_py = e + z                 #Distance to p-y curve of interest from the point of load application (inches)
    
    y = zeros((x,32))            #Initialize array to hold displacement vs frames for each row of sensors

    #Force equilibrium            
    F = sum(p, axis=1)              #Sum pressure along pile length
    F = F - F[0]                    #Remove zero error from the force vector

    #Moment equilibrium
    M = zeros(frames)
    
    for t in range(0,frames):
        for i in range(0,32):
            M[t]   += A*(p[t,i] - p[0,i])*(i+delta_z)   #Resultant moment about pile tip due to mobilized soil resistance,
                                                      #distance from pile tip to middle of sensel row 0 = delta_z
    
    #rcParams['figure.figsize'] = 18,8

    #Plot load-displacement curve
    #subplot2grid((2,3),(0,0), colspan=2)
    subplot(3,2,1)
    plot(disp[0:x],abs(F[0:x]), label='Pressure Sensor, Force Eqbm.')
    plot(disp[0:x], M[0:x]/(L+e), label='Pressure Sensor, Moment Eqbm.')
    xlabel('Displacement (in)'), ylabel('Load (lb)'), grid(True)  
    
    #Plot p vs z along pile
    #subplot2grid((2,3),(0,2), rowspan=2)
    subplot(3,2,2)
    for j in range(0,p.shape[0]):
        plot(p[j,:], h, 'b.')
        
    xlabel('p (lb/in)'), ylabel('Depth (in)'), grid(True)
    ax = gca()
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    #Plot p vs y at specified depths
    #subplot2grid((2,3),(1,0), colspan=2)
    j=2
    for i in depths:
        j += 1
        subplot(3,2,j)
        y[:,i] = disp[0:x]*(z_rot - z[i])/(z_rot + e)             #Convert displacement at load application point to displacement at depth z
        plot(abs(y[0:x,i]),abs(p[0:x,31-i]), label='z = %.1f-in' %z[i])
        
    xlabel('Displacement (in)'), ylabel('p (lb/in)'), grid(True), legend(loc='best', fontsize='small')
    
    


###################################################
###     ANALYZE CYCLIC TESTS                    ###
###################################################
    
def plot_cyclic_test(filename, cal_load=16.285, cal_disp=0.1313, method='LVDT', DAQ='OLD',
                     T = 10, n1=0, n2=1000, n3 = 100, dn=400, repeat=3):
    """Plot results from the cyclic lateral load tests on piles.

    Input:
    ------
    filename    - Location of raw data file 
    cal_load    - Load cell calibration factor (16.285 by default)
    cal_disp    - LVDT/LMT calibration factor (0.1313 by default)
    method      - Method of displacement measurement ('LVDT', 'LMT', 'MOTOR')
    DAQ         - Labview DAQ, Datafile format ('OLD', 'NEW', 'NEW_2')
    T           - Period of load/disp cycle (s)
    n1, n2, n3 dn, repeat - start, end for time history, end for load-disp, gap load-disp, and number of families of loops

    Output:
    -------
    Plots the following graphs in a 3x1 grid:
    Load vs Time, Displacement vs Time, Load vs Displacement"""
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename, 'r')
    
    #Read Data
    if DAQ=="OLD":
        [data_time, data_load, data_disp, zero_load, loc] = read_data_file_1(filename=filename, cal_load=cal_load, cal_disp=cal_disp, method=method)
    elif DAQ=="NEW":
        [data_time, data_load, data_disp, zero_load, loc] = read_data_file_2(filename=filename, cal_load=cal_load, cal_disp=cal_disp, method=method)
    elif DAQ=="NEW_2":
        [data_time, data_load, data_disp, zero_load, loc] = read_data_file_3(filename=filename, cal_load=cal_load, cal_disp=cal_disp, method=method)

    
    print "Number of data points = %d" %loc
    
    x1 = n1*(5*T); x2 = n2*(5*T); x3 = n3*(5*T); dx = dn*(5*T)          #5 readings/sec, 10 sec/cycle

    #rcParams['figure.figsize']=18,12
    #rcParams['figure.figsize'] = 8, 6  #Change the size of the inline figures/graphs
    rcParams['lines.linewidth'] = 1
    
    figure()
    #Load time history
    subplot(3,1,1)
    plot(data_time[x1:min(x2,loc)],data_load[x1:min(x2,loc)], label = 'Cycles %d to %d' %(n1,n2))
    #plot(data_time[x1+dx:x2+dx], data_load[x1+dx:x2+dx], label = 'Cycles %d to %d' %(n1+dn, n2+dn))
    xlabel('Time (seconds)'), ylabel('Load (lb)'), grid(True), legend(loc='best')
    
    #Displacement time history
    subplot(3,1,2)
    plot(data_time[x1:min(x2,loc)],data_disp[x1:min(x2,loc)], label = 'Cycles %d to %d' %(n1,n2))
    xlabel('Time (seconds)'), ylabel('Displacement (in)'),grid(True), legend(loc='best')
    
    #Hysterisis loops
    subplot(3,1,3)
    for i in range(0,repeat):
        dx = i*50*dn
        plot(data_disp[x1+dx:min(x3+dx, loc)], data_load[x1+dx:min(x3+dx,loc)], label = 'Cycles %d to %d' %(n1+i*dn, n3+i*dn))
        xlabel('Displacement (in)'), ylabel('Load (lb)'), grid(True), legend(loc='best')
        

def plot_stiffness(filename, cal_load=16.285, cal_disp=0.1313, method='LVDT', DAQ='OLD',
                   n1=0, n2=1000, name=r'D = ?, $\theta$ = ?', disp_zero=0.1, load_zero = 0):
    """Plot results from the cyclic lateral load tests on piles.

    Input:
    ------
    filename    - Location of raw data file 
    cal_load    - Load cell calibration factor (16.285 by default)
    cal_disp    - LVDT/LMT calibration factor (0.1313 by default)
    method      - Method of displacement measurement ('LVDT', 'LMT', 'MOTOR')
    DAQ         - Labview DAQ, Datafile format ('OLD', 'NEW', 'NEW_2')
    n1, n2      - start, end for time history
    disp_zero, load_zero' - The point of reference used to separate the peaks and troughs in the sub-function 'extract_ampl'

    Output:
    -------
    Plots the following graphs in a 3x1 grid:
    Load Amplitude vs Cycle, Displacement Amplitude vs Cycle, Secant Stiffness vs Cycle
    
    Note: If an 'out of bounds' error results:
    1. Decrease n2
    2. Check if 'zero' is within the bounds of the disp/load values. 
    3. Also make sure that the displacements from one-way load tests are positive since the default point of reference is +0.1"""
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename, 'r')
    
    #Read Data
    if DAQ=="OLD":
        [data_time, data_load, data_disp, zero_load, loc] = read_data_file_1(filename=filename, cal_load=cal_load, cal_disp=cal_disp, method=method)
    elif DAQ=="NEW":
        [data_time, data_load, data_disp, zero_load, loc] = read_data_file_2(filename=filename, cal_load=cal_load, cal_disp=cal_disp, method=method)
    elif DAQ=="NEW_2":
        [data_time, data_load, data_disp, zero_load, loc] = read_data_file_3(filename=filename, cal_load=cal_load, cal_disp=cal_disp, method=method)
    
    print "Number of data points = %d" %loc
    
    #Plot Data
    
    rcParams['figure.figsize'] = 18, 12  #Change the size of the inline figures/graphs
    rcParams['lines.linewidth'] = 2
    
    def extract_ampl(data,zero):
        """Input array with cyclic data in order to extract the amplitude of each cycle. The algorithm isolates each peak and trough
        in order to calculate the amplitude."""
        i=0
        max_data = zeros(1000)
        min_data = zeros(1000)
        
        for m in range(0,n2):
            temp1 = [zero]    #Initilize with a single value to avoid error if max() is called on an empty list
            temp2 = [zero]    #Initilize with a single value to avoid error if min() is called on an empty list
            
            #It is more intuitive to use 0 as the point of reference in the while loops 
            #(and it works well with load cell data as it oscillates around 0 during the test).
            #However, the displacement in one-way cyclic testing does not oscillate around zero, rather it is the start or end point. 
            
            while data[i] >= zero and i<100000-1:
                temp1.append(data[i])
                i += 1
            
            while data[i] < zero and i<100000-1:
                temp2.append(data[i])
                i += 1
            
            max_data[m] = max(array(temp1))
            min_data[m] = min(array(temp2))
            
            #plot(array(temp1))     #Use these plots to check if peaks and troughs are being isolated correctly
            #plot(array(temp2))     #-do-
            
        ampl = abs(max_data - min_data)
            
        return ampl
        
    load_ampl = extract_ampl(data_load, load_zero)
    disp_ampl = extract_ampl(data_disp, disp_zero)
    stiffness = load_ampl/disp_ampl
    cycles    = range(0,n2)
    
    subplot(3,1,1), plot(cycles[n1:n2], load_ampl[n1:n2], label=name)
    ylabel('Load Amplitude (lb)'), grid(True), #legend(), #ylim(ymin=0)
    
    subplot(3,1,2), plot(cycles[n1:n2], disp_ampl[n1:n2], label=name)
    ylabel('Displacement Amplitude (in)'), grid(True), #legend()
    
    subplot(3,1,3), 
    plot(cycles[n1:n2], stiffness[n1:n2], label=name)
    xlabel('Number of Cycles'), ylabel('Secant Stiffness (lb/in)'), grid(True), legend(), #ylim(ymin=0), #, xscale('log')
        

###################################################
###     ANALYZE T-BAR TESTS                     ###
###################################################

def plot_tbar_test_1(filename1, filename2, method='LMT', start_depth=1, test_name='', x1=0, x2=180,
                     plot_OCR='YES', gamma_sub=30, m=0.67, lamda=0.19, depth_max=40, Su_max=70, A=1, B=1):
    """Plots results from T-bar test results using datafiles from the old OTC-SCF and OTC-MPILE Labiew VIs.
    
    filename1   - location of raw data from T-bar test
    filename2   - location of raw data from Rod test
    method      - Displacement measurement using LMT only.
    start_depth - depth at which t-bar test was started (i.e. embedment depth) in inches  (default=1)
    test_name   - Label for plotting
    x1, x2      - start point and end point for plot
    plot_OCR    - 'YES' by default but the calculation is only valid for undisturbed Su,
                  'NO' for remolded Su
    gamma_sub   - submerged unit weight of the soil, used in the OCR calculations
    m           - Exponent used in the equation (c/p)_NC = (c/p)_OC*OCR^m,
                  m=.67 for kaolin (Jeanjean, 2009), m=0.8 in general
    lamda       - c/p ratio, 0.19 for kaolin, 0.25 for Gulf of Mexico clay
    depth_max   - Maximum depth in plots
    Su_max      - Maximum Su in plots
    A, B        - Coefficients to switch the sign of load and displacement respectively (Default A=1,B=1)
        
    Output:
    ------   
    
    Plots Displacement vs Time, Undrained Shear Strength vs. Depth, and OCR vs Depth
    """
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename1, 'r') #T-bar Data File
    file2 = open(filename2, 'r') #Rod Data File
    
    if method == 'LMT':     column = 4; cal_disp = 7.806*B
    elif method == 'Motor': column = 22; cal_disp = -0.000062*3.25*B
    
    #Calibration Factors
    cal_load = -16.285*A
    
    #Initialize arrays in which to store extracted data
    time_stamp = array([[datetime.timedelta(hours=i) for i in xrange(100000)],[datetime.timedelta(hours=i) for i in xrange(100000)]])
    data_load  = zeros((2,10000))
    data_disp  = zeros((2,10000))
    data_time  = zeros((2,10000))
    #data_motor = zeros((2,10000))
    net_load   = zeros(10000)
    
    loc1  = 0               
    for line in file1:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[0,loc1] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell Voltage, LMT/LVDT Voltage, and Motor Location
            data_load[0,loc1] = float(temp[3])
            data_disp[0,loc1] = float(temp[column])
            #data_motor[0,loc1] = float(temp[22])   #Change column number to 5 if using cleaned up input file
            loc1 += 1   #Counter to keep track of the last data entry into the data arrays
    
           
    loc2 = 0
    for line in file2:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[1,loc2] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell Voltage, LVDT/LMT Voltage, and Motor Location
            data_load[1,loc2] = float(temp[3])
            data_disp[1,loc2] = float(temp[column])
            #data_motor[0,loc2] = float(temp[22])  #Change column number to 5 if using cleaned up input file
            loc2 += 1   #Counter to keep track of the last data entry into the data arrays
          
    #Identify zero errors 
    zero_time = array(time_stamp[:,0])
    zero_load = array(data_load[:,0])
    zero_disp = array(data_disp[:,0])
    #zero_motor = array(data_motor[:,0])
    
    #Convert Timestamps in to Seconds, Remove zero errors in Load & Disp, Apply Calibration Factors
    for i in range(0,2):
        for j in range(0,min(loc1,loc2)):
            data_time[i,j] = (time_stamp[i,j]- zero_time[i]).total_seconds()
            data_load[i,j] = (data_load[i,j] - zero_load[i])*cal_load
            data_disp[i,j] = (data_disp[i,j] - zero_disp[i])*cal_disp - start_depth
            #data_motor[i,j] = (data_motor[i,j] - zero_motor[i])*cal_motor
            
    #Subtract rod friction from T-bar load
    i = 0
    for z0 in data_disp[0,0:loc1]:
        j = 0
        for z1 in data_disp[1,j:loc2]:
            if abs(z0-z1)<0.2:
                net_load[i] = data_load[0,i] - data_load[1,j]
                break
            j += 1
        i += 1
    
    print "Number of data points = %d" %min(loc1,loc2)
    
    #rcParams['figure.figsize'] = 18, 10  #Change the size of the inline figures/graphs
    rcParams['figure.figsize'] = 8, 6
    rcParams['lines.linewidth'] = 1
    
    subplot(1,3,1)
    plot(data_time[0,x1:x2], abs(data_disp[0,x1:x2]), label= test_name)
    xlabel('Time (s)'), ylabel('Depth (in)'), yticks(arange(0,depth_max,5))
    grid(True), legend()
    
    subplot(1,3,2)
    #Su Calculation, N=10.5, Area of tbar = 0.028 sq.ft.
    Su = net_load/10.5/0.028 
    
    plot(Su[x1:x2], abs(data_disp[0,x1:x2]), label= test_name)   
    xlabel('Undrained Shear Strength (psf)'), #ylabel('Depth (in)')
    yticks(arange(0,depth_max,5)), xlim([0,Su_max]) 
    legend(), grid(True)
    
    if plot_OCR=='YES':
        subplot(1,3,3)
    
        #OCR Calculation
        sigma_v_eff = gamma_sub*abs(data_disp[0,:])/12.0 #See draft thesis chapter on Test System for calculations
        #plot(sigma_v_eff[x1:x2], abs(data_disp[0,x1:x2]), label=r'$\sigma_v"$')
        OCR = (Su/sigma_v_eff/lamda)**(1.0/m)       #lamda = 0.19 or Kaolin based on 0.11+0.0037PI Madabhushi & Haiderali (2013)
        
        #Start OCR plot from 2-in below the mudline
        for i in range(0,loc1):
            x1 = i
            if abs(data_disp[0,i]) - 2 > 0: break
        
        plot(OCR[x1:x2], abs(data_disp[0,x1:x2]), label= test_name)
        xlabel('Overconsolidation Ratio (OCR)'), #ylabel('Depth (in)')
        yticks(arange(0,depth_max,5)), xlim([0.1,100]), xscale('log')
        grid(True), legend(loc='lower right')


def plot_tbar_test_2(filename1, filename2, method='LMT', start_depth=1, test_name='', x1=0, x2=180,
                     plot_OCR='YES', gamma_sub=30, m=0.67, lamda=0.19, depth_max=25, Su_max=30, cal_load=1, cal_disp=1):
    """Plots results from T-bar tests using data from Yunhan's Labview VI.
    
    Input Parameters:
    ----------------
    
    filename1   - location of raw data from T-bar test
    filename2   - location of raw data from Rod test
    method      - Displacement measurement using LMT only.
    start_depth - depth at which t-bar test was started (i.e. embedment depth) in inches  (default=1)
    test_name   - Label for plotting
    x1, x2      - start point and end point for plot
    plot_OCR    - 'YES' by default but the calculation is only valid for undisturbed Su,
                  'NO' for remolded Su
    gamma_sub   - submerged unit weight of the soil, used in the OCR calculations
    m           - Exponent used in the equation (c/p)_NC = (c/p)_OC*OCR^m,
                  m=.67 for kaolin (Jeanjean, 2009), m=0.8 in general
    depth_max   - Maximum depth in plots
    Su_max      - Maximum Su in plots
    lamda       - c/p ratio, 0.19 for kaolin, 0.25 for Gulf of Mexico clay
    cal_load    - Load calibration factor
    cal_disp    - Displacement calibration factor
    
    Output:
    ------ 
    
    Plots Displacement vs Time, Undrained Shear Strength vs. Depth, and OCR vs Depth
    """
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename1, 'r') #T-bar Data File
    file2 = open(filename2, 'r') #Rod Data File

    
    #Initialize arrays in which to store extracted data
    time_stamp = array([[datetime.timedelta(hours=i) for i in xrange(100000)],[datetime.timedelta(hours=i) for i in xrange(100000)]])
    data_load  = zeros((2,10000))
    data_disp  = zeros((2,10000))
    data_time  = zeros((2,10000))
    net_load   = zeros(10000)
    
    loc1  = 0               
    for line in file1:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[0,loc1] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell Voltage, LMT/LVDT Voltage, and Motor Location
            data_load[0,loc1] = float(temp[4])
            data_disp[0,loc1] = float(temp[5])
            loc1 += 1   #Counter to keep track of the last data entry into the data arrays
    
           
    loc2 = 0
    for line in file2:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[1,loc2] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell Voltage, LVDT/LMT Voltage, and Motor Location
            data_load[1,loc2] = float(temp[4])
            data_disp[1,loc2] = float(temp[5])
            loc2 += 1   #Counter to keep track of the last data entry into the data arrays
          
    #Identify zero errors 
    zero_time = array(time_stamp[:,0])
    zero_load = array(data_load[:,0])
    zero_disp = array(data_disp[:,0])
    
    #Convert Timestamps in to Seconds, Remove zero errors in Load & Disp, Apply Calibration Factors
    for i in range(0,2):
        for j in range(0,min(loc1,loc2)):
            data_time[i,j] = (time_stamp[i,j]- zero_time[i]).total_seconds()
            data_load[i,j] = (data_load[i,j] - zero_load[i])*cal_load
            data_disp[i,j] = (data_disp[i,j] - zero_disp[i])*cal_disp + start_depth
            
    #Subtract rod friction from T-bar load
    i = 0
    for z0 in data_disp[0,0:loc1]:
        j = 0
        for z1 in data_disp[1,j:loc2]:
            if abs(z0-z1)<0.3:
                net_load[i] = data_load[0,i] - data_load[1,j]
                break
            j += 1
        i += 1
    
    print "Number of data points = %d" %min(loc1,loc2)
    
    #rcParams['figure.figsize'] = 18, 10  #Change the size of the inline figures/graphs
    rcParams['figure.figsize'] = 8, 6
    rcParams['lines.linewidth'] = 2
    
    #Su Calculation, N=10.5, Area of tbar = 0.028 sq.ft.
    Su = -net_load/10.5/0.028 
    
    subplot(1,3,1)
    plot(data_time[0,x1:x2], data_disp[0,x1:x2], label= test_name)
    xlabel('Time (s)'), ylabel('Depth (in)')
    yticks(arange(0,depth_max,5)), grid(True), legend()
    
    subplot(1,3,2)
    plot(Su[x1:x2], data_disp[0,x1:x2], label= test_name)   #N=10.5, Area of tbar = 0.028 sq.ft.
    xlabel('Undrained Shear Strength (psf)'), #ylabel('Depth (in)') 
    yticks(arange(0,depth_max,5)), xlim([0,Su_max]) 
    legend(), grid(True)

    if plot_OCR=='YES':
        subplot(1,3,3)
    
        #OCR Calculation
        sigma_v_eff = gamma_sub*data_disp[0,:]/12.0 #See draft thesis chapter on Test System for calculations
        OCR = (Su/sigma_v_eff/lamda)**(1.0/m)       #lamda = 0.19 or Kaolin based on 0.11+0.0037PI Madabhushi & Haiderali (2013)
        
        #Start OCR plot from 2-in below the mudline
        for i in range(0,loc1):
            x1 = i
            if data_disp[0,i] - 2 > 0: break
        
        plot(OCR[x1:x2], data_disp[0,x1:x2], label= test_name)
        xlabel('Overconsolidation Ratio (OCR)'), #ylabel('Depth (in)'), 
        yticks(arange(0,depth_max,5)), xlim([0.1,100]), xscale('log')
        grid(True), legend(loc='lower right')
    
def plot_tbar_test_3(filename1, filename2, method='LMT', start_depth=1, test_name='', x1=0, x2=180,
                     gamma_sub=30, plot_OCR='YES', m=0.67, lamda=0.19, depth_max=20, Su_max=70,
                     tbar_diameter=1.0, tbar_length=4.0, cal_load=1, cal_disp=1):
    """Plots results from T-bar tests using data from Yunhan's Labview VI (latest output format).
    
    Input Parameters:
    ----------------
    
    filename1   - location of raw data from T-bar test
    filename2   - location of raw data from Rod test
    method      - Displacement measurement using LMT only.
    start_depth - depth at which t-bar test was started (i.e. embedment depth) in inches  (default=1)
    test_name   - Label for plotting
    x1, x2      - start point and end point for plot
    plot_OCR    - 'YES' by default but the calculation is only valid for undisturbed Su,
                  'NO' for remolded Su
    gamma_sub   - submerged unit weight of the soil, used in the OCR calculations
    m           - Exponent used in the equation (c/p)_NC = (c/p)_OC*OCR^m,
                  m=.67 for kaolin (Jeanjean, 2009), m=0.8 in general
    lamda       - c/p ratio, 0.19 for kaolin, 0.25 for Gulf of Mexico clay
    depth_max   - Maximum depth in plots
    Su_max      - Maximum Su in plots
    tbar_diameter - Diameter of tbar that is used (default 1.0-inch)
    tbar_length - Length of tbar that is used (default 4.0-inches)
    cal_load    - Load calibration factor
    cal_disp    - Displacement calibration factor
    
    Output:
    ------ 
    
    Plots Displacement vs Time, Undrained Shear Strength vs. Depth, and OCR vs Depth
    """
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename1, 'r') #T-bar Data File
    file2 = open(filename2, 'r') #Rod Data File
    
    #Initialize arrays in which to store extracted data
    time_stamp = array([[datetime.timedelta(hours=i) for i in xrange(100000)],[datetime.timedelta(hours=i) for i in xrange(100000)]])
    data_load  = zeros((2,10000))
    data_disp  = zeros((2,10000))
    data_time  = zeros((2,10000))
    net_load   = zeros(10000)
    
    loc1  = 0               
    for line in file1:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[0,loc1] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell Voltage, LMT/LVDT Voltage, and Motor Location
            data_load[0,loc1] = float(temp[3])
            data_disp[0,loc1] = float(temp[7])
            loc1 += 1   #Counter to keep track of the last data entry into the data arrays
    
           
    loc2 = 0
    for line in file2:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[1,loc2] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell Voltage, LVDT/LMT Voltage, and Motor Location
            data_load[1,loc2] = float(temp[3])
            data_disp[1,loc2] = float(temp[7])
            loc2 += 1   #Counter to keep track of the last data entry into the data arrays
          
    #Identify zero errors 
    zero_time = array(time_stamp[:,0])
    zero_load = array(data_load[:,0])
    zero_disp = array(data_disp[:,0])
    
    #Convert Timestamps in to Seconds, Remove zero errors in Load & Disp, Apply Calibration Factors
    for i in range(0,2):
        for j in range(0,min(loc1,loc2)):
            data_time[i,j] = (time_stamp[i,j]- zero_time[i]).total_seconds()
            data_load[i,j] = (data_load[i,j] - zero_load[i])*cal_load
            data_disp[i,j] = (data_disp[i,j] - zero_disp[i])*cal_disp + start_depth
            
    #Subtract rod friction from T-bar load
    i = 0
    for z0 in data_disp[0,0:loc1]:
        j = 0
        for z1 in data_disp[1,j:loc2]:
            if abs(z0-z1)<0.3:
                net_load[i] = data_load[0,i] - data_load[1,j]
                break
            j += 1
        i += 1
    
    print "Number of data points = %d" %min(loc1,loc2)
    
    #rcParams['figure.figsize'] = 18, 10  #Change the size of the inline figures/graphs
    rcParams['figure.figsize'] = 8, 6
    rcParams['lines.linewidth'] = 2
    
    #Su Calculation, N=10.5
    A = tbar_diameter*tbar_length/144.0 #Projected area of tbar, Area = 0.028 sq.ft. (default)
    
    Su = -net_load/10.5/A 
    
    subplot(1,3,1)
    plot(data_time[0,x1:x2], data_disp[0,x1:x2], label= test_name)
    xlabel('Time (s)'), ylabel('Depth (in)'),
    yticks(arange(0,depth_max,5)), grid(True), legend()
    
    subplot(1,3,2)
    plot(Su[x1:x2], data_disp[0,x1:x2], label= test_name)   #N=10.5
    xlabel('Undrained Shear Strength (psf)'), #ylabel('Depth (in)'), 
    yticks(arange(0,depth_max,5)),
    xlim([0,Su_max]), legend(), grid(True)
    
    if plot_OCR=='YES':
        subplot(1,3,3)
        
        #OCR Calculation
        sigma_v_eff = gamma_sub*data_disp[0,:]/12.0 #See draft thesis chapter on Test System for calculations
        OCR = (Su/sigma_v_eff/lamda)**(1.0/m)       #lamda = 0.19 or Kaolin based on 0.11+0.0037PI Madabhushi & Haiderali (2013)
        
        #Start OCR plot from 2-in below the mudline
        for i in range(0,loc1):
            x1 = i
            if data_disp[0,i] - 2 > 0: break
        
        plot(OCR[x1:x2], data_disp[0,x1:x2], label= test_name)
        xlabel('Overconsolidation Ratio (OCR)'), #ylabel('Depth (in)'), 
        yticks(arange(0,depth_max,5)), xlim([0.1,100]), xscale('log')
        grid(True), legend(loc='lower right')

def plot_tbar_test_without_rod_friction(filename1, method='LMT', start_depth=1, test_name='', x1=0, x2=180,
                                        gamma_sub=30, plot_OCR='YES', m=0.67, lamda=0.19, depth_max=20, Su_max=70,
                                        tbar_diameter=1.0, tbar_length=4.0):
    """Plots results from T-bar tests using data from Yunhan's Labview VI.
    
    Input Parameters:
    ----------------
    
    filename1   - location of raw data from T-bar test
    method      - Displacement measurement using LMT only.
    start_depth - depth at which t-bar test was started (i.e. embedment depth) in inches  (default=1)
    test_name   - Label for plotting
    x1, x2      - start point and end point for plot
    plot_OCR    - 'YES' by default but the calculation is only valid for undisturbed Su,
                  'NO' for remolded Su
    gamma_sub   - submerged unit weight of the soil, used in the OCR calculations
    m           - Exponent used in the equation (c/p)_NC = (c/p)_OC*OCR^m,
                  m=.67 for kaolin (Jeanjean, 2009), m=0.8 in general
    lamda       - c/p ratio, 0.19 for kaolin, 0.25 for Gulf of Mexico clay
    depth_max   - Maximum depth in plots
    Su_max      - Maximum Su in plots
    tbar_diameter - Diameter of tbar that is used (default 1.0-inch)
    tbar_length - Length of tbar that is used (default 4.0-inches)
    
    Output:
    ------ 
    
    Plots Displacement vs Time, Undrained Shear Strength vs. Depth, and OCR vs Depth
    """
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename1, 'r') #T-bar Data File
    
    #Calibration Factors are directly applied in the Labview VI
    cal_load = 1
    cal_disp = 1
    
    #Initialize arrays in which to store extracted data
    time_stamp = array([datetime.timedelta(hours=i) for i in xrange(100000)])
    data_load  = zeros((10000))
    data_disp  = zeros((10000))
    data_time  = zeros((10000))
    
    loc1  = 0               
    for line in file1:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[loc1] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell Voltage, LMT/LVDT Voltage, and Motor Location
            data_load[loc1] = float(temp[3])
            data_disp[loc1] = float(temp[7])
            loc1 += 1   #Counter to keep track of the last data entry into the data arrays
    
                     
    #Identify zero errors 
    zero_time = array(time_stamp[0])
    zero_load = array(data_load[0])
    zero_disp = array(data_disp[0])
    
    #Convert Timestamps in to Seconds, Remove zero errors in Load & Disp, Apply Calibration Factors
    for j in range(0,loc1):
        data_time[j] = (time_stamp[j]- zero_time).total_seconds()
        data_load[j] = (data_load[j] - zero_load)*cal_load
        data_disp[j] = (data_disp[j] - zero_disp)*cal_disp + start_depth
            
    net_load = data_load
    
    print "Number of data points = %d" %loc1
    
    rcParams['figure.figsize'] = 18, 10  #Change the size of the inline figures/graphs
    rcParams['lines.linewidth'] = 2
    
    #Su Calculation, N=10.5
    A = tbar_diameter*tbar_length/144.0 #Projected area of tbar, Area = 0.028 sq.ft. (default)
    Su = -net_load/10.5/A 
    
    subplot(1,3,1)
    plot(data_time[x1:x2], data_disp[x1:x2], label= test_name)
    xlabel('Time (s)'), ylabel('Depth (in)'),
    yticks(arange(0,depth_max,5)), grid(True), legend()
    
    subplot(1,3,2)
    plot(Su[x1:x2], data_disp[x1:x2], label= test_name)   #N=10.5, Area of tbar = 0.028 sq.ft.
    xlabel('Undrained Shear Strength (psf)'), #ylabel('Depth (in)'),
    yticks(arange(0,depth_max,5)),
    xlim([0,Su_max]), legend(), grid(True)
    
    if plot_OCR=='YES':
        subplot(1,3,3)
        
        #OCR Calculation
        sigma_v_eff = gamma_sub*data_disp[:]/12.0 #See draft thesis chapter on Test System for calculations
        OCR = (Su/sigma_v_eff/lamda)**(1.0/m)       #lamda = 0.19 or Kaolin based on 0.11+0.0037PI Madabhushi & Haiderali (2013)
        
        #Start OCR plot from 2-in below the mudline
        for i in range(0,loc1):
            x1 = i
            if data_disp[i] - 2 > 0: break
        
        plot(OCR[x1:x2], data_disp[x1:x2], label= test_name)
        xlabel('Overconsolidation Ratio (OCR)'), #ylabel('Depth (in)'),
        yticks(arange(0,depth_max,5)), xlim([0.1,100]), xscale('log')
        grid(True), legend(loc='lower right')

def plot_tbar_test_4(filename1, filename2, method='LMT', start_depth=0, test_name='', x1=0, x2=180,
                     gamma_sub=30, plot_OCR='YES', m=0.67, lamda=0.19, depth_max=20, Su_max=70,
                     tbar_diameter=1.0, tbar_length=4.0, critical_depth=2.5, cal_load=1, cal_disp=1):
    """Plots results from T-bar tests using data from Yunhan's Labview VI (latest output format). 
    Applies a linearly increasing bearing capacity factor (N_c) from 5.14 at the mudline to 10.5 at depth. 
    
    This function assumes that the T-bar starts just above the mudline and assigns N_c = 5.14. If the T-bar was 
    embedded in the soil at the start of the test then, the embedment depth can be specified. However, if the T-bar
    was too far above the mudline at the start of the test then, the linearly interpolation will not be done correctly.
    
    This function should only be used to analyze penetration data. The extraction portion of the test will not be 
    correctly analyzed due to limitations in the algorithms used to find the rod friction and the N_c values.
    
    Input Parameters:
    ----------------
    
    filename1   - location of raw data from T-bar test
    filename2   - location of raw data from Rod test
    method      - Displacement measurement using LMT only.
    start_depth - depth at which t-bar test was started (i.e. embedment depth) in inches  (default=0)
    test_name   - Label for plotting
    x1, x2      - start point and end point for plot
    plot_OCR    - 'YES' by default but the calculation is only valid for undisturbed Su,
                  'NO' for remolded Su
    gamma_sub   - submerged unit weight of the soil, used in the OCR calculations
    m           - Exponent used in the equation (c/p)_NC = (c/p)_OC*OCR^m,
                  m=.67 for kaolin (Jeanjean, 2009), m=0.8 in general
    lamda       - c/p ratio, 0.19 for kaolin, 0.25 for Gulf of Mexico clay
    depth_max   - Maximum depth in plots
    Su_max      - Maximum Su in plots
    tbar_diameter - Diameter of tbar that is used (default 1.0-inch)
    tbar_length - Length of tbar that is used (default 4.0-inches)
    critical_depth  - Depth at which full plastic flow is achieved i.e. N_c = 10.5 (Unit = T-bar diameter)
                  2.5 recommended for OC clay beds (default value), 5.0 is recommended for NC clay beds
    cal_load    - Load calibration factor
    cal_disp    - Displacement calibration factor
    
    Output:
    ------ 
    
    Plots Displacement vs Time, Undrained Shear Strength vs. Depth, and OCR vs Depth
    """
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename1, 'r') #T-bar Data File
    file2 = open(filename2, 'r') #Rod Data File
    
    #Initialize arrays in which to store extracted data
    time_stamp = array([[datetime.timedelta(hours=i) for i in xrange(100000)],[datetime.timedelta(hours=i) for i in xrange(100000)]])
    data_load  = zeros((2,10000))
    data_disp  = zeros((2,10000))
    data_time  = zeros((2,10000))
    net_load   = zeros(10000)
    
    loc1  = 0               
    for line in file1:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[0,loc1] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell Voltage, LMT/LVDT Voltage, and Motor Location
            data_load[0,loc1] = float(temp[3])
            data_disp[0,loc1] = float(temp[7])
            loc1 += 1   #Counter to keep track of the last data entry into the data arrays
    
           
    loc2 = 0
    for line in file2:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[1,loc2] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*100000)
            
            #Extract Load Cell Voltage, LVDT/LMT Voltage, and Motor Location
            data_load[1,loc2] = float(temp[3])
            data_disp[1,loc2] = float(temp[7])
            loc2 += 1   #Counter to keep track of the last data entry into the data arrays
          
    #Identify zero errors 
    zero_time = array(time_stamp[:,0])
    zero_load = array(data_load[:,0])
    zero_disp = array(data_disp[:,0])
    
    #Convert Timestamps in to Seconds, Remove zero errors in Load & Disp, Apply Calibration Factors
    for i in range(0,2):
        for j in range(0,min(loc1,loc2)):
            data_time[i,j] = (time_stamp[i,j]- zero_time[i]).total_seconds()
            data_load[i,j] = (data_load[i,j] - zero_load[i])*cal_load
            data_disp[i,j] = (data_disp[i,j] - zero_disp[i])*cal_disp + start_depth
            
    #Subtract rod friction from T-bar load
    i = 0
    for z0 in data_disp[0,0:loc1]:
        j = 0
        for z1 in data_disp[1,j:loc2]:
            if abs(z0-z1)<0.3:
                net_load[i] = data_load[0,i] - data_load[1,j]
                break
            j += 1
        i += 1
    
    print "Number of data points = %d" %min(loc1,loc2)
    
    #rcParams['figure.figsize'] = 18, 10  #Change the size of the inline figures/graphs
    rcParams['figure.figsize'] = 8, 6
    rcParams['lines.linewidth'] = 2
    
    
    #Su and N_c Calculation
    Su  = zeros(loc1)
    N_c = zeros(loc1)
    
    A = tbar_diameter*tbar_length/144.0 #Projected area of tbar, Area = 0.028 sq.ft. (default)
    
    for i in range(loc1):
        if data_disp[0,i] < critical_depth*tbar_diameter:
            N_c[i] = 5.14 + (10.5 - 5.14)/(critical_depth*tbar_diameter)*data_disp[0,i]      
        else:
            N_c[i] = 10.5
            
        Su[i] = -net_load[i]/N_c[i]/A 
            
    #plot(N_c[x1:x2],data_disp[0,x1:x2]), xlabel('$N_c$'), ylabel('Displacement (in)')
    #grid(True)
    #ax = gca()
    #ax.invert_yaxis(), ax.xaxis.tick_top(), ax.xaxis.set_label_position('top')
    #figure()
    
    subplot(1,3,1)
    plot(data_time[0,x1:x2], data_disp[0,x1:x2], label= test_name)
    xlabel('Time (s)'), ylabel('Depth (in)'),
    yticks(arange(0,depth_max,5)), grid(True), legend()
    
    subplot(1,3,2)
    plot(Su[x1:x2], data_disp[0,x1:x2], label= test_name) 
    xlabel('Undrained Shear Strength (psf)'), #ylabel('Depth (in)'), 
    yticks(arange(0,depth_max,5)),
    xlim([0,Su_max]), legend(), grid(True)
    
    if plot_OCR=='YES':
        subplot(1,3,3)
        
        #OCR Calculation
        sigma_v_eff = gamma_sub*data_disp[0,:]/12.0 #See draft thesis chapter on Test System for calculations
        OCR = (Su/sigma_v_eff/lamda)**(1.0/m)       #lamda = 0.19 or Kaolin based on 0.11+0.0037PI Madabhushi & Haiderali (2013)
        
        #Start OCR plot from 2-in below the mudline
        for i in range(0,loc1):
            x1 = i
            if data_disp[0,i] - 2 > 0: break
        
        plot(OCR[x1:x2], data_disp[0,x1:x2], label= test_name)
        xlabel('Overconsolidation Ratio (OCR)'), ylabel('Depth (in)'), 
        yticks(arange(0,depth_max,5)), xlim([0.1,100]), xscale('log')
        grid(True), legend(loc='lower right')


###################################################
###           PUBLISHED P-Y MODELS              ###
###################################################
    
def plot_matlock_py(D, S_u0, k, z, gamma_eff=40, epsilon_50=0.020, graph='static', ls='-'):
    """
    Plots Matlock (1970) p-y curves.

    Input:
    ------
    D           - Pile diameter (inches)
    S_u0        - Undrained shear strength at the mudline (psf)
    k           - Rate of increase of shear strength versus depth (psf/ft)
    z           - Depth of p-y curve
    gamma_eff   - Effective unit weight (pcf)
    epsilon_50  - Strain factor
    graph       - Type of curve 'static', 'cyclic', or 'both'
    ls          - Linestyle following the PyPlot format eg: 'b-', 'g:'

    Output:
    -------
    Plot of Matlock (1970) p-y curves.
    """
    
    #rcParams['figure.figsize'] = 10,8
    #rcParams['font.size'] = 14
    
    #Normalized p-y curve
    Y = linspace(0,15,1000)
    P = 0.5*Y**(1.0/3.0)
    
    y_c = 2.5*D*epsilon_50
    S_u = S_u0 + k*(z/12.0)
    
    for i in range(0,len(Y)): 
        if P[i] > 1.0: P[i] = 1.0
    
    #Normalized depth
    Z = z/D
    
    J = 0.5
    sigma_v = gamma_eff*(z/12.0)     #z has to be converted to feet since it was input in inches
    N_p = 3 + sigma_v/S_u + J*Z
    
    if N_p >= 9.0:
        N_p = 9.0
    
    p_ult = N_p*(S_u/144.0)*D
    p = P*p_ult
    y = Y*y_c
    
    #Normalized cyclic p-y curve
    P_cyc = P
    z_r   = D/J*(9-3- k/gamma_eff)
    
    for i in range(0,len(Y)): 
        if Y[i] > 8 and z<z_r: 
            P_cyc[i] = 0.72*z/z_r
        elif 3 < Y[i] < 8 and z<z_r: 
            P_cyc[i] = 0.72 - (Y[i]-3)/(8.0-3.0)*(0.72 - 0.72*z/z_r)
        elif Y[i] > 8 and z>=z_r: 
            P_cyc[i] = 0.72
        elif 3 < Y[i] < 8 and z>=z_r: 
            P_cyc[i] = 0.72
            
    p_cyc = P_cyc*p_ult
    
    print 'S_u = %.1f psf, p_ult = %.2f lb/in' %(S_u, p_ult)
    print 'N_p = %.1f' %N_p

    if graph=='static':   
        plot(y,p, ls, label='Matlock (1970), Static')
    elif graph=='cyclic':
        plot(y,p_cyc, ls, label='Matlock (1970), Cyclic')
    elif graph=='both':
        plot(y,p, ls, label='Matlock (1970), Static')
        plot(y,p_cyc, ls, label='Matlock (1970), Cyclic')
        
    xlabel(r'y (in)'), ylabel(r'p (lb/in)'), grid(True), legend(loc='upper right')


def plot_matlock_py_normalized(Z,Z_r, J=0.5, graph='static', ls='-'):
    '''Plots nomalized Matlock (1970) p-y curves.

    Input:
    -----
    Z   - Normalized depth
    Z_r - Critical normalized depth (for cyclic loading)
    
    Output:
    ------
    Plot of Matlock (1970) p-y curve.
    '''
    #Normalized p-y curve
    Y = linspace(0,20,1000)
    P = 0.5*Y**(1.0/3.0)
    
    for i in range(0,len(Y)): 
        if P[i] > 1.0: P[i] = 1.0
    
    
    #Normalized cyclic p-y curve
    P_cyc = 0.5*Y**(1.0/3.0)
    
    for i in range(0,len(Y)): 
        if Y[i] >= 15 and Z<Z_r: 
            P_cyc[i] = 0.72*Z/Z_r
        elif 3 < Y[i] < 15 and Z<Z_r: 
            P_cyc[i] = 0.72 - (Y[i]-3)/(15.0-3.0)*(0.72 - 0.72*Z/Z_r)
        elif Y[i] > 3 and Z>=Z_r: 
            P_cyc[i] = 0.72
    
    if graph=='static':   
        plot(Y,P, ls, label='Matlock (1970), Static', lw=3)
    elif graph=='cyclic':
        plot(Y,P_cyc, ls, label='Matlock (1970), Cyclic', lw=3)
    elif graph=='both':
        plot(Y,P, ls, label='Matlock (1970), Static')
        plot(Y,P_cyc, ls, label='Matlock (1970), Cyclic')
        
    xlabel(r'$y/y_{50}$'), ylabel(r'$p/p_{ult}$'), ylim(ymax=1.1)
    grid(True), legend(loc='lower right')



def plot_jeanjean_py(D=4, S_u0=0, k=0, A=550, Y_max=1, z=0, graph='static', ls='-'):
    """ Plot Jeanjean (2009) p-y curves.

    Input:
    ------
    D       - Pile diameter (inches)
    S_u0    - Undrained shear strength at the mudline (psf)
    k       - Rate of increase of shear strength versus depth (psf/ft)
    Y_max   - Maximum lateral displacement in terms of the pile diameters (to fit plots)
    A       - Maximum shear modulus, G_max/S_u ratio
    graph   - Type of curve 'static', 'cyclic', or 'both'
    ls      - Linestyle following the PyPlot format eg: 'b-', 'g:'
    
    Note: I have calculated G_max as 550*S_u following Jeanjean (2009).
          G_max should be explicity defined if data is availale or if a better empirical model is found.

    Output:
    -------
    Plot of Jeanjean (2009) p-y curves.
    """
    
    #rcParams['figure.figsize'] = 10,8
    
    #Undrained shear stregnth profile
    S_u = S_u0 + k*(z/12.0)  #psf
    
    #Let G_max = 100*S_u (approximate empirical correlation)
    G_max = A*S_u  #psf
    
    #Normalized p-y curve
    Y = linspace(0,Y_max,100)   #Y_max can be defined in terms of the pile diameter so the curves can be made to fit plots
    P = tanh(G_max/100.0/S_u*Y**(0.5))
    
    #Calculation of N_p

    lamda = S_u0/k/(D/12.0)
    
    if lamda < 6:
        xi = 0.25 + 0.05*lamda
    else:
        xi = 0.55
    
    N_p = 12 - 4*exp(-xi*z/D)
    
    
    p_ult = N_p*(S_u/144.0)*D
    p = P*p_ult
    y = Y*D
    
    print 'N_p = %.2f' %N_p
    #plot(Y,P)
    #xlabel(r'$y/y_c$'), ylabel(r'$P/P_u$'), ylim(0,1.1), grid(True), legend(loc='lower right')
       
    if graph=='static':   
        plot(y,p, ls, label='Jeanjean (2009), Static')
    elif graph=='cyclic':
        plot(y,p_cyc, ls, label='Jeanjean (2009), Cyclic')
    elif graph=='both':
        plot(y,p, ls, label='Jeanjean (2009), Static')
        plot(y,p_cyc, ls, label='Jeanjean (2009), Cyclic')
        
    #plot(-disp[(275):(loc+275)],-p_net[0:loc]*4)
    xlabel(r'y (in)'), ylabel(r'p (lb/in)'), grid(True), legend(loc='lower right')


def plot_jeanjean_py_normalized(A=550, Y_max=1.0, S_u=1.0, graph='static', ls='-'):
    '''
    Input:
    -----
    A       - G_max/Su (default = 550)
    Y_max   - Maximum value of normalized displacement for plot
    S_u     - Undrained shear strength (only a dummy value needed)
    '''
    #Let G_max = 100*S_u (approximate empirical correlation)
    G_max = A*S_u  #psf
    
    #Normalized p-y curve
    Y = linspace(0,Y_max,100)   #Y_max can be defined in terms of the pile diameter so the curves can be made to fit plots
    P = tanh(G_max/100.0/S_u*Y**(0.5))

    if graph=='static':   
        plot(Y,P, ls, label='Jeanjean (2009), Static', lw=3)
    #elif graph=='cyclic':
    #    plot(Y,P_cyc, ls, label='Jeanjean (2009), Cyclic')
    #elif graph=='both':
    #    plot(Y,P, ls, label='Jeanjean (2009), Static')
    #    plot(Y,P_cyc, ls, label='Jeanjean (2009), Cyclic')
        
    #plot(-disp[(275):(loc+275)],-p_net[0:loc]*4)
    xlabel(r'$y/y_{50}$'), ylabel(r'$p/p_{ult}$'), ylim(ymax=1.1)
    grid(True), legend(loc='lower right')


def plot_Np(Su_0 = 0.0, k = 10.0, D=4.0):
    """ Compares Np versus depth for Matlock (1970) and Jeanjean (2009) models.

    Input:
    -----
    Su_0    - Shear strength at mudline
    k       - psf/ft
    D       - pile diameter (in)

    Output:
    ------
    Plot 1  - N_p verus z/D
    Plot 2  - p-multipliers verus depth to convert Matlock (1970) to Jeanjean (2009)"""
    
    D=float(D)/12.0 # Convert in to ft
    
    #Let Z =z/D
    Z= arange(0,20,.1)
    Su = Su_0 + k*Z*D  #Shear strength profile
    
    #Jeanjean (2009)
    if k == 0:
        lamda = 6
    else:
        lamda = Su_0/(k*D)
        
    if lamda < 6:
        xi = 0.25 + 0.05*lamda
    else:
        xi = 0.55
        
    Np_jeanjean = 12 - 4*exp(-xi*Z)
    
    #Matlock (1970)
    J = 0.5
    sigma_v = (100-62.4)*Z*D
    
    Np_matlock = 3 + sigma_v/Su + J*Z
    
    for i in range(len(Np_matlock)):
        if Np_matlock[i] > 9:
            Np_matlock[i] = 9
    
    #Plot, invert axes and reposition axis labels
    rcParams['figure.figsize'] = 12,6
    figure()
    
    subplot(1,3,1)
    plot(Np_matlock, Z, label='Matlock (1970)')
    plot(Np_jeanjean,Z, label='Jeanjean (2009)')
    xlim(xmin=0, xmax=14), xlabel('$N_p$'), ylabel('$z/D$'), grid(True)
    
    ax = gca()
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    legend(loc='lower left')
    
    subplot(1,3,2)
    plot(Np_jeanjean/Np_matlock, Z, label='Calculated'), xlim(xmin=0)
    xlabel('$N_{p,Jeanjean}/N_{p,Matlock}$'), grid(True)
    
    ax = gca()
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    subplot(1,3,3)
    plot(Su, Z)
    xlabel(r'$S_u$ (psf)'), grid(True)
    
    ax = gca()
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    legend(loc='lower left')


    
###################################################
###     CURVE FITTING TO EXPERIMENTAL DATA      ###
###################################################
    
def fit_hyperbolic_curve(filename,method='LVDT', ylim = 4,DAQ='OLD',start=0,end=0, P_0 = 1):
    """Fits a hyperblic curve to load-displacement curve from a monotonic test. 
    The fitted curve is of the form 'y = a*tanh(b*x^c)'
    
    Input:
    -----
    
    filename - Name of load-disp data file
    method   - Displacment measurement method (LVDT, LVDT_2, or LMT)
    ylim     - Range of curve, maximum lateral displacement (inches)
    DAQ      - DAQ system used old Labview program or new Labview program (the output 
               files are slight different so they need different function to read the
               data. Enter 'NEW' or 'OLD'.
    P_0      - Initial guess for pile capacity i.e. 'a'. Try guessing the 
               capacity of the pile if convergence error is encountered. Initial guesses
               for parameter 'b' and 'c' can be added from within the code.
               
    Output:
    ------
    Plot with actual data and the fitted hyperbolic curve.
    The equation of the fitted hyperbolic curve is also printed.

    Example:
    -------
    fit_hyperbolic_curve(r'D:\Box Sync\Lab Tests\2014_07_29\M4_monotonic_1rpm.N.txt',
                         method='LVDT', ylim=4, DAQ='OLD')
    """
    
    from scipy import optimize
    from Monopile_Project import mpile 
    
    
    def fun(x,a,b,c):
        #Hyperbolic function
        y = a*tanh(b*x**c)
        return y
    
    if DAQ=='OLD':
        [t, load, disp, zero_load, loc2] = mpile.read_data_file_1(filename, method=method)
    elif DAQ=='NEW':
        [t, load, disp, zero_load, loc2] = mpile.read_data_file_2(filename, method=method)
    
    loc1 = start
    
    if end != 0: loc2 = end
        
    rcParams['figure.figsize'] = 12,8
    
    plot(abs(disp[0:loc2]),abs(load[0:loc2]),':', label='Measured')
    
    x1 = linspace(0,ylim,100)

    fitpars, covmat = optimize.curve_fit(fun, abs(disp[loc1:loc2]),abs(load[loc1:loc2]), p0=[P_0,.1,.1])
    y1 = fun(x1,fitpars[0],fitpars[1],fitpars[2])
    
    plot(x1,y1)
    print 'y = %.2ftanh(%.2fx**%.2f)' %(fitpars[0],fitpars[1],fitpars[2])
    #print 'sigma_yy = %.2f,%.2f,%.2f' %(sqrt(diag(covmat)))
    print 'Standard deviation of each parameter,' + str(sqrt(diag(covmat))) + '\n'                      
    xlabel('Displacement (in)'), ylabel('Load (lb)'), grid(True)



##########################
#### Free-fall Tests #####
##########################

def free_fall(filename, cal_load=-1, cal_disp=1, x1=1, x2=100, ls='-'):
    """Read data from Yunhan's new DAQ VI output files with both engineering units and raw voltages. 

    Input:
    -----
    filename - Location of raw data file 
    cal_load - Load cell calibration factor 
    cal_disp - LVDT1/LMT calibration factor
    x1, x2   - Start/End points
    ls       - Linestyle for Matplotlib plots

    Output:
    ------
    Displacement versus time plot of a free-fall test
    """
    
    import datetime as datetime
    import time as time
    
    file1 = open(filename, 'r')

    method = 'LMT'
    
    if   method == 'LVDT':   column = 5;
    elif method == 'LMT':    column = 7;
    elif method == 'LVDT_2': column = 9;
        
    #Initialize arrays in which to store extracted data
    time_stamp = array([datetime.timedelta(hours=i) for i in xrange(100000)])
    data_load  = zeros((10000))
    data_disp  = zeros((10000))
    data_time  = zeros((10000))
    
    loc  = 0               
    for line in file1:
        if line.startswith('D')==True:
            temp = line.split()
            
            #Extract Timestamps
            [mo,dy,yr] = [int(x) for x in temp[1].split('/')]
            [hh,mm,ss] = temp[2].split(':')
            [ss,micro] = [int(x) for x in ss.split('.')]
            time_stamp[loc] = datetime.datetime(yr,mo,dy,int(hh),int(mm),ss,micro*10000) #Works with timestamp with a precision of 2 decimals
            
            #Extract Load Cell and LVDT/LMT Voltage
            data_load[loc] = float(temp[3])
            data_disp[loc] = float(temp[column])
            loc += 1                           #Counter to keep track of the last data entry into the data arrays
    
            
    #Identify zero errors 
    zero_time = time_stamp[0]
    zero_load = data_load[0]
    zero_disp = data_disp[0]
    
    #Convert Timestamps in to Seconds, Remove zero errors in Load & Disp, Apply Calibration Factors
    for i in range(loc):
        data_time[i] = (time_stamp[i] - time_stamp[0]).total_seconds()
        data_load[i] = (data_load[i] - zero_load)*cal_load #- 1.6
        data_disp[i] = (data_disp[i] - zero_disp)*cal_disp
        
    file1.close()
    
    print 'Number of data points = %d' %loc
    
    plot(data_time[x1:x2], data_disp[x1:x2], 'o')
    xlabel('Time (s)'), ylabel('Displacement (in)')
    
    inst_vel = zeros(1000)
    inst_acc = zeros(1000)
    
    for i in range(0,loc):
        dt = data_time[i+1] - data_time[i]
        
        ds_1 = data_disp[i+1] - data_disp[i]
        ds_2 = data_disp[i+2] - data_disp[i+1]
        ds_3 = data_disp[i+3] - data_disp[i+2]
        
        inst_vel[i] = ds_1/dt
        inst_acc[i] = (ds_1 - 2*ds_2 + ds_3)/dt**2
    
    subplot(3,1,1)
    plot(data_time[x1:x2], data_disp[x1:x2], ls)
    xlabel('Time (s)'), ylabel('Displacement (in)')
    
    subplot(3,1,2)
    plot(data_time[x1:x2], inst_vel[x1:x2], ls)
    xlabel('Time (s)'), ylabel('Velocity (in/s)')
    
    subplot(3,1,3)
    plot(data_time[x1:x2], inst_acc[x1:x2], ls)
    xlabel('Time (s)'), ylabel('Acceleration (in/s/s)')
