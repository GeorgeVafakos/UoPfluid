import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from numpy import array
import sys 
from matplotlib.ticker import AutoMinorLocator

# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------

# Clear Console (only when executed with Spyder IDE)
#print(chr(27) + "[2J")

# Fonts Configuration
font_title = {'family':'sans-serif', 'weight':'bold', 'size':14}
font_labels = {'family':'sans-serif', 'weight':'normal', 'size':14}

# Scientific notation in colorbar
def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

# Round to the nearest integer (and not only enarest even integer)/
def iround(x):
   return int(round(x) - .5) + (x > 0)

# Function that converts the 1D data into 2D array
def convert1Dto2D(u_vec,Ni,Nj):
    u = np.empty((Ni,Nj))
    for j in range(0,Nj+1):
        for i in range(0,Ni+1):
            u[i-1,j-1]=u_vec[i+(j-1)*Ni-1]
    return u;

# Function that converts the 1D data into 3D array
def convert1Dto3D(u_vec,Ni,Nj,Nk):
    u = np.empty((Ni,Nj,Nk))
    for k in range(0,Nk+1):
        for j in range(0,Nj+1):
            for i in range(0,Ni+1):
                u[i-1,j-1,k-1]=u_vec[i+(j-1)*Ni+(k-1)*Ni*Nj-1]
    return u;

# Function to read 1D data where every line contains a vector's element
def readData(filename):
    return np.array(open(filename).read().splitlines()).astype(np.float)

# Load Data
x = readData('x_nodes.csv')
y = readData('y_nodes.csv')
psi_vec = readData('psi.csv')
psi = convert1Dto2D(psi_vec,len(x),len(y))
psi = psi.transpose()

# Save csv data files
#np.savetxt('p_drop.csv', p[iround((p.shape[0]-1)/2.0),:], delimiter=',')
#np.savetxt('u_axial.csv', u[iround(u.shape[0]/2.0)-1,:], delimiter=',')
#np.savetxt('u_profile.csv', u[:,iround(u.shape[1])-1], delimiter=',')

# Define Streamlines
Strm_a=-1.0e-10
Strm_b=-1.0e-7
Strm_c=-1.0e-5
Strm_d=-1.0e-4
Strm_e=-0.01
Strm_f=-0.03
Strm_g=-0.05
Strm_h=-0.07
Strm_i=-0.09
Strm_j=-0.1
Strm_k=-0.11
Strm_l=-0.115
Strm_m=-0.1175
Strm_0=1.0e-8
Strm_1=1.0e-7
Strm_2=1.0e-6
Strm_3=1.0e-5
Strm_4=5.0e-5
Strm_5=1.0e-4
Strm_6=2.5e-4
Strm_7=5.0e-4
Strm_8=1.0e-3
Strm_9=1.5e-3
Strm_10=3.0e-3

# Srteamfunction levels (streamlines)
streamlevels100=  [Strm_j, Strm_i, Strm_h, Strm_g, Strm_f, Strm_e, Strm_d, Strm_c, Strm_b, Strm_1, Strm_2, Strm_3]
streamlevels1000= [-0.12, -0.116, Strm_k, Strm_j, Strm_i, Strm_h, Strm_g, Strm_f, Strm_e, -0.003, Strm_d, Strm_b, Strm_2, Strm_3, 7.0e-5, 2.0e-4, Strm_7, Strm_8, Strm_9]
streamlevels3200= [-0.121, -0.116, Strm_k, Strm_j, Strm_i, Strm_h, Strm_g, Strm_f, Strm_e, -0.003, Strm_d, Strm_c, -5.2e-8,  Strm_2, Strm_3, Strm_5, Strm_6, Strm_7, Strm_8, Strm_9, 2.3e-3, 2.8e-3]
streamlevels5000= [-0.121, -0.116, Strm_k, Strm_j, Strm_i, Strm_h, Strm_g, Strm_f, Strm_e, -0.003, Strm_d, Strm_c, -1.5e-6, -2.0e-7, Strm_2,  Strm_3, Strm_5, Strm_6, Strm_7, Strm_8, Strm_9, 2.3e-3, 2.8e-3]

# Velocity Contour Plot
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1,  aspect='equal')
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
CF =  ax1.contour(x, y, psi, streamlevels5000, colors='k',linewidths=0.8)
plt.xticks(np.arange(min(x), max(x)+0.1, 0.5))
plt.yticks(np.arange(min(y), max(y)+0.1, 0.5))
ax1.set_xticks([0.25, 0.75], minor=True)
ax1.set_yticks([0.25, 0.75], minor=True)
#ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=3)
ax1.tick_params(axis='both',direction='out', length=5)
plt.title('Re=5000', **font_title)
plt.xlabel('x', **font_labels)
plt.ylabel('y', **font_labels)
divider = make_axes_locatable(ax1)
fig.tight_layout()
plt.savefig('CavityStreamRe5000.eps', format='eps', dpi=1000, transparent=True, bbox_inches = 'tight', pad_inches = 0) # print graph


#plt.show()
