import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from numpy import array
import sys 

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

def obtain_variable(filename, var_name, var_type):
    f = open(filename, "r")
    lines = f.read().splitlines()
    for i in range(0,len(lines)):
        if lines[i] == var_name:
            if var_type == 'float':
                var = float(lines[i+1])
            elif var_type == 'int':
                var = int(lines[i+1])
            elif var_type == 'string':
                var = lines[i+1]
            else:
                sys.exit('Error: Invalid var_type in function obtain_variable')
    f.close()
    return var;

# Read variables from info file
flow_2D = obtain_variable('.simulation_info','flow_2D','string')

# Load Data
x = readData('x_nodes.csv')
y = readData('y_nodes.csv')

u_vec = readData('u.csv')
v_vec = readData('v.csv')
p_vec = readData('p.csv')

Iter = readData('Iterations.csv')
Err_u = readData('Err_u.csv')
Err_v = readData('Err_v.csv')
Err_p = readData('Err_p.csv')

if flow_2D == 'false':
    z = readData('z_nodes.csv')
    w_vec = readData('w.csv')
    Err_w = readData('Err_w.csv')

# Transpose variable size for ploting
if flow_2D == 'false':
    u = convert1Dto3D(u_vec,len(x),len(y),len(z))
    v = convert1Dto3D(v_vec,len(x),len(y),len(z))
    w = convert1Dto3D(w_vec,len(x),len(y),len(z))
    p = convert1Dto3D(p_vec,len(x),len(y),len(z))
elif flow_2D == 'true':
    u = convert1Dto2D(u_vec,len(x),len(y))
    v = convert1Dto2D(v_vec,len(x),len(y))
    p = convert1Dto2D(p_vec,len(x),len(y))
    u = u.transpose()
    v = v.transpose()
    p = p.transpose()

# Save csv data files
#np.savetxt('p_drop.csv', p[iround((p.shape[0]-1)/2.0),:], delimiter=',')
#np.savetxt('u_axial.csv', u[iround(u.shape[0]/2.0)-1,:], delimiter=',')
#np.savetxt('u_profile.csv', u[:,iround(u.shape[1])-1], delimiter=',')

# --------------------------------------------------------------------------------------
# Main Plots
# --------------------------------------------------------------------------------------
if flow_2D == 'true':
    # Velocity Contour Plot
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(x, y, u, 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(x, y, u, 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(x), max(x)+0.1))
    plt.yticks(np.arange(min(y), max(y)+0.1))
    plt.title('Velocity u', **font_title)
    plt.xlabel('x', **font_labels)
    plt.ylabel('y', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    ax1 = fig.add_subplot(1, 2, 2,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(x, y, v, 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(x, y, v, 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(x), max(x)+0.1))
    plt.yticks(np.arange(min(y), max(y)+0.1))
    plt.title('Velocity v', **font_title)
    plt.xlabel('x', **font_labels)
    plt.ylabel('y', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Pressure Contour Plot
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(x, y, p, 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(x, y, p, 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(x), max(x)+0.1))
    plt.yticks(np.arange(min(y), max(y)+0.1))
    plt.title('Pressure p', **font_title)
    plt.xlabel('x', **font_labels)
    plt.ylabel('y', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Exit velocity profile
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    plt.plot(u[:,iround(u.shape[1]/2.0)-1],y)
    ax1.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    # ax1.set_xlim(min(u[:,iround(u.shape[1]/2.0)-1]), max(u[:,iround(u.shape[1]/2.0)-1]))
    ax1.set_xlim(-1.0, 1.0)
    ax1.set_ylim(min(y), max(y))
    plt.title('Velocity profile - Vertical', **font_title)
    plt.xlabel('u', **font_labels)
    plt.ylabel('y', **font_labels)
    fig.tight_layout()
    # Axial velocity profile
    ax1 = fig.add_subplot(1, 2, 2)
    plt.plot(x, v[iround(v.shape[0]/2.0)-1,:])
    ax1.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax1.set_xlim(min(x), max(x))
    # ax1.set_ylim(min(v[iround(v.shape[0]/2.0)-1,:]), max(v[iround(v.shape[0]/2.0)-1,:]))
    ax1.set_ylim(-1.0,1.0)
    plt.title('Velocity profile - horizontal', **font_title)
    plt.xlabel('x', **font_labels)
    plt.ylabel('v', **font_labels)
    fig.tight_layout()
elif flow_2D == 'false':
    # Velocity u Contour Plot in xz plane
    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(z, x, u[:,iround(u.shape[1]/2.0)-1,:], 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(z, x, u[:,iround(u.shape[1]/2.0)-1,:], 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(z), max(z)+0.1))
    plt.yticks(np.arange(min(x), max(x)+0.1))
    plt.title('Velocity u - xz plane', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('x', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Velocity v Contour Plot in xz plane
    ax1 = fig.add_subplot(3, 1, 2,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(z, x, v[:,iround(v.shape[1]/2.0)-1,:], 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(z, x, v[:,iround(v.shape[1]/2.0)-1,:], 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(z), max(z)+0.1))
    plt.yticks(np.arange(min(x), max(x)+0.1))
    plt.title('Velocity v - xz plane', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('x', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Velocity w Contour Plot in xz plane
    ax1 = fig.add_subplot(3, 1, 3,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(z, x, w[:,iround(w.shape[1]/2.0)-1,:], 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(z, x, w[:,iround(w.shape[1]/2.0)-1,:], 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(z), max(z)+0.1))
    plt.yticks(np.arange(min(x), max(x)+0.1))
    plt.title('Velocity w - xz plane', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('x', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Velocity u Contour Plot in yz plane
    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(z, y, u[iround(u.shape[0]/2.0)-1,:,:], 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(z, y, u[iround(u.shape[0]/2.0)-1,:,:], 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(z), max(z)+0.1))
    plt.yticks(np.arange(min(y), max(y)+0.1))
    plt.title('Velocity u - yz plane', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('y', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Velocity v Contour Plot in yz plane
    ax1 = fig.add_subplot(3, 1, 2,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(z, y, v[iround(v.shape[0]/2.0)-1,:,:], 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(z, y, v[iround(v.shape[0]/2.0)-1,:,:], 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(z), max(z)+0.1))
    plt.yticks(np.arange(min(y), max(y)+0.1))
    plt.title('Velocity v - yz plane', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('y', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Velocity w Contour Plot in yz plane
    ax1 = fig.add_subplot(3, 1, 3,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(z, y, w[iround(w.shape[0]/2.0)-1,:,:], 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(z, y, w[iround(w.shape[0]/2.0)-1,:,:], 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(z), max(z)+0.1))
    plt.yticks(np.arange(min(y), max(y)+0.1))
    plt.title('Velocity w - yz plane', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('y', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Pressure Contour Plot in xz plane
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(z, x, p[:,iround(p.shape[1]/2.0)-1,:], 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(z, x, p[:,iround(p.shape[1]/2.0)-1,:], 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(z), max(z)+0.1))
    plt.yticks(np.arange(min(x), max(x)+0.1))
    plt.title('Pressure p - xz plane', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('x', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Pressure Contour Plot in yz plane
    ax1 = fig.add_subplot(2, 1, 2,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(z, y, p[iround(p.shape[0]/2.0)-1,:,:], 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(z, y, p[iround(p.shape[0]/2.0)-1,:,:], 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(z), max(z)+0.1))
    plt.yticks(np.arange(min(y), max(y)+0.1))
    plt.title('Pressure p - yz plane', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('y', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Exit Velocity Profiles
    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1)
    plt.plot(y, w[iround(w.shape[0]/2.0)-1,:,w.shape[2]-1])
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax.set_xlim(min(y), max(y))
    ax.set_ylim(min(w[iround(w.shape[0]/2.0)-1,:,w.shape[2]-1]), max(w[iround(w.shape[0]/2.0)-1,:,w.shape[2]-1])+0.001)
    plt.xticks(np.arange(min(x), max(x)))
    plt.title('Exit Velocity - y axis', **font_title)
    plt.xlabel('y', **font_labels)
    plt.ylabel('w', **font_labels)
    fig.tight_layout()
    ax = fig.add_subplot(1, 2, 2)
    plt.plot(x, w[:,iround(w.shape[1]/2.0)-1,w.shape[2]-1])
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax.set_xlim(min(x), max(x))
    ax.set_ylim(min(w[:,iround(w.shape[1]/2.0)-1,w.shape[2]-1]), max(w[:,iround(w.shape[1]/2.0)-1,w.shape[2]-1])+0.001)
    plt.xticks(np.arange(min(y), max(y)))
    plt.title('Exit Velocity - x axis', **font_title)
    plt.xlabel('x', **font_labels)
    plt.ylabel('w', **font_labels)
    fig.tight_layout()
    # Exit Velocity and Magnetic Contour
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(x, y, np.transpose(w[:,:,w.shape[2]-1]), 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(x, y, np.transpose(w[:,:,w.shape[2]-1]), 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(x), max(x)))
    plt.yticks(np.arange(min(y), max(y)))
    plt.title('Exit Velocity', **font_title)
    plt.xlabel('x', **font_labels)
    plt.ylabel('y', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Velocity axial and pressure drop plot
    fig = plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    plt.plot(z, w[iround(w.shape[0]/2.0)-1,iround(w.shape[1]/2.0)-1,:])
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax.set_xlim(min(z), max(z))
    ax.set_ylim(min(w[iround(w.shape[0]/2.0)-1,iround(w.shape[1]/2.0)-1,:]), max(w[iround(w.shape[0]/2.0)-1,iround(w.shape[1]/2.0)-1,:])+0.001)
    plt.xticks(np.arange(min(z), max(z)))
    plt.title('Axial Velocity', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('w', **font_labels)
    fig.tight_layout()
    ax = fig.add_subplot(2, 1, 2)
    plt.plot(z, p[iround(p.shape[0]/2.0)-1,iround(p.shape[1]/2.0)-1,:])
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax.set_xlim(min(z), max(z))
    ax.set_ylim(min(p[iround(p.shape[0]/2.0)-1,iround(p.shape[1]/2.0)-1,:]), max(p[iround(p.shape[0]/2.0)-1,iround(p.shape[1]/2.0)-1,:])+0.001)
    plt.xticks(np.arange(min(z), max(z)))
    plt.title('Pressure Drop', **font_title)
    plt.xlabel('z', **font_labels)
    plt.ylabel('p', **font_labels)
    fig.tight_layout()
    # Inlet Velocity a
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1,  aspect='equal')
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    PC = ax1.contourf(x, y, np.transpose(w[:,:,0]), 10, cmap='jet',antialiased=True)
    CF =  ax1.contour(x, y, np.transpose(w[:,:,0]), 10, colors='k',linewidths=0.6)
    plt.xticks(np.arange(min(x), max(x)+0.1))
    plt.yticks(np.arange(min(y), max(y)+0.1))
    plt.title('Inlet Velocity', **font_title)
    plt.xlabel('x', **font_labels)
    plt.ylabel('y', **font_labels)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="2%", pad=0.2)
    cbar = plt.colorbar(PC, cax = cax1)
    fig.tight_layout()
    # Inlet Velocity Profiles
    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1)
    plt.plot(y, w[iround(w.shape[0]/2.0)-1,:,0])
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax.set_xlim(min(y), max(y))
    ax.set_ylim(min(w[iround(w.shape[0]/2.0)-1,:,0]), max(w[iround(w.shape[0]/2.0)-1,:,0])+0.001)
    plt.xticks(np.arange(min(x), max(x)))
    plt.title('Inlet Velocity - y axis', **font_title)
    plt.xlabel('y', **font_labels)
    plt.ylabel('w', **font_labels)
    fig.tight_layout()
    ax = fig.add_subplot(1, 2, 2)
    plt.plot(x, w[:,iround(w.shape[1]/2.0)-1,0])
    ax.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax.set_xlim(min(x), max(x))
    ax.set_ylim(min(w[:,iround(w.shape[1]/2.0)-1,0]), max(w[:,iround(w.shape[1]/2.0)-1,0])+0.001)
    plt.xticks(np.arange(min(y), max(y)))
    plt.title('Inlet Velocity - x axis', **font_title)
    plt.xlabel('x', **font_labels)
    plt.ylabel('w', **font_labels)
    fig.tight_layout()

# --------------------------------------------------------------------------------------
# Momentum error plots
# --------------------------------------------------------------------------------------
if flow_2D == 'true':
    # Residual u
    fig = plt.figure()
    plt.subplots_adjust(left=0.2, hspace=1.3)
    ax1 = fig.add_subplot(3, 1, 1)
    plt.semilogy(Iter, Err_u)
    ax1.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax1.set_xlim(0, max(Iter))
    ax1.set_ylim(min(Err_u), max(Err_u))
    plt.title('Residual u', **font_title)
    plt.xlabel('Time Steps', **font_labels)
    plt.ylabel('Residual', **font_labels)
    # Residual v
    ax1 = fig.add_subplot(3, 1, 2)
    plt.semilogy(Iter, Err_v)
    ax1.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax1.set_xlim(0, max(Iter))
    ax1.set_ylim(min(Err_v), max(Err_v))
    plt.title('Residual v', **font_title)
    plt.xlabel('Time Steps', **font_labels)
    plt.ylabel('Residual', **font_labels)
    # Residual p
    ax1 = fig.add_subplot(3, 1, 3)
    plt.semilogy(Iter, Err_p)
    ax1.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax1.set_xlim(0, max(Iter))
    ax1.set_ylim(min(Err_p), max(Err_p))
    plt.title('Residual p', **font_title)
    plt.xlabel('Time Steps', **font_labels)
    plt.ylabel('Residual', **font_labels)
elif flow_2D == 'false':
    # Residual u
    fig = plt.figure()
    plt.subplots_adjust(left=0.2, hspace=1.3)
    ax1 = fig.add_subplot(3, 1, 1)
    plt.semilogy(Iter, Err_u)
    ax1.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax1.set_xlim(0, max(Iter))
    ax1.set_ylim(min(Err_u), max(Err_u))
    plt.title('Residual u', **font_title)
    plt.xlabel('Time Steps', **font_labels)
    plt.ylabel('Residual', **font_labels)
    # Residual v
    ax1 = fig.add_subplot(3, 1, 2)
    plt.semilogy(Iter, Err_v)
    ax1.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax1.set_xlim(0, max(Iter))
    ax1.set_ylim(min(Err_v), max(Err_v))
    plt.title('Residual v', **font_title)
    plt.xlabel('Time Steps', **font_labels)
    plt.ylabel('Residual', **font_labels)
    # Residual w
    ax1 = fig.add_subplot(3, 1, 3)
    plt.semilogy(Iter, Err_w)
    ax1.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax1.set_xlim(0, max(Iter))
    ax1.set_ylim(min(Err_w), max(Err_w))
    plt.title('Residual w', **font_title)
    plt.xlabel('Time Steps', **font_labels)
    plt.ylabel('Residual', **font_labels)
    # Residual p
    fig = plt.figure()
    plt.semilogy(Iter, Err_p)
    ax1.grid(color='lightgrey', linestyle='-', linewidth=0.5)
    ax1.set_xlim(0, max(Iter))
    ax1.set_ylim(min(Err_p), max(Err_p))
    plt.title('Residual p', **font_title)
    plt.xlabel('Time Steps', **font_labels)
    plt.ylabel('Residual', **font_labels)

plt.show()
