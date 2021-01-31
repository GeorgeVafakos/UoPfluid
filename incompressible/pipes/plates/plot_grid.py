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

# Round to the nearest integer (and not on;y enarest even integer)
def iround(x):
   return int(round(x) - .5) + (x > 0)

# Function that converts the 1D data into 3D array
def convert1Dto2D(u_vec,Ni,Nj):
    u = np.empty((Nj,Ni))
    for j in range(0,Nj+1):
        for i in range(0,Ni+1):
            u[j-1,i-1]=u_vec[i+(j-1)*Ni-1]
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

# Load Data
x_nodes = readData('x_nodes.csv')
y_nodes = readData('y_nodes.csv')
z_nodes = readData('z_nodes.csv')
x_faces = readData('x_faces.csv')
y_faces = readData('y_faces.csv')
z_faces = readData('z_faces.csv')


# Grid Plot - Cross Section
fig = plt.figure(tight_layout='true')
ax1 = fig.add_subplot(1, 1, 1,  aspect='equal')

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['bottom'].set_visible(False)

for i in range(0, len(y_faces)):
    plt.plot(x_faces, y_faces[i]*np.ones(len(x_faces)), 'k-', linewidth=0.1)
    
for i in range(0, len(x_faces)):
    plt.plot(x_faces[i]*np.ones(len(y_faces)), y_faces, 'k-', linewidth=0.1)

plt.xticks(np.arange(min(x_faces), max(x_faces)+0.1, 0.5))
plt.yticks(np.arange(min(y_faces), max(y_faces)+0.1, 0.5))
plt.title('Grid - Cross Section', **font_title)
plt.xlabel('x', **font_labels)
plt.ylabel('y', **font_labels)


# Grid Plot - Axial Plane
fig = plt.figure(tight_layout='true')
ax1 = fig.add_subplot(1, 1, 1,  aspect='equal')

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['bottom'].set_visible(False)

for i in range(0, len(y_faces)):
    plt.plot(z_faces, y_faces[i]*np.ones(len(z_faces)), 'k-', linewidth=0.1)
    
for i in range(0, len(z_faces)):
    plt.plot(z_faces[i]*np.ones(len(y_faces)), y_faces, 'k-', linewidth=0.1)

plt.xticks(np.arange(min(z_faces), max(z_faces)+0.1, 0.5))
plt.yticks(np.arange(min(y_faces), max(y_faces)+0.1, 0.5))
plt.title('Grid - Axial', **font_title)
plt.xlabel('z', **font_labels)
plt.ylabel('y', **font_labels)



plt.show()
