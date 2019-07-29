#!/usr/bin/env python
'''
Program for finding the steady state voltage distribution and capacitance
per unit length of a shielded microstrip. Solves laplace equation using a
Finite Difference Method (FDM). Microstrip has outer shielding at 0V, and
a conducting strip (at V_0) on top of a dielectric material.
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import epsilon_0
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

def plot_3d(ms):
    '''
    Create a 3D plot of the voltage distribution inside a microstrip line.
    '''
    fig1 = plt.figure()
    ax1 = fig1.gca(projection='3d')
    X, Y = np.meshgrid(np.arange(ms.W), np.arange(ms.H))

    # Add plot to figure and change axis labels
    ax1.plot_surface(X, Y, ms.V, cmap=plt.get_cmap('jet'), rstride=1,
                     cstride=1, linewidth=0, antialiased=False)
    ax1.set_zlim(0, ms.V_0)
    ax1.zaxis.set_major_locator(LinearLocator(5))
    ax1.zaxis.set_major_formatter(FormatStrFormatter("%.2f"))

def plot_contour(ms):
    '''
    Create a contour plot for the voltage distribution.
    '''
    __, ax2 = plt.subplots()
    X, Y = np.meshgrid(np.arange(ms.W), np.arange(ms.H))

    # Create the axis tick labels to display physical units
    xlabels = ["{0:.1f}".format(x) for x in
               np.linspace(0, ms.W/ms.g_ratio, num=8, endpoint=True)]
    ylabels = ["{0:.1f}".format(y) for y in
               np.linspace(0, ms.H/ms.g_ratio, num=8, endpoint=True)]

    # Create contour plot and add reference line for dielectric
    plt.contourf(X, Y, ms.V, levels=40, cmap=plt.get_cmap('jet'))
    plt.hlines(ms.h, 0, ms.W-1, linestyles='dashed', alpha=0.5)

    # Change tick labels and add info text to plot
    ax2.set_xticks(np.linspace(0, ms.W-1, num=8, endpoint=True))
    ax2.set_xticklabels(xlabels)
    ax2.set_yticks(np.linspace(0, ms.H-1, num=8, endpoint=True))
    ax2.set_yticklabels(ylabels)
    plt.xlabel("Width ({})".format(ms.units))
    plt.ylabel("Height ({})".format(ms.units))
    plt.axis('scaled')
    plt.text(ms.W-12, ms.H-10, r"$\epsilon_0$", fontsize=16, color='w')
    plt.text(ms.W-16, 6, r"$\epsilon_r\epsilon_0$", fontsize=16, color='w')
    plt.colorbar().set_label("Voltage")

class Microstrip:
    '''
    Class for storing Microstrip properties and FDM function.
    '''
    def __init__(self, V_0=1, eps_r=1, max_res=0.001, **kwargs):
        # Parameters for FDM calculation
        self.V_0 = V_0
        self.eps_r = eps_r
        self.max_res = max_res

        # Add microstrip dimensions to object
        self.__dict__.update(kwargs)

        # Convert from physical units to grid units
        self.g_ratio = 100/self.W
        self.W = 100
        self.H = round(self.H * self.g_ratio)
        self.h = round(self.h * self.g_ratio)
        self.w = round(self.w * self.g_ratio)
        self.t = round(self.t * self.g_ratio)
        self.V = np.zeros((self.H, self.W))
        self.cpul = None

        # Conducting strip is usually centered
        if self.x is None:
            self.x = round((self.W - self.w)/2)
        else:
            self.x = round(self.x * self.g_ratio)

    def fdm(self):
        '''
        Funcion for finding the voltage distribution of a shielded microstrip.
        Uses an iterative process to find the voltage values, and stops when
        the residual between iterations falls below the desired tolerance.
        '''
        V = self.V
        inc_node = np.ones((self.H, self.W), dtype=bool)
        slice_rows = slice(self.h, self.h+self.t)
        slice_cols = slice(self.x, self.x+self.w)

        # Determine which nodes are part of the conducting strip
        # and can be ignored in the calculation
        V[slice_rows, slice_cols] = self.V_0
        inc_node[slice_rows, slice_cols] = False

        while True:
            largest_res = 0
            for i in range(1, self.H-1):
                for j in range(1, self.W-1):
                    if inc_node[i, j]:
                        V_old = V[i, j]

                        # Average the voltages of adjacent nodes:
                        # left, right, above, and below
                        V_new = 0.25 * (V[i][j+1] + V[i][j-1])
                        if i == self.h:
                            # Change calculation for nodes on the
                            # air-dielectric boundary
                            V_new += (1/(2+2*self.eps_r)) * V[i+1][j]
                            V_new += (self.eps_r/(2+2*self.eps_r)) * V[i-1][j]
                        else:
                            V_new += 0.25 * V[i+1][j]
                            V_new += 0.25 * V[i-1][j]

                        V[i, j] = V_new
                        residual = abs(V_new - V_old)
                        if residual > largest_res:
                            # Find the largest deviation for this iteration
                            largest_res = residual
            if largest_res < self.max_res:
                # Finish when allowed tolerance is met
                self.V = V
                break

    def compute_cpul(s):
        '''
        Function for computing the capacitance per unit length of
        the shielded microstrip. Find the charge per unit length
        of the conducting strip using Gauss's Law (in 2D). Then use
        the equation C=Q*V to find capacitance (per unit length).
        '''
        cpul = 0
        V = s.V

        # Process nodes at the air-dielectric boundary first.
        cpul += 0.5 * (s.eps_r+1) * (V[s.h, s.x] - V[s.h, s.x-1])
        cpul += 0.5 * (s.eps_r+1) * (V[s.h, s.x+s.w-1] - V[s.h, s.x+s.w])
        for i in range(s.h+1, s.h+s.t):
            cpul += V[i, s.x] - V[i, s.x-1]
            cpul += V[i, s.x+s.w-1] - V[i, s.x+s.w]
        for j in range(s.x, s.x+s.w):
            cpul += s.eps_r * (V[s.h, j] - V[s.h-1, j])
            cpul += V[s.h+s.t-1, j] - V[s.h+s.t, j]
        cpul *= epsilon_0/s.V_0

        # Store value and print to screen.
        s.cpul = cpul
        print("C/L = {:.4e} F/m".format(s.cpul))

DIMS = {
    'units': 'cm',
    'W': 4,
    'H': 3,
    'h': 1,
    'w': 1.5,
    't': 0.4,
    'x': None
}

if __name__ == '__main__':
    M1 = Microstrip(V_0=5, eps_r=1, max_res=0.0001, **DIMS)

    M1.fdm()
    M1.compute_cpul()
    
    plot_contour(M1)
    plt.show()
