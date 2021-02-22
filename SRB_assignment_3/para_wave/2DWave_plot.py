import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

M = 2306
Z1 = np.loadtxt("U_horiz_t0.csv")
Z2 = np.loadtxt("U_horiz_t0.33.csv")
Z3 = np.loadtxt("U_horiz_t0.67.csv")
Z4 = np.loadtxt("U_horiz_t1.00.csv")

Z5 = np.loadtxt("U_square_t0.csv")
Z6 = np.loadtxt("U_square_t0.33.csv")
Z7 = np.loadtxt("U_square_t0.67.csv")
Z8 = np.loadtxt("U_square_t1.00.csv")

Z9 = np.loadtxt("U_vert_t0.csv")
Z10 = np.loadtxt("U_vert_t0.33.csv")
Z11= np.loadtxt("U_vert_t0.67.csv")
Z12 = np.loadtxt("U_vert_t1.00.csv")

dx = 2/(M-1)
X = np.linspace(-1, 1, M)
Y = np.linspace(-1, 1, M)
X, Y = np.meshgrid(X, Y)

fig1 = plt.figure()
ax1 = fig1.gca(projection = '3d')
surf = ax1.plot_surface(X, Y, Z1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax1.set_title('Horizontal strips, t=0')

fig2 = plt.figure()
ax2 = fig2.gca(projection='3d')
surf = ax2.plot_surface(X, Y, Z2, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax2.set_title('Horizontal strips, t=0.33')

fig3 = plt.figure()
ax3 = fig3.gca(projection='3d')
surf = ax3.plot_surface(X, Y, Z3, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax3.set_title('Horizontal strips, t=0.67')

fig4 = plt.figure()
ax4 = fig4.gca(projection='3d')
surf = ax4.plot_surface(X, Y, Z4, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax4.set_title('Horizontal strips, t=1')

fig5 = plt.figure()
ax5 = fig5.gca(projection='3d')
surf = ax5.plot_surface(X, Y, Z5, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax5.set_title('Squares, t=0')

fig6 = plt.figure()
ax6 = fig6.gca(projection='3d')
surf = ax6.plot_surface(X, Y, Z6, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax6.set_title('Squares, t=0.33')

fig7 = plt.figure()
ax7 = fig7.gca(projection='3d')
surf = ax7.plot_surface(X, Y, Z7, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax3.set_title('Squares, t=0.67')

fig8 = plt.figure()
ax8 = fig8.gca(projection='3d')
surf = ax8.plot_surface(X, Y, Z8, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax8.set_title('Squares, t=1')

fig9 = plt.figure()
ax9 = fig9.gca(projection='3d')
surf = ax9.plot_surface(X, Y, Z9, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax9.set_title('Vertical strips, t=0')

fig10 = plt.figure()
ax10 = fig10.gca(projection='3d')
surf = ax10.plot_surface(X, Y, Z10, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax10.set_title('Vertical strips, t=0.33')

fig11 = plt.figure()
ax11 = fig11.gca(projection='3d')
surf = ax11.plot_surface(X, Y, Z11, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax1.set_title('Vertical strips, t=0.67')

fig12 = plt.figure()
ax12 = fig12.gca(projection='3d')
surf = ax12.plot_surface(X, Y, Z12, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax12.set_title('Vertical strips, t=1')