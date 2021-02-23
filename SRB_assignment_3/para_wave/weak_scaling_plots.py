import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def speed_up(a):
    speed_up =  a[0]/a
    return speed_up
    
P = np.array([1,4,16,24,36])
P_square = np.array([1,4,16,25,36])

run_time_square = np.array([1143.29, 298.748,93.9222, 84.9656, 46.224])
run_time_vert = np.array([1140.44,300.632,95.7239,91.0576,48.556])
run_time_horiz = np.array([1028.21,270.011,90.9404,88.1855,46.192])


run_time_no_io_square =np.array([294.812,85.6856,76.6362,37.9617])
run_time_no_io_vert =np.array([285.008,87.3275,81.764,39.3733])
run_time_no_io_horiz =np.array([254.679,82.0836,78.5395,36.826])

io_frac_square = (run_time_square[1:5]-run_time_no_io_square)/run_time_square[1:5]
io_frac_vert = (run_time_vert[1:5]-run_time_no_io_vert)/run_time_vert[1:5]
io_frac_horiz = (run_time_horiz[1:5]-run_time_no_io_horiz)/run_time_horiz[1:5]


speed_up_square =  speed_up(run_time_square)
speed_up_horiz =  speed_up(run_time_horiz)
speed_up_vert =  speed_up(run_time_vert)

parallel_efficiency_square = speed_up_square/P_square
parallel_efficiency_vert = speed_up_vert/P
parallel_efficiency_horiz = speed_up_horiz/P

fig1,ax1 = plt.subplots(1,1)
ax1.plot(P,speed_up_vert,marker = 'x', label = 'Vertical strips')
ax1.plot(P_square,speed_up_square,marker = 'x', label = 'Squares')
ax1.plot(P,speed_up_horiz,marker = 'x', label = 'Horizontal strips')
ax1.set_xlabel('Number of processors')
ax1.set_ylabel('Speed up')
ax1.legend(loc='upper left')
ax1.set_title('Speed-up')

fig2,ax2 = plt.subplots(1,1)
ax2.plot(P,run_time_vert,marker = 'x', label = 'Vertical strips')
ax2.plot(P_square,run_time_square,marker = 'x', label = 'Squares')
ax2.plot(P,run_time_horiz,marker = 'x', label = 'Horizontal strips')
ax2.set_xlabel('Number of processors')
ax2.set_ylabel('Run time')
ax2.legend(loc = 'upper right')
ax2.set_title('Run time')

fig3,ax3 = plt.subplots(1,1)
ax3.plot(P,parallel_efficiency_vert,marker = 'x', label = 'Vertical strips')
ax3.plot(P_square,parallel_efficiency_square,marker = 'x', label = 'Squares')
ax3.plot(P,parallel_efficiency_horiz,marker = 'x', label = 'Horizontal strips')
ax3.set_xlabel('Number of processors')
ax3.set_ylabel('Parallel efficiency')
ax3.legend(loc = 'upper right')
ax3.set_title('Parallel efficiency')

fig4,ax4 = plt.subplots(1,1)
ax4.plot(P[1:5],io_frac_vert,marker = 'x', label = 'Vertical strips')
ax4.plot(P_square[1:5],io_frac_square,marker = 'x', label = 'Squares')
ax4.plot(P[1:5],io_frac_horiz,marker = 'x', label = 'Horizontal strips')
ax4.set_xlabel('Number of processors')
ax4.set_ylabel('Fraction of time I/O takes')
ax4.legend(loc = 'upper left')
ax4.set_title('Fraction of time used for I/O')


