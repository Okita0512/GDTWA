import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots() # 创建图实例

fs_to_au = 41.341

# ==============================================================================================
#                                         data reading     
# ==============================================================================================

fig = plt.figure(figsize=(8, 4), dpi = 512)

data = np.loadtxt("MCTDH.txt", dtype = float)
data = data[np.lexsort(data[:,::-1].T)]
plt.plot(data[:,0], data[:,1], "k-", color = 'black', label = "MCTDH")

# data2 = np.loadtxt("DM_Q_focused.txt", dtype = float)
# plt.plot(data2[:,0] / fs_to_au, data2[:,7], "k--", color = 'green', label = "SPINLSC_Q_focused")

# data3 = np.loadtxt("DM_P_focused.txt", dtype = float)
# plt.plot(data3[:,0] / fs_to_au, data3[:,7], "k-", color = 'purple', label = "SPINLSC_P_focused")

data4 = np.loadtxt("Ehrenfest.txt", dtype = float)
plt.plot(data4[:,0], data4[:,1], "k--", color = 'green', label = "Ehrenfest")

data4 = np.loadtxt("DM_W_focused.txt", dtype = float)
plt.plot(data4[:,0] / fs_to_au, data4[:,7], "k--", color = 'red', label = "SPINLSC_W_focused")

# data5 = np.loadtxt("DM_W_sampled.txt", dtype = float)
# plt.plot(data5[:,0] / fs_to_au, data5[:,7], "k-", color = 'red', label = "SPINLSC_W_sampled")

data4 = np.loadtxt("DM_GDTWA.txt", dtype = float)
plt.plot(data4[:,0] / fs_to_au, data4[:,7], "k-", color = 'blue', label = "GDTWA")



# ==============================================================================================
#                                      plotting set up     
# ==============================================================================================

# x and y range of plotting 
time = 500             # x-axis range: (0, time)
y1, y2 = 0.0, 1.0     # y-axis range: (y1, y2)

# ==============================================================================================

from matplotlib.pyplot import MultipleLocator, tick_params
font = {'family':'Times New Roman', 'weight': 'roman', 'size':12}

# scale for major and minor locator
x_major_locator = MultipleLocator(100)
x_minor_locator = MultipleLocator(20)
y_major_locator = MultipleLocator(0.2)
y_minor_locator = MultipleLocator(0.04)

# x-axis and LHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 8, labelsize = 10)
ax.tick_params(which = 'minor', length = 4)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(labelsize = 20, which = 'both', direction = 'in')
plt.xlim(0.0, time)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 8)
ax2.tick_params(which = 'minor', length = 4)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

# name of x, y axis and the panel
ax.set_xlabel('time / fs', font = 'Times New Roman', size = 18)
ax.set_ylabel('P2', font = 'Times New Roman', size = 18)
ax.set_title('Pyrazine', font = 'Times New Roman', size = 24)

# legend location, font & markersize
ax.legend(loc = 'upper right', prop = font, markerscale = 1)
plt.legend()

plt.savefig("figure_1.png", bbox_inches='tight')

#plt.show()