import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots()

fs_to_au = 41.341

# x and y range of plotting 
time = 300             # x-axis range: (0, time)
y1, y2 = 0.0, 1.0     # y-axis range: (y1, y2)

# ==============================================================================================
#                                      Fig 1: P1     
# ==============================================================================================

plt.subplot(3, 1, 1)        # 3 lines + 1 column, here plot the first line

data = np.loadtxt("MCTDH_P1.txt", dtype = float)
plt.plot(data[:,0], data[:,1], "k-", color = 'black', label = "MCTDH")

data2 = np.loadtxt("DM_Q_focused.txt", dtype = float)
plt.plot(data2[:,0] / fs_to_au, data2[:,1], "k--", color = 'green', label = "Ehrenfest")

# data3 = np.loadtxt("DM_P_focused.txt", dtype = float)
# plt.plot(data3[:,0] / fs_to_au, data3[:,1], "k-", color = 'purple', label = "SPINLSC_P_focused")

data4 = np.loadtxt("DM_W_focused.txt", dtype = float)
plt.plot(data4[:,0] / fs_to_au, data4[:,1], "k--", color = 'red', label = "SPINLSC_W_focused")

# data14 = np.loadtxt("1_0.25.txt", dtype = float)
# plt.plot(data14[:,0] / fs_to_au, data14[:,1], "k--", color = 'green', label = "SPINLSC_W_focused_1+0.25")
# 
# data14 = np.loadtxt("1_0.txt", dtype = float)
# plt.plot(data14[:,0] / fs_to_au, data14[:,1], "k--", color = 'purple', label = "SPINLSC_W_focused_1+0")

# data4 = np.loadtxt("Cnn_1_benzene_W_fcs.txt", dtype = float)
# plt.plot(data4[:,0] / fs_to_au, data4[:,1], "k--", color = 'orange', label = "SPINLSC_W_focused_Duncan")


# data5 = np.loadtxt("DM_W_sampled.txt", dtype = float)
# plt.plot(data5[:,0] / fs_to_au, data5[:,1] / 0.6686, "k-", color = 'red', label = "SPINLSC_W_sampled")

# there're some issues with the sampling method, probably a constant that contains Nstates is missed. Here the data is rescaled

data4 = np.loadtxt("DM_GDTWA.txt", dtype = float)
plt.plot(data4[:,0] / fs_to_au, data4[:,1], "k-", color = 'blue', label = "GDTWA")

# plt.plot(data5[:,0]  / fs_to_au, ( data2[:,1] ) * 0.75 + ( data3[:,1] ) * 0.25, "k-", color = 'yellow', label = "Mixed_Q3:P1")

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
# ax.set_xlabel('time / fs', font = 'Times New Roman', size = 18)
ax.set_ylabel('P1', font = 'Times New Roman', size = 18)
ax.set_title('Benzene', font = 'Times New Roman', size = 24)

# legend location, font & markersize
ax.legend(loc = 'upper left', prop = font, markerscale = 1)
plt.legend()

# ==============================================================================================
#                                      Fig 2: P2     
# ==============================================================================================

plt.subplot(3, 1, 2)

data6 = np.loadtxt("MCTDH_P2.txt", dtype = float)
plt.plot(data6[:,0], data6[:,1], "k-", color = 'black', label = "MCTDH")

data7 = np.loadtxt("DM_Q_focused.txt", dtype = float)
plt.plot(data7[:,0] / fs_to_au, data7[:,9], "k--", color = 'green', label = "Ehrenfest")

# data8 = np.loadtxt("DM_P_focused.txt", dtype = float)
# plt.plot(data8[:,0] / fs_to_au, data8[:,9], "k-", color = 'purple', label = "SPINLSC_P_focused")

data9 = np.loadtxt("DM_W_focused.txt", dtype = float)
plt.plot(data9[:,0] / fs_to_au, data9[:,9], "k--", color = 'red', label = "SPINLSC_W_focused")

# data14 = np.loadtxt("1_0.25.txt", dtype = float)
# plt.plot(data14[:,0] / fs_to_au, data14[:,9], "k--", color = 'green', label = "SPINLSC_W_focused_1+0.25")
# 
# data14 = np.loadtxt("1_0.txt", dtype = float)
# plt.plot(data14[:,0] / fs_to_au, data14[:,9], "k--", color = 'purple', label = "SPINLSC_W_focused_1+0")

# data4 = np.loadtxt("Cnn_1_benzene_W_fcs.txt", dtype = float)
# plt.plot(data4[:,0] / fs_to_au, data4[:,2], "k--", color = 'orange', label = "SPINLSC_W_focused_Duncan")

#data10 = np.loadtxt("DM_W_sampled.txt", dtype = float)
#plt.plot(data10[:,0] / fs_to_au, data10[:,9] / 0.6686, "k-", color = 'red', label = "SPINLSC_W_sampled")

data4 = np.loadtxt("DM_GDTWA.txt", dtype = float)
plt.plot(data4[:,0] / fs_to_au, data4[:,9], "k-", color = 'blue', label = "GDTWA")

# plt.plot(data7[:,0]  / fs_to_au, ( data7[:,9] ) * 0.75 + ( data8[:,9] ) * 0.25, "k-", color = 'yellow', label = "Mixed_Q3:P1")

# x and y range of plotting 
# time = 300             # x-axis range: (0, time)
# y1, y2 = 0.0, 1.0     # y-axis range: (y1, y2)

# ==============================================================================================

# from matplotlib.pyplot import MultipleLocator, tick_params
# font = {'family':'Times New Roman', 'weight': 'roman', 'size':12}

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
# ax.set_xlabel('time / fs', font = 'Times New Roman', size = 18)
ax.set_ylabel('P2', font = 'Times New Roman', size = 18)
# ax.set_title('Benzene', font = 'Times New Roman', size = 24)

# legend location, font & markersize
# ax.legend(loc = 'upper right', prop = font, markerscale = 1)
# plt.legend()

# ==============================================================================================
#                                      Fig 3: P3     
# ==============================================================================================

plt.subplot(3, 1, 3)

data11 = np.loadtxt("MCTDH_P3.txt", dtype = float)
plt.plot(data11[:,0], data11[:,1], "k-", color = 'black', label = "MCTDH")

data12 = np.loadtxt("DM_Q_focused.txt", dtype = float)
plt.plot(data12[:,0] / fs_to_au, data12[:,17], "k--", color = 'green', label = "Ehrenfest")

# data13 = np.loadtxt("DM_P_focused.txt", dtype = float)
# plt.plot(data13[:,0] / fs_to_au, data13[:,17], "k-", color = 'purple', label = "SPINLSC_P_focused")
# data4 = np.loadtxt("Cnn_1_benzene_W_fcs.txt", dtype = float)
# plt.plot(data4[:,0] / fs_to_au, data4[:,3], "k--", color = 'orange', label = "SPINLSC_W_focused_Duncan")

data14 = np.loadtxt("DM_W_focused.txt", dtype = float)
plt.plot(data14[:,0] / fs_to_au, data14[:,17], "k--", color = 'red', label = "SPINLSC_W_focused") # _0.25+1

# data14 = np.loadtxt("1_0.25.txt", dtype = float)
# plt.plot(data14[:,0] / fs_to_au, data14[:,17], "k--", color = 'green', label = "SPINLSC_W_focused_1+0.25")
# 
# data14 = np.loadtxt("1_0.txt", dtype = float)
# plt.plot(data14[:,0] / fs_to_au, data14[:,17], "k--", color = 'purple', label = "SPINLSC_W_focused_1+0")

# data15 = np.loadtxt("DM_W_sampled.txt", dtype = float)
# plt.plot(data15[:,0] / fs_to_au, data15[:,17] / 0.6686, "k-", color = 'red', label = "SPINLSC_W_sampled")

data4 = np.loadtxt("DM_GDTWA.txt", dtype = float)
plt.plot(data4[:,0] / fs_to_au, data4[:,17], "k-", color = 'blue', label = "GDTWA")

# plt.plot(data12[:,0]  / fs_to_au, ( data12[:,17] ) * 0.75 + ( data13[:,17] ) * 0.25, "k-", color = 'yellow', label = "Mixed_Q3:P1")

# ==============================================================================================

# from matplotlib.pyplot import MultipleLocator, tick_params
# font = {'family':'Times New Roman', 'weight': 'roman', 'size':12}

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
ax.set_ylabel('P3', font = 'Times New Roman', size = 18)
# ax.set_title('Benzene', font = 'Times New Roman', size = 24)

# legend location, font & markersize
# ax.legend(loc = 'upper right', prop = font, markerscale = 1)
# plt.legend()

plt.show()

