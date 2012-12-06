#! /usr/bin/python

import os, sys
from math import *
import numpy as np
import pylab
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d import Axes3D

##### GRAB DATA ######################

results_mio = file("results/mio/runtimes.txt", 'r')
results_sayers = file("results/sayers/runtimes.txt", 'r')

mio = []
sayers = []

while 1:
	temp_mio = results_mio.readline().split()
	temp_sayers = results_sayers.readline().split()
	if len(temp_mio) == 0 and len(temp_sayers) == 0:
		break

	if len(temp_mio) != 0:
		mio.append(float(temp_mio[2]))
		mio.append(float(temp_mio[3]))
	if len(temp_sayers) != 0:
		sayers.append(float(temp_sayers[2]))
		sayers.append(float(temp_sayers[3]))

results_mio.close()
results_sayers.close()

mio = np.array(mio).reshape(len(mio)/2, 2)
mio_i = np.lexsort((mio[:,1], mio[:,0]))
sayers = np.array(sayers).reshape(len(sayers)/2, 2)
sayers_i = np.lexsort((sayers[:,1], sayers[:,0]))

mio = mio[mio_i]
sayers = sayers[sayers_i]

##### RUN TIMES #######################

temp_mio = np.log2(mio)
temp_sayers = np.log2(sayers)

# plt.plot(temp_mio[:,0], temp_mio[:,1], "r*-", label="Mio")
# plt.plot(temp_sayers[:,0], temp_sayers[:,1], "bo-", label="Sayers")

# plt.title("Number of Cores Vs. Runtime for n = 50,000 and m = 250,000")
# plt.xlabel("lg(Number of Cores)")
# plt.ylabel("lg(Time in seconds)")

# plt.legend(loc=1)

# plt.show()


### Table ###
runtimes = file("results/tables/runtimes_table.txt", 'w')
runtimes.write("\\begin{tabular}{rrr}\n")
runtimes.write("Cores".rjust(7))
runtimes.write("&Sayers Lab".rjust(25))
runtimes.write("&Mio".rjust(25))
runtimes.write(" \\\\ \n\\hline\n")
for i in range(len(mio[:,0])):
	runtimes.write(repr(mio[i,0]).rjust(7))
	if i >= len(sayers[:,1]):
		runtimes.write("&N/A".rjust(25))
	else:
		temp = "&" + repr(sayers[i,1])
		runtimes.write(temp.rjust(25))
	temp = "&" + repr(mio[i,1])
	runtimes.write(temp.rjust(25))
	runtimes.write(" \\\\ \n")
runtimes.write("\\end{tabular}")

runtimes.close()


##### SPEEDUP #########################

speedup_mio = mio[0,1]/mio[:,1]
speedup_sayers = sayers[0,1]/sayers[:,1]

# temp_mio = np.log2(speedup_mio)
# temp_sayers = np.log2(speedup_sayers)

# plt.plot(np.log2(mio[:,0]), temp_mio, "r*-", label="Mio")
# plt.plot(np.log2(sayers[:,0]), temp_sayers, "bo-", label="Sayers")

# plt.title("Number of Cores Vs. Speedup for n = 50,000 and m = 250,000")
# plt.xlabel("lg(Number of Cores)")
# plt.ylabel("lg(Speedup)")

# plt.legend(loc=2)
# plt.show()


### Table ###
runtimes = file("results/tables/speedup_table.txt", 'w')
runtimes.write("\\begin{tabular}{rrr}\n")
runtimes.write("Cores".rjust(7))
runtimes.write("&Sayers Lab".rjust(25))
runtimes.write("&Mio".rjust(25))
runtimes.write(" \\\\ \n\\hline\n")
for i in range(len(mio[:,0])):
	runtimes.write(repr(mio[i,0]).rjust(7))
	if i >= len(sayers[:,1]):
		runtimes.write("&N/A".rjust(25))
	else:
		temp = "&" + repr(speedup_sayers[i])
		runtimes.write(temp.rjust(25))
	temp = "&" + repr(speedup_mio[i])
	runtimes.write(temp.rjust(25))
	runtimes.write(" \\\\ \n")
runtimes.write("\\end{tabular}")

runtimes.close()


##### EFFICIENCY ####################

temp_mio = mio[0,1] / (mio[:,1] * mio[:,0])
temp_sayers = sayers[0,1] / (sayers[:,1] * sayers[:,0])

# plt.plot(np.log2(mio[:,0]), temp_mio, "r*-", label="Mio")
# plt.plot(np.log2(sayers[:,0]), temp_sayers, "bo-", label="Sayers")

# plt.title("Number of Cores Vs. Efficiency for n = 50,000 and m = 250,000")
# plt.xlabel("lg(Number of Cores)")
# plt.ylabel("Efficiency")

# plt.legend(loc=1)
# plt.show()


### Table ###
runtimes = file("results/tables/efficiency_table.txt", 'w')
runtimes.write("\\begin{tabular}{rrr}\n")
runtimes.write("Cores".rjust(7))
runtimes.write("&Sayers Lab".rjust(25))
runtimes.write("&Mio".rjust(25))
runtimes.write(" \\\\ \n\\hline\n")
for i in range(len(mio[:,0])):
	runtimes.write(repr(mio[i,0]).rjust(7))
	if i >= len(sayers[:,1]):
		runtimes.write("&N/A".rjust(25))
	else:
		temp = "&" + repr(temp_sayers[i])
		runtimes.write(temp.rjust(25))
	temp = "&" + repr(temp_mio[i])
	runtimes.write(temp.rjust(25))
	runtimes.write(" \\\\ \n")
runtimes.write("\\end{tabular}")

runtimes.close()


##### EXPERIMENTAL SERIAL FRACTION ###

temp_mio = (1/speedup_mio[1:(len(speedup_mio))] - 1/mio[1:(len(mio)),0])/(1 - 1/mio[1:(len(mio)),0])
temp_sayers = (1/speedup_sayers[1:(len(speedup_sayers))] - 1/sayers[1:(len(sayers)),0])/(1 - 1/sayers[1:(len(sayers)),0])

# plt.plot(np.log2(mio[1:len(mio),0]), temp_mio, "r*-", label="Mio")
# plt.plot(np.log2(sayers[1:len(sayers):,0]), temp_sayers, "bo-", label="Sayers")

# plt.title("Number of Cores Vs. Experimental Serial Fraction \nfor n = 50,000 and m = 250,000")
# plt.xlabel("lg(Number of Cores)")
# plt.ylabel("Experimental Serial Fraction")

# plt.legend(loc=1)
# plt.show()


### Table ###
runtimes = file("results/tables/experimental_serial_fraction_table.txt", 'w')
runtimes.write("\\begin{tabular}{rrr}\n")
runtimes.write("Cores".rjust(7))
runtimes.write("&Sayers Lab".rjust(25))
runtimes.write("&Mio".rjust(25))
runtimes.write(" \\\\ \n\\hline\n")
for i in range(len(mio[:,0])-1):
	runtimes.write(repr(mio[i+1,0]).rjust(7))
	if i >= len(temp_sayers):
		runtimes.write("&N/A".rjust(25))
	else:
		temp = "&" + repr(temp_sayers[i])
		runtimes.write(temp.rjust(25))
	temp = "&" + repr(temp_mio[i])
	runtimes.write(temp.rjust(25))
	runtimes.write(" \\\\ \n")
runtimes.write("\\end{tabular}")

runtimes.close()


#### FUNCTION PLOTS #################

# def u(x, t):
# 	return (1+t)*np.exp(-t)*np.sin(x) + np.cos(t)*np.exp(-2*t)*np.sin(2*x)

# # def u(x, t):
# # return (1+t)*np.exp(-t)*np.sin(x) + np.cos(t)*np.exp(-2*t)*np.sin(2*x)

# big_u_values_file = file("results/u_values.txt", 'r')
# t_vect = []
# x_vect = []
# u_vect = []

# while 1:
# 	t = big_u_values_file.readline()
# 	if len(t) == 0:
# 		break

# 	t_vect.append(float(t))

# 	t = big_u_values_file.readline().split()

# 	x_vect2 = []
# 	for x in t:
# 		x_vect2.append(float(x))
# 	x_vect.append(x_vect2)

# 	t = big_u_values_file.readline().split()

# 	u_vect2 = []
# 	for x in t:
# 		u_vect2.append(float(x))
# 	u_vect.append(u_vect2)

# big_u_values_file.close()

# # fig = plt.figure()
# # ax = fig.gca(projection='3d')
# fig = pylab.figure()

# #ax = Axes3D(fig)
# ax = fig.add_subplot(1, 2, 1, projection='3d')
# ax.set_autoscale_on(False)
# ax.set_xlim3d(0,pi)
# ax.set_ylim3d(0,4)
# ax.set_zlim3d(-.15,1.5)

# x_values = np.arange(0, pi+pi/380, pi/380)
# t_values = np.arange(.2, 4, .01)

# x_values, t_values = np.meshgrid(x_values, t_values)
# u_values = u(x_values, t_values)
# ax.plot_surface(x_values, t_values, u_values, linewidth=0.3, cmap=cm.jet)

# ax.set_title("Actual u")
# ax.set_xlabel("x")
# ax.set_ylabel("t")
# ax.set_zlabel("u(x,t)")

# ax = fig.add_subplot(1, 2, 2, projection='3d')
# ax.set_autoscale_on(False)
# ax.set_xlim3d(0,pi)
# ax.set_ylim3d(0,4)
# ax.set_zlim3d(-.15,1.5)

# x_vect2, t_vect = np.meshgrid(x_vect2, t_vect)
# u_vect = np.array(u_vect)

# ax.plot_surface(x_vect2, t_vect, u_vect, linewidth=0.3, cstride=1, cmap=cm.jet)

# ax.set_title("Approximate u with n=20, m=24")
# ax.set_xlabel("x")
# ax.set_ylabel("t")
# ax.set_zlabel("u(x,t)")

# plt.show()


##### GENERATE MOVIE #####################################

# files = []
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# for i in range(380):  # 50 frames
#     ax.cla()
#     ax.plot(x_vect2[i], u_vect[i], 'k-')
#     ax.axis([0, pi, np.amin(u_vect), np.amax(u_vect)])
#     fname = '_tmp%03d.png'%i
#     print 'Saving frame', fname
#     fig.savefig(fname)
#     files.app\\end(fname)

# print 'Making movie animation.mpg - this may take a while'
# os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=19 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")


##### MAKE TABLES ####################

data = np.genfromtxt("results/table5.txt", dtype=str)
table5 = file("results/tables/table_5_table.txt", 'w')
table5.write("\\begin{tabular}{rrrrrrrr}\n")
for i in data:
	table5.write(i[0].rjust(5))
	for j in range(1,len(i)):
		temp = ('&' + i[j]).rjust(14)
		table5.write(temp)

	table5.write(" \\\\ \n")
table5.write("\\end{tabular}")

table5.close()

data = np.genfromtxt("results/table6.txt", dtype=str, skip_header=1)
table6 = file("results/tables/table_6_table.txt", 'w')
table6.write("\\begin{tabular}{rrrrrrr}\n")
for i in data:
	table6.write(i[0].rjust(5))
	for j in range(1,len(i)):
		temp = ('&' + i[j]).rjust(14)
		table6.write(temp)

	table6.write(" \\\\ \n")
table6.write("\\end{tabular}")

table6.close()

data = np.genfromtxt("results/table8.txt", dtype=str)
table8 = file("results/tables/table_8_table.txt", 'w')
table8.write("\\begin{tabular}{rrrrrrrr}\n")
for i in data:
	table8.write(i[0].rjust(5))
	for j in range(1,len(i)):
		temp = ('&' + i[j]).rjust(14)
		table8.write(temp)

	table8.write(" \\\\ \n")
table8.write("\\end{tabular}")

table8.close()