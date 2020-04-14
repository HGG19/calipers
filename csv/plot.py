import matplotlib.pyplot as plt
import csv


x_c = []
y_c = []
with open('cones.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_c.append(float(row[0]))
        y_c.append(float(row[1]))

x_p = []
y_p = []
with open('points.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_p.append(float(row[0]))
        y_p.append(float(row[1]))

x_j = []
y_j = []
with open('jarvis.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_j.append(float(row[0]))
        y_j.append(float(row[1]))

x_r = []
y_r = []
with open('rotating_calipers.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_r.append(float(row[0]))
        y_r.append(float(row[1]))

#plt.plot(x_c, y_c, 'ro', color='g')
plt.plot(x_p, y_p, 'ro', color='r')
#plt.plot(x_j, y_j, '-ro', color='y')
plt.plot(x_r, y_r, '-ro', color='b')
plt.show()
