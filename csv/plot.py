import matplotlib.pyplot as plt
import csv

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

plt.plot(x_j, y_j, '-ro', color='y')
plt.plot(x_r, y_r, 'ro', color='b')


plt.show()
