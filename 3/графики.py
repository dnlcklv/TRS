import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Считываем данные из файла
data1 = np.loadtxt('G:/Programming/Repos/TRS3/TRS3/task1.txt')
data2 = np.loadtxt('G:/Programming/Repos/TRS3/TRS3/task2.txt')
data3= np.loadtxt('G:/Programming/Repos/TRS3/TRS3/task3.txt')
data4 = np.loadtxt('G:/Programming/Repos/TRS3/TRS3/task4.txt')
data5 = np.loadtxt('G:/Programming/Repos/TRS3/TRS3/analitic.txt')


#size = data1.shape
#size = data2.shape
#size = data3.shape
# size = data4.shape
size = data5.shape

# Получаем размеры сетки
y_dim = size[1]
x_dim = size[0]

# Получаем значения z
#u1_values = data1[0:].reshape(x_dim, y_dim)
#u2_values = data2[0:].reshape(x_dim, y_dim)
#u3_values = data3[0:].reshape(x_dim, y_dim)
#u4_values = data3[0:].reshape(x_dim, y_dim)
u5_values = data5[0:].reshape(x_dim, y_dim)

# Создаем сетку x, y
x = np.linspace(0, 1, x_dim)
y = np.linspace(0, 1, y_dim)
X, Y = np.meshgrid(x, y)

fig = plt.figure(figsize = (14,9))
ax = plt.axes(projection = '3d')

ax.plot_wireframe(X,Y,u5_values, cmap = 'viridis')

# Отображаем
plt.show()