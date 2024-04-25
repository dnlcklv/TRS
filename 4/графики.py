#1 задание график
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Считываем данные из файла
#1
#data1 = np.loadtxt('G:\Programming\Repos\TRS4\TRS4\z1_0.txt')

#2
#data1 = np.loadtxt('G:\Programming\Repos\TRS4\TRS4\z2_0.txt')

#3
#data1 = np.loadtxt('G:\Programming\Repos\TRS4\TRS4\z2_0.txt')

#4
#data1 = np.loadtxt('G:\Programming\Repos\TRS4\TRS4\z4_0.txt')

#5
data1 = np.loadtxt('G:\Programming\Repos\TRS4\TRS4\z5_0.txt')


size1 = data1.shape

# Получаем размеры сетки
t_dim1 = size1[1]
x_dim1 = size1[0]

# Получаем значения z
u1_values = data1[0:].reshape(x_dim1, t_dim1)
# Создаем сетку x, y
x1 = np.linspace(0, 1, x_dim1)
t1 = np.linspace(0, 10, t_dim1)
X1, T1 = np.meshgrid(t1, x1)

fig = plt.figure(figsize = (14,9))
ax = plt.axes(projection = '3d')

ax.plot_wireframe(T1,X1,u1_values, cmap = 'viridis')

# Отображаем
plt.show()