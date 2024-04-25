#1 задание график
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Считываем данные из файла
#1
# data1 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z1_1.txt')
# data2 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z1_2.txt')
# data3 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z1_3.txt')
#2.1
# data1 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z2_p11.txt')
# data2 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z2_p12.txt')
# data3 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z2_p13.txt')
#2.2
# data1 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z2_p21.txt')
# data2 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z2_p22.txt')
# data3 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z2_p23.txt')
#3
data1 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z3_u0.txt')
data2 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z3_u1.txt')
data3 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z3_u2.txt')
#4
# data1 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z4_u0.txt')
# data2 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z4_u0.txt')
# data3 = np.loadtxt('G:\Programming\Repos\TRS2\TRS2\z4_u0.txt')

size1 = data1.shape
size2 = data2.shape
size3 = data3.shape

# Получаем размеры сетки
t_dim1 = size1[1]
x_dim1 = size1[0]

t_dim2 = size2[1]
x_dim2 = size2[0]

t_dim3 = size3[1]
x_dim3 = size3[0]


# Получаем значения z
u1_values = data1[0:].reshape(x_dim1, t_dim1)
u2_values = data2[0:].reshape(x_dim2, t_dim2)
u3_values = data3[0:].reshape(x_dim3, t_dim3)
# Создаем сетку x, y
x1 = np.linspace(0, 1, x_dim1)
t1 = np.linspace(0, 1, t_dim1)
X1, T1 = np.meshgrid(t1, x1)

x2 = np.linspace(0, 1, x_dim2)
t2 = np.linspace(0, 1, t_dim2)
X2, T2 = np.meshgrid(t2, x2)

x3 = np.linspace(0, 1, x_dim3)
t3 = np.linspace(0, 1, t_dim3)
X3, T3 = np.meshgrid(t3, x3)

# Строим график поверхности
# fig = plt.figure()
# ax = []
# ax.append(fig.add_subplot(221, projection='3d'))
# ax.append(fig.add_subplot(222,projection="3d"))
# ax.append( fig.add_subplot(223,projection='3d'))
# ax[0].plot_wireframe(T1, X1, u1_values, cmap='viridis')
# ax[1].plot_wireframe(T2, X2, u2_values, cmap='viridis')
# ax[2].plot_wireframe(T3, X3, u3_values, cmap='viridis')

fig = plt.figure(figsize = (14,9))
ax = plt.axes(projection = '3d')

ax.plot_wireframe(T3,X3,u3_values, cmap = 'viridis')

# Отображаем
plt.show()