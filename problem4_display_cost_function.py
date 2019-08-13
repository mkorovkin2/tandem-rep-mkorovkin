import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def func(x, y, z):
    return np.log(x * y) / z + z

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

lbound = 0.05
nbum = 0.05
ubound = 1.

steps = np.int32((ubound - lbound) / nbum)
# z_static_list = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]

z_static = 0.7 # z_static_list[0]

for i in range(steps):
    x = np.arange(lbound, ubound, nbum)
    y = np.array([i / steps + 1 / steps for x in range(len(x))])
    z = np.array([z_static for x in range(len(x))])
    cs = [] # np.arange(lbound, ubound, nbum)
    for index in range(len(x)):
        cs.append(func(x[index], y[index], z[index]))

    cs = (np.tanh(np.array(cs)) + 1) / 2
    img = ax.scatter(x, y, z, c=cs, cmap=plt.hot())

ax.set_xlabel('Precision')
ax.set_ylabel('Recall')
ax.set_zlabel('Specificity')
plt.title("[Precision, Recall] against Specificity={}\nScore indicated by lightness of color (lighter color -> higher score)".format(z_static))
plt.show()

fig.savefig("specificity07.png")