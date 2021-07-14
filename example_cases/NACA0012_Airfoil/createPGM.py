from operator import ge
import numpy as np
from numpy.core.numeric import zeros_like
import pandas as pd
from matplotlib import pyplot as plt


data = pd.read_csv("geom.csv", header=None)

N = 100
theta = 3.62 / 25.75 # define the radians for rotation  

c, s = np.cos(theta), np.sin(theta)
R = np.array(((c, -s), (s, c)))

x_range = data.iloc[:, 0].to_numpy()
y_top = data.iloc[:, 1].to_numpy()
top = np.stack((x_range, y_top))
top = np.matmul(R, top)
y_top = top[1, :]
y_top = np.flip(y_top, axis=0)

x_range = np.flip(x_range, axis=0)
y_bot = data.iloc[:, 3].to_numpy()
bot = np.stack((x_range, y_bot))
bot = np.matmul(R, bot)
x_range = bot[0, :]
y_bot = bot[1, :]


geo = np.zeros((N + 1, N + 1))

dx = 100 / N
dy = 100 / N

def checkPoint(x, y):
    result = 0
    left = np.where(x_range >= x)[0][0] - 1
    interval = x_range[left+1] - x_range[left]
    interp1 = y_top[left] * (x_range[left+1] - x) / interval \
            + y_top[left+1] * (x - x_range[left]) / interval

    interp2 = y_bot[left] * (x_range[left+1] - x) / interval \
            + y_bot[left+1] * (x - x_range[left]) / interval

    if interp1 >= y and interp2 <= y:
        result = 10



    return result


for i in range(N):
    for j in range(N):
        x = i * dx
        y = j * dy - 50
        
        geo[i, j] = checkPoint(x, y)


geo = np.swapaxes(geo, 0, 1)

geo_padded = np.zeros((N+1, 2*N + N+1))
geo_padded[:geo.shape[0], N : N + geo.shape[1]] = geo
geo_padded[:, 0] = 2
geo_padded[-1, :] = 2
geo_padded[0, :] = 1
geo_padded[:, -1] = 1


geo_frame = pd.DataFrame(geo_padded).astype('int32')
geo_frame.to_csv("data.csv", sep=" ", header=False, index=False)

print("Size", geo_padded.shape[1], " ", geo_padded.shape[0])

plt.imshow(geo_padded)
plt.show()