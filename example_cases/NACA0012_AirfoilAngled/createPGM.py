from operator import ge
import numpy as np
from numpy.core.numeric import zeros_like
import pandas as pd
from matplotlib import pyplot as plt
import sys


data = pd.read_csv("geom.csv", header=None)

N = 100
theta = np.radians(12)
print(np.degrees(theta))

c, s = np.cos(theta), np.sin(theta)
R = np.array(((c, -s), (s, c)))

x_range = data.iloc[:, 0].to_numpy()
num = x_range.size
y_top = data.iloc[:, 1].to_numpy()
top = np.vstack((x_range, y_top))

x_range_flip = np.flip(x_range, axis=0)
y_bot = data.iloc[:, 3].to_numpy()
bot = np.vstack((x_range_flip, y_bot))

# combine bot and top then rotate
all = np.hstack((top, bot))
mean = np.mean(all, axis=1)
mean = np.expand_dims(mean, axis=1)
all = np.matmul(R, all - mean) + mean

# get the value back
y_top = all[1, :num]
y_top = np.flip(y_top, axis=0)
x_range = all[0, num:]
y_bot = all[1, num:]


geo = np.zeros((N + 1, N + 1))

dx = 100 / N
dy = 100 / N

def checkPoint(x, y):
    result = 0
    left = np.where(x_range >= x)[0]
    if(left.size == 0):
        return 0
    left = left[0] - 1

    interval = x_range[left+1] - x_range[left]
    interp1 = y_top[left] * (x_range[left+1] - x) / interval \
            + y_top[left+1] * (x - x_range[left]) / interval

    interp2 = y_bot[left] * (x_range[left+1] - x) / interval \
            + y_bot[left+1] * (x - x_range[left]) / interval

    if interp1 >= y and interp2 <= y:
        result = 9



    return result


for i in range(N):
    for j in range(N):
        x = i * dx
        y = j * dy - 50
        
        geo[i, j] = checkPoint(x, y)


# check invalid and modifiy the cell number
def check_validity(geo):
    offset = [1, 0, -1, 0, 1]
    valid = 1
    for i in range(N):
        for j in range(N):
            if(geo[i][j] == 9):
                zeros = []
                for k in range(4):
                    x = i + offset[k]
                    y = j + offset[k+1]
                    if(x < N+1 and x > 0 and y < N+1 and y > 0):
                        if(geo[x][y] == 0):
                            zeros.append(k)
                if(len(zeros)>2):
                    print("invalid", i, j)
                    geo[i][j] = 5
                    valid = 0
                    for k in range(len(zeros)):
                        x = i + offset[zeros[k]]
                        y = j + offset[zeros[k]]
                        count = 0
                        for kk in range(4):
                            xx = x + offset[k]
                            yy = y + offset[k+1]
                            if(xx < N+1 and xx > 0 and yy < N+1 and yy > 0):
                                if(geo[xx][yy] == 0):
                                    count = count + 1
                        if(count <= 2):
                            geo[x][y] = 9
    if(valid):
        print("valid")
        return True
    return False


geo = np.swapaxes(geo, 0, 1)

max_iter = 5
iter = 0
while(iter < max_iter and not check_validity(geo)):
    iter = iter + 1 
    print(iter)

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