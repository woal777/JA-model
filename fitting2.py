#%%
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate

class Hysteresis:
    def __init__(self):
        Nfirst = 125
        Ndown = 250
        Nup = 250
        self.DeltaH = 10000

        self.Ms = 337147.847 # A/m
        self.a = 470
        self.alpha = 9.38e-4
        self.k = 483
        self.c = 0.0889
        self.H = [0]
        for i in range(Nfirst):
            self.H.append(self.H[i] + self.DeltaH)

        for i in range(Ndown):
            self.H.append(self.H[-1] - self.DeltaH)

        for i in range(Nup):
            self.H.append(self.H[-1] + self.DeltaH)
        
    def plot(self):

        Nfirst = 125
        Ndown = 250
        Nup = 250

        delta = [0]
        Man = [0]
        dMirrdH = [0]
        Mirr = [0]
        M = [0]

        for i in range(len(self.H) - 1):
            if self.H[i + 1] > self.H[i]:
                delta.append(1)
            else:
                delta.append(-1)
            
        def L(x):
            return 1 / np.tanh(x) - (1 / x)

        for i in range(Nfirst + Ndown + Nup):
            Man.append(self.Ms * L((self.H[i + 1] + self.alpha * M[i]) / self.a))
            if self.H[i + 1] > self.H[i] and Man[i] > Mirr[i]:
                delta_m = 1
            elif self.H[i + 1] < self.H[i] and Man[i] < Mirr[i]:
                delta_m = 1
            else:
                delta_m = 0
            dMirrdH.append(delta_m * (Man[i+1] - Mirr[i]) / (self.k * delta[i+1] - self.alpha * (Man[i + 1] - Mirr[i])))
            Mirr.append(Mirr[i] + dMirrdH[i + 1] * (self.H[i+1] - self.H[i]))
            M.append(self.c * Man[i + 1] + (1 - self.c) * Mirr[i + 1])

        return np.array(M[Nfirst:])

arr = np.genfromtxt('exp.dat')
# ind = np.where(np.abs(arr[:,1] - arr[:,1][::-1])> 2e+3)[0]
# arr = arr[ind, :]

# plt.plot(arr[:,0], arr[:,1])
m = Hysteresis()
m.plot()
#%%
ynew = []
f = interpolate.interp1d(arr[:len(arr)//2 + 1,0], arr[:len(arr)//2 + 1,1])
ynew.extend(f(m.H[125:125 + 250]))
f = interpolate.interp1d(arr[len(arr)//2:,0], arr[len(arr)//2:,1])
ynew.extend(f(m.H[125 + 250:]))
# plt.plot(m.H[125:], ynew)
np.savetxt('inter_exp.dat', np.array([m.H[125:], ynew]).T)

m = Hysteresis()
#%%
msd = []
xarr = np.linspace(0.9, 1.1, 20) # b
yarr = np.linspace(0.9, 1.1, 20) # c
# plt.plot(m.H[125:], ynew, color='r')

m.Ms = 337147.847 * 1.3

for b in xarr:
    for c in yarr:
        m.c = 0.0414261 * 1
        m.a = 4000 * 132.6 * c # slope
        m.alpha = 0.09131052 * 34 * b   # slope
        m.k = 157993.396 # coercivity
        y = m.plot()
        # plt.plot(m.H[125:], y, '--')
        # plt.plot(m.H[125:], (m.plot() - ynew) ** 2, label=a)
        msd.append(np.sqrt(np.sum((ynew - y) ** 2) / len(ynew)))

msd = np.array(msd)
msd = np.where(msd < 30000, msd, 0)
fig = plt.figure()
xi, yi = np.meshgrid(xarr, yarr)
cs = plt.contourf(xi, yi, np.reshape(msd, (20, -1)), levels=50)
cbar = fig.colorbar(cs)
plt.scatter(1, 1, color='r')

# %%
m.Ms = 337147.847 * 1.3
m.c = 0.0414261 *  1
m.a = 4000 * 132.6 # mag of H
m.alpha = 0.09131052 * 34 # neigligiable
m.k = 157993.396 # coercivity
plt.plot(m.H[125:], ynew)
plt.plot(m.H[125:], m.plot())

#%%
# xarr[np.argmin(msd)]
# plt.legend()
# plt.show()

# plt.plot(xarr, msd)
# xarr[np.argmin(msd)]
# plt.legend()
# plt.show()

# plt.plot(arr[:,1] - arr[:,1][::-1])


# %%
