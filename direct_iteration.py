#%%
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
import math

class Hysteresis:
    def __init__(self):
        Nfirst = 125
        Ndown = 250
        Nup = 250
        self.DeltaH = 10000

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
        Ms = 337147.847 # A/m

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
            Man.append(Ms * L((self.H[i + 1] + self.alpha * M[i]) / self.a))
            if self.H[i + 1] > self.H[i] and Man[i] > Mirr[i]:
                delta_m = 1
            elif self.H[i + 1] < self.H[i] and Man[i] < Mirr[i]:
                delta_m = 1
            else:
                delta_m = 0
            dMirrdH.append(delta_m * (Man[i+1] - Mirr[i]) / (self.k * delta[i+1] - self.alpha * (Man[i + 1] - Mirr[i])))
            Mirr.append(Mirr[i] + dMirrdH[i + 1] * (self.H[i+1] - self.H[i]))
            M.append(self.c * Man[i + 1] + (1 - self.c) * Mirr[i + 1])

        r_point = Nfirst + Ndown // 2
        hi_r = (M[r_point - 1] - M[r_point]) / self.DeltaH
        hi_an =  (Man[r_point - 1] - Man[r_point]) / self.DeltaH
        for i in range(len(self.H) - 1):
            if M[i + 1] * M[i] < 0:
                if delta[i] == 1:
                    hi_c = (M[i + 1] - M[i]) / self.DeltaH
                    hi_can =  (Man[i + 1] - Man[i]) / self.DeltaH
                    self.k = Man[i] / (1 - self.c) * (self.alpha + (1 - self.c) / (1 * hi_c - self.c * hi_can))
        r_point = Nfirst + Ndown // 2
        hi_td = (M[r_point - 1] - M[r_point]) / self.DeltaH
        hi_tdan =  (Man[r_point - 1] - Man[r_point]) / self.DeltaH
        '''
        Remanence point:
        '''
        # print(Man[r_point] + self.k / (self.alpha / (1 - self.c) \
                                                #    + 1 / (hi_r - self.c * hi_an)))
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
xarr = np.linspace(0.3, 1.5, 50)
for c in xarr:
    m.c = 0.0414261 * c
    m.a = 4000 * 16 * 1.01724 * 1.11 # mag of H
    m.alpha = 0.09131052 * .806  # neigligiable
    m.k = 68993.396 # coercivity
    y = m.plot()
    # plt.plot(m.H[125:], y)
    # plt.plot(m.H[125:], (m.plot() - ynew) ** 2, label=a)
    msd.append(np.sqrt(np.sum((ynew - y) ** 2) / len(ynew)))

plt.plot(xarr, msd)
xarr[np.argmin(msd)]
# plt.legend()
# plt.show()

# plt.plot(arr[:,1] - arr[:,1][::-1])


# %%
m.c = 0.0414261 
m.a = 4000 * 16 * 1.01724 # mag of H
m.alpha = 0.09131052  # neigligiable
m.k = 168993.396 # coercivity
plt.plot(m.H[125:], ynew)
plt.plot(m.H[125:], m.plot())
# %%
