#%%
import numpy as np
from matplotlib import pyplot as plt

class Hysteresis:
    def __init__(self):
        self.a = 40
        self.alpha = 8e-5
        self.k = 60
        self.c = 0.9
        
    def plot(self):
        Ms = 337147.847 # A/m
        # Ms = 1470000
        H = [0]
        delta = [0]
        Man = [0]
        dMirrdH = [0]
        Mirr = [0]
        M = [0]

        DeltaH = 8000
        Nfirst = 125
        Ndown = 250
        Nup = 250

        for i in range(Nfirst):
            H.append(H[i] + DeltaH)

        for i in range(Ndown):
            H.append(H[-1] - DeltaH)

        for i in range(Nup):
            H.append(H[-1] + DeltaH)

        for i in range(len(H) - 1):
            if H[i + 1] > H[i]:
                delta.append(1)
            else:
                delta.append(-1)
            
        def L(x):
            return 1 / np.tanh(x) - (1 / x)

        for i in range(Nfirst + Ndown + Nup):
            Man.append(Ms * L((H[i + 1] + self.alpha * M[i]) / self.a))
            if H[i + 1] > H[i] and Man[i] > Mirr[i]:
                delta_m = 1
            elif H[i + 1] < H[i] and Man[i] < Mirr[i]:
                delta_m = 1
            else:
                delta_m = 0
            dMirrdH.append(delta_m * (Man[i+1] - Mirr[i]) / (self.k * delta[i+1] - self.alpha * (Man[i + 1] - Mirr[i])))
            Mirr.append(Mirr[i] + dMirrdH[i + 1] * (H[i+1] - H[i]))
            M.append(self.c * Man[i + 1] + (1 - self.c) * Mirr[i + 1])

        plt.plot(H, M, label=self.c)
        r_point = Nfirst + Ndown // 2
        hi_r = (M[r_point - 1] - M[r_point]) / DeltaH
        hi_an =  (Man[r_point - 1] - Man[r_point]) / DeltaH
        for i in range(len(H) - 1):
            if M[i + 1] * M[i] < 0:
                if delta[i] == 1:
                    hi_c = (M[i + 1] - M[i]) / DeltaH
                    hi_can =  (Man[i + 1] - Man[i]) / DeltaH
                    self.k = Man[i] / (1 - self.c) * (self.alpha + (1 - self.c) / (1 * hi_c - self.c * hi_can))
        r_point = Nfirst + Ndown // 2
        hi_td = (M[r_point - 1] - M[r_point]) / DeltaH
        hi_tdan =  (Man[r_point - 1] - Man[r_point]) / DeltaH
        '''
        Remanence point:
        '''
        print(Man[r_point] + self.k / (self.alpha / (1 - self.c) \
                                                   + 1 / (hi_r - self.c * hi_an)))
        return M[Nfirst + Ndown // 2] #, M[-1]

#%%
arr = np.genfromtxt('exp.dat')
# ind = np.where(np.abs(arr[:,1] - arr[:,1][::-1])> 2e+3)[0]
# arr = arr[ind, :]
plt.plot(arr[:,0], arr[:,1])
plt.show()
# plt.plot(arr[:,0][::-1])
for i in range(1, 5):
    m = Hysteresis()
    m.c = 1e-4
    m.a = 4000 * 2 ** i # mag of H
    m.alpha = 1e-3 # neigligiable
    m.k = 60 * 3000 # coercivity
    print(m.plot())
plt.legend()

# for i in range(5):
#     m.alpha *= .5
    
#     print(m.plot())


# plt.plot(arr[:,1] - arr[:,1][::-1])

# %%
