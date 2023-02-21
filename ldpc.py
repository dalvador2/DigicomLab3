import numpy as np
class ldpc:
    def __init__(self,bits,pbits,links,sead=10) -> None:
        np.random.seed = sead
        self.bits = bits
        self.pbits = pbits
        self.code_ray = np.zeros(shape=(bits,bits+pbits),dtype=int)
        for i,row in enumerate(self.code_ray):
            taps = np.random.choice(self.bits+self.pbits,links)
            self.code_ray[i,taps] = np.ones(links,dtype=int)

        

    def encode(self):
        pass

    def decode(self):
        pass

x = ldpc(20,20,8)
print(x.code_ray)
print("done")