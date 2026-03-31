from numpy import *
from matplotlib.pyplot import *
import mrphantom as pht
from time import time

tScan = 10
tRes = 10e-3
nT = int(tScan/tRes)
nPix = 256

cycRes = pi/2
cycCar = 1
arrAmpRes = pht.genResAmp(tScan, tRes, cycRes)
arrAmpCar = pht.genCarAmp(tScan, tRes, cycCar)

fig = figure()
ax = fig.add_subplot(211)
ax.plot(arrAmpRes, ".-")
ax.set_title("Respiratory")
ax = fig.add_subplot(212)
ax.plot(arrAmpCar, ".-")
ax.set_title("Cardiac")

fig = figure(figsize=(6,6), dpi=120)
ax = fig.add_subplot(111)
axim = ax.imshow(zeros([nPix,nPix]), cmap='gray', vmin=0, vmax=1)

while 1:
    for iT in range(nT):
        t = time()
        arrPhant = pht.genPhant(2, nPix, arrAmpRes[iT], arrAmpCar[iT])
        arrM0 = pht.Enum2M0(arrPhant, arrAmpCar[iT])
        if iT%10==0: print(f"{(time()-t)*1e3:.3f} ms / frame")
        
        axim.set_data(arrM0)
        ax.set_title(f"time: {iT*tRes:.2f}s")
        draw()
        pause(tRes/10)
        # pause(1e-2)