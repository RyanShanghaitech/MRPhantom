from numpy import *
from matplotlib.pyplot import *
from matplotlib.colors import ListedColormap
import mrphantom as pht

nAx = 3
nPix = 128
B0 = 5

random.seed(0)
mapPh = pht.genPhMap((nPix,)*nAx, std=pi/16)
mapB0 = pht.genB0Map((nPix,)*nAx, std=1) # unit: ppm

arrPhant = pht.genPhant((nPix,)*nAx)
mapPD = pht.Enum2PD(arrPhant, B0)
mapT1 = pht.Enum2T1(arrPhant, B0)
mapT2 = pht.Enum2T2(arrPhant, B0)
mapOm = pht.Enum2Om(arrPhant, B0)
mapC = pht.genCsm((nPix,)*nAx, mean=0, std=pi/16)

# plot
cmT1 = ListedColormap(loadtxt("./Resource/lipari.csv"), name="T1")
cmT2 = ListedColormap(loadtxt("./Resource/navia.csv"), name="T2")

# phase map, B0 map
figure(figsize=(9,5), dpi=120)
subplot(121)
imshow(angle(mapPh[:,nPix//2,:]), cmap="hsv", vmin=-pi, vmax=pi)
colorbar()
title("phase map")
subplot(122)
imshow(mapB0[:,nPix//2,:], vmin=-3, vmax=3)
colorbar().set_label("ppm")
title("B0 map")

# PD, T1, T2, B0 map
figure(figsize=(9,9), dpi=120)
subplot(221)
imshow(mapPD[:,nPix//2,:], cmap="gray"); colorbar(); title("PD map")
subplot(222)
imshow(mapT1[:,nPix//2,:]*1000, cmap=cmT1); colorbar().set_label("ms"); title("T1 map")
subplot(223)
imshow(mapT2[:,nPix//2,:]*1000, cmap=cmT2); colorbar().set_label("ms"); title("T2 map")
subplot(224)
imshow(mapOm[:,nPix//2,:]); colorbar(); title("B0 map")
tight_layout()

# coil sensitivity map
mapCsmAbs = abs(mapC)
mapCsmAng = angle(mapC)
fig = figure(figsize=(9,9), dpi=120)
gs = fig.add_gridspec(2,1)

subfig = fig.add_subfigure(gs[0])
subfig.suptitle("Coil Sensitivity Maps, Magnitude")
for iFig in range(3*4):
    ax = subfig.add_subplot(3,4,iFig+1)
    axim = ax.imshow(mapCsmAbs[iFig,nPix-2,:,:], cmap="gray", vmin=0, vmax=1)
    subfig.colorbar(axim)

subfig = fig.add_subfigure(gs[1])
subfig.suptitle("Coil Sensitivity Maps, Phase")
for iFig in range(3*4):
    ax = subfig.add_subplot(3,4,iFig+1)
    axim = ax.imshow(mapCsmAng[iFig,nPix-2,:,:], cmap="hsv", vmin=-pi, vmax=pi)
    subfig.colorbar(axim)

show()