from numpy import *
from matplotlib.pyplot import *
import mrphantom as pht

# 2D
nPix = 256
arrPhant = pht.genPhant(nAx=2, nPix=nPix)
mapPh=[pht.genPhMap(2,nPix,std=pi/4) for _ in range(pht.Tissue.NTissue.value)]
arrM0 = pht.Enum2SS(arrPhant, mapPh=mapPh)

figure(figsize=(6,3), dpi=150)
subplot(121)
imshow(abs(arrM0), cmap="gray"); colorbar()
subplot(122)
imshow(angle(arrM0), cmap="hsv", vmin=-pi, vmax=pi); colorbar()

# 3D
nPix = 256
arrPhant = pht.genPhant(nAx=3, nPix=nPix)
arrM0 = pht.Enum2SS(arrPhant)

arrM0Abs = abs(arrM0)
figure(figsize=(9,3), dpi=150)
subplot(131)
imshow(arrM0Abs[nPix//2,:,:], cmap="gray"); colorbar()
subplot(132)
imshow(arrM0Abs[:,nPix//2,:], cmap="gray"); colorbar()
subplot(133)
imshow(arrM0Abs[:,:,nPix//2], cmap="gray"); colorbar()

show()

# --- NEW FIGURE: animated 3D slice viewer (does not modify the static figure) ---
V = arrM0Abs  # already computed above (shape: [nPix,nPix,nPix])

vmin = V.min()
vmax = V.max()

figA = figure(figsize=(9, 3), dpi=150)
ax1 = figA.add_subplot(131)
ax2 = figA.add_subplot(132)
ax3 = figA.add_subplot(133)

k0 = nPix // 2
im1 = ax1.imshow(V[k0, :, :], cmap="gray", vmin=vmin, vmax=vmax, origin="lower")
im2 = ax2.imshow(V[:, k0, :], cmap="gray", vmin=vmin, vmax=vmax, origin="lower")
im3 = ax3.imshow(V[:, :, k0], cmap="gray", vmin=vmin, vmax=vmax, origin="lower")

ax1.set_title(f"axial z={k0}")
ax2.set_title(f"coronal y={k0}")
ax3.set_title(f"sagittal x={k0}")

tight_layout()
draw()
pause(0.25)

pause_dt = 1e-6  # seconds; increase if too fast
for k in list(range(nPix)) + list(range(nPix-2, 0, -1)):
    im1.set_data(V[k, :, :]); ax1.set_title(f"axial z={k}")
    im2.set_data(V[:, k, :]); ax2.set_title(f"coronal y={k}")
    im3.set_data(V[:, :, k]); ax3.set_title(f"sagittal x={k}")
    figA.canvas.draw_idle()
    pause(pause_dt)
# --- end block ---

show()