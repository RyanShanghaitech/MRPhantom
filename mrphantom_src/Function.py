from . import ext
from numpy import *
from numpy.fft import fftn, ifftn, fftshift, ifftshift
from numpy.typing import NDArray
from typing import *
from scipy.signal.windows import gaussian
from scipy.ndimage import gaussian_filter
from . import dictNmrPara, lstTissue

def genPhant(shape:Tuple, ampRes:float=0, ampCar:float=0) -> NDArray: # call C++ backend to generate a phantom
    """
    generate a phantom in Enum type
    
    Args:
        shape: shape of the phantom (independent from FOV)
        ampRes: respiratory motion amplitude
        ampCar: cardiac motion amplitude
        
    Returns:
        NDArray contains elements in `Tissue` enum type
    """
    nAx = len(shape)
    shape = (1,)*(3-len(shape)) + shape
    return ext.genPhant(nAx, *shape, ampRes, ampCar)

def LPF(arr:NDArray, std:float) -> NDArray:
    lstWind = []
    for l in arr.shape:
       lstWind.append(gaussian(l, 1/(2*pi*std)))
    wind = ones(arr.shape)
    for a, w in enumerate(lstWind):
        shape = [1] * arr.ndim
        shape[a] = len(w)
        wind *= w.reshape(shape)

    arr = fftshift(fftn(ifftshift(arr)))
    arr *= wind
    arr = fftshift(ifftn(ifftshift(arr)))
    return arr

def genPhMap(shape:Tuple, mean:int|float|None=None, std:int|float=pi/16) -> NDArray:
    """
    generate random phase map
    
    Args:
        shape: shape of the phantom (independent from FOV)
        mean: mean of the noise
        std: std of the noise
        
    Returns:
        smooth complex noise with unity magnitude
    """
    if mean is None: mean = random.uniform(-pi,pi)
    mapPh = random.uniform(-pi, pi, shape)
    mapPh = LPF(mapPh, 1/4).real
    # normalize
    mapPh -= mapPh.mean(); mapPh = asarray(mapPh)
    mapPh /= mapPh.std()
    mapPh *= std
    mapPh += mean
    # convert to rotation factor
    mapPh = exp(1j*mapPh)
    mapPh = mapPh/abs(mapPh)
    return mapPh

def genB0Map(shape:Tuple, mean:int|float=0, std:int|float=1e-6*(2*pi*42.58e6*3)) -> NDArray:
    """
    generate random B0 map
    
    Args:
        shape: shape of the phantom (independent from FOV)
        mean: mean of the noise
        std: std of the noise
        
    Returns:
        smooth random noise in `rad/s`
    """
    mapB0 = random.uniform(-1, 1, shape)
    mapB0 = LPF(mapB0, 1/4).real
    # normalize
    mapB0 -= mapB0.mean(); mapB0 = asarray(mapB0)
    mapB0 /= mapB0.std()
    mapB0 *= std
    mapB0 += mean
    return mapB0

def genCsm(shape:Tuple, nCh:int=12, mean:int|float|None=None, std:int|float=pi/16) -> NDArray:
    """
    generate random coil sensitivity map
    
    Args:
        shape: shape of the phantom (independent from FOV)
        nCh: number of coils
        mean: mean of the noise
        std: std of the noise
        
    Returns:
        complex smooth and inhomogeneous map
    """
    if mean is None: mean = random.uniform(-pi,pi)
    nAx = len(shape)
    mapC = zeros([nCh,*shape], dtype=complex128)
    arrCoor = meshgrid\
    (
        *(linspace(-0.5,0.5,s,0) for s in shape),
        indexing="ij"
    ); arrCoor = array(arrCoor).transpose(*arange(1,nAx+1), 0)
    arrTht = linspace(0,2*pi,nCh,0)
    arrCoorCoil = zeros([nCh,nAx], dtype=float64)
    arrCoorCoil[:,-2:] = 1*array([sin(arrTht), cos(arrTht)]).T
    if nAx == 3:
        arrCoorCoil[0::2,0] = 0.2
        arrCoorCoil[1::2,0] = -0.2
    for iCh in range(nCh):
        mapC[iCh] = genPhMap(shape, mean=mean, std=std)
        dist = sqrt(sum((arrCoor - arrCoorCoil[iCh])**2, axis=-1))
        mapC[iCh] *= exp(-dist)
    return mapC

def genAmp(tScan:int|float, tRes:int|float, cyc:int|float, isRand:bool=True) -> NDArray:
    """
    generate amplitude curve
    
    Args:
        tScan: length of the signal in `s`
        tRes: temporal resolution in `s`
        cyc: period of the signal in `s`
        isRand: make the signal have irregular period
        
    Returns:
        generated amplitude
    """
    nT = around(tScan/tRes).astype(int)

    if isRand:
        arrT = sort(random.rand(nT)*tScan)
        arrAmp = sin(2*pi/cyc*arrT)

        sigma = cyc/tRes/8
        arrAmp = gaussian_filter(arrAmp, sigma)
    else:
        arrT = linspace(0, tScan, nT)
        arrAmp = sin(2*pi/cyc*arrT)

    return arrAmp

def genResAmp(tScan:int|float, tRes:int|float, cyc:int|float=pi/2) -> NDArray:
    """
    generate respiratory amplitude curve
    
    Args:
        tScan: length of the signal in `s`
        tRes: temporal resolution in `s`
        cyc: period of the signal in `s`
        
    Returns:
        generated amplitude, approx. -0.02~0.02
    """
    return 20e-3*genAmp(tScan, tRes, cyc, 1)

def genCarAmp(tScan:int|float, tRes:int|float, cyc:int|float=1) -> NDArray:
    """
    generate cardiac amplitude curve
    
    Args:
        tScan: length of the signal in `s`
        tRes: temporal resolution in `s`
        cyc: period of the signal in `s`
        
    Returns:
        generated amplitude, approx. -0.01~0.01
    """
    return 10e-3*genAmp(tScan, tRes, cyc, 0)

def fB02strB0(B0:int|float) -> str:
    """
    convert B0 data type from float/int to string.

    Args:
        B0: B0 in number format
    
    Returns:
        B0 in string format
    """
    if isclose(B0,0.55): return "B0_0T55"
    if isclose(B0,1.5): return "B0_1T5"
    if isclose(B0,3.0): return "B0_3T"
    if isclose(B0,5.0): return "B0_5T0"
    if isclose(B0,9.4): return "B0_9T4"
    raise RuntimeError("unsupported B0")

def initSS_bSSFP(
    B0: float,
    TR: float = 5e-3,
    FA_deg: float = 60.0,
) -> None:
    """
    Precalculate and store the bSSFP steady-state signal Mss for every tissue.
    
    Args:
        B0: field strength
        TR: repetition time
        FA_deg: flip angle in degree
    """
    strB0 = fB02strB0(B0)
    FA = deg2rad(FA_deg)

    for strTissue in lstTissue:
        dictTissue = dictNmrPara[strTissue]

        PD = dictTissue[strB0]["PD"]
        T1 = dictTissue[strB0]["T1"]
        T2 = dictTissue[strB0]["T2"]

        E1 = exp(-TR / T1)
        E2 = exp(-TR / T2)

        dictTissue["Mss"] = (
            PD
            * (1 - E1)
            * sqrt(E2)
            * sin(FA)
            / (1 - (E1 - E2) * cos(FA) - E1 * E2)
        )

def initSS_FLASH(
    B0: float,
    TE: float = 1e-3,
    TR: float = 10e-3,
    FA_deg: float = 10.0,
) -> None:
    """
    Precalculate and store the FLASH steady-state signal Mss for every tissue.
    
    Args:
        B0: field strength
        TE: echo time
        TR: repetition time
        FA_deg: flip angle in degree
    """
    strB0 = fB02strB0(B0)
    FA = deg2rad(FA_deg)

    for strTissue in lstTissue:
        dictTissue = dictNmrPara[strTissue]

        PD = dictTissue[strB0]["PD"]
        T1 = dictTissue[strB0]["T1"]
        T2s = dictTissue[strB0]["T2s"]

        E1 = exp(-TR / T1)
        E2 = exp(-TE / T2s)

        dictTissue["Mss"] = (
            PD
            * sin(FA)
            * (1 - E1)
            / (1 - cos(FA) * E1)
            * E2
        )

def Enum2SS(arrPht:NDArray) -> NDArray:
    """
    get steady-state signal map of a phantom generated by `genPhant()`
    
    Args:
        arrPht: phantom
        
    Returns:
        steady-state signal map of the given phantom
    """
    mapSS = zeros_like(arrPht, dtype=float64)
    for strTissue in lstTissue:
        try: mapSS[arrPht==dictNmrPara[strTissue]["enum"]] = dictNmrPara[strTissue]["Mss"]
        except KeyError: raise RuntimeError("Please call `initS_FLASH()` or `initS_bSSFP()` before `Enum2SS()`.")
    return mapSS

def Enum2Para(arrPht:NDArray, B0:str|float="B0_1T5", strPara:str="PD") -> NDArray:
    """
    get PD map of a phantom generated by `genPhant()`
    
    Args:
        arrPht: phantom
        B0: "B0_0T55" / "B0_1T5" / "B0_3T" / "B0_5T" / "B0_9T4" / 0.55 / 1.5 / 3.0 / 5.0 / 9.4
        strPara: "PD" / "T1" / "T2" / "T2s" / "ADC" / "Om"
        
    Returns:
        Proton density map of the given phantom, relevant to water
    """
    mapPara = zeros_like(arrPht, dtype=float64)
    if not isinstance(B0, str): B0 = fB02strB0(B0)
    for strTissue in lstTissue:
        mapPara[arrPht==dictNmrPara[strTissue]["enum"]] = dictNmrPara[strTissue][B0][strPara]
    return mapPara

def Enum2PD(arrPht:NDArray, B0:str|int|float) -> NDArray:
    """
    get PD map of a phantom generated by `genPhant()`
    
    Args:
        arrPht: phantom
        B0: field strength in Tesla
        
    Returns:
        Proton density map of the given phantom, relevant to water
    """
    return Enum2Para(arrPht, B0, "PD")

def Enum2T1(arrPht:NDArray, B0:str|int|float) -> NDArray:
    """
    get T1 map of a phantom generated by `genPhant()`
    
    Args:
        arrPht: phantom
        B0: field strength in Tesla
        
    Returns:
        T1 map of the given phantom
    """
    return Enum2Para(arrPht, B0, "T1")

def Enum2T2(arrPht:NDArray, B0:str|int|float) -> NDArray:
    """
    get T2 map of a phantom generated by `genPhant()`
    
    Args:
        arrPht: phantom
        B0: field strength in Tesla
        
    Returns:
        T2 map of the given phantom
    """
    return Enum2Para(arrPht, B0, "T2")

def Enum2T2s(arrPht:NDArray, B0:str|int|float) -> NDArray:
    """
    get T2* map of a phantom generated by `genPhant()`
    
    Args:
        arrPht: phantom
        B0: field strength in Tesla
        
    Returns:
        T2* map of the given phantom
    """
    return Enum2Para(arrPht, B0, "T2s")

def Enum2Adc(arrPht:NDArray, B0:str|int|float) -> NDArray:
    """
    get Apparent Diffusion Coefficient (ADC) map (in `m^2/s`) of a phantom generated by `genPhant()`
    
    Args:
        arrPht: phantom
        B0: field strength in Tesla
        
    Returns:
        ADC map of the given phantom
    """
    return Enum2Para(arrPht, B0, "ADC")

def Enum2Om(arrPht:NDArray, B0:str|int|float) -> NDArray:
    """
    get off-resonance map (in `rad/s`) of a phantom generated by `genPhant()`
    
    Args:
        arrPht: phantom
        B0: field strength in Tesla
        
    Returns:
        off-resonance map of the given phantom
    """
    return Enum2Para(arrPht, B0, "Om")
