from . import utility
from importlib.resources import files

with files(__package__).joinpath("nmr_para.xml").open("r") as f:
    dictNmrPara = utility.xml2dict(f)
lstTissue = list(dictNmrPara.keys())

from .Function import genPhant, genPhMap, genB0Map, genCsm, genAmp, genResAmp, genCarAmp, Enum2Para, Enum2PD, Enum2T1, Enum2T2, Enum2Adc, Enum2Om, initSS_bSSFP, initSS_FLASH, Enum2SS, Enum2T2s