import math

import numpy as np
from scipy.io import loadmat
from scipy.signal import hilbert
from typing import Tuple

from .objects import DataOutputStruct, InfoStruct
from ..transforms import scanConvert

def readFileInfo(path: str) -> Tuple[InfoStruct, np.ndarray]:
    """Reads the metadata and RF data from a GE .mat file.
    
    Args:
        path (str): Path to the .mat file.
        
    Returns:
        Tuple[InfoStruct, np.ndarray]: Metadata and RF data.
    """
    input = loadmat(path)
    rfData = input["frameRF"][0][0][0]

    Info = InfoStruct()
    Info.minFrequency = 1000000
    Info.maxFrequency = 6000000
    Info.lowBandFreq = 2000000
    Info.upBandFreq = 5000000
    Info.centerFrequency = 3500000
    Info.clipFact = 1
    Info.dynRange = 55

    # Scan Convert Settings
    Info.tilt1 = 0
    Info.width_rad = input["frameRF"][0][0][6][0][0]
    Info.width1 = np.rad2deg(Info.width_rad)

    Info.lines = input["frameRF"][0][0][2][0][0]  # number of scan lines (beams)
    Info.samples = input["frameRF"][0][0][3][0][0]  # samples per beam
    Info.depth = input["frameRF"][0][0][7][0][0] * 10  # in mm
    Info.startDepth = input['frameRF'][0][0][5][0][0] # m
    Info.endDepth = input['frameRF'][0][0][7][0][0]/100 # m
    Info.samplingFrequency = input["frameRF"][0][0][10][0][0]  # in Hz

    return Info, rfData

def readFileImg(rfData: np.ndarray, Info: InfoStruct) -> Tuple[DataOutputStruct, InfoStruct]:
    """Reads the RF data and adds to scan metadata.
    
    Args:
        rfData (np.ndarray): RF data.
        Info (InfoStruct): RF Scan Metadata.
        
    Returns:
        Tuple[DataOutputStruct, InfoStruct]: RF data and metadata.
    """   
    # Apply Hilbert transform to each beam to get envelope-detected signal
    bmode = np.zeros_like(rfData)
    for i in range(rfData.shape[1]):
        envelope = np.abs(hilbert(rfData[:, i]))
        bmode[:, i] = 20 * np.log10(envelope + 1e-8)  # avoid log(0) with small constant

    # Clip and normalize B-mode image
    clippedMax = Info.clipFact*np.amax(bmode)
    bmode = np.clip(bmode, clippedMax-Info.dynRange, clippedMax)
    bmode -= np.amin(bmode)
    bmode *= (255/np.amax(bmode))

    scBmodeStruct, hCm1, wCm1 = scanConvert(bmode, Info.width1, Info.tilt1, Info.startDepth, Info.startDepth+Info.endDepth)
    Info.depth = hCm1*10
    Info.width = wCm1*10
    Info.lateralRes = Info.width/scBmodeStruct.scArr.shape[1]
    Info.axialRes = Info.depth/scBmodeStruct.scArr.shape[0]

    Data = DataOutputStruct()
    Data.scBmodeStruct = scBmodeStruct
    Data.rf = rfData
    
    Data.bMode = bmode
    Data.scBmode = Data.scBmodeStruct.scArr
    
    return Data, Info


def geRfParser(filePath: str, phantomPath: str) \
    -> Tuple[DataOutputStruct, InfoStruct, DataOutputStruct, InfoStruct]:
    """Parses GE RF data and metadata.
    
    Args:
        filePath (str): Path to the RF data.
        phantomPath (str): Path to the phantom data.
    
    Returns:
        Tuple: RF data and metadata for image and phantom.
    """

    imgInfoStruct, rfData = readFileInfo(filePath)
    imgDataStruct, imgInfoStruct = readFileImg(rfData, imgInfoStruct)

    refInfoStruct, refRfData = readFileInfo(phantomPath)
    refDataStruct, refInfoStruct = readFileImg(refRfData, refInfoStruct)

    return imgDataStruct, imgInfoStruct, refDataStruct, refInfoStruct
