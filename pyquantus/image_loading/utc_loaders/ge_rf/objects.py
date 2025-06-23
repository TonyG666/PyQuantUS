import numpy as np

class OutImStruct():
    """Output image structure for scan converted images, converted to proper x-y layout"""
    scArr: np.ndarray # The scan-converted image (e.g. B-mode)
    xmap: np.ndarray # sc (y,x) --> preSC x
    ymap: np.ndarray # sc (y,x) --> preSC y

class DataOutputStruct():
    """Data output structure for general RF/IQ data.
    holding raw and processed data together"""
    def __init__(self):
        self.scBmodeStruct: OutImStruct
        self.scBmode: np.ndarray
        self.rf: np.ndarray # directly holds frameRF.data
        self.bMode: np.ndarray # B-mode image reconstructed from rf
        self.widthPixels: int # Corresponds to nBeams
        self.depthPixels: int # Corresponds to nSamplesPerBeam
        
class InfoStruct():
    """Metadata structure for RF/IQ data.
    help interpret or reconstruct the image from RF data"""
    def __init__(self):
        self.minFrequency: int # transducer freq (Hz)
        self.maxFrequency: int # transducer freq (Hz)
        self.lowBandFreq: int # analysis freq (Hz)
        self.upBandFreq: int # analysis freq (Hz)
        self.centerFrequency: int #Hz
        self.clipFact: float # %
        self.dynRange: int
        self.samples: int # nSamplesPerBeam = 13216
        self.lines: int # nBeams = 272
        self.depthOffset: float
        self.depth: float # DispDepth_cm = 20.3 cm
        self.width: float 
        self.samplingFrequency: int # SamplingRateHz = 50 MHz
        self.lineDensity: int
        self.lateralRes: float
        self.axialRes: float
        self.yResRF: float
        self.xResRF: float
        self.quad2x: float
        self.numSonoCTAngles: int
        self.numSamplesDrOut: int
        
        # For scan conversion
        self.tilt1: float
        self.width_rad: float
        self.width1: float
        self.startDepth: float
        self.endDepth: float
        self.endHeight: float