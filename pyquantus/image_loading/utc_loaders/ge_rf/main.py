from pathlib import Path

from ....data_objs.image import UltrasoundRfImage
from .parser import geRfParser

class EntryClass(UltrasoundRfImage):
    """
    Class for GE RF image data.
    This class is used to parse RF data from GE ultrasound machines.
    """
    extensions = [".mat"]
    spatial_dims = 2
    
    def __init__(self, scan_path: str, phantom_path: str, **kwargs):
        super().__init__(scan_path, phantom_path)
        
        # Supported file extensions for this loader
        assert Path(scan_path).suffix in self.extensions, f"File must end with {self.extensions}"
        
        imgData, imgInfo, refData, refInfo = geRfParser(scan_path, phantom_path)
        self.rf_data = imgData.rf
        self.phantom_rf_data = refData.rf
        self.bmode = imgData.bMode
        self.axial_res = imgInfo.depth / imgData.rf.shape[0]
        self.lateral_res = self.axial_res * (
            imgData.rf.shape[0] / imgData.rf.shape[1]
        ) # placeholder
        self.sc_axial_res = imgInfo.axialRes
        self.sc_lateral_res = imgInfo.lateralRes
        self.sc_bmode = imgData.scBmode
        self.xmap = imgData.scBmodeStruct.xmap
        self.ymap = imgData.scBmodeStruct.ymap
        self.tilt = imgInfo.tilt1
        self.width = imgInfo.width1
        self.start_depth = imgInfo.startDepth
        self.end_depth = imgInfo.endDepth