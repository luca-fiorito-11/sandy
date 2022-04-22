import pdb

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_mf2",
        "write_mf2"
        ]
import sandy

allowed_mt = ( 151 )

mf = 2
mt=151

def read_mf2(tape, mat):
    out = _read_respar(tape, mat)
    if mt in allowed_mt:
        raise sandy.Error("'MF={mf}/MT={mt}' not yet implemented")
    else:
        raise ValueError("'MF={mf}/MT={mt}' not allowed")
    return out   
        
        
def _read_respar(tape, mat):
    
    return out 




def write_mf2(sec):
    out = _write_respar(sec)
    return out    
        
        
def _write_respar(sec):
    
    
    return out 
    
