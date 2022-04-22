import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_mf2",
        "write_mf2"
        ]


allowed_mt = ( 151 )

mf = 2
mt=151

def read_mf2(tape, mat):
   
    if mt in allowed_mt:
        raise sandy.Error("'MF={mf}/MT={mt}' not yet implemented")
    else:
        raise ValueError("'MF={mf}/MT={mt}' not allowed")
       
        

    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    add = {
            "ZAI": C.C1,
            "ABN": C.C2,
            "LFW": C.L2,
            "NER": C.N1,
            }
    out.update(add)
    C, i = sandy.read_cont(df, i)
    add = {
            "EL": C.C1,
            "EH": C.C2,
            "LRU": C.L1,
            "LRF": C.L2,
            "NRO": C.N1,
            "NAPS":C.N2,
            }
    out.update(add)
    if "LRU" == 0:
        C, i = sandy.read_cont(df, i)
        add = {
                "SPI": C.C1,
                "AP": C.C2,
                "NLS": C.N1,
                }
        out.update(add)
    if "LRU" == 1:
        if "LRF" == 1 and  "LRF" == 2 :
            T, i = sandy.read_tab1(df, i)
            add = {
                    "NR": T.NBT,
                    "NP": T.INT,
                    "E_int": T.x,
                    "AP": T.y,
                    }
            out.update(add)
            
            C, i = sandy.read_cont(df, i)
            add = {
                    "SPI": C.C1,
                    "AP": C.C2,
                    "NLS": C.N1,
                    }
            out.update(add)
            
            """L, i = sandy.read_list(df, i)
            add = {
                    "AWRI": L.C1,
                    "QX": L.C2,
                    "L": L.L1,
                    "LRX": L.L2,
                    "6*NRS": L.NPL,
                    "NRS": L.N2
                    }
            List = np.array(L.B).reshape(NRS/6,6)
            out.update(add)"""
            
        if "LRF" == 3:
           
       
        if "LRF" == 4:
            
        
            
        if "LRF" == 7:
            
        
    return out 




def write_mf2(sec):
    
    
    return out 

