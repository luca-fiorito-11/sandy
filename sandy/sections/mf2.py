import sandy
import numpy as np 

__author__ = "Rayan HADDAD"
__all__ = [
        "read_mf2",
        ]


mf = 2
mt = 151

def read_mf2(tape, mat):
        

    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,         
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    add = {
            "ZA": C.C1,             
            "AWR": C.C2,            #AWR is defines as the ratio of the mass of the material to that of the neutron
            "NIS": C.N1,            #NIS is the number of isotopes in the material
            }
    out.update(add)
    C, i = sandy.read_cont(df, i)
    add = {
            "ZAI": C.C1,
            "ABN": C.C2,           #Abundance of an isotope in the material 
            "LFW": C.L2,           #indication whether average fission wifths are given in the unresolbed resonance region
            "NER": C.N1,           #Number of resonance energy ranges for this isotope. 
            }    
    LFW = int(C.L2) 
    out.update(add)
    
    C, i = sandy.read_cont(df, i)
    add = {
            "EL": C.C1,             #Lower limit for an energy range
            "EH": C.C2,             #Upper limit for an energy range
            "LRU": C.L1,            #Flag indicating whether this energy range contains data for resolved or unresolved resonance parameters:
            "LRF": C.L2,            #Flag indicating which representation has been used for the energy range. 
            "NRO": C.N1,            #Flag designating possible energy dependence of the scattering radiu
            "NAPS":C.N2,            #Flag controlling the use of the two radii
            }
    LRU = int(C.L1)
    LRF = int(C.L2)
    NRO = int(C.N1)
    
    out.update(add)
    
    print (LRU, LRF)
    
    if LRU == 0:
        C, i = sandy.read_cont(df, i)
        add = {
                "SPI": C.C1,        #Spin, I, of the target nucleus.
                "AP": C.C2,         #Scattering radius in units of 10e-12cm.
                "NLS": C.N1,        #Number of l-values (neutron orbital angular momentum) in this energy region.
                }
        out.update(add)
        
    if LRU == 1:
        if LRF == 1 or LRF == 2 :
            
            if NRO == 0 :
                C, i = sandy.read_cont(df, i)
                add = {
                        "SPI": C.C1,
                        "AP": C.C2,
                        "NLS": C.N1,
                        }
                out.update(add)
                L, i = sandy.read_list(df, i)
                NRS= int(L.N2)
                add = {
                        "AWRI": L.C1,       #Ratio of the mass of a particular isotope to that of a neutron
                        "QX": L.C2,         #Q-value to be added to the incident particle’s center-of-mass energy to determine the channel energy for use in the penetrability factor.
                        "L": L.L1,          #Value of l.
                        "LRX": L.L2,        #Flag indicating whether this energy range contains a competitive width
                        "6*NRS": L.NPL,     
                        "NRS": L.N2,        #Number of resolved resonances for a given l-value.
                        "List" : np.array(L.B).reshape(NRS,6),
                        }
                out.update(add)
                print (NRS)
            
            else :                
                T, i = sandy.read_tab1(df, i)
                add = {
                        "NR": T.NBT,    
                        "NP": T.INT,
                        "E_int": T.x,
                        "AP": T.y,
                        }
                out.update(add)
        
        if LRF == 3 :
            
            if NRO == 0 :
               C, i = sandy.read_cont(df, i)
               add = {
                       "SPI": C.C1,
                       "AP": C.C2,
                       "LAD": C.L1,     #Flag indicating whether these parameters can be used to compute angular distributions
                       "NLS": C.N1,     
                       "NLSC": C.N2,    #Number of l-values
                       }
               NLSC = int(C.N2)
               out.update(add)
               for j in range(NLSC) :
                   L, i = sandy.read_list(df, i)
                   NRS= int(L.N2)
                   add = {
                            "AWRI": L.C1,
                            "QX": L.C2,
                            "L": L.L1,
                            "LRX": L.L2,
                            "6*NRS": L.NPL,
                            "NRS": L.N2,
                            "List": L.B,
                            }
                   print (NRS)
                   out.update(add)
               """ the code here isn't complete"""
               
            else :
                T, i = sandy.read_tab1(df, i)
                add = {
                        "NR": T.NBT,    
                        "NP": T.INT,
                        "E_int": T.x,
                        "AP": T.y,
                        }
                out.update(add)
        
        if LRF == 4 :
            
            
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
            L, i = sandy.read_list(df, i)
            NPP= int(L.N2)
            add = { 
                    "AWRI": L.C1,       #Number of sets of background constants given.
                    "LI": L.L1,         #Flag to indicate the kind of parameters given
                    "6*NX": L.NPL,
                    "NX": L.N2,         #Number of sets of background constants given.
                    "AT1": L.B[0],      #Background constants for the total cross section
                    "AT2": L.B[1],      #Background constants for the total cross section
                    "AT3": L.B[2],      #Background constants for the total cross section
                    "AT4": L.B[3],      #Background constants for the total cross section
                    "BT1": L.B[4],      #Background constants for the total cross section
                    "BT2": L.B[5],      #Background constants for the total cross section
                    "AF1": L.B[6],      #Background constants for the fission cross section.
                    "AF2": L.B[7],      #Background constants for the fission cross section.
                    "AF3": L.B[8],      #Background constants for the fission cross section.
                    "AF4": L.B[9],      #Background constants for the fission cross section.
                    "BF1": L.B[10],     #Background constants for the fission cross section.
                    "BF2": L.B[11],     #Background constants for the fission cross section.
                    "AC1": L.B[12],     #Background constants for the radiative capture cross
                    "AC2": L.B[13],     #Background constants for the radiative capture cross
                    "AC3": L.B[14],     #Background constants for the radiative capture cross
                    "AC4": L.B[15],     #Background constants for the radiative capture cross
                    "BC1": L.B[16],     #Background constants for the radiative capture cross
                    "BC2": L.B[17],     #Background constants for the radiative capture cross
                    }
            out.update(add)
            C, i = sandy.read_cont(df, i)
            add = {
                    "L": C.L1,
                    "NJS": C.N1,        #Number of sets of resolved resonance parameters
                    }
            NJS = int(C.N1)
            out.update(add)
            for j in range (NJS):
                C, i = sandy.read_cont(df, i)
                add = {
                        "L": C.L1,
                        "NJS": C.N1,
                        }
                out.update(add)
                L, i = sandy.read_list(df, i)
                NLJ = int(L.N2)
                add = { 
                        "AJ": L.C1,
                        "12*NLJ": L.NPL,        #Number of resonances for which parameters are given, for a specified AJ and L.
                        "NLJ": L.N2,
                        "List": np.array(L.B).reshape(NLJ,12),
                        }
                out.update(add)
            
        
        if LRF == 7 : 
            
            C, i = sandy.read_cont(df, i)
            add = {
                    "IFG": C.L1,
                    "KMR": C.L2,
                    "NJS": C.N1,
                    "KRL": C.N2,
                    }
            KRM = int(C.L2)
            NJS = int(C.N1)
            out.update(add)
            
            if KRM == 1 or KRM == 2 or KRM == 3:
                L, i = sandy.read_list(df, i)
                NPP= int(L.N2)
                add = { 
                        "NPP": L.L1,
                        "12*NPP": L.NPL,
                        "2*NPP": L.N2,                  #Total number of particle-pairs.
                        "List": np.array(L.B).reshape(NPP,12),
                        }
                out.update(add)
                
                for j in range(NJS):
                    L, i = sandy.read_list(df, i)
                    NCH = int(L.L2)
                    add = {
                            "AJ": L.C1,             #Floating point value of J (spin); sign indicates parity
                            "PJ": L.C2,             #Parity (used only if AJ = 0.0).
                            "KBK": L.L1,            #Non-zero if background R-matrix exists
                            "KPS": L.L2,            #Non-zero if non-hard-sphere phase shift are to be specified.
                            "6*NCH": L.NPL,
                            "NCH": L.N2,            #Number of channels for the given J pi.   
                            "List":np.array(L.B).reshape(NCH,6),
                            }                        
                    out.update(add)
                    L, i = sandy.read_list(df, i)
                    NRS = int(L.L2) 
                    add = {
                            "NRS": L.L2,            #Number of resonances for the given J pi 
                            "6*NX": L.NPL,
                            "NX": L.N2,
                            "List": np.array(L.B).reshape(NCH,NRS),
                            }
                    out.update(add)      

            
    if LRU == 2 :
        if LFW == 0 and LRF == 1 : 
            C, i = sandy.read_cont(df, i)
            add = {
                    "SPI": C.C1,
                    "AP": C.C2,
                    "LSSF": C.L1,                   #Flag governing the interpretation of the File 3 cross sections
                    "NLS": C.N1,                    #Number of l-values.
                    }
            out.update(add)
            L, i = sandy.read_list(df, i)
            NJS= int(L.N2)
            add = {
                    "AWRI": L.C1,
                    "L": L.L1,
                    "6*NJS": L.NPL,
                    "NJS": L.N2,                    #Number of J-states for a particular l-state
                    "List" : np.array(L.B).reshape(NJS,6),
                    }
            out.update(add)
            
        if LFW == 1 and LRF == 1 :
            C, i = sandy.read_cont(df, i)
            add = {
                    "SPI": C.C1,
                    "AP": C.C2,
                    "LSSF": C.L1,
                    "NE": C.N1,
                    "NLS": C.N2,
                    }
            out.update(add)
            
            L, i = sandy.read_list(df, i)
            add = {
                    "AWRI": L.C1,
                    "L": L.L1,
                    "6*NjS": L.NPL,
                    "NJS": L.N2,
                    "ES" : np.array(L.B),
                    }
            out.update(add)
            
            C, i = sandy.read_cont(df, i)
            add = {
                    "AWRI": C.C1,
                    "L": C.L1,
                    "NJS": C.N1,
                    }
            out.update(add)
            
            L, i = sandy.read_list(df, i)
            add = {
                    "L": L.L1,
                    "MUF": L.L2,            #Integer value of the number of degrees of freedom for fission widths
                    "NE+6": L.NPL,          #Number of energy points at which energy-dependent widths are tabulated
                    "D": L.B[0],            #Average level spacing for resonances with spin J
                    "AJ": L.B[1],           #Floating-point value of J
                    "AMUN": L.B[2],         #Number of degrees of freedom in the neutron width distribution.
                    "GNO": L.B[3],          #Average reduced neutron width.
                    "GG": L.B[4],           #Average radiation width.
                    "ES" : np.array(L.B[6::]),          #Energy of the i th point used to tabulate energy-dependent widths
                    }
            out.update(add)
            
        if LRF == 2 : 
            
            C, i = sandy.read_cont(df, i)
            add = {
                    "SPI": C.C1,
                    "AP": C.C2,
                    "LSSF": C.L1,
                    "NLS": C.N1,
                    }
            out.update(add)
            
            C, i = sandy.read_cont(df, i)
            add = {
                    "AWRI": C.C1,
                    "L": C.L1,
                    "NJS": C.N1,
                    }
            out.update(add)
            
            L, i = sandy.read_list(df, i)
            NE= int(L.N2)
            add = {
                    "AJ": L.C1,
                    "INT": L.L2,                #Interpolation scheme to be used for interpolating between the cross sections obtained from average resonance parameters.
                    "6*NE+6": L.NPL,
                    "NE": L.N2,
                    "AMUX": L.B[2],             #Number of degrees of freedom used in the competitive width distribution.
                    "AMUN": L.B[3],             #Number of degrees of freedom in the neutron width distribution.
                    "AMUG": L.B[4],             #Number of degrees of freedom in the radiation width distribution.
                    "AMUF": L.B[5],             #Integer value of the number of degrees of freedom for fission widths.
                    "List" : np.array(L.B[6::]).reshape(NE,6),
                    }
            out.update(add)
    
    return out 





      
            

