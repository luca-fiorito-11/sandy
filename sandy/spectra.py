import pandas as pd

__author__ = "Luca Fiorito"
__all__ = [
        "custom_spectra",
        "get_custom_spectrum",
        ]


pd.options.display.float_format = '{:.5e}'.format


custom_spectra = {
    "PWR_UO2_0_1102": "https://fispact.ukaea.uk/wiki/images/2/29/1102_PWR-UO2-0.txt",
    "PWR_UO2_15_1102": "https://fispact.ukaea.uk/wiki/images/3/33/1102_PWR-UO2-15.txt",
    "BIGTEN_407": "https://fispact.ukaea.uk/wiki/images/a/a6/407_Bigten.txt",
    }


def get_custom_spectrum(key):
    file = custom_spectra[key]
    data = pd.read_csv(file, header=None).iloc[:-1].squeeze().astype(float)
    split = data.size // 2
    e = data.iloc[:split].values[::-1]
    # Set lower bin to 1e-5 eV
    e[0] = 1e-5
    e = pd.IntervalIndex.from_breaks(e, name="E")
    f = data.iloc[split:-1].values[::-1]
    spe = pd.Series(f, index=e).rename("SPE")
    return spe


# Allocate all spectra in "custom_spectra" as module attributes
for k in custom_spectra:
    exec(f"{k} = get_custom_spectrum(k)")
