import math

#FFMC
def calculate_ffmc(temp, rh, wind, rain, prev_ffmc):
    mo = 147.2 * (101.0 - prev_ffmc) / (59.5 + prev_ffmc)
    if rain > 0.5:
        rf = rain - 0.5
        if mo > 150.0:
            mo = mo + 42.5 * rf * math.exp(-100.0 / (251.0 - mo)) * (1 - math.exp(-6.93 / rf))
        else:
            mo = mo + 42.5 * rf * math.exp(-100.0 / (251.0 - mo)) * (1 - math.exp(-6.93 / rf))
        if mo > 250.0:
            mo = 250.0

    ed = 0.942 * rh**0.679 + (11 * math.exp((rh - 100) / 10)) + 0.18 * (21.1 - temp) * (1 - math.exp(-0.115 * rh))
    if mo < ed:
        ew = 0.618 * rh**0.753 + (10 * math.exp((rh - 100) / 10)) + 0.18 * (21.1 - temp) * (1 - math.exp(-0.115 * rh))
        kl = 0.424 * (1.0 - ((100.0 - rh) / 100.0)**1.7) + 0.0694 * math.sqrt(wind) * (1.0 - ((100.0 - rh) / 100.0)**8)
        kw = kl * 0.581 * math.exp(0.0365 * temp)
        mo = ew - (ew - mo) * math.exp(-kw)
    elif mo > ed:
        kl = 0.424 * (1.0 - (rh / 100.0)**1.7) + 0.0694 * math.sqrt(wind) * (1.0 - (rh / 100.0)**8)
        kw = kl * 0.581 * math.exp(0.0365 * temp)
        mo = ed + (mo - ed) * math.exp(-kw)

    ffmc = (59.5 * (250.0 - mo)) / (147.2 + mo)
    if ffmc > 101.0:
        ffmc = 101.0
    if ffmc < 0.0:
        ffmc = 0.0
    return ffmc

# DMC
def calculate_dmc(temp, rh, rain, month, prev_dmc):
    el = [6.5, 7.5, 9.0, 12.8, 15.6, 16.4, 16.0, 14.0, 12.0, 10.8, 9.0, 7.0]
    rk = 1.894 * (temp + 1.1) * (100.0 - rh) * el[month - 1] * 0.0001

    if rain > 1.5:
        ra = rain
        rw = 0.92 * ra - 1.27
        smi = 800.0 * math.exp(-prev_dmc / 43.43)
        dr = prev_dmc - 43.43 * math.log(1.0 + ((3.937 * rw) / smi))
        if dr < 0.0:
            dr = 0.0
    else:
        dr = prev_dmc

    dmc = dr + rk
    if dmc < 0.0:
        dmc = 0.0
    return dmc

# DC
def calculate_dc(temp, rain, month, prev_dc):
    fl = [0.75, 1.0, 1.4, 2.0, 2.7, 3.2, 3.1, 2.7, 2.0, 1.5, 1.0, 0.75]
    pe = (0.36 * (temp + 2.8) + fl[month - 1]) / 2

    if rain > 2.8:
        ra = rain
        rw = 0.83 * ra - 1.27
        smi = 800.0 * math.exp(-prev_dc / 400.0)
        dr = prev_dc - 400.0 * math.log(1.0 + (3.937 * rw / smi))
        if dr < 0.0:
            dr = 0.0
    else:
        dr = prev_dc

    dc = dr + pe
    return dc

def calculate_isi(ffmc, wind):
    mo = 147.2 * (101.0 - ffmc) / (59.5 + ffmc)
    ff = 19.115 * math.exp(-0.1386 * mo) * (1.0 + (mo**5.31) / 4.93e7)
    isi = ff * math.exp(0.05039 * wind)
    return isi

def calculate_bui(dmc, dc):
    if dmc <= 0.4 * dc:
        bui = (0.8 * dc * dmc) / (dmc + 0.4 * dc)
    else:
        bui = dmc - (1.0 - (0.8 * dc) / (dmc + 0.4 * dc)) * (0.92 + (0.0114 * dmc)**1.7)
    if bui < 0:
        bui = 0
    return bui

def calculate_fwi(isi, bui):
    if bui <= 80.0:
        fD = 0.626 * bui**0.809 + 2.0
    else:
        fD = 1000.0 / (25.0 + 108.64 * math.exp(-0.023 * bui))
    b = 0.1 * isi * fD
    if b > 1.0:
        fwi = math.exp(2.72 * (0.434 * math.log(b))**0.647)
    else:
        fwi = b
    return fwi

#  Example usage with weather inputs:
# You can plug in daily values here
temperature = -4  # °C
humidity = 40     # %
wind_speed = 2.6   # km/h
rainfall = 0      # mm
month = 1         # May
day =1
yr = 2003

# Previous day's moisture codes (can start with typical defaults)
ffmc_prev = 85.0
dmc_prev = 6.0
dc_prev = 150.0

# 🔁 Calculate each index
ffmc = calculate_ffmc(temperature, humidity, wind_speed, rainfall, ffmc_prev)
dmc = calculate_dmc(temperature, humidity, rainfall, month, dmc_prev)
dc = calculate_dc(temperature, rainfall, month, dc_prev)
isi = calculate_isi(ffmc, wind_speed)
bui = calculate_bui(dmc, dc)
fwi = calculate_fwi(isi, bui)

# ✅ Output results
print("FFMC:", round(ffmc, 2))
print("DMC:", round(dmc, 2))
print("DC:", round(dc, 2))
print("ISI:", round(isi, 2))
print("BUI:", round(bui, 2))
print("FWI:", round(fwi, 2))
