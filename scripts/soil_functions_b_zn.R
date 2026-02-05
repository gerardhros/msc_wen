b_conc  = 1000*b/KD[nut %in% 'b']

# estimate B-water (mg/kg) from B-mehlig3 (mg/kg)
.cfunB <- function(bm3){exp(-1.958732 + 1.295274 * log(bm3+0.0001))}

# calculate B-water (mg/kg) estimated from B-mehlig3 (100* mg/kg) using funciton .cfunB
inp[,b := .cfunB(m3.B)]

# calculate C concentraiotn in soil soluion (ug / L)
b_conc  = 1000*b/KD[nut %in% 'b']

# the OCP model calculates supply as concentration x transpirationflux x stream concentraiton factor
# transpiration => from Nikos
# scf is in table modinp
dt[,scf_b := 0.894]
dt[,tr := 750]
dt[,bs := b_conc * scf_b * tr * 0.001]


# similar for zinc

# soil solution concentration (all in mg/L, except for B, Zn and Cu = ug/L)
dt[,zn_conc = .zn_cf(zn, soc, clay, ph)]
dt[,scf_zn := 1.113]
dt[,zns := zn_conc * scf_zn * tr* 0.001] # check the 0.001 for unit correction, just a placeholer here for now

calc_zn_supply(zn_m3, soc, clay, ph){
  
  # output gives Zn supply in kg S/ha
}



# function to calculate zn concentration in soil solution 
.zn_cf <- function(sct, soc, clay, ph, a0 = -4.51, a1 = 0.39, a2 = 0.35,
                   a3 = 0.45, aN = 0.74, b0 = 0.428, b1 = 1.235, b2 = 0.183,
                   b3 = -0.298, e = 65.4){ 
  # reactieve soil content, mmol/ kg
  zn_sctr <- 10^(b0 + b1*log10(sct/(e*1000)) + b2*log10(soc) + b3*log10(clay))
  # Kf, -
  zn_kf   <- 10^(a0 + a1*log10(soc) + a2*log10(clay) + a3*ph)
  # soil solution concentration, ug/ L (= mmol/L * mg/mmol * ug/mg)
  zn_conc <- (zn_sctr/zn_kf)^(1/aN)*e*1000
  # return soil solution concentration
  return(zn_conc)
}

# irrigation, precipitation, evapotranspiration, transpiration rate (mm/year) and precipitation surplus sp
inp[,ir := ifelse(country == 'Mali', ifelse(season == 'wet', 1100, 1500),ifelse(season == 'wet', 1000, 1400))]

# general model input for running soil model OCP
modinp = data.frame(
  # name of the nutrients used in the calculations
  nut = c('n', 'pox', 'k', 'ca', 'mg', 's', 'b', 'zn', 'cu'),
  # name of the fertilizer nutrient
  nutF = c('N', 'P', 'K', 'Ca', 'Mg', 'S', 'B', 'Zn', 'Cu'),
  # unit correction for calculating availability in kg/ ha
  ucor = c(0.01, rep(0.01, 5), rep(1e-5, 3)),
  # unit correction to calculate soil status (e = 31 for P)
  sucor = c(0.1, 100 / 31, rep(0.1, 4), rep(100, 3)),
  # optimum crop content, g/kg product
  Xcavopt = c(20, 4.5, 12, 0.2, 0.2, 2, 0.3, 1.5, 0.4),
  # ration between carbon and nutrient
  XCratio = pmax(1 / c(15, 405, 300, 110, 540, 119, NA, NA, NA),0, na.rm = TRUE),
  # kD, linear relation between solid phase and soil solution, L/kg
  KD = c(NA, NA, 13.75, 11.33, 24, NA, 17, NA, NA),
  # fraction of applied fertilizer that will be taken up by the crop
  #frXav = c(0.4, 0.4, 0.6, 0.5, 0.25, 0.6, 0.5, 0.5, 0.5),
  frXav = c(0.4, 0.4, 0.6, 0.5, 0.50, 0.6, 0.25,0.2,0.4),
  # stream concentration factors (-)
  Xscf	= c(0.5937, 0.46,1.15, 0.007,0.025,0.788,0.894,1.113,0.047),
  # Crop content, derived from Reuter et al, 1986 (g/kg), whole plant
  ctXcr = c(24.7,4.05,25.6,4.1,3.35,2.75,0.005215,0.0475,0.0058),
  # manual correction possible to adjust soil supply
  ssc = rep(1,9),
  # nutrient input via deposition
  Xdpst = c(2.5,0.3,5.4,rep(0,6)),
  # nurient input via irrigation
  Xirri = c(0.25,0.2,3.5,rep(0,6)),
  # prevent strings to be converted to factors
  stringsAsFactors = FALSE
) 
