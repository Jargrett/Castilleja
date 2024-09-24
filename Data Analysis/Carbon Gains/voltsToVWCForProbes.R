
# write a function for mv to vwc

# lets get it

VoltsToVWC <- function(V) {
  # first, convert volts to mV
  mV <- V*1000
  # now, run equation from the Taros guide
  VWC <- (5.439*(10^-10)*(mV^3)) - (2.731*(10^-6)*(mV^2)) + (4.868*(10^-3)*mV) - 2.683
  print(VWC)
}

VoltsToVWC()

#save(VoltsToVWC, file = "/Users/rileythoen/Dropbox/105_dissertation/13_usefulRFunctions/voltsToVWC.Rdata")

#load("/Users/rileythoen/Dropbox/105_dissertation/13_usefulRFunctions/voltsToVWC.Rdata")
VoltsToVWC(1.1)

# test

VoltsToVWC(1.478)
#30.2%

VoltsToVWC(1.464)
#29.7%

VoltsToVWC(1.585)
#33.8%

VoltsToVWC(1.840)
#41.6%

VoltsToVWC(1.667)
#36.2%


