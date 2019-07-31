# Predict fitness from AA comp.
library(stringr)
# Function for calculating fitness given an amino acid sequence.
# Requires stringr. Can only handle standard amino acids, but can be
# set to ignore all non-standard/unknown amino acids. To break when a
# non-standard/unknown amino acid is encountered, set "non.standard = 'stop'".
# To continue anyway, set "non.standard = 'continue'".
fitness.calculator <- function(aa.sequence, non.standard = "continue"){
  L.count <- str_count(aa.sequence, pattern = "L")
  V.count <- str_count(aa.sequence, pattern = "V")
  M.count <- str_count(aa.sequence, pattern = "M")
  F.count <- str_count(aa.sequence, pattern = "F")
  I.count <- str_count(aa.sequence, pattern = "I")
  H.count <- str_count(aa.sequence, pattern = "H")
  K.count <- str_count(aa.sequence, pattern = "K")
  R.count <- str_count(aa.sequence, pattern = "R")
  D.count <- str_count(aa.sequence, pattern = "D")
  E.count <- str_count(aa.sequence, pattern = "E")
  N.count <- str_count(aa.sequence, pattern = "N")
  Q.count <- str_count(aa.sequence, pattern = "Q")
  Y.count <- str_count(aa.sequence, pattern = "Y")
  W.count <- str_count(aa.sequence, pattern = "W")
  T.count <- str_count(aa.sequence, pattern = "T")
  S.count <- str_count(aa.sequence, pattern = "S")
  C.count <- str_count(aa.sequence, pattern = "C")
  P.count <- str_count(aa.sequence, pattern = "P")
  G.count <- str_count(aa.sequence, pattern = "G")
  A.count <- str_count(aa.sequence, pattern = "A")
  seq.length <- L.count + V.count + M.count + F.count + I.count +
    H.count + K.count + R.count + D.count + E.count +
    N.count + Q.count + Y.count + W.count + T.count +
    S.count + C.count + P.count + G.count + A.count
  L.freq <- str_count(aa.sequence, pattern = "L") / seq.length
  V.freq <- str_count(aa.sequence, pattern = "V") / seq.length
  M.freq <- str_count(aa.sequence, pattern = "M") / seq.length
  F.freq <- str_count(aa.sequence, pattern = "F") / seq.length
  I.freq <- str_count(aa.sequence, pattern = "I") / seq.length
  H.freq <- str_count(aa.sequence, pattern = "H") / seq.length
  K.freq <- str_count(aa.sequence, pattern = "K") / seq.length
  R.freq <- str_count(aa.sequence, pattern = "R") / seq.length
  D.freq <- str_count(aa.sequence, pattern = "D") / seq.length
  E.freq <- str_count(aa.sequence, pattern = "E") / seq.length
  N.freq <- str_count(aa.sequence, pattern = "N") / seq.length
  Q.freq <- str_count(aa.sequence, pattern = "Q") / seq.length
  Y.freq <- str_count(aa.sequence, pattern = "Y") / seq.length
  W.freq <- str_count(aa.sequence, pattern = "W") / seq.length
  T.freq <- str_count(aa.sequence, pattern = "T") / seq.length
  S.freq <- str_count(aa.sequence, pattern = "S") / seq.length
  C.freq <- str_count(aa.sequence, pattern = "C") / seq.length
  P.freq <- str_count(aa.sequence, pattern = "P") / seq.length
  G.freq <- str_count(aa.sequence, pattern = "G") / seq.length
  A.freq <- str_count(aa.sequence, pattern = "A") / seq.length
  #print(seq.length)
  #print(str_length(aa.sequence))
  if (isFALSE(identical(seq.length, str_length(aa.sequence)))){
    print("Non standard amino acids present.")
    if (non.standard == "stop"){
      print("Exiting.")
      stopifnot(isTRUE(identical(seq.length, str_length(aa.sequence))))
    } else if (non.standard == "continue") {
      print("Ignoring non-standard amino acids.")
    } else {
      print("Unknown entry for \"non.standard\". Please set equal to \"stop\" or \"continue\".")
      stopifnot(isTRUE(identical(seq.length, str_length(aa.sequence))))
    }
  }
  fitness.ln <- - 2.24810*L.freq + 3.39860*P.freq - 4.19921*M.freq - 0.06214*W.freq + 2.92215*A.freq - 
    2.53068*V.freq - 5.31797*F.freq - 7.89332*I.freq + 0.84818*G.freq + 1.74630*S.freq -
    0.47038*T.freq - 1.21066*C.freq - 2.34759*N.freq + 0.26526*Q.freq - 3.83773*Y.freq -
    3.77098*H.freq + 0.49807*D.freq - 1.43710*E.freq - 4.48956*K.freq - 2.15506*R.freq
  fitness <- exp(fitness.ln)
  return(fitness)
}

