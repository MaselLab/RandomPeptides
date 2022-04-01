# Predict fitness from AA comp.
library(stringr)
# Function for calculating fitness given an amino acid sequence.
# Requires stringr. Can only handle standard amino acids, but can be
# set to ignore all non-standard/unknown amino acids. To break when a
# non-standard/unknown amino acid is encountered, set "non.standard = 'stop'".
# To continue anyway, set "non.standard = 'continue'".
#
# Several models with different transforms are listed. It is highly recommended
# to use the untransformed model. Simply comment out whichever one is being
# used, and made sure the others are comments.
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
  # Old log2 fitness model with the wrong weights.
  # fitness.ln <- - 2.24810*L.freq + 3.39860*P.freq - 4.19921*M.freq - 0.06214*W.freq + 2.92215*A.freq - 
  #   2.53068*V.freq - 5.31797*F.freq - 7.89332*I.freq + 0.84818*G.freq + 1.74630*S.freq -
  #   0.47038*T.freq - 1.21066*C.freq - 2.34759*N.freq + 0.26526*Q.freq - 3.83773*Y.freq -
  #   3.77098*H.freq + 0.49807*D.freq - 1.43710*E.freq - 4.48956*K.freq - 2.15506*R.freq
  # New log2 fitness model with the right weights.
  # fitness.ln <- - 2.9001*L.freq + 3.8826*P.freq - 4.9114*M.freq - 0.4731*W.freq + 2.7177*A.freq -
  #   3.0156*V.freq - 5.6437*F.freq - 9.6877*I.freq + 1.3503*G.freq + 0.8549*S.freq -
  #   0.8635*T.freq - 1.7250*C.freq - 3.2173*N.freq + 0.1350*Q.freq - 4.0855*Y.freq -
  #   4.4536*H.freq + 2.2677*D.freq - 1.9224*E.freq - 4.1495*K.freq - 1.7283*R.freq
  # Untransformed fitness model.
  fitness.ln <- - 0.3703*L.freq + 2.9803*P.freq - 1.3868*M.freq + 1.8869*W.freq + 3.2314*A.freq -
    0.7650*V.freq - 0.6700*F.freq - 3.7534*I.freq + 1.9088*G.freq + 2.1824*S.freq +
    0.3213*T.freq + 0.5560*C.freq - 2.4857*N.freq + 2.2325*Q.freq - 1.9062*Y.freq -
    2.2881*H.freq + 3.3342*D.freq - 1.1222*E.freq + 0.3742*K.freq + 0.3866*R.freq
  # Fixed-effects only log2 fitness transformed.
  # fitness.ln <- - 0.44510*L.freq + 0.65312*P.freq - 3.21922*M.freq + 2.07656*W.freq + 0.67450*A.freq -
  #   0.74113*V.freq + 0.83128*F.freq + 0.63621*I.freq - 0.09014*G.freq + 2.93294*S.freq -
  #   0.22391*T.freq - 0.56780*C.freq - 3.64141*N.freq - 2.35929*Q.freq + 1.11864*Y.freq +
  #   0.21332*H.freq - 0.18255*D.freq - 4.82163*E.freq + 5.17770*K.freq + 1.13105*R.freq
  # Fixed-effect only non-transformed.
  # fitness.ln <- - 0.3669*L.freq + 2.8726*P.freq - 1.6441*M.freq + 1.7699*W.freq + 2.9968*A.freq -
  #   0.7867*V.freq - 0.8554*F.freq - 3.5042*I.freq + 1.8184*G.freq + 2.3805*S.freq +
  #   0.5263*T.freq + 0.4163*C.freq - 2.3516*N.freq + 1.9419*Q.freq - 1.5624*Y.freq -
  #   1.8172*H.freq + 3.2542*D.freq - 1.3465*E.freq + 0.7604*K.freq + 0.4455*R.freq
  fitness <- exp(fitness.ln)
  return(fitness)
}

