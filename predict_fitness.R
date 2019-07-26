# Predict fitness from AA comp.
library(stringr)
# Function for calculating fitness given an amino acid sequence.
# Requires stringr. Can only handle standard amino acids, but can be
# set to ignore all non-standard/unknown amino acids. To break when a
# non-standard/unknown amino acid is encountered, set "non.standard = 'stop'".
# To continue anyway, set "non.standard = 'continue'".
fitness.calculator <- function(aa.sequence, non.standard = "continue"){
  seq.length <- str_length(aa.sequence)
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
  known.aa.freq <- sum(L.freq, V.freq, M.freq, F.freq, I.freq,
                       H.freq, K.freq, R.freq, D.freq, E.freq,
                       N.freq, Q.freq, Y.freq, W.freq, T.freq,
                       S.freq, C.freq, P.freq, G.freq, A.freq)
  if (isFALSE(all.equal(known.aa.freq, 1))){
    cat("Non standard amino acids present.")
    if (non.standard == "stop"){
      cat("Exiting.")
      stopifnot(isTRUE(all.equal(known.aa.freq, 1)))
    } else if (non.standard == "continue") {
      cat("Ignoring non-standard amino acids.")
    } else {
      cat("Unknown entry for \"non.standard\". Please set equal to \"stop\" or \"continue\".")
      stopifnot(known.aa.freq == 1)
    }
  }
  fitness.ln <- - 2.24810*L.freq + 3.39860*P.freq - 4.19921*M.freq - 0.06214*W.freq + 2.92215*A.freq - 
    2.53068*V.freq - 5.31797*F.freq - 7.89332*I.freq + 0.84818*G.freq + 1.74630*S.freq -
    0.47038*T.freq - 1.21066*C.freq - 2.34759*N.freq + 0.26526*Q.freq - 3.83773*Y.freq -
    3.77098*H.freq + 0.49807*D.freq - 1.43710*E.freq - 4.48956*K.freq - 2.15506*R.freq
  fitness <- exp(fitness.ln)
  return(fitness)
}
