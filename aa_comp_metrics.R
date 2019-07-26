# Amino acid composition-based predictors.
mean.metric.calculator <- function(aa.sequence, metric){
  # Metrics available:  "disorder" for disoder propensity
  #                     "stickiness" for interaction interface propensity
  #                     "RSA-hydrophobicity" for 1 - mean RSA, a measure of hydrophobicity
  #                     "RSA-hydrophilicity" for RSA, a measure of hydrophilicity
  #                     "Buried-100" for fraction of 100% buried residues with RSA = 0.
  #                     "Buried-95" for fraction of 95% buried residues with RSA < 0.05.
  #                     "fitness" for fitness effects calculated from random peptides in E. coli.
  seq.length <- str_length(aa.sequence)
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
  if (metric == "disorder"){
    disorder.propensity <- 
      0.000*C.count + 0.004*W.count + 0.090*I.count + 0.113*Y.count + 0.117*F.count +
      0.195*L.count + 0.259*H.count + 0.263*V.count + 0.285*N.count + 0.291*M.count +
      0.394*R.count + 0.401*T.count + 0.407*D.count + 0.437*G.count + 0.450*A.count +
      0.588*K.count + 0.665*Q.count + 0.713*S.count + 0.781*E.count + 1.000*P.count
    mean.disorder <- disorder.propensity / seq.length
    return(mean.disorder)
  } else if (metric == "stickiness"){
    stickiness <- 
      1.0372*C.count + 0.7925*W.count + 1.1109*I.count + 0.8806*Y.count + 1.2727*F.count +
      0.9138*L.count + 0.1204*H.count + 0.7599*V.count - 0.2693*N.count + 1.0124*M.count -
      0.0876*R.count + 0.1031*T.count - 0.7485*D.count - 0.1771*G.count + 0.0062*A.count -
      1.1806*K.count - 0.4114*Q.count + 0.1376*S.count - 0.7893*E.count - 0.1799*P.count
    mean.stickiness <- stickiness / seq.length
    return(mean.stickiness)
  } else if (metric == "RSA-hydrophobicity") {
    rsa.total <-
      0.899*C.count + 0.837*W.count + 0.875*I.count + 0.813*Y.count + 0.864*F.count +
      0.853*L.count + 0.731*H.count + 0.857*V.count + 0.658*N.count + 0.841*M.count +
      0.639*R.count + 0.728*T.count + 0.634*D.count + 0.731*G.count + 0.782*A.count +
      0.554*K.count + 0.636*Q.count + 0.722*S.count + 0.589*E.count + 0.658*P.count
    mean.rsa <- rsa.total / seq.length
    return(mean.rsa)
  } else if (metric == "RSA-hydrophilicity") {
    rsa.total <-
      0.101*C.count + 0.163*W.count + 0.125*I.count + 0.187*Y.count + 0.136*F.count +
      0.147*L.count + 0.269*H.count + 0.143*V.count + 0.342*N.count + 0.159*M.count +
      0.361*R.count + 0.272*T.count + 0.366*D.count + 0.269*G.count + 0.218*A.count +
      0.446*K.count + 0.364*Q.count + 0.278*S.count + 0.411*E.count + 0.342*P.count
    mean.rsa <- rsa.total / seq.length
    return(mean.rsa)
  } else if (metric == "Buried-100") {
    buried100 <-
      0.287*C.count + 0.0979*W.count + 0.247*I.count + 0.0797*Y.count + 0.186*F.count +
      0.213*L.count + 0.0532*H.count + 0.25*V.count + 0.0451*N.count + 0.217*M.count +
      0.0121*R.count + 0.0987*T.count + 0.0276*D.count + 0.166*G.count + 0.228*A.count +
      0.00597*K.count + 0.0289*Q.count + 0.105*S.count + 0.0183*E.count + 0.0607*P.count
    mean.buried100 <- buried100 / seq.length
    return(mean.buried100)
  } else if (metric == "Buried-95") {
    buried95 <-
      0.576*C.count + 0.368*W.count + 0.516*I.count + 0.306*Y.count + 0.483*F.count +
      0.486*L.count + 0.198*H.count + 0.494*V.count + 0.146*N.count + 0.484*M.count +
      0.0749*R.count + 0.237*T.count + 0.104*D.count + 0.291*G.count + 0.399*A.count +
      0.0283*K.count + 0.109*Q.count + 0.241*S.count + 0.0717*E.count + 0.162*P.count
    mean.buried95 <- buried95 / seq.length
    return(mean.buried95)
  } else if (metric == "fitness") {
    L.freq <- L.count / seq.length
    V.freq <- V.count / seq.length
    M.freq <- M.count / seq.length
    F.freq <- F.count / seq.length
    I.freq <- I.count / seq.length
    H.freq <- H.count / seq.length
    K.freq <- K.count / seq.length
    R.freq <- R.count / seq.length
    D.freq <- D.count / seq.length
    E.freq <- E.count / seq.length
    N.freq <- N.count / seq.length
    Q.freq <- Q.count / seq.length
    Y.freq <- Y.count / seq.length
    W.freq <- W.count / seq.length
    T.freq <- T.count / seq.length
    S.freq <- S.count / seq.length
    C.freq <- C.count / seq.length
    P.freq <- P.count / seq.length
    G.freq <- G.count / seq.length
    A.freq <- A.count / seq.length
    fitness.ln <-
      2.24810*L.freq + 3.39860*P.freq - 4.19921*M.freq - 0.06214*W.freq + 2.92215*A.freq - 
      2.53068*V.freq - 5.31797*F.freq - 7.89332*I.freq + 0.84818*G.freq + 1.74630*S.freq -
      0.47038*T.freq - 1.21066*C.freq - 2.34759*N.freq + 0.26526*Q.freq - 3.83773*Y.freq -
      3.77098*H.freq + 0.49807*D.freq - 1.43710*E.freq - 4.48956*K.freq - 2.15506*R.freq
    fitness <- exp(fitness.ln)
    return(fitness)
  } else {
    print("Error -- unkown or missing metric. Please select a metric from the list provided.")
  }
}