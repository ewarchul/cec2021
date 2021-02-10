library(cecs)
source("alg.R")
print("start")

results = list()
for (D in c(10, 20)) {
  startTimeT0 = base::Sys.time()
  x = 0.55
  for (i in 1:1000000) {
    x = x + x;
    x = x / 2;
    x = x * x;
    x = sqrt(x);
    x = log(x);
    x = exp(x);
    x = x / (x + 2);
  }
  endTimeT0 = base::Sys.time()
  T0 = endTimeT0 - startTimeT0
  
  starTimeT1 = base::Sys.time()
  rb_ipop_cma_esr_ppmf(
    par = runif(D, -100, 100),
    fn = function(x) { cecs::cec2021(1, x, "basic") }
  )
  endTimeT1 = base::Sys.time()
  T1 = endTimeT1 - startTimeT1

  vecT2 = c()
  for (k in 1:5) {
    starTimeT2 = base::Sys.time()
    rb_ipop_cma_esr_ppmf(
      par = runif(D, -100, 100),
      fn = function(x) { cecs::cec2021(1, x, "basic") }
    )
    endTimeT2 = base::Sys.time()
    T2 = endTimeT2 - startTimeT2
    vecT2 = c(vecT2, T2)
  }
  results[[D]] = list(
    T0 = T0
    T1 = T1,
    vecT2 = vecT2
  )
}
