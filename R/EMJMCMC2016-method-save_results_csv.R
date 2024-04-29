EMJMCMC2016$methods(
  # save big.data results, if the latter are available
  save_results_csv = function(statistics1, filename = tmpfile()) {
    bigmemory::write.big.matrix(x = statistics1, filename = filename)
  }
)
