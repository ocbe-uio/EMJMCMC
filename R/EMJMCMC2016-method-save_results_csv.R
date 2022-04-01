EMJMCMC2016$methods(
  # save big.data results, if the latter are available
  save_results_csv = function(statistics1, filename) {
    write.bigmemory::big.matrix(x = statistics1, filename = filename)
  }
)
