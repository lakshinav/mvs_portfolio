require(rmgarchMy)
uniGarchSpec <- ugarchspec(mean.model=list(armaOrder=c(1,0)),
                           variance.model=list(model="fGARCH",submodel="GARCH",
                                               garchOrder=c(1,1)))
multiGarchSpec <- multispec(replicate(2, uniGarchSpec))
copSpec <- cgarchspec(multiGarchSpec,
                      dccOrder = c(1, 1), asymmetric = T,
                      distribution.model = list(copula = "mvt",
                                                method = "ML", 
                                                time.varying = T,
                                                transformation = "empirical")) # parametric
copFit <- cgarchfit(copSpec, ret[[1]][,2:3], out.sample = 0, fit.control = list(eval.se=FALSE))
