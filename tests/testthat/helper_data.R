data(GvHD)
library(flowStats)
morphGate <- norm2Filter("FSC-H", "SSC-H", filterId = "MorphologyGate",scale = 2)
rectGate <- rectangleGate(filterId="nonDebris","FSC-H"=c(200,Inf))
