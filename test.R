response <- read.csv("../synergyfinder_issue/2019.09.25_Viral/cellline1.csv", 
                     fileEncoding = "latin1", stringsAsFactors = FALSE)

matrix <- ReshapeData(response)
PlotDoseResponse(matrix)
