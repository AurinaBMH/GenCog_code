
library(R.matlab)
library(tnet)
#data = readMat("/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/FinalCorrections/iFOD2_cust100AVERAGEcountCortex.mat") 

data = readMat("/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/FinalCorrections/FACT_cust100_commonConnections_edges.mat") 
# Adj_cust100=data$Average.matrix.50
#g  <- graph.adjacency(Adj,weighted=TRUE)
df <- data$edges;
head(df)
# define prominent factor (degree (k) or strength (s) and resufle "links" ar "weights")
out <- weighted_richclub_w(df, rich="k", reshuffle="links", NR=150, nbins = 100; seed=NULL)
# Plot output
plot(out[,c("x","y")], type="b", log="x", xlim=c( min(out[,"x"]), max(out[,"x"]+1)), ylim=c(0,max(out[,"y"])+0.5), xaxs = "i", yaxs = "i", ylab="Weighted Rich-club Effect", xlab="Prominence (strength greater than)")
lines(out[,"x"], rep(1, nrow(out)))
View(out)

