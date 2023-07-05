#comparison with the Oracle version of loggle
rm(list = ls())

library(devtools)
#install_github(repo="jlyang1990/loggle")

library(loggle)

library(BDgraph)
library(igraph)
library(GGMselect)
library(ggplot2)
library(ggthemes)
library(crayon)
library(Matrix)
library(dplyr)
library(reshape2)
library(RColorBrewer); mycol = brewer.pal(12, "Paired")
source("DGGM.R")

#SIMULATIONS SCENARIO3 

set.seed(0)

sim_data = simulate_data(scenario = 3, replica = 1)

Y = sim_data$Y

Gtruth = sim_data$G

#thres = seq(0,100)/100 
thres = 0.5
TN = rep(0, length(thres))
TP = rep(0, length(thres))
FP = rep(0, length(thres))
FN = rep(0, length(thres))
FPR = rep(0, length(thres))
TPR = rep(0, length(thres))

for (r in 1:20){
  Y = read.csv(paste("SimData_scenario3/Ys3r",r,".csv", sep = ""), sep=" "); p = 10
  pos = c(35,135)
  result = loggle.cv(t(Y), pos = c(35,135))
  par(mfrow = c(1, 2))
  for(k in 1:length(pos)) {
    adj.matrix <- result$cv.select.result$adj.mat.opt[[k]] != 0
    net <- graph.adjacency(adj.matrix, mode =
                             "undirected", diag = FALSE)
    set.seed(0)
    plot(net, vertex.size = 10, vertex.color =
           "lightblue", vertex.label = NULL, edge.color =
           "black", layout = layout.circle)
    title(main = paste("t =",
                       round(pos[k]/(156-1), 2)), cex.main = 0.8)
  }
  i = 1
  for (tt in thres){
    print(tt)
    #1
    Gtemp = result$cv.select.result$adj.mat.opt[[1]]
    Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0
    
    TP[i] = TP[i] + sum(Gtemp * Gtruth[1,,]) / 2  #0 in independence 
    FN[i] = FN[i] + sum((1 - Gtemp) * Gtruth[1,,]) / 2 #0 in independence
    TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth[1,,])) - p) / 2 
    FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth[1,,])) / 2
    #2
    Gtemp =  result$cv.select.result$adj.mat.opt[[2]]
    Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0
    
    TP[i] = TP[i] + sum(Gtemp * Gtruth[2,,]) / 2  #0 in independence 
    FN[i] = FN[i] + sum((1 - Gtemp) * Gtruth[2,,]) / 2 #0 in independence
    TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth[2,,])) - p) / 2 
    FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth[2,,])) / 2
    i = i + 1
  }
}

#TN[1] = 0
i = 1
for (tt in thres){
  TPR[i] = TP[i] / (TP[i] + FN[i])
  FPR[i] = FP[i] / (FP[i] + TN[i]) #fall-out: among negative who turns out positive
  i = i + 1
}
Performance_scenario3_loggle = FPR
ROC3_loggle = cbind(thres, TPR, FPR)

#SIMULATIONS SCENARIO4

set.seed(0)

sim_data = simulate_data(scenario = 4, replica = 1, p = 20)

Y = sim_data$Y

Gtruth = sim_data$G

#thres = seq(0,100)/100 
thres = 0.5
TN = rep(0, length(thres))
TP = rep(0, length(thres))
FP = rep(0, length(thres))
FN = rep(0, length(thres))
FPR = rep(0, length(thres))
TPR = rep(0, length(thres))

for (r in 1:1){
  pos = c(30,80,125,175)
  result = loggle.cv(t(Y), pos = pos)
  par(mfrow = c(1, 4))
  for(k in 1:length(pos)) {
    adj.matrix <- result$cv.select.result$adj.mat.opt[[k]] != 0
    net <- graph.adjacency(adj.matrix, mode =
                             "undirected", diag = FALSE)
    set.seed(0)
    plot(net, vertex.size = 20, vertex.color =
           "lightblue", vertex.label = NULL, edge.color =
           "black", layout = layout.circle)
    title(main = paste("t =",
                       round(pos[k]/(200-1), 2)), cex.main = 1.8)
  }
  p = 20
  i = 1
  for (tt in thres){
    print(tt)
    #1
    Gtemp = result$cv.select.result$adj.mat.opt[[1]]
    Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0
    
    TP[i] = TP[i] + sum(Gtemp * Gtruth[1,,]) / 2  #0 in independence 
    FN[i] = FN[i] + sum((1 - Gtemp) * Gtruth[1,,]) / 2 #0 in independence
    TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth[1,,])) - p) / 2 
    FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth[1,,])) / 2
    #2
    Gtemp =  result$cv.select.result$adj.mat.opt[[2]]
    Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0
    
    TP[i] = TP[i] + sum(Gtemp * Gtruth[2,,]) / 2  #0 in independence 
    FN[i] = FN[i] + sum((1 - Gtemp) * Gtruth[2,,]) / 2 #0 in independence
    TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth[2,,])) - p) / 2 
    FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth[2,,])) / 2
    
    #3
    Gtemp = result$cv.select.result$adj.mat.opt[[3]]
    Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0
    
    TP[i] = TP[i] + sum(Gtemp * Gtruth[3,,]) / 2  #0 in independence 
    FN[i] = FN[i] + sum((1 - Gtemp) * Gtruth[3,,]) / 2 #0 in independence
    TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth[3,,])) - p) / 2 
    FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth[3,,])) / 2
    
    #1
    Gtemp = result$cv.select.result$adj.mat.opt[[4]]
    Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0
    
    TP[i] = TP[i] + sum(Gtemp * Gtruth[4,,]) / 2  #0 in independence 
    FN[i] = FN[i] + sum((1 - Gtemp) * Gtruth[4,,]) / 2 #0 in independence
    TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth[4,,])) - p) / 2 
    FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth[4,,])) / 2
    i = i + 1
  }
}

#TN[1] = 0
i = 1
for (tt in thres){
  TPR[i] = TP[i] / (TP[i] + FN[i])
  FPR[i] = FP[i] / (FP[i] + TN[i]) #fall-out: among negative who turns out positive
  i = i + 1
}
Performance_scenario3_loggle = FPR
ROC3_loggle = cbind(thres, TPR, FPR)



#REAL DATA 
data = read.csv("data/daily_returns.csv"); p = 9
data["date"] = as.Date(data$date, format = "%d/%m/%Y")
data[,2:(p+1)] = data[,2:(p+1)] / 100 + 1

weekly_returns = data %>% group_by(date) %>% summarise_all(prod) 
weekly_returns[3:(p+2)] = as.data.frame(weekly_returns[3:(p+2)])

mean = colMeans(log(weekly_returns[3:(p+2)]))
var = var(log(weekly_returns[3:(p+2)]))
weekly_returns[3:(p+2)] = scale( log(weekly_returns[3:(p+2)]) )
weekly_returns = as.data.frame(weekly_returns[,-c(2,12,13)])

Y = ( as.matrix( weekly_returns[,2:10] ) )
data_long = melt(weekly_returns, id.vars = "date") 
names(data_long)[2] = "Portfolio"

par(mfrow = c(1,1))
p = ggplot(data_long,                          
           aes(x = date,
               y = value,
               col = Portfolio)) +
  geom_line() + xlab("Time") + 
  ylab("Returns ") + 
  theme_classic()+ theme(text = element_text(size = 20, family = "Times"))
print(p)

ptf = c("NoDur", "Durbl", "Manuf", "Enrgy",
        "HiTec", "Telcm", "Shops", 
        "Hlth", "Utils")

window_size = 80
par(mfrow = c(2,3))
for(i in c(30,50,70,90,110,130)){
  #dates compute 
  if (i ==30){ #1 to #70
    start = "December 31,2018"
    end = "May 1, 2020"
  }else if (i==50){ #20 to #90
    start = "March 4,2019"
    end = "Sepetmber 18, 2020"
  }
  else if (i==70){ #30 to #110
    start = "July 22,2019"
    end = "February 5, 2021"
  }
  else if (i==90){ #50 to #130
    start = "December 9,2019"
    end = "June 25, 2021"
  }
  else if (i==110){ #70 to #150
    start = "April 27, 2020"
    end = "November 12, 2021"
  }
  else if (i==130){ #90 to #157
    start = "Sepetmber 14, 2020"
    end = "December 31, 2021"
  }
  #random estimate by selectFast 
  set.seed(NULL)
  G_window = selectFast(Y[max((i-window_size/2),1):min((i+window_size/2),dim(Y)[1]),])$EW$G
  #G_window = selectQE(Y[max((i-window_size/2),1):min((i+window_size/2),dim(Y)[1]),])$G
  
  plot(graph_from_adjacency_matrix(G_window, weighted = TRUE, mode = c("undirected")),
       edge.width=3, edge.color="black",
       vertex.label = ptf,
       vertex.color=mycol[11], vertex.size=20,
       vertex.label.font=2, vertex.label.color="black", vertex.label.cex=1.5,
       edge.label.font=2, edge.label.cex=1.5, 
       main = paste("From", start, "\n", "to", end ))
}

result = loggle.cv(t(Y), pos = c(31,70,118))


ptf = c("NoDur", "Durbl", "Manuf", "Enrgy",
        "HiTec", "Telcm", "Shops", 
        "Hlth", "Utils")
pos = c(31,70,118)
par(mfrow = c(1, 3))
for(k in 1:length(pos)) {
  adj.matrix <- result$cv.select.result$adj.mat.opt[[k]] != 0
  net <- graph.adjacency(adj.matrix, mode =
                           "undirected", diag = FALSE)
  set.seed(0)
  plot(net, vertex.size = 20, vertex.color =
         "lightblue", vertex.label = ptf, edge.color =
         "black", layout = layout.circle)
  title(main = paste("t =",
                     round(pos[k]/(156-1), 2)), cex.main = 1.8)
}
