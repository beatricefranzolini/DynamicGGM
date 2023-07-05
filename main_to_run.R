#rm(list = ls())

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
windowsFonts(Times=windowsFont("TT Times New Roman"))
source("DGGM.R")

set.seed(0)
#A. select scenario 
scenario = 0 #1,2,3,4,5 (corresponding to a,b, 1,2,3) -> for simulated data; 
#0 -> weekly real data; -1 -> monthly real data 

#B. select replica for simulated data
replica = 1 #1,2,3,4,5,6,7,8,9,10 (not use for real data)

if(scenario > 0){
  
  #to replicate results in Franzolini et al. (2022)
  sim_data = simulate_data(scenario, replica)
  Y = sim_data$Y
  
  #alternatively:
  #sim_data = simulate_data(scenario, replica, seed = 1, Tmax = 200, p = 10) 
  #change seed to change graph structures, change replica to change data.
  
}else if(scenario == 0){
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
}else if (scenario == -1){
  data_m = read.csv("data/monthly_returns.csv"); p = 9
  data_m["yyyymm"] = as.Date(data_m$yyyymm, format = "%m/%d/%Y")
  data_m[,2:10] = scale( log( data_m[,2:10] / 100 + 1 ) )
  
  
  Y = ( as.matrix( data_m[,2:10] ) )
  data_long = melt((data_m[,1:10]), id.vars = "yyyymm") 
  names(data_long)[2] = "Portfolio"
  
  p = ggplot(data_long,                          
         aes(x = yyyymm,
             y = value,
             col = Portfolio)) +
    geom_line() + xlab("Time") + 
    ylab("Returns ") + 
    theme_classic()+ theme(text = element_text(size = 20, family = "Times"))
  print(p)
}

p = dim(Y)[2]; Tmax = dim(Y)[1]

#model hyperparam
w = (sum(selectFast(Y)$EW$G)) / (4 * p); if (w == 0){w = 0.1} #Bernoulli param 1
z = 0.1 #Bernoulli param 2
delta = 3.0; D = diag(1, p) #G-Wishart params
l = 2 + p #minimum-span constrain
p0 = 0.2 #geometric parameter (sparsity of change points)

#algorithm parameters
N = 20; M = 5 #particles and mutation steps

qdprime = 1/3  #MH proposal probabilities
qb = 1/6 #MH proposal probabilities
qd = 1/6 #MH proposal probabilities
lambda = 0.2 #param of local move

if (scenario > 0){ #simulation studies
  burnin = 2000 #burnin-period
  totiter = 5000 #tot-iterations (including burnin)
}else{
  burnin = 2000 #burnin-period
  totiter = 32000 #tot-iterations (including burnin)
}

#if print_change_point == TRUE,
#print the number of the iteration and change points 
#at every outer iteration:
print_change_point = TRUE 

thres = 1e-8  #threshod for Laplace approx (not used by default)

row = matrix(rep(c(1:p), p), ncol = p)
col = matrix(rep(c(1:p), p), ncol = p, byrow = TRUE)

#initialization
set.seed(0)
prior = compute_prior(l, Tmax)

if(scenario %in% c(1,2)){c = seq(51,Tmax,by=50)}else{c = NULL}

freep = detect_free_points(c, l, Tmax)
S = temperatures_tuning(Y, c, N, delta, D, w, z)
#S = matrix(rep(seq(0,1,by=0.2),length(c)+1), nrow = length(c)+1, byrow= TRUE)

IC = particle_filter(Y, c, N, S, delta, D, w, z, M)
log_lik = IC$lik
G = IC$G

if (length(c)>0){
  log_prior = length(c)*log(p0) + (length(c)>0) * log(prior[length(c)]) 
}else{
  log_prior = 0
}

for (i in 1:totiter){
  out = sample_ChangePoints(Y, c, freep, log_lik, log_prior, G, l, Tmax,
                            qb, qd, qdprime, lambda, S_old = S, reest = rbinom(1,1,0),
                            w, z, prior, N, M, p0)
  c = out$c
  freep = out$freep 
  log_lik = out$log_lik
  log_prior = out$log_prior
  S = out$temp
  if(print_change_point == TRUE){print(c(i, c))}
  if (i == 1){
    cat( paste( magenta( bold( "**Program has started**","\n" ) ) ) )
    cat(paste("|"))
    for (j in 1:(48)){
      cat(paste("-"))
    }
    cat(paste("|","\n"))
    write.table(t(seq(1,length(prior))), file = "cp_results.csv",append = TRUE, row.names = FALSE,
                col.names = FALSE)
  }
  if(floor((i-1)/(totiter/50))==((i-1)/(totiter/50))){cat(paste("*"))}
  if (length(c)!=0){
    write.table(t(c), file = "cp_results.csv",append = TRUE, row.names = FALSE,
                col.names = FALSE)
  }else{
    write.table(NA, file = "cp_results.csv",append = TRUE, row.names = FALSE,
                col.names = FALSE)
  }
    
}
cat( paste( magenta( bold( "**First run completed - change points estimated**","\n" ) ) ) )
windowsFonts(Times=windowsFont("TT Times New Roman"))


#plot the posterior of the number of change points
cp = read.csv("cp_results.csv", sep = " ")
cp = cp[seq(dim(cp)[1]- (totiter - burnin),dim(cp)[1],10),]
howmany = apply(cp, 1 , function(x) sum(!is.na(x)))

table = table(seq(0, length(prior)))
table[0:(length(prior) + 1)] = 0
table[sort(unique(howmany+1))] = table(howmany) / (dim(cp)[1])
t = as.data.frame(t(table))
colnames(t) = c("A", "Number_of_change_points", "Posterior_probability")

nofcp<-ggplot(data=t, aes(x=Number_of_change_points, y=Posterior_probability)) +
  geom_bar(stat="identity", color = mycol[1], fill = mycol[1]) + scale_y_continuous(limit = c(0, 1)) +
 xlab("Number of Change Points") + 
  ylab("Posterior probability ")
nofcp + theme_hc(base_size = 30) + theme(text = element_text(family = "serif")) 

#plot the marginal of each change point
where = table(as.vector(cp[!is.na(cp)]))
table = table(seq(1, Tmax))
table[1:Tmax] = 0 
table[sort(unique(as.vector(cp[!is.na(cp)])))] = where/  (dim(cp)[1])
t = as.data.frame(t(table))
colnames(t) = c("A", "Change_points", "Posterior_probability") 

library(extrafont)
loadfonts(device = "win")
labelx = rep("",Tmax)
if (scenario > 0){
  labelx[table > 0] = (seq(1,Tmax)[table > 0])
  labelx[where==max(where)] = (seq(1,Tmax)[where==max(where)])
}else{
date_labels = unique(data["date"])
labelx[seq(1, Tmax, 10)] =  as.character(date_labels[seq(1, Tmax, 10), ])
labelx[Tmax] = as.character(date_labels[Tmax, ])
}
nofcp<-ggplot(data=t, aes(x=Change_points, y=Posterior_probability)) +
  geom_bar(stat="identity",color = mycol[2], fill = mycol[1]) 
nofcp  + theme_hc(base_size = 30)+ scale_colour_hc()+ scale_x_discrete(label = labelx) + 
  theme(axis.text.x=element_text(size=rel(0.7), angle=90, vjust=-0.001, hjust = 0.5), text = element_text(family = "serif")) +
  xlab("Change Points") + 
  ylab("Posterior probability ") 

#print the joint of all change points
library(data.table)
setDT(cp)[,list(Count=.N),names(cp)]

################################################################################################################


#Estimate the graphs 
if(scenario == 0){
  #c_est must be the MAP of the first run
  c_est = c(61,79)
  temp = 0 
  set.seed(0)
  S_last = temperatures_tuning(Y, c_est, 100, delta, D, w, z, precision = "high")
  est = particle_filter(Y, c_est, 1000, S_last, delta, D, w, z, M = 20)
  temp = est$G

  temp1 = temp[1,,]
  temp2 = temp[2,,]
  temp3 = temp[3,,]
  colnames(temp1) = ptf; colnames(temp2) = ptf; colnames(temp3) = ptf
  rownames(temp1) = ptf; rownames(temp2) = ptf; rownames(temp3) = ptf
  
  heatmap(temp1,Rowv=NA,Colv=NA,col=paste("gray",1:99,sep=""))
  #determine threshold to have specificity equal 0.95
  #t = 0.5 ---> 0.84
  1 - (sum(1 - temp1[temp1>0.5],1 - temp2[temp2>0.5],1 - temp3[temp3>0.5]) / 
    length(c(temp1[temp1>0.5],temp2[temp2>0.5],temp3[temp3>0.5])))
  #t = 0.7 ---> 0.92
  1 - sum(1 - temp1[temp1>0.7],1 - temp2[temp2>0.7],1 - temp3[temp3>0.7]) / 
    length(c(temp1[temp1>0.7],temp2[temp2>0.7],temp3[temp3>0.7]))
  #t = 0.8 ---> 0.96
  1 - sum(1 - temp1[temp1>0.8],1 - temp2[temp2>0.8],1 - temp3[temp3>0.8]) / 
    length(c(temp1[temp1>0.8],temp2[temp2>0.8],temp3[temp3>0.8]))
  
  temp1 = (temp1 + t(temp1))
  temp2 = (temp2 + t(temp2))
  temp3 = (temp3 + t(temp3))
  write.table(temp1, file = "temp1.csv",append = TRUE, row.names = FALSE,
            col.names = FALSE)
  write.table(temp2, file = "temp2.csv",append = TRUE, row.names = FALSE,
            col.names = FALSE)
  write.table(temp3, file = "temp3.csv",append = TRUE, row.names = FALSE,
            col.names = FALSE)
  temp1[temp1<=0.80] = 0 
  temp2[temp2<=0.80] = 0 
  temp3[temp3<=0.80] = 0 

  label1 = round(as.vector(temp1),2)
  label2 = round(as.vector(temp2),2)
  label3 = round(as.vector(temp3),2)
   ptf = c("NoDur", "Durbl", "Manuf", "Enrgy",
            "HiTec", "Telcm", "Shops", 
            "Hlth", "Utils")
  #label = round(temp[row<col],2)
  par(mfrow = c(1,3))
  plot(graph_from_adjacency_matrix(temp1, weighted = TRUE, mode = c("undirected")),
       edge.label = label1[label1>0], edge.width=3, edge.color="green",
       vertex.color=mycol[11], vertex.size=20,
       vertex.label = ptf,
       vertex.label.font=2, vertex.label.color="black", vertex.label.cex=1.5,
       edge.label.font=2, edge.label.cex=1.5,
       main = "Graph from December 31, 2018 \n to February 21, 2020")
  
  plot(graph_from_adjacency_matrix(temp2, weighted = TRUE, mode = c("undirected")),
       edge.label = label2[label2>0], edge.width=3, edge.color="green",
       vertex.color=mycol[11], vertex.size=20,
       vertex.label = ptf,
       vertex.label.font=2, vertex.label.color="black", vertex.label.cex=1.5,
       edge.label.font=2, edge.label.cex=1.5,
       main = "Graph from February 24, 2020 \n to June 26, 2020")
  
  plot(graph_from_adjacency_matrix(temp3, weighted = TRUE, mode = c("undirected")),
       edge.label = label3[label3>0], edge.width=3, edge.color="green",
       vertex.color=mycol[11], vertex.size=20,
       vertex.label = ptf,
       vertex.label.font=2, vertex.label.color="black", vertex.label.cex=1.5,
       edge.label.font=2, edge.label.cex=1.5,
       main = "Graph from June 29, 2020 \n to December 31, 2021")
  
  
  set.seed(0)
  #Var-cov estimates 
  temp1[temp1>0] = 1
  temp1[row>col] = 0 
  G = temp1
  Y1 = Y[1:(c_est[1] - 1),]
  post1 = rgwish( n = 1000, adj = temp1, b = (c_est[1] - 1 + delta),
                  D = D + t(Y1)%*%Y1 , threshold = 1e-8 )
  Sigma1 = solve(apply(post1,c(1,2),mean))
  
  temp2[temp2>0] = 1
  temp2[row>col] = 0 
  G = temp2
  Y2 = Y[c_est[1]:(c_est[2]-1),]
  post2 = rgwish( n = 1000, adj = temp2, b = (c_est[2] - c_est[1] + delta),
                  D = D + t(Y2)%*%Y2 , threshold = 1e-8 )
  Sigma2 = solve(apply(post2,c(1,2),mean))
  
  
  temp3[temp3>0] = 1
  temp3[row>col] = 0 
  G = temp3
  Y3 = Y[c_est[2]:Tmax,]
  post3 = rgwish( n = 1000, adj = temp3, b = (Tmax - c_est[2] + 1 + delta),
                  D = D + t(Y3)%*%Y3 , threshold = 1e-8 )
  Sigma3 = solve(apply(post3,c(1,2),mean))
  
  
  colnames(Sigma1) = ptf; colnames(Sigma2) = ptf; colnames(Sigma3) = ptf
  rownames(Sigma1) = ptf; rownames(Sigma2) = ptf; rownames(Sigma3) = ptf
  
  par(mfrow=c(1,3))
  library(reshape2)
  temp = Sigma1 * ((diag(var)) %*% t(diag(var)))^(1/2)
  corr1 = temp / ((diag(temp)) %*% t(diag(temp)))^(1/2)
  temp = temp[,9:1]
  melted_Sigma1 = melt(temp)
  ggplot(data = melted_Sigma1, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + scale_fill_continuous(low = "white", 
                                       high = "blue",
                                       limits = c(0,0.009)) + 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position="none")+
    geom_text(aes(Var1, Var2, label = round(value * 10000, 2)), color = "black", size = 4)
  
  temp = Sigma2 * ((diag(var)) %*% t(diag(var)))^(1/2)
  temp = temp[,9:1]
  corr2 = temp / ((diag(temp)) %*% t(diag(temp)))^(1/2)
  melted_Sigma2 = melt(temp)
  ggplot(data = melted_Sigma2, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + scale_fill_continuous(low = "white", 
                                        high = "blue",
                                        limits = c(0,0.009))+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position="none")+
    geom_text(aes(Var1, Var2, label = round(value * 10000, 2)), color = "black", size = 4)
  
  temp = Sigma3 * ((diag(var)) %*% t(diag(var)))^(1/2)
  temp = temp[,9:1]
  corr3 = temp / ((diag(temp)) %*% t(diag(temp)))^(1/2)
  melted_Sigma3 = melt(temp)
  ggplot(data = melted_Sigma3, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + scale_fill_continuous(low = "white", 
                                        high = "blue",
                                        limits = c(0,0.009))+ 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position="none")+
    geom_text(aes(Var1, Var2, label = round(value * 10000, 2)), color = "black", size = 4)
  
  corr1[row==col] = NA
  corr1 = corr1[,9:1]
  melted_corr1 = melt(corr1)
  ggplot(data = melted_corr1, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + scale_fill_continuous(low = "white", 
                                        high = "blue",
                                        limits = c(0,1)) + 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position="none")+
    geom_text(aes(Var1, Var2, label = round(value, 2)), color = "black", size = 4)
  
  corr2[row==col] = NA
  corr2 = corr2[,9:1]
  melted_corr2 = melt(corr2)
  ggplot(data = melted_corr2, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + scale_fill_continuous(low = "white", 
                                        high = "blue",
                                        limits = c(0,1)) + 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position="none")+
    geom_text(aes(Var1, Var2, label = round(value, 2)), color = "black", size = 4)
  
  corr3[row==col] = NA
  corr3 = corr3[,9:1]
  melted_corr3 = melt(corr3)
  ggplot(data = melted_corr3, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + scale_fill_continuous(low = "white", 
                                        high = "blue",
                                        limits = c(0,1)) + 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position="none")+
    geom_text(aes(Var1, Var2, label = round(value, 2)), color = "black", size = 4)
  
}else if(scenario > 0){
  thres = seq(0,100)/100 
  TN = rep(0, length(thres))
  TP = rep(0, length(thres))
  FP = rep(0, length(thres))
  FN = rep(0, length(thres))
  FPR = rep(0, length(thres))
  TPR = rep(0, length(thres))
  Performance12 = data.frame(thres = thres)
  Performance_scenario3 = data.frame(thres = thres)
  
  if(scenario == 1){
    for (replica in 1:10){
      print(replica)
      #Scenario1 independence and no change points
      #simulate data
      p = 10; Tmax = 200
      G = array(0, dim = c(1, p, p))
      temp = G[1, , ]
      temp = (temp + t(temp))
      delta <- 3.0; D <- diag(1, p); 
      Omega = rgwish(1, G[1,,], delta, D, threshold=1e-8)
      Omega = round(Omega, 3); Sigma <- solve(Omega)
      
      set.seed(replica) 
      Y = rmvnorm(Tmax, rep(0,p), Sigma)
      Y = scale(Y)
      
      Gtruth = temp
      
      #hyperparam 
      w <- (sum(selectFast(Y)$EW$G)) / (4 * p); if (w == 0){w = 0.1}
      z <- 0.1
      delta <- 3.0; D <- diag(1, p) 
      l <- 2 * p;
      p = dim(Y)[2]; Tmax = dim(Y)[1]
      
      #conditional particle filter
      c_est = NULL
      temp = 0 
      set.seed(0)
      S_last = temperatures_tuning(Y, c_est, 100, delta, D, w, z)
      for (iter in 1:10){
        est = particle_filter(Y, c_est, 100, S_last, delta, D, w, z, M = 10)
        temp = temp + est$G 
      }
      temp = temp/(10)
      
      Gest = temp[1,,]; Gest = Gest + t(Gest)
      
      i = 1
      for (tt in thres){
        print(tt)
        Gtemp = Gest
        Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0; 
        
        TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth)) - p) / 2 
        FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth)) / 2
        i = i + 1
      }
    }
    TN[1] = 0;
    i = 1
    for (tt in thres){
      FPR[i] = FP[i] / (FP[i] + TN[i]) #fall-out: among negative who turns out positive
      i = i + 1
    }
    Performance12$scenario1 = FPR
  }
  if(scenario == 2){
    for (replica in 1:10){
      print(replica)
      #Scenario2 dependence and no change points
      #simulate data
      set.seed(1)
      p = 10; Tmax = 200
      G = array(0,dim=c(1,p,p))
      G[1,2,5] = 1; 
      G[1,4,5] = 1; 
      G[1,3,7] = 1; 
      G[1,4,9] = 1; 
      G[1,1,10] = 1; 
      G[1,4,10] = 1; 
      G[1,6,10] = 1; 
      G[1,7,10] = 1; 
      G[1,8,10] = 1; 
      
      temp = G[1, , ]
      temp = (temp + t(temp))
      
      delta <- 3.0; D <- diag(1, p) 
      Omega = rgwish(1, G[1,,], delta, D, threshold=1e-8)
      Omega = round(Omega, 3); Sigma <- solve(Omega)
      
      set.seed(replica) 
      Y = rmvnorm(Tmax, rep(0,p), Sigma)
      Y = scale(Y)
      
      Gtruth = temp
      
      #hyperparam 
      w <- (sum(selectFast(Y)$EW$G)) / (4 * p); if (w == 0){w = 0.1}
      z <- 0.1
      delta <- 3.0; D <- diag(1, p) 
      l <- 2 * p;
      p = dim(Y)[2]; Tmax = dim(Y)[1]
      
      #conditional particle filter
      c_est = NULL
      temp = 0 
      set.seed(0)
      S_last = temperatures_tuning(Y, c_est, 100, delta, D, w, z)
      for (iter in 1:10){
        est = particle_filter(Y, c_est, 100, S_last, delta, D, w, z, M = 10)
        temp = temp + est$G 
      }
      temp = temp/(10)
      
      Gest = temp[1,,]; Gest = Gest + t(Gest)
      
      i = 1
      for (tt in thres){
        print(tt)
        Gtemp = Gest
        Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0
        
        TP[i] = TP[i] + sum(Gtemp * Gtruth) / 2  #0 in independence 
        FN[i] = FN[i] + sum((1 - Gtemp) * Gtruth) / 2 #0 in independence
        TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth)) - p) / 2 
        FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth)) / 2
        i = i + 1
      }
    }
    TN[1] = 0
    i = 1
    for (tt in thres){
      TPR[i] = TP[i] / (TP[i] + FN[i])
      FPR[i] = FP[i] / (FP[i] + TN[i]) #fall-out: among negative who turns out positive
      i = i + 1
    }
    Performance12$scenario2 = FPR
    ROC2 = cbind(thres, TPR, FPR)
  write.table(Performance12, file = "Performance12_FPR.csv", row.names = FALSE,
              col.names = TRUE)
  write.table(ROC2, file = "ROC2.csv", append = TRUE, row.names = FALSE,
              col.names = TRUE)
  
  melt_performace12 = data.frame( thres = rep(Performance12$thres, 2), 
                                  FPR = c(Performance12$scenario1,Performance12$scenario2), 
                                  scenario = factor(c(rep(1,length(thres)), rep(2,length(thres)))) )
  p<-ggplot(melt_performace12, aes(x=thres, y=FPR, group = scenario)) +
    geom_line(aes(color = scenario))+
    geom_point(aes(color = scenario))
  p + scale_color_brewer(palette="Dark2") + theme_minimal() + 
    theme(text = element_text(size = 20, family = "Times")) +
    labs(x = "PPI Threshold", colour = "Scenario")
  auc = sum(diff(FPR[(i-1):1]) * TPR[(i-2):1])
  #FPR plot 
  temp = data.frame("TPR" = c(TPR,0), "FPR" = c(FPR,0))
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 4, 0, 0)) 
  plot(temp$FPR, temp$TPR, type = "l",
       xlab = "FPR (fall-out)", ylab = "TPR (sensitivity)", 
       cex.main = 4, cex.lab = 2.4,  
       cex.axis = 1.8)
  par(xpd=FALSE)
  abline(a=0, b=1, col=c("red"),lwd=3, lty=2)
  text(0.5, 0.1, "AUC = " , cex = 1.7)
  text(0.7, 0.1, round(auc, 5), cex = 1.7 )
  }
  if (scenario==3){
    #B. select replica for simulated data
    #replica = 1 #1,2,3,4,5,6,7,8,9,10
    MAPest = c(69,70,70,72,70,73,70,69,72,70)
    
    #algorithm parameters
    #N = 100; M = 5 #particles and mutation steps
    qdprime = 1/3  #MH proposal probabilities
    qb = 1/6 #MH proposal probabilities
    qd = 1/6 #MH proposal probabilities
    
    #burnin = 2000 #burnin-period
    #totiter = 5000 #tot-iterations (including burnin)
    
    
    for (replica in 1:10){
      p = 10; Tmax = 200 #dependence and one change points
      delta <- 3.0; D <- diag(1, p)
      set.seed(0)
      
      sim_data = simulate_data(scenario, replica)
      
      Y = sim_data$Y
      
      Gtruth = sim_data$G
      #Gtruth[1,,] first graph
      
      #hyperparam 
      w = 1
      z <- 0.1
      delta <- 3.0; D <- diag(1, p) 
      l <- 2 * p;
      Y = scale(Y)
      p = dim(Y)[2]; Tmax = dim(Y)[1]
      set.seed(0)
      
      #conditional particle filter
      c_est = MAPest[replica]
      temp = 0 
      set.seed(0)
      S_last = temperatures_tuning(Y, c_est, 100, delta, D, w, z)
      for (iter in 1:10){
        est = particle_filter(Y, c_est, 100, S_last, delta, D, w, z, M = 10)
        temp = temp + est$G  #est$G[1,,]
      }
      temp = temp/(10)
      
      G1 = temp[1,,]; G1 = G1  + t( G1 )
      G2 = temp[2,,]; G2 = G2  + t( G2 )
      i = 1
      for (tt in thres){
        print(tt)
        #1
        Gtemp = G1
        Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0
        
        TP[i] = TP[i] + sum(Gtemp * Gtruth[1,,]) / 2  #0 in independence 
        FN[i] = FN[i] + sum((1 - Gtemp) * Gtruth[1,,]) / 2 #0 in independence
        TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth[1,,])) - p) / 2 
        FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth[1,,])) / 2
        #2
        Gtemp = G2
        Gtemp[Gtemp >= tt] = 1; Gtemp[Gtemp < tt] = 0
        
        TP[i] = TP[i] + sum(Gtemp * Gtruth[2,,]) / 2  #0 in independence 
        FN[i] = FN[i] + sum((1 - Gtemp) * Gtruth[2,,]) / 2 #0 in independence
        TN[i] = TN[i] + (sum((1 - Gtemp) * (1 - Gtruth[2,,])) - p) / 2 
        FP[i] = FP[i] + sum(Gtemp *  (1 - Gtruth[2,,])) / 2
        i = i + 1
      }
    }
    TN[1] = 0
    i = 1
    for (tt in thres){
      TPR[i] = TP[i] / (TP[i] + FN[i])
      FPR[i] = FP[i] / (FP[i] + TN[i]) #fall-out: among negative who turns out positive
      i = i + 1
    }
    Performance_scenario3 = FPR
    ROC3 = cbind(thres, TPR, FPR)
    
    write.table(Performance3, file = "Performance3_FPR.csv", row.names = FALSE,
                col.names = TRUE)
    write.table(ROC3, file = "ROC3.csv", append = TRUE, row.names = FALSE,
                col.names = TRUE)
    
    performace3 = data.frame( thres = thres, 
                              FPR = Performance_scenario3) 
    
    p<-ggplot(performace3, aes(x=thres, y=FPR)) +
      geom_line()+
      geom_point()
    p + scale_color_brewer(palette="Dark2") + theme_minimal() + 
      theme(text = element_text(size = 20, family = "Times")) +
      labs(x = "PPI Threshold")
    
    
    par(mfrow = c(1, 1))
    auc = sum(diff(FPR[(i-1):1]) * TPR[(i-2):1])
    #FPR plot 
    temp = data.frame("TPR" = c(TPR,0), "FPR" = c(FPR,0))
    mar.default <- c(5,4,4,2) + 0.1
    par(mar = mar.default + c(0, 4, 0, 0)) 
    plot(temp$FPR, temp$TPR, type = "l",
         xlab = "FPR (fall-out)", ylab = "TPR (sensitivity)", 
         cex.main = 4, cex.lab = 2.4,  
         cex.axis = 1.8)
    par(xpd=FALSE)
    abline(a=0, b=1, col=c("red"),lwd=3, lty=2)
    text(0.5, 0.1, "AUC = " , cex = 1.7)
    text(0.7, 0.1, round(auc, 5), cex = 1.7 )
  }
}