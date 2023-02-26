#DGGM.r, this file contains all functions used to estimate the dynamic Gaussian graphical model 
#in Franzolini et al. 2022, 
#"Change point detection in dynamic Gaussian graphical models: the impact of COVID-19 pandemic on the US stock market"

#functions description
# - compute_prior: returns the inverse of the normalizing constant of uniform prior for each possible number of change point
# - detect_free_points: given a sequence of change points find available change points
# - detect_free_points_btw2: given a sequence of change points find available change points between two change points
# - ess: compute the effective sample size
# - lik_block: compute the likelihood of a block given the corresponding graph
# - mutate_all_G: perform a MH step for the whole sequence of graphs and all particles to sample for their posterior (not used)
# - mutate_G: perform MH step for one graph and all particles to sample for their temperated-posterior
# - particle_filter: inner component of the algorithm, compute the marg lik and sample the graphs
# - sample_ChangePoints: outer component, MH for change points
# - sim.data: simulate data from the model (not used)
# - simulate_data: simulate data from scenarios described in the Franzolini et al. (2022)
# - sim.G: sample a graph at time 0 from the prior
# - sim.N.G: sample N graphs at time 0 from the prior
# - sim.GG: sample a graph from the prior given the previous graph
# - sim.N.GG: sample N graphs from the prior given N previous graphs
# - temperatures_tuning: compute temperatures adaptively 

compute_prior <- function(l, Tmax){
  prior = rep(0, floor(Tmax/l))
  combination = rep(0, Tmax)
  s_last = seq(l, Tmax - l)#
  combination[s_last] = 1
  prior[1] = 1 / sum(combination)
  for (k in 2:floor(Tmax/l)){
    combination_old = combination
    combination = rep(0, Tmax)
    for (c_last in ((k*l):(Tmax-l))){
      combination[c_last] = sum(combination_old[1:(c_last-l)])
    }
    prior[k] = 1 / sum(combination)
  }
  return(prior)
}


detect_free_points <- function(c, l, Tmax){
  #works only for l integer(!)
  if(length(c)>=1){
    taken = matrix(rep(c,2*l-1)*(1-is.null(c)),
                   nrow=(2*l-1)*(1-is.null(c)), byrow=TRUE) +
      matrix(rep(c((-l+1):(l-1)), length(c)), ncol=length(c))
  }else{
    taken = NULL
  }
  return(setdiff(c((l+1):(Tmax-l)), taken))
}

detect_free_points_btw2 <- function(c, l, Tmax, i){
  #works only for l integer(!)
  cold = c; c = c[-i]
  if(length(c)>=1){
    taken = matrix(rep(c,2*l-1)*(1-is.null(c)),
                   nrow=(2*l-1)*(1-is.null(c)), byrow=TRUE) +
      matrix(rep(c((-l+1):(l-1)), length(c)), ncol=length(c))
  }else{
    taken = NULL
  }
  c = c(l+1, cold, Tmax-l)
  return(setdiff(c(c[i]:c[i+2]), taken))
}

ess <-function(w, phi1, phi2){
  ess = sum(w^(phi2 - phi1))^2 / sum(w^(2*phi2 - 2*phi1))
  return(ess)
}

lik_block <- function(Y, G, delta, D, s2, s1, log = FALSE){
  S = t(Y)%*%Y 
  out = gnorm( G, b = delta + dim(Y)[1], D = D + S, iter = 100) -   gnorm( G, b = delta, D = D, iter = 100) -
  (dim(Y)[1] * dim(Y)[2] / 2) * log(2*pi)
  if(log){
    return((s2-s1)*out)
  }else{
    return(exp(out)^(s2-s1))
  }
}


mutate_all_G <-function(G, Y, z, c, delta, D){
  p = dim(Y)[2]
  N = dim(G)[4]
  G_new = G
  for (j in c( 1:(length(c)-1) ) ){
    G_new[j,,,] = sim.N.GG(N, G[j,,,], p, 0.5)
  }
  for (n in 1:N){
    lik_new = 1; lik_old = 1
    for (j in c( 1:(length(c)-1) ) ){
      lik_new = lik_new * lik_block(Y[c[j]:(c[j+1]-1),], G_new[j, , ,n], delta, D, 1, 0)
      lik_old = lik_old * lik_block(Y[c[j]:(c[j+1]-1),], G[j, , ,n], delta, D, 1, 0)
    }
    if (lik_old == lik_new){diff = 0}else{diff = lik_new - lik_old }
    if (runif(1) < diff){
      G[,,,n] = G_new[, , ,n]
    }
  }
  return(G)
}

mutate_G <-function(G, Y, z, c, delta, D, phi, previous_graph =  array(0, dim=c(1,p,p,N))){
  p = dim(Y)[2]
  N = dim(G)[4]
  G_new = G
  G_new[1, , , ] = sim.N.GG(N, G[1, , ,], p, 0.1)
  row = matrix(rep(c(1:p), p), ncol = p)
  col = matrix(rep(c(1:p), p), ncol = p, byrow = TRUE)
  for (n in 1:N){
    temp1 = G_new[1, , , n]; temp2 = G[1, , , n]
    temp1_prev = previous_graph[1,,,n]
    lik_new = lik_block(Y, G_new[1, , , n], delta, D, phi, 0, log = TRUE)
    lik_old = lik_block(Y, G[1, , , n], delta, D, phi, 0, log = TRUE)
    prior_new = sum(dbinom(temp1[col>row],1, abs( temp1_prev[col>row]-(2 * z/ (p - 1.0))), log = TRUE ))
    prior_old = sum(dbinom(temp2[col>row],1, abs( temp1_prev[col>row]-(2 * z/ (p - 1.0))), log = TRUE ))
    if (lik_old  + prior_old == lik_new + prior_new){diff = 0}else{diff = lik_new - lik_old + prior_new - prior_old }
    if (log(runif(1)) < diff ){
      G[1 , , ,n] = temp1
    }
  }
  return(G)
}

particle_filter <-function(Y, c, N, S, delta, D, w, z, M = 1, M2 = 0){
  Tmax = dim(Y)[1]
  p = dim(Y)[2]
  c = c(1, c, Tmax+1)
  k = length(c)-1
  
  G = array(rep(NA, k * p * p * N), dim=c(k, p, p, N))
  logww = array(rep(NA, N ), dim=c(N))
  weight = array(rep(NA, N * k ), dim=c(N, k))
  lik_temp = array(rep(NA, N * k * dim(S)[2] ), dim=c(N, k, dim(S)[2]))
  norm_weight = array(rep(NA, N * k ), dim=c(N, k))
  
  G[1, , , ] = sim.N.G(N, p, w)
  if((length(S[1,!is.na(S[1,])])-1)>0){temp_ss = 1:(length(S[1,!is.na(S[1,])])-1)}else{temp_ss = NULL}
  if(is.null(temp_ss)){
    for (n in (1:N)){
      logww[n] = lik_block(Y[c[1]:(c[2]-1),], G[1, , ,n], delta, D, 1, 0 , log = TRUE)
    }
    lik_temp[, 1, 1] = exp(logww)
    if (length(logww[logww!= -Inf]) == 0){ 
      logww[1:N] = 0
    }
    weight[,1] = exp(logww)
    if ((sum(weight[,1])) == 0){ 
      logww = logww - min(logww[logww!= -Inf])
      weight[,1] = exp (logww)
    }else if ((sum(weight[,1])) == Inf){ 
      logww = logww - max(logww)
      weight[,1] = exp (logww)
    }
    norm_weight[,1] = weight[,1] / sum(weight[,1])
    which = sample(N, prob = norm_weight[,1], replace = TRUE)
    G[1, , , ] = G[1, , , which]
    if ( (sum(norm_weight[,1])^2 / sum(norm_weight[,1]^(2))) < (0.5 * N) ){
      if (M>0){
        for(mcmc in (1:M)){
          G[1, , , ] = mutate_G(G[1, , , ,drop = FALSE], Y[c[1]:(c[2]-1), ,drop = FALSE], w, c, delta, D, 1)
        }}
    }
  }else{
    for (s in temp_ss){
      for (n in (1:N)){
        logww[n] = lik_block(Y[c[1]:(c[2]-1),], G[1, , ,n], delta, D, S[1, s+1], S[1, s] , log = TRUE)
      }
      lik_temp[, 1, s] = exp(logww)
      if (length(logww[logww!= -Inf]) == 0){ 
        logww[1:N] = 0
      }
      weight[,1] = exp(logww) 
      if ((sum(weight[,1])) == 0){ 
        logww = logww - min(logww[logww!= -Inf])
        weight[,1] = exp (logww)
      }else if ((sum(weight[,1])) == Inf){ 
        logww = logww - max(logww)
        weight[,1] = exp (logww)
      }
      norm_weight[,1] = weight[,1] / sum(weight[,1])
      which = sample(N, prob = norm_weight[,1], replace = TRUE)
      G[1, , , ] = G[1, , , which]
      if ( sum(norm_weight[,1])^2 / sum(norm_weight[,1]^(2)) < (0.5 * N) ){
        if (M>0){
          #print(s)
          for(mcmc in (1:M)){
            G[1, , , ] = mutate_G(G[1, , , ,drop = FALSE], Y[c[1]:(c[2]-1), ,drop = FALSE], w, c, delta, D, S[1,s+1])
          }}
      }
    } }
  
  if (length(c)>2){
    for (j in c( 2:(length(c)-1) ) ){
      G[j, , , ] = sim.N.GG(N, G[j-1, , , ], p, z)   
      #logww = array(rep(NA, N ), dim=c(N))
      if((length(S[j,!is.na(S[j,])])-1)>0){temp_ss = 1:(length(S[j,!is.na(S[j,])])-1)}else{temp_ss = NULL}
      if(is.null(temp_ss)){
        for (n in (1:N)){
          logww[n] = lik_block(Y[c[j]:(c[j+1]-1),], G[j, , ,n], delta, D, 1, 0 , log = TRUE)
        }
        lik_temp[, j, 1] = exp(logww)
        if (length(logww[logww!= -Inf]) == 0){ 
          logww[1:N] = 0
        }
        weight[,j] = exp(logww) 
        if ((sum(weight[,j])) == 0){ 
          logww = logww - min(logww[logww!= -Inf])
          weight[, j]= exp (logww)
        }else if ((sum(weight[,j])) == Inf){ 
          logww = logww - max(logww)
          weight[,j] = exp (logww)
        }
        norm_weight[,j] = weight[,j] / sum(weight[,j])
        which = sample(N, prob=norm_weight[,j], replace = TRUE)
        G[j, , ,] = G[j, , , which]
        if ( sum(norm_weight[,j])^2 / sum(norm_weight[,j]^(2)) < (0.5 * N) ){
          if (M>0){
            for(mcmc in (1:M)){
              G[j, , , ] = mutate_G(G[j, , , ,drop = FALSE], Y[c[j]:(c[j+1]-1), ,drop = FALSE], z, c, delta, D, 1, previous_graph = G[j-1, , , ,drop = FALSE])
            }}
        }
      }else{
        for (s in temp_ss){
          for (n in (1:N)){
            logww[n] = lik_block(Y[c[j]:(c[j+1]-1),], G[j, , ,n], delta, D, S[j, s+1], S[j, s] , log = TRUE)
          }
          lik_temp[, j, s] = exp(logww)
          if (length(logww[logww!= -Inf]) == 0){ 
            logww[1:N] = 0
          }
          weight[,j] = exp(logww) 
          if ((sum(weight[,j])) == 0){ 
            logww = logww - min(logww[logww!= -Inf])
            weight[,j] = exp (logww)
          }else if ((sum(weight[,j])) == Inf){ 
            logww = logww - max(logww)
            weight[,j] = exp (logww)
          }
          norm_weight[,j] = weight[,j] / sum(weight[,j])
          which = sample(N, prob=norm_weight[,j], replace = TRUE)
          G[j, , ,] = G[j, , , which]
          if ( sum(norm_weight[,j])^2 / sum(norm_weight[,j]^(2)) < (0.5 * N) ){
            if (M>0){
              for(mcmc in (1:M)){
                G[j, , , ] = mutate_G(G[j, , , ,drop = FALSE], Y[c[j]:(c[j+1]-1), ,drop = FALSE], z, c, delta, D, S[j,s+1], previous_graph = G[j-1, , , ,drop = FALSE])
              }}
          }
        }}
    } }
  lik = sum(log(colMeans(lik_temp, na.rm = TRUE)), na.rm = TRUE)
  if(M2>0){for (mcmc in 1:M2){G = mutate_all_G(G, Y, z, c, delta, D)}}
  G = apply(G, c(1,2,3), mean)
  return(list(lik=lik, G = G))
}



sample_ChangePoints<- function (Y, c, freep, log_lik, log_prior, G_old, l, Tmax,
                                qb, qd, qdprime, lambda = 0.2, S_old = NULL, reest = FALSE , 
                                w, z, prior, N, M, p0 = 0.2){
  n = length(freep) #number of available points
  k = length(c) #number of change points
  #sample the BDM event (B=1, D=2, M=3)
  E = (k == 0) + 
    sample(c(2, 3), size = 1, prob = c(qdprime, 1-qdprime)) * (n==0) + 
    sample(c(1, 2, 3), size = 1, prob = c(qb, qd, 1-qb-qd)) * !(n==0 || k==0) 
  #sample the new vector
  if (E==1){
    #birth
    c_star = sample(freep, size = 1)
    cnew = append(c, c_star, after = length(which(c < c_star)))
    qnew = 1/n * ((k == 0) + qb * (k!=0))
  }else{
    #death
    i = sample(k, size = 1)
    cnew = c[-i]
    qnew = 1/k * ((n==0) * (qdprime * (E==2) + (1-qdprime)/2 * (E==3)) +
                    (n!=0) * (qd * (E==2) + (1-qd-qb)/2 * (E==3)))
    if(E==3){
      which = sample(2, 1)
      freetemploc = detect_free_points_btw2(c, l, Tmax, i) 
      ploc = exp(- lambda * abs(freetemploc - c[i]))
      freetempglob = detect_free_points(cnew, l, Tmax)
      if (which == 1){
        c_star = sample(freetemploc, size = 1, prob = ploc)
      }else{
        c_star = sample(freetempglob, size = 1)
      }
      cnew = append(cnew, c_star, after = length(which(cnew < c_star)))
      c_temp = c(l+1, c, Tmax-l)
      qnew = identical(c, cnew)  +  (!identical(c, cnew)) * 
        (qnew / length(freetempglob) + (c_star<c_temp[i+2])*(c_star>c_temp[i])*
        qnew * exp(- lambda * abs(c_star - c[i])) / sum(ploc))
    }
  }
  #compute qold
  #inverse event B=2, D=1, M=3
  freep_new = detect_free_points(cnew, l, Tmax)
  n_new = length(freep_new)
  k_new = length(cnew)
  if (E == 2){
    qold =  1/n_new * ( (k_new == 0) + qb * (k_new!=0))
  }else if (E==1){
    qold = 1/k_new * (n_new==0) * (qdprime) + (n_new!=0) * (qd)
  }else{
    if (identical(c, cnew)){ qold = 1 }else{
    qold =  1/k_new * (n_new==0) * (1 - qdprime) + (n_new!=0) * (1 - qd - qb)
    c_to_rem = setdiff(cnew,c); i = which(cnew == c_to_rem)
    c_star_new = setdiff(c, cnew)
    freetemploc_new = detect_free_points_btw2(cnew, l, Tmax, i) 
    ploc_new = exp(- lambda * abs(freetemploc_new - c_to_rem))
    freetempglob_new = detect_free_points(cnew[-i], l, Tmax)
    c_temp = c(l+1, cnew , Tmax-l)
    qold = (qold / length(freetempglob_new) + (c_star_new<c_temp[i+2])*(c_star_new>c_temp[i])*
         qold * exp(- lambda * abs(c_star_new - cnew[i])) / sum(ploc_new))
    }
  }
  #S = temperatures_tuning(Y, cnew, N, delta, D, w, z)
  S = matrix(rep(seq(0,1,by=0.2),length(cnew)+1), nrow = length(cnew)+1, byrow= TRUE)
  IC = particle_filter(Y, cnew, N, S, delta, D, w, z, M)
  log_lik_new = IC$lik
  if (length(cnew)>0){
    log_prior_new = length(cnew)*log(p0) + log(prior[length(cnew)]) 
  }else{
    log_prior_new = 0
  }
  if(reest){
    if(length(S_old)==0){S_old = temperatures_tuning(Y, c, N, delta, D, w, z)}
    ICold = particle_filter(Y, c, N, S_old, delta, D, w, z, M )
    log_lik = ICold$lik
    G_old = ICold$G
  }

  A =  log_lik_new + log_prior_new + log(qold) -
    log_lik - log_prior - log(qnew) 
  if((log_lik_new + log_prior_new + log(qold)) == (log_lik + log_prior + log(qnew) ) ){A = 0}
  if (log(runif(1))<=A){
    output = list("c" = cnew, "freep" = freep_new, "log_lik" = log_lik_new, 
                  "log_prior" = log_prior_new, "G" = IC$G, "temp" = S)
  }else{
    output = list("c" = c, "freep" = freep, "log_lik" = log_lik, 
                  "log_prior" = log_prior, "G" = G_old, "temp" = S_old)
  }
}


sim.data <- function(p, Tau, w, z, delta, D, thres, S=NULL, ind=TRUE) {
  Y = matrix(ncol=p, nrow=Tau)
  l = length(S) # may be 0
  G = array(dim=c(l+1,p,p))
  Omega = array(dim=c(l+1,p,p))
  if(!ind){G[1,,] <- sim.G(p, w)}else{G[1,,]=matrix(rep(0, p * p), nrow = p)}
  Omega[1,,] = rgwish(1, G[1,,], delta, D, threshold=thres)
  Omega[1,,]  = round(Omega[1,,] , 3); Sigma <- solve(Omega[1,,] )
  sL = 1
  if (l>0) sR = S[1]-1 else sR = Tau
  Y[sL:sR,] = rmvnorm(sR-sL+1, rep(0,p), Sigma)
  if (l>0){ 
    for (i in 1:l) {
      G[i+1,,] = sim.GG(G[i,,], p, z)
      Omega[i+1,,]  = rgwish(1, G[i+1,,], delta, D, threshold = thres)
      Omega[i+1,,] = round(Omega[i+1,,], 3); Sigma <- solve(Omega[i+1,,])
      sL = S[i]
      if (i<l) sR = S[i+1]-1 else sR = Tau 
      Y[sL:sR,] = rmvnorm(sR-sL+1, rep(0,p), Sigma)
    }  
  }
  return(list(Y = Y, c = S, G = G, Omega = Omega))
}


simulate_data <- function(scenario, replica, seed = 1, Tmax = 200, p = 10){
  if(scenario == 1){
    #Scenario1 independence and no change points
    #p = 10; Tmax = 200
    G = array(0, dim = c(1, p, p))
    
    temp = G[1, , ]
    temp = (temp + t(temp))
    par(mfrow=c(1,1))
    plot(graph_from_adjacency_matrix(temp, weighted = TRUE, mode = c("undirected")),
         edge.width=3, edge.color="black",
         vertex.color=mycol[11], vertex.size=20,
         vertex.label.font=2, vertex.label.color="black", vertex.label.cex=1.5,
         edge.label.font=2, edge.label.cex=1.5,
         main = "True graph from t=1 to t=200")
    
    delta <- 3.0; D <- diag(1, p); thres <- 1e-8 
    Omega = rgwish(1, G[1,,], delta, D, threshold=thres)
    Omega = round(Omega, 3); Sigma <- solve(Omega)
    
    set.seed(replica) 
    Y = rmvnorm(Tmax, rep(0,p), Sigma)
    Y = scale(Y)
    
  }else if(scenario == 2){ 
    set.seed(seed)
    #Scenario 2 dependence and no change points
    #p = 10; Tmax = 200
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
    par(mfrow=c(1,1))
    plot(graph_from_adjacency_matrix(temp, weighted = TRUE, mode = c("undirected")),
         edge.width=3, edge.color="black",
         vertex.color=mycol[11], vertex.size=20,
         vertex.label.font=2, vertex.label.color="black", vertex.label.cex=1.5,
         edge.label.font=2, edge.label.cex=1.5,
         main = "True graph from t=1 to t=200")
    
    
    delta <- 3.0; D <- diag(1, p); thres <- 1e-8 
    Omega = rgwish(1, G[1,,], delta, D, threshold=thres)
    Omega = round(Omega, 3); Sigma <- solve(Omega)
    
    set.seed(replica) 
    Y = rmvnorm(Tmax, rep(0,p), Sigma)
    Y = scale(Y)
    
  }else if(scenario == 3){
    set.seed(seed)
    #p = 10; Tmax = 200 #dependence and one change points
    #graphs dependency as in Molinari et al.(2020) and Peterson et al.(2015)
    c = 70
    
    Omega1 = array(0, dim=c(p,p))
    row = matrix(rep(c(1:p), p), ncol = p)
    col = matrix(rep(c(1:p), p), ncol = p, byrow = TRUE)
    
    Omega1[col == row + 1] = 0.5
    #Omega1[col == row + 2] = 0.4
    G = array(0, dim=c(2,p,p))
    temp = Omega1; temp[Omega1!=0] = 1; temp[row == col] = 0
    G[1, ,] = temp
    Omega1 = Omega1 + t(Omega1)
    Omega1[row == col] = 1
    
    Sigma1 <- solve(Omega1)
    
    Omega2 = Omega1
    Omega2[,] = 0 
    temp = Omega2
    upper1 = Omega1[col>row]
    upper2 = upper1
    
    #remove 5 edges
    active = upper1[upper1!=0]
    active[sample(length(active), 5)] = 0
    upper2[upper1!=0] = active
    
    #add 5 edges
    inactive = upper1[upper1==0]
    #inactive[sample(length(inactive), 5)] = sample(c(-.4, -.6, .4, .6), 5, replace = TRUE)
    inactive[sample(length(inactive), 5)] = 0.2
    upper2[upper1==0] = inactive
    
    Omega2[col>row] = upper2
    Omega2 = Omega2 + t(Omega2)
    Omega2[row == col] = Omega1[row == col]
    temp = Omega2; temp[round(Omega2, 2)!=0] = 1; temp[round(Omega2, 2)==0] = 0; temp [row == col] = 0
    temp[col<row] = 0
    G[2,,] = as.matrix(temp)
    
    #for (i in 1:p){Omega2[i,] = Omega2[i,] / sum(Omega2[i,] - 1)  }
    #Omega2[row == col] = Omega1[row == col]
    #Omega2 = (Omega2 + t(Omega2))
    Omega2 = as.matrix(nearPD(Omega2, )$mat)
    Sigma2 <- solve(Omega2)
    
    
    par(mfrow=c(1,2))
    temp = G[1, , ]
    temp = (temp + t(temp))
    plot(graph_from_adjacency_matrix(temp, weighted = TRUE, mode = c("undirected")),
         edge.width=3, edge.color="black",
         vertex.color=mycol[11], vertex.size=20,
         vertex.label.font=2, vertex.label.color="black", vertex.label.cex=1.5,
         edge.label.font=2, edge.label.cex=1.5,
         main = "True graph from t=1 to t=69")
    temp = G[2, , ]
    temp = (temp + t(temp))
    plot(graph_from_adjacency_matrix(temp, weighted = TRUE, mode = c("undirected")),
         edge.width=3, edge.color="black",
         vertex.color=mycol[11], vertex.size=20,
         vertex.label.font=2, vertex.label.color="black", vertex.label.cex=1.5,
         edge.label.font=2, edge.label.cex=1.5,
         main = "True graph from t=70 to t=200")
    
    
    set.seed(replica)
    Y = matrix(0, ncol = p, nrow = Tmax)
    Y[1:(c-1),] =  rmvnorm(c-1, rep(0,p), Sigma1)
    Y[c:Tmax,] =  rmvnorm(Tmax - c + 1, rep(0,p), Sigma2)
    Y = scale(Y)
    Sigma = list("Sigma1" = Sigma1, "Sigma2" = Sigma2)
    Omega = list("Omega1" = Omega1, "Omega2" = Omega2)
  }
  return(list("Y"= Y, "G" = G, "Omega" = Omega, "Sigma" = Sigma))
}

sim.G <- function(p, w) {
  G = matrix(rep(0, p*p), ncol=p, nrow=p)
  pB = 2.0 * w / (p - 1.0) # mean no of edges is w*p
  lv = round(p * (p - 1.0)/2.0)
  row = matrix(rep(c(1:p), p), ncol = p)
  col = matrix(rep(c(1:p), p), ncol = p, byrow = TRUE)
  G[col > row] = rbinom(lv, 1, pB)
  return(G)
}

sim.N.G <- function(N, p, w) {
  G = array(rep(0, N*p*p), dim = c( p, p, N))
  pB = 2.0 * w / (p - 1.0) # mean no of edges is w*p
  lv = round(p * (p - 1.0)/2.0)
  row = matrix(rep(c(1:p), p), ncol = p)
  col = matrix(rep(c(1:p), p), ncol = p, byrow = TRUE)
  G[col > row] = rbinom(N * lv, 1, pB)
  return(G)
}

sim.GG <- function(Gp, p, z) {
  G = matrix(rep(0,p*p), ncol=p, nrow=p)
  pB = 2.0*z / (p - 1.0) # mean no of edges selected is z*p
  lv = round(p*(p-1.0) / 2.0)
  row = matrix(rep(c(1:p), p), ncol = p)
  col = matrix(rep(c(1:p), p), ncol = p, byrow = TRUE)
  G[col > row] = abs(Gp[col > row] - rbinom(lv, 1, pB))
  return(G)
}

sim.N.GG <- function(N, Gp, p, z) {
  G = array(rep(0, p * p * N), dim = c(p, p, N))
  pB = 2.0*z / (p - 1.0) # mean no of edges selected is z*p
  lv = round(p*(p-1.0) / 2.0)
  row = matrix(rep(c(1:p), p), ncol = p)
  col = matrix(rep(c(1:p), p), ncol = p, byrow = TRUE)
  G[col > row] = abs(Gp[col > row] - rbinom(N * lv, 1, pB))
  return(G)
}

temperatures_tuning <-function(Y, c, N, delta, D, w, z, precision = "normal"){
  
  minESS = 0.5 * N
  if(precision == "high"){eps = 0.05}else{eps = 0.10}
  S = matrix(NA, nrow = length(c)+1, ncol = 1 / eps + 1)
  S[,1] = 0
  Tmax = dim(Y)[1]
  p = dim(Y)[2]
  c = c(1, c, Tmax+1)
  k = length(c)
  
  G = array(rep(NA, N * k * p * p), dim=c(k, p, p, N))
  logww = rep(NA, N )
  weight = array(rep(NA, N * k ), dim=c(N, k))
  norm_weight = array(rep(NA, N * k ), dim=c(N, k))
  
  G[1, , , ] = sim.N.G(N, p, w)
  phi = 0; ss = 2
  Sm = t(Y[c[1]:(c[2]-1),])%*%Y[c[1]:(c[2]-1),]
  while (phi < (1-eps)){
    for (n in (1:N)){
      logww[n] = gnorm( G[1, , , n] , b = delta - c[1] + c[2], D = D + Sm, iter = 100 ) -
        gnorm( G[1, , , n] , b = delta, D = D, iter = 100 )
    }
    if (length(logww[logww!= -Inf]) == 0){ 
      logww[1:N] = 0
      }
    ww = exp(logww)
    candidate_temp = seq(phi+eps,1,by=eps) 
    ESS = rep(0, length(candidate_temp))
    i = 1 
    for (phi_cand in candidate_temp){
      ESS[i] = ess(ww, phi, phi_cand)
      i = i + 1
    }
    ESS[is.na(ESS)] = 0
    if(sum(ESS>=minESS)==0){phi = candidate_temp[1] }else{phi = candidate_temp[ESS>=minESS][sum(ESS>=minESS)]}
    S[1, ss] = phi
    weight[,1] = ww ^ (phi - S[1, ss-1])
    if ((sum(weight[,1])) == 0){ 
      logww = logww - min(logww[logww!= -Inf])
      ww = exp (logww)
      weight[,1] = ww ^ (phi - S[1, ss-1])
    }else if ((sum(weight[,1])) == Inf){ 
      logww = logww - max(logww)
      ww = exp (logww)
      weight[,1] = ww ^ (phi - S[1, ss-1])
    }
      norm_weight[,1] = weight[,1] / sum(weight[,1])
      which = sample(N, prob = norm_weight[,1], replace = TRUE)
      G[1, , , ] = G[1, , , which]
    ss = ss + 1 
  }
  if(!isTRUE(all.equal(phi,1))){S[1, ss] = 1}
  
  if (length(c)>2){
    for (j in c( 2:(length(c)-1) ) ){
      G[j, , , ] = sim.N.GG(N, G[j-1, , , ], p, z)
      phi = 0; ss = 2
      Sm = t(Y[c[j]:(c[j+1]-1),])%*%Y[c[j]:(c[j+1]-1),]
      while (phi < (1-eps)){
        for (n in (1:N)){
          logww[n] = gnorm( G[j, , , n] , b = delta - c[j] + c[j+1], D = D + Sm, iter = 50 ) -
            gnorm( G[j, , , n] , b = delta, D = D, iter = 50 )
        }
        if (length(logww[logww!= -Inf]) == 0){ 
          logww[1:N] = 0
        }
        ww = exp(logww)
        candidate_temp = seq(phi+eps,1,by=eps) 
        ESS = rep(0, length(candidate_temp))
        i = 1
        for (phi_cand in candidate_temp){
          ESS[i] = ess(ww, phi, phi_cand)
          i = i + 1
        }
        ESS[is.na(ESS)] = 0
        if(sum(ESS>=minESS)==0){phi = candidate_temp[1] }else{phi = candidate_temp[ESS>=minESS][sum(ESS>=minESS)]}
        S[j, ss] = phi
          weight[,j] = ww ^ (phi - S[j, ss-1])
          if ((sum(weight[,j])) == 0){ 
            logww = logww - min(logww[logww!= -Inf])
            ww = exp (logww)
            weight[,j] = ww ^ (phi - S[j, ss-1])
          }else if ((sum(weight[,j])) == Inf){ 
            logww = logww - max(logww)
            ww = exp (logww)
            weight[,j] = ww ^ (phi - S[1, ss-1])
          }
          norm_weight[,j] = weight[,j] / sum(weight[,j])
          which = sample(N, prob = norm_weight[,j], replace = TRUE)
          G[j, , , ] = G[j, , , which]
        ss = ss + 1 
      }
      if(!isTRUE(all.equal(phi,1))){S[j, ss] = 1}
    }
  }
  S = S[,colSums(!is.na(S))>0, drop = FALSE]
  return(S)
}



