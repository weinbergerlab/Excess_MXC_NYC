model_string_pois2<-
  "model
{
  # Likelihood
  for( t in 1:n.dates ){
     
  for(d in 1:(D+1)){
  
  n[t,d] ~ dpois(lambda[t,d])
  n.train[t,d] ~dpois(baseline[t,d])

  log(lambda[t,d]) <- (int +     
            beta.logged2[d]*step(D-d) +
            sum.beta.logged2[t]*(1-step(D-d)) +
                    sin26[t]*delta1[1] +
                    cos26[t]*delta1[2] +
                    sin52[t]*delta1[3] +
                    cos52[t]*delta1[4] +
                    epsilon1[epiyr.index[t]] +
                    phi[t] +disp[t]
                    ) 
      #baseline is just standard harmonic reg, fit to train period     
    log(baseline[t,d])<- 
   (int +     
            beta.logged2[d]*step(D-d) +
            sum.beta.logged2[t]*(1-step(D-d)) +
                    sin26[t]*delta1[1] +
                    cos26[t]*delta1[2] +
                    sin52[t]*delta1[3] +
                    cos52[t]*delta1[4] +
                    epsilon1[epiyr.index[t]]
                    +disp[t]
                    ) 
  }       
  
  #Baseline has harmonics and iid prediction noise
    log(baseline.unobs[t])<-
        (int + sin26[t]*delta1[1] +
                    cos26[t]*delta1[2] +
                    sin52[t]*delta1[3] +
                    cos52[t]*delta1[4] +
                    epsilon1[epiyr.index[t]]+
                    disp2[t]
           )
  sum.n[t]  <-    sum(n[t,1:D])
  sum.lambda[t] <- sum(lambda[t,1:D])

  baseline.n[t] ~ dpois(baseline.unobs[t])

  excessN[t] <- sum.n[t] - baseline.n[t]
  }

    #Time Series Random Effect
    disp[1]  <-0
    disp2[1] <- 0

    for( t in 2:n.dates){
    disp[t] ~dnorm(0, tau.disp)
    disp2[t] ~ dnorm(0, tau.disp)
    }
    
    for(t in 1:(n.dates-D-1)){
    phi[t] <- 0
    }
  for(t in (n.dates-D):n.dates){
     phi[t] ~ dnorm((phi[t-1]), tau2.alpha)
     }

   ##  
  int~dnorm(0, 1e-4)

  ## Prior for beta
  beta.logged <- log(beta)
  beta.logged2 <- c(beta.logged,0)
  beta ~ ddirch(beta.priors)
  
  tau.disp <- 1/sd.disp^2
  sd.disp~dunif(0,100)
  
  for( t in 1:n.dates ){
    sum.beta[t] <- sum(beta[1:N.first.obs[t]])
    sum.beta.logged2[t] <- log(sum.beta[t])
  }

  # Prior for variance
  tau2.alpha ~ dgamma(alphat.shape.prior,alphat.rate.prior)
  
  for(i in 1:4){
    delta1[i] ~ dnorm(0,1e-4)
  }

  epsilon1[1]~dnorm(0,1e-4)
  for(k in 2:n.epiyr){
   epsilon1[k]~dnorm(epsilon1[k-1],prec.epsilon1)
  }
  
  prec.epsilon1 <- 1/sd.epsilon1^2
  sd.epsilon1~dunif(0,100)

  
}
"
