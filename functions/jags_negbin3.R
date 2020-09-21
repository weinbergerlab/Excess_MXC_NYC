model_string_negbin3<-
  "model
{
  # Likelihood
  for( t in 1:n.dates ){
     
  for(d in 1:(D+1)){
  
  #First likelihood is for full data
  n[t,d] ~ dnegbin(p1[t,d],r1)
  # Conversions
  p1[t,d] <- r1/(r1+lambda[t,d])
  
  log(lambda[t,d]) <- (int +     
            beta.logged2[d]*step(D-d) +
            sum.beta.logged2[t]*(1-step(D-d)) +
                    sin26[t]*delta1[1] +
                    cos26[t]*delta1[2] +
                    sin52[t]*delta1[3] +
                    cos52[t]*delta1[4] +
                    epsilon1[epiyr.index[t]] +
                    phi[t] 
                    ) 
  }
    log(baseline[t]) <- (int +     
                    sin26[t]*delta1[1] +
                    cos26[t]*delta1[2] +
                    sin52[t]*delta1[3] +
                    cos52[t]*delta1[4] +
                    epsilon1[epiyr.index[t]]
    )

  sum.lambda[t] <- sum(lambda[t,1:D])
    
    # Conversions
  p3[t] <- r1/(r1+baseline[t])
  baseline.n[t] ~dnegbin(p3[t],r1)
  
  sum.n[t]  <- sum(n[t,1:D])

  excessN[t] <- sum.n[t] - baseline.n[t]
  }

    #Time Series Random Effect--starts when extrapolation starts
    for(t in 1:(extrap.index-1)){
    phi[t] <- 0
    }
  for(t in extrap.index:n.dates){
     phi[t] ~ dnorm((phi[t-1]), tau2.alpha)
     }

   ##  
  int~dnorm(0, 1e-4)
  
  # Prior for neg binomial rate
  r1 ~ dgamma(0.001,0.001)

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
