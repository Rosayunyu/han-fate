
sink("model.txt")
cat("
    
model 
{  
    for (p in 1:num_pH_Cl2)
    {
      pH[p] ~ dnorm(pH_levels[p],1/((0.06/3)^2))
      H[p] <- 10^(-pH[p])
      OH[p] <- 10^(-14+pH[p])
      alpha_0[p] <- H[p]/(H[p]+10^(-pKa[p]))
      alpha_1[p] <- 10^(-pKa[p])/(H[p]+10^(-pKa[p]))
    }

      for (com in 1:num_cmpd) 
      {		
        k_H2O[com] ~ dunif(-1e+07,1e+07)
        k_OH[com] ~ dunif(-1e+07,1e+07)
        k_HOCl[com] ~ dunif(-1e+07,1e+07)
        k_OCl[com] ~ dunif(-1e+07,1e+07)
    
        for (p in 1:num_pH_Cl2) 
        {
        
          for (r in 1:replicate_array[p+1]) 
          {
            lnC0[com,(r+sum(replicate_array[1:p]))] ~ dnorm(ln_dose[com],1/((log(0.5)/3)^2))
            precission_e[com,(r + sum(replicate_array[1:p]))] ~ dgamma(0.1,0.1)
    
            for (t in 1:num_Time) 
            {     
              HAN_Loss[t,(r + sum(replicate_array[1:p])),com] <- (k_H2O[com] + k_OH[com]*OH[p])*Time[t] + (k_HOCl[com]*alpha_0[p]+k_OCl[com]*alpha_1[p])*CT_INTEGRAL[t,(r + sum(replicate_array[1:p]))]
              FINAL_DATA[t,(r + sum(replicate_array[1:p])),com] ~ dnorm(lnC0[com,(r+sum(replicate_array[1:p]))] - HAN_Loss[t,(r + sum(replicate_array[1:p])),com], precission_e[com,(r + sum(replicate_array[1:p]))])
            }		
          }
        }
      }
    
}
",fill=TRUE)
sink()


