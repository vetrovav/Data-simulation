library(moments)
SimData <- function(X,S,m,lowerbound, upperbound)
# Code in support of the paper: 
# "R Code for Data Simulation with Moment Matching" 
# 
# by Varvara Vetrova and Earl Bardsley.  Open Water Journal, 2017.

# Simulates data from a finite mixture of generalised beta distributions,
# with a bootstrap simulation being used if a beta random variable cannot be
# simulated (see the paper for conditions preventing a beta random variable
             # being generated).

# The simulated data (for a sufficiently large number of
                      # simulations) will have approximately the same mean, variance, skewness,
# and kurtosis as the original data set.

# Input:
  #   X     array of N real numbers X(N,1).
#   S     number of data simulations to be carried out.
#   m     integer value set to 4, 5, or 6 (sample size).
#   lowerbound          lower bound to the simulation range.
#   upperbound          upper bound to the simulation range.

# The lower bound should be <= the smallest data value
#    (set the lower bound to a value considerably less than the smallest data value if it
      #    is desired that no lower bound effect should be present).

# The upper bound should be >= the largest data value
#    (set the upper bound to a value considerably greater than the largest data value if it
      #    is desired that no upper bound effect should be present).


# Output:
  #   sims(S,1) contains the S simulated values.
#   Boot    number of bootstrap simulations.

{
 
  d=0.0004        # Factor of safety for the beta distribution region (may need to
                   # be increased if numerical issues arise.
  sims=matrix(0,S,1) # Returns the simulated values.
  rootb1=matrix(0,S,1) # Skewness values of the generated samples in the simulation process.
  b1=matrix(0,S,1)
  b2=matrix(0,S,1) # Kurtosis values of the generated samples in the simulation process. 
  
  # All generated samples are of size m.
  
  p=matrix(0,S,1) # Beta distribution first shape parameter.
  q=matrix(0,S,1) # Beta distribution second shape parameter.
  a=matrix(0,S,1) # Beta distribution lower location parameter.
  b=matrix(0,S,1) # Beta distribution upper location parameter.
  
  
  B=matrix(0,S,m)    # Holds rescaled values of the S samples generated.
  
  m2=matrix(0,S,1) # Variances of standardised samples.
  av=matrix(0,S,1) # Means of standardised samples.
  
  # Work arrays.
  D=matrix(0,S,1)
  r=matrix(0,S,1)
  b_a=matrix(0,S,1)
  
  N=dim(X)[1] # Find the number of data points passed to the function.
  
  # Generate S samples (each time without replacement), of size m, from the N
  # data values.Samples are in rows.
  #X needs to be a row vector
  A=replicate(S,sample(X,size=m,replace=FALSE))
  A=t(A)
 
  # Set all the simulation default output as bootstrap values (one from each
                                                               # of the generated S samples).
  # To be replaced later (where permitted) with simulated values from the finite mixture of
  # generalised beta distributions.
  sims[,1]=A[,1]
  
  # Scale all sample values with nonzero variance to the range 0-1:
    
  # Find the maximum and minimum value of each of the S samples.
  max_value=apply(A,1,max)
  min_value=apply(A,1,min)
  
  # Find the indices for nonzero variances.
  II=which(max_value>min_value)
  
  # Scale the values within the S samples and store in the array B.
  for (j in 1:m)
  {
  B[II,j]=(A[II,j]-min_value[II])/(max_value[II]-min_value[II]);
  }
  
  # Find skewness and kurtosis values for all samples with nonzero variance.
  # b1 and b2 are respectively the horizontal and vertical axes on the
  # Pearson plot of systems of distributions (beta1 and beta2 in the usual
                                              # Pearson plot notation).
  rootb1=matrix(0,S,1)
  b1=matrix(0,S,1)
  b2=matrix(0,S,1)
  rootb1[II]=apply(B[II,],1,skewness)
  b1[II]=rootb1[II]^2
  b2[II]=apply(B[II,],1,kurtosis)

  # Find the indices of samples permitting a beta simulation from within the beta 
  # distribution region as specified by the value of d.
  JJ=which(b2[,1]>0 &  b2[,1]> b1[,1] + 1.0 + d  & b2[,1] < 1.5*b1[,1] + 3.0 - d)
  
  # Find the average value and variance of the scaled samples permitting a
  # beta simulation.
  av[JJ]=apply(B[JJ,],1,mean)#mean(B[JJ,],2);
  m2[JJ]=apply(B[JJ,],1,var)#moment(B(JJ,:),2,2); ## Might need to check it
  
  # The symbols used in the code here largely follow those of the paper.
  r[JJ]=6.0*(b2[JJ]-b1[JJ]-1)/(6.0+3.0*b1[JJ]-2.0*b2[JJ])
  D[JJ]=b1[JJ]*(r[JJ]+2)^2+ 16*(r[JJ]+1)
  b_a[JJ]=m2[JJ]^0.5* D[JJ]^0.5/2 # b_a is "b-a" in the paper.
  
  # Obtain the beta distribution two shape parameters (moment-based
                                                       # estimates).
  p[JJ]=(r[JJ]/2)*(1 - (r[JJ] + 2)*rootb1[JJ]/ D[JJ]^0.5)
  q[JJ]=(r[JJ]/2)*(1 + (r[JJ] + 2)*rootb1[JJ]/D[JJ]^0.5)
  
  # Obtain the beta distribution two location parameters (moment-based
                                                          # estimates for scaled sample values).
  a[JJ] = av[JJ] - b_a[JJ]*p[JJ]/r[JJ] # lower location parameter
  b[JJ] = a[JJ] + b_a[JJ]                    # upper location parameter
  
  
  # Convert the location parameter values back to the original scale.
  a[JJ]=min_value[JJ]+a[JJ]*(max_value[JJ]-min_value[JJ])
  b[JJ]=min_value[JJ]+(max_value[JJ]-min_value[JJ])*b[JJ]
  
  # Exclude beta distributions with a location parameter outside
  # the permissible data simulation range as defined by the specified values
  # for lowerbound and upperbound.
  KK = which(a[JJ]>=lowerbound & b[JJ]<=upperbound)
  KK = JJ[KK] # This might be written more efficiently but is left for ease of reading.
  
  # Simulate the beta random variables (where beta random variable generation is permissible)
  # from the standard beta distribution.
  sims[KK]=mapply(rbeta,p[KK],q[KK],n=1)#rbeta(n, shape1, shape2, ncp = 0)#betarnd(p(KK),q(KK));
  
  # Rescale the beta random variables to allow for different location
  # parameters. This gives the required KK data simulations.
  sims[KK] =a[KK]+sims[KK]*(b[KK]-a[KK])
  
  # Find the number of bootstrap simulations that were used in the KK
  # simulations (the remainder of the KK values being random variables from
                 # the finite mixture of generalised beta distributions).
  Boot= S - nrow(KK)
  
  res_list<-list(sims,Boot)
  return(res_list)
}

