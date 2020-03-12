function jsr_sos(v)

  # ---------------------------------------------------------------- #
  # this function computes bounds on the joint spectral radius using #
  # the SOS approximation method, the code is based on the paper     #
  #                                                                  #
  # Pablo A. Parillo, Ali Jadbabaise; Approximation of the joint     #
  # spectral radius using sum of squares; Linear Algebra and its     #
  # Applications, 428, 10 (2008), pp. 2385-2402                      #
  # ---------------------------------------------------------------- #

#  m = size(v)[1] # number of matrices
#  n = size(v[1])[1] # size of matrices
#  d = 1 # ???
#  eta = min(m, binomial(n+d-1,d))
#  
#  rho_sos_2d = 1
  
#  lower_bound = eta^(-1/(2*d))*rho_sos_2d
#  upper_bound = rho_sos_2d

#  return lower_bound, upper_bound
  return NaN, NaN

end

