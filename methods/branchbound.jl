using LinearAlgebra
using IterTools

function jsr_branchbound(v; delta=1e-4, max_iter=10)

  # ---------------------------------------------------------------- #
  # this function computes bounds on the joint spectral radius using #
  # the branch and bound method, the code is based on the paper      #
  #                                                                  #
  # Gustaf Gripenberg; Computing the Joint Spectral Radius; Linear   #
  # Algebra and its Applications, 234 (1996), pp. 43-60              #
  # ---------------------------------------------------------------- #

  spectral_radiuses(T) = map(x -> maximum(abs.(eigvals(x))), T)
  norms(T) = map(x -> norm(x), T)
  product(T) = reduce((x,y) -> x*y, T, init=I)

  # iteration 1
  sp1 = spectral_radiuses(v)
  alpha1 = maximum(sp1)
  nm1 = norms(v)
  beta1 = maximum(nm1)
  
  # branch and bound
  To = copy(v) # initializing the set of the previous iteration
  alphan = alpha1 # initializing the bounds for the previous iteration
  betan = beta1
  
  # advance only for a maximum number of iterations, may just finish
  # before the maximum iterations are reached
  for num_iter = 2:max_iter

    # we are creating the set Tn for the n-th iteration using the
    # set for the (n-1)-th iteration in To, and the initial set v
    # note: we need to multiply an entry in v with an entry in To
    Tn = Matrix{Float64}[]
    for x in Iterators.product(v, To)
      # we are going to compute p(candidate) and add it to the list
      # for the next iteration only if it is greater than the current
      # alpha + delta, if Tn is empty then we can't improve with
      # branch and bound and we should exit and return the bounds we
      # found - alphan and betan
      p = min(norm(x[1]), norm(x[1]*x[2])^(1/2))
      if (p > alphan + delta) push!(Tn, x[1]*x[2]) end
    end
    
    if isempty(Tn) break end # get out if we can't improve
    
    sp_pow = map(x -> x^(1/num_iter), spectral_radiuses(Tn))
    p = [minimum(norm(product(Tn[1:idx]))^(1/idx)) for idx = 1:length(Tn)]
    betan = min(betan, max(alphan + delta, maximum(p)))
    alphan = max(alphan, maximum(sp_pow))
    To = Tn # preparing for the next iteration
    
  end
  
  lower_bound = alphan
  upper_bound = betan

  return lower_bound, upper_bound

end

