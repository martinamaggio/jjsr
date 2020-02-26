using LinearAlgebra
using JuMP, MosekTools

function jsr_ellipsoid(v)

  # ---------------------------------------------------------------- #
  # this function computes bounds on the joint spectral radius using #
  # the ELLIPSOID approximation method, the code is based on the     #
  # paper                                                            #
  #                                                                  #
  # Vincent D. Blondel, Yurii Nesterov, Jacquest Theys; On the       #
  # accuracy of the ellipsoid norm approximation of the joint        #
  # spectral radius; Linear Algebra and its Applications, 394 (2005) #
  # pp. 91-107                                                       #
  # ---------------------------------------------------------------- #
  
  # function to check if the entries in a matrix are non-negative
  check_all_nonnegative(m) = all(x -> x>=0, m)
  
  # if all the elements of all the matrices in the set are
  # non-negative, it is possible to approximate the joint spectral
  # radius using the spectral radius of a single matrix that is
  # constructed using the elementwise maximum elements of the
  # original matrices
  lower_bound_nn = NaN
  upper_bound_nn = NaN
  if all(x -> check_all_nonnegative(x), v)

    # generate the matrix of elementwise maximum elements
    S = reduce((x,y) -> max.(x,y), v, init=zeros(size(v[1])))
    m = size(v)[1]
    # computing the eigenvalues of S
    eigval, eigvec = eigen(S)
    # computing the spectral radius of S)
    r = maximum(abs.(eigval))
    # computing lower bound and upper bound according to equation 4
    lower_bound_nn = r/m
    upper_bound_nn = r
  end

  # in any case, we now compute the other bound, obtained by solving
  # the linear matrix inequalities with convex optimization
  # define tolerance for positive definite matrix (JuMP handles only
  # positive semidefinite matrices)
  tol = 1e-4
  n = size(v[1])[1]
  eye = Matrix{Float64}(I, n, n)
  zer = zeros(n, n)

  best_gamma = NaN # initialization
  gamma = 10 # starting value for gamma
  max_iter = 1000
  num_iter = 0
  choice = 2 # initially we will reduce of a half of gamma

  # we are going to use bisection to find the minimum value of gamma
  # which is going to be our upper bound for the joint spectral radius
  # we are going to iterate max_iter times and set up the optimization
  # problem - if the problem is feasible then the termination status
  # is OPTIMAL, otherwise it is INFEASIBLE - we use the termination
  # status to determine the next attempted value of gamma
  while(num_iter < max_iter)
    num_iter += 1

    m = Model(Mosek.Optimizer)
    set_optimizer_attribute(m, "QUIET", true)
    @objective(m, Min, 0) # just checking feasibility
    @variable(m, P[1:2, 1:2], Symmetric) # hunting P
    @SDconstraint(m, P - tol*eye >= zer) # P is positive definite
    for A in v # for all the matrices in v, adding constraint
      @SDconstraint(m, transpose(A)*P*A - gamma*gamma*P + tol*eye <= zer);
    end

    optimize!(m)
    if termination_status(m) == MOI.OPTIMAL
      best_gamma = gamma
      gamma -= gamma/choice # reduce of gamma/choice
    else
      gamma += gamma/choice # go back to last good iteration
      choice += 1 # and next time you are going to reduce less
    end
  end

  lower_bound = best_gamma/sqrt(n)
  upper_bound = best_gamma

  # mix the two results returning the tightest bound
  if ~isnan(lower_bound_nn)
    lower_bound = max(lower_bound, lower_bound_nn)
    upper_bound = min(upper_bound, upper_bound_nn)
  end

  return lower_bound, upper_bound;

end

