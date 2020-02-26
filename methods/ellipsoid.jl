using LinearAlgebra
using Convex, SCS

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
  if all(x -> check_all_nonnegative(x), v)

    # generate the matrix of elementwise maximum elements
    S = reduce((x,y) -> max.(x,y), v, init=zeros(size(v[1])))
    m = size(v)[1]
    # computing the eigenvalues of S
    eigval, eigvec = eigen(S)
    # computing the spectral radius of S)
    r = maximum(abs.(eigval))
    # computing lower bound and upper bound according to equation 4
    lower_bound = r/m
    upper_bound = r
  
  # in this case not all the elements are non-negative
  else
      
    lower_bound = NaN
    upper_bound = NaN
  
  end

  return lower_bound, upper_bound;

end

