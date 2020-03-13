include("methods/ellipsoid.jl")
include("methods/sos.jl")
include("methods/branchbound.jl")

function jsr(v; verbose = false)

  # this function computes different approximations of the joint
  # spectral radius and then returns the tightest bounds
  
  # matrices should all have the same size
  if ~all(x -> size(x)==size(v[1]), v)
    println("jsr::error, matrices do not have the same size")
    return -1;
  end

  if (verbose)
    println("called jsr function with ", size(v)[1], " matrices")
    println("-------------------------------------------------------")
  end

  # computing bounds using different methods
  lb_bnb, ub_bnb = jsr_branchbound(v)
  lb_ellipsoid, ub_ellipsoid = jsr_ellipsoid(v)
  lb_sos, ub_sos = jsr_sos(v)
  if (verbose)
    println("branch&bound: [", lb_bnb, ", ", ub_bnb, "]")
    println("   ellipsoid: [", lb_ellipsoid, ", ", ub_ellipsoid, "]")
    println("         sos: [", lb_sos, ", ", ub_sos, "]")
  end
  
  lb = max(lb_bnb, lb_ellipsoid)#, lb_sos)
  ub = min(ub_bnb, ub_ellipsoid)#, ub_sos)
  if (verbose)
    println("-------------------------------------------------------")
    println("[", lb, ", ", ub, "]")
    println("#######################################################")
  end
  
  return lb, ub

end
