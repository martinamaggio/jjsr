include("methods/ellipsoid.jl")

function jsr(v)

  # this function computes different approximations of the joint
  # spectral radius and then returns the tightest bounds
  
  # matrices should all have the same size
  if ~all(x -> size(x)==size(v[1]), v)
  	println("jsr::error, matrices do not have the same size");
  	return -1;
  end

  println("called jsr function with ", size(v)[1], " matrices");

  # computing bounds using different methods
	lb_ellipsoid, up_ellipsoid = jsr_ellipsoid(v)
	
	# selecting the best bounds
	lb, ub = lb_ellipsoid, up_ellipsoid

end
