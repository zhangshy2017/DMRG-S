"""
   dmrgs(H::MPO,H0::MPO,H1::MPO,H2::MPO,psi0::MPS, sweeps::Sweeps; kwargs...)
                    
Use the DMRG-S algorithm to extract quantum many-body scarred eigenstates for generic Hamiltonians.
The MPS `psi0` is used to initialize the MPS to be optimized
The MPO `H` is the shift MPO (Ĥ-ξ), ξ is the target energy
The MPO `H0` is the Identity MPO Î
The MPO `H1` is the Hamiltonian MPO Ĥ used to evaluate the energy E and variance σ during the optimization
The MPO `H2` is any other obsevable MPO which need to evaluate during the optimization	
and the `sweeps` object determines the parameters used to 
control the DMRG-S algorithm.

Returns:
* `psi::MPS` - optimized MPS
"""



function dmrgs(H::MPO,H0::MPO,H1::MPO,H2::MPO,psi0::MPS, sweeps::Sweeps; kwargs...)
  check_hascommoninds(siteinds, H, psi0)
  check_hascommoninds(siteinds, H, psi0')
  # Permute the indices to have a better memory layout
  # and minimize permutations
  H = permute(H, (linkind, siteinds, linkind))
  H0 = permute(H0, (linkind, siteinds, linkind))
  H1 = permute(H1, (linkind, siteinds, linkind))
  PHH0 = ProjMPO(H1)
  PH0 = ProjMPO(H1)
  PHH = ProjMPO(H)
  PHH.product_label=2
  Z2H = ProjMPO(H0)
  NN = ProjMPO(H2)
  return dmrgs(PHH,PHH0,PH0,Z2H,NN,psi0,sweeps;kwargs...)
end



"""
   simps(H::MPO,H1::MPO,H2::MPO,psi0::MPS, sweeps::Sweeps; kwargs...)

                    
Use the SIMPS algorithm to obtaining quantum many-body localized eigenstates for MBL phases.
The MPS `psi0` is used to initialize the MPS to be optimized,
and the `sweeps` object determines the parameters used to 
control the DMRG algorithm.

Returns:
* `psi::MPS` - optimized MPS
"""





function simps(H::MPO,H1::MPO,H2::MPO,psi0::MPS, sweeps::Sweeps; kwargs...)
  check_hascommoninds(siteinds, H, psi0)
  check_hascommoninds(siteinds, H, psi0')
  # Permute the indices to have a better memory layout
  # and minimize permutations
  H = permute(H, (linkind, siteinds, linkind))
  H1 = permute(H1, (linkind, siteinds, linkind))
  H2 = permute(H2, (linkind, siteinds, linkind))
  PHH = ProjMPO(H)
  PHH0 = ProjMPO(H1)
  PH0 = ProjMPO(H1)
  PHH.product_label=2
  Z2H = ProjMPO(H)
  NN = ProjMPO(H2)
  return simps(PHH,PHH0,PH0,Z2H,NN,psi0,sweeps;kwargs...)
end




function dmrgs(PHH,PHH0, PH0, Z2H, NN, psi0::MPS, sweeps::Sweeps; kwargs...)

  which_decomp::Union{String, Nothing} = get(kwargs, :which_decomp, nothing)
  svd_alg::String = get(kwargs, :svd_alg, "recursive")
  obs = get(kwargs, :observer, NoObserver())
  outputlevel::Int = get(kwargs, :outputlevel, 1)

  eigsolve_tol::Float64 = get(kwargs, :eigsolve_tol, 1e-14) #can be adjusted
  eigsolve_krylovdim::Int = get(kwargs, :eigsolve_krylovdim, 20) #can be adjusted
  eigsolve_maxiter::Int = get(kwargs, :eigsolve_maxiter, 50) #can be adjusted
  eigsolve_verbosity::Int = get(kwargs, :eigsolve_verbosity, 0)

  ishermitian::Bool = get(kwargs, :ishermitian, true)

  sites = siteinds(psi0)

  isdone = false

  if haskey(kwargs, :errgoal)
    error("errgoal keyword has been replaced by eigsolve_tol.")
  end

  if haskey(kwargs, :quiet)
    error("quiet keyword has been replaced by outputlevel")
  end


  N = length(psi0)
  psi = copy(psi0)
  psi2 = copy(psi0)




  if !isortho(psi) || orthocenter(psi) != 1
    orthogonalize!(psi,1)
  end
  @assert isortho(psi) && orthocenter(psi) == 1

  if !isortho(psi2) || orthocenter(psi2) != 1
    orthogonalize!(psi2,1)
  end
  @assert isortho(psi2) && orthocenter(psi2) == 1

  position2!(PHH, psi, 1)
  position2!(PHH0, psi, 1)
  position!(PH0, psi, 1)
  position!(NN, psi, 1)
  position3!(Z2H, psi,psi2, 1)

  step=0
 
  for sw=1:nsweep(sweeps)
       sw_time = @elapsed begin
       for (b, ha) in sweepnext(N)
            step += 1     
            if isfile("STOP_DMRGS")
                 rm("STOP_DMRGS")
                 isdone = true
            end
            @timeit_debug timer "dmrgs: position!" begin
            position2!(PHH, psi, b)
            position2!(PHH0, psi, b)
            position!(PH0, psi, b)
            position!(NN, psi, b)
            position3!(Z2H, psi, psi2, b)
            end


            @timeit_debug timer "dmrgs: psi[b]*psi[b+1]" begin
            phi = psi[b] * psi[b+1]
            phiZ2 = psi2[b] * psi2[b+1]
            end


            bb = product2(Z2H,dag(prime(phiZ2)))
            
            @timeit_debug timer "dmrgs: linsolve" begin
            vec, info =  linsolve(PHH, bb, phi, ishermitian = ishermitian,
                                   maxiter = eigsolve_maxiter,
                                   krylovdim = eigsolve_krylovdim,
                                   issymmetric = true,
                                   isposdef = true,
                                   tol = eigsolve_tol)  
            end
           



            phi = vec
            normalize!(phi)

            vec2=copy(vec)

            temp_E = scalar(product2(PH0,vec2)*dag(vec2))
            variance = scalar(product3(PHH0,vec2)*dag(vec2)) - temp_E^2
            NN_int = scalar(product2(NN,vec2)*dag(vec2))



           
            @printf("step %d, E= %.12f, variance= %.12f,NN= %.8f\n" ,step, temp_E, variance, NN_int)
            ortho = ha == 1 ? "left" : "right"

            drho = nothing
 
            
            
            @timeit_debug timer "dmrgs: replacebond!" begin
            spec = replacebond!(psi, b, phi; maxdim = maxdim(sweeps, sw),
                                             mindim = mindim(sweeps, sw),
                                             cutoff = cutoff(sweeps, sw),
                                             eigen_perturbation = drho,
                                             ortho = ortho,
                                             normalize = true,
                                             which_decomp = which_decomp,
                                             svd_alg = svd_alg)
            end

            isdone && break

            if outputlevel >= 2
              @printf("(Truncated using cutoff=%.1E maxdim=%d mindim=%d)\n",
                      cutoff(sweeps, sw),maxdim(sweeps, sw),mindim(sweeps, sw))
              @printf("Trunc. err=%.1E, bond dimension %d\n\n",spec.truncerr,dim(linkind(psi,b)))
            end

            measure!(obs; psi = psi,
                          bond = b,
                          sweep = sw,
                          half_sweep = ha,
                          spec = spec,
                          outputlevel = outputlevel)



          end

       end
    if outputlevel >= 1
      @printf("After sweep %d  maxlinkdim=%d time=%.3f\n",
              sw, maxlinkdim(psi), sw_time)
    end
    
    isdone && break

  end
  
  return  psi
end




function simps(PHH, PHH0,PH0, Z2H, NN, psi0::MPS, sweeps::Sweeps; kwargs...)

  which_decomp::Union{String, Nothing} = get(kwargs, :which_decomp, nothing)
  svd_alg::String = get(kwargs, :svd_alg, "recursive")
  obs = get(kwargs, :observer, NoObserver())
  outputlevel::Int = get(kwargs, :outputlevel, 1)

  eigsolve_tol::Float64 = get(kwargs, :eigsolve_tol, 1e-14)
  eigsolve_krylovdim::Int = get(kwargs, :eigsolve_krylovdim, 20)
  eigsolve_maxiter::Int = get(kwargs, :eigsolve_maxiter, 50)
  eigsolve_verbosity::Int = get(kwargs, :eigsolve_verbosity, 0)

  ishermitian::Bool = get(kwargs, :ishermitian, true)

  isdone = false

  if haskey(kwargs, :errgoal)
    error("errgoal keyword has been replaced by eigsolve_tol.")
  end

  if haskey(kwargs, :quiet)
    error("quiet keyword has been replaced by outputlevel")
  end


  N = length(psi0)
  psi = copy(psi0)
  psi2 = copy(psi0)



  if !isortho(psi) || orthocenter(psi) != 1
    orthogonalize!(psi,1)
  end
  @assert isortho(psi) && orthocenter(psi) == 1

  if !isortho(psi2) || orthocenter(psi2) != 1
    orthogonalize!(psi2,1)
  end
  @assert isortho(psi2) && orthocenter(psi2) == 1

  position2!(PHH0, psi, 1)
  position!(PH0, psi, 1)
  position2!(PHH, psi, 1)
  position3!(Z2H, psi,psi2, 1)
  position!(NN, psi, 1)

  step=0
 
  for sw=1:nsweep(sweeps)
       sw_time = @elapsed begin
       for (b, ha) in sweepnext(N)
            step += 1
            if isfile("STOP_DMRGS")
                 rm("STOP_DMRGS")
                 isdone = true
            end     
            @timeit_debug timer "dmrg: position!" begin
            position2!(PHH0, psi, b)
            position!(PH0, psi, b)
            position2!(PHH, psi, b)
            position3!(Z2H, psi, psi2, b)
            position!(NN, psi, b)
            end


            @timeit_debug timer "simps: psi[b]*psi[b+1]" begin
            phi = psi[b] * psi[b+1]
            phiZ2 = psi2[b] * psi2[b+1]
            end


            bb = product2(Z2H,dag(prime(phiZ2)))
            
            @timeit_debug timer "simps: eigsolve" begin
            vec, info =  linsolve(PHH, bb, phi, ishermitian = ishermitian,
                                   maxiter = eigsolve_maxiter,
                                   krylovdim = eigsolve_krylovdim,
                                   issymmetric = true,
                                   isposdef = true,
                                   tol = eigsolve_tol)  
            end
            vec2 = copy(vec)
            normalize!(vec2)
            temp_E = scalar(product2(PH0,vec2)*dag(vec2))
            NN_int = scalar(product2(NN,vec2)*dag(vec2))
            variance = scalar(product3(PHH0,vec2)*dag(vec2)) - temp_E^2
            error = 1 - abs(scalar(vec*bb))


            phi = vec

            @printf("step %d, 1-|<ψ|ϕ>| = %.8f, E = %.8f, variance = %.12f, NN = %.8f\n",step, error, temp_E, variance, NN_int)
            ortho = ha == 1 ? "left" : "right"

            drho = nothing
 
            
            
            @timeit_debug timer "simps: replacebond!" begin
            spec = replacebond!(psi, b, phi; maxdim = maxdim(sweeps, sw),
                                             mindim = mindim(sweeps, sw),
                                             cutoff = cutoff(sweeps, sw),
                                             eigen_perturbation = drho,
                                             ortho = ortho,
                                             normalize = true,
                                             which_decomp = which_decomp,
                                             svd_alg = svd_alg)
            end


            if outputlevel >= 2
              @printf("(Truncated using cutoff=%.1E maxdim=%d mindim=%d)\n",
                      cutoff(sweeps, sw),maxdim(sweeps, sw),mindim(sweeps, sw))
              @printf("Trunc. err=%.1E, bond dimension %d\n\n",spec.truncerr,dim(linkind(psi,b)))
            end

            measure!(obs; psi = psi,
                          bond = b,
                          sweep = sw,
                          half_sweep = ha,
                          spec = spec,
                          outputlevel = outputlevel)

            isdone && break 

          end

       end
    if outputlevel >= 1
      @printf("After sweep %d  maxlinkdim=%d time=%.3f\n",
              sw, maxlinkdim(psi), sw_time)
    end
    

    isdone && break

  end
  
  return  psi
end



