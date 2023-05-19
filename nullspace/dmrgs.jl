using ITensors
using Printf
using Random
using LinearAlgebra
using ITensors.HDF5


let
   N = 20
   maxD = 10
   U = 1000
   nk=200

   nth=1
   initial_energy = 0.01*(2*rand()-1.0)


   minvalue = 0.1
   stop_value = 0.00000001

   sites = siteinds("S=1/2",N,conserve_qns=false)

   ampo = AutoMPO()

   for j=1:N
      if j==1
        ampo += 1.0,"X",1,"ProjDn",2
      elseif j==N
        ampo += 1.0,"ProjDn",N-1,"X",N
      else
        ampo += 1.0,"ProjDn",j-1,"X",j,"ProjDn",j+1
      end
   end

   H = MPO(ampo, sites) 


   ampo1 = AutoMPO()

   for j=1:N-1
        ampo1 += U,"ProjUp",j,"ProjUp",j+1
   end

   H1 = MPO(ampo1, sites) 


   ampo3 = AutoMPO()

   for j=1:N-1
        ampo3 += U,"ProjUp",j,"ProjUp",j+1
   end
   ampo3 += U,"ProjUp",N,"ProjUp",1
   ampo3 += 1.0,"Id",1
   H3 = MPO(ampo3, sites) 


   ampo0 = AutoMPO()
   ampo0 += 1.0,"Id",1
   H0 = MPO(ampo0, sites)

   numlist = rand(1:2,N)
   states = [isodd(numlist[n]) ? "Up" : "Dn" for n=1:N]
   psi0 = randomMPS(sites,states,2)
  
   H2 = H  - (initial_energy)*H0

   sweeps = Sweeps(1)
   maxdim!(sweeps, maxD,maxD)
   cutoff!(sweeps, 1E-16)
   
   @show sweeps

   @printf("------------------------------------\n")   
   psi = copy(psi0)   

   sw_time = @elapsed begin
       for k=1:nk
           @show k
           psi = dmrgs(H2,H0,H,H1,psi0, sweeps)
           psi0 = copy(psi)


           orthogonalize!(psi0,1)
           for j=1:N-1
                 @show j
                 phi = psi0[j] * psi0[j+1]
                 phi0 = copy(phi)
                 phi1 = mapprime(phi*op("ProjUp",sites[j]),1 => 0)
                 phi2 = mapprime(phi1*op("ProjUp",sites[j+1]),1 => 0)
                 phi3= phi0 - phi2

                 spec = replacebond!(psi0, j, phi3; maxdim = maxD,
                                                  mindim = 1,
                                                  cutoff = 1e-15,
                                                  eigen_perturbation = nothing,
                                                  ortho = "left",
                                                  normalize = true,
                                                  which_decomp = nothing,
                                                  svd_alg = "recursive")

           end

           psi =copy(psi0)


           Elist = inner(psi,H,psi)
           Elist1 = inner(psi,H1,psi)
           varlist = inner(H,psi,H,psi) - Elist^2
           @printf("PXP E variance=%.12f,    E=%.12f\n",varlist ,Elist)

           if (varlist < 0.001)&&(varlist<=minvalue)
              H2 = H - (Elist-(-1)^k*3*varlist)*H0
           end


           if varlist<minvalue
              f = h5open(string("N",string(N),"_",string(nth),".h5"),"w")
              write(f,"psi",psi)
              close(f)
              minvalue = varlist
              if varlist<stop_value
                  break
              end
           end
       end
   end
 
   Elist = inner(psi,H,psi)
   Elist1 = inner(psi,H1,psi)
   varlist = inner(H,psi,H,psi) - Elist^2


   @printf("------------------------\n")
   @printf("PXP E variance=%.12f,    E=%.12f\n",varlist ,Elist)
   @printf("E=%.12f\n" ,Elist1)


   @printf("------------------------------------\n")
   @printf("%.12f\n" ,sw_time)


   
 end


