using ITensors
using Printf
using Random
using LinearAlgebra
using ITensors.HDF5


let
   N = 30
   maxD = 200
   initial_energy = -1.6
   U = 1000
   nk=100

   varray = zeros(1,nk)

   minvalue = 0.1
   stop_value = 0.0000000001

   sites = siteinds("S=1/2",N,conserve_qns=false)

   ampo = AutoMPO()

   for j=1:N
      if j==1
        ampo += 1.0,"ProjDn",N,"X",1,"ProjDn",2
      elseif j==N
        ampo += 1.0,"ProjDn",N-1,"X",N,"ProjDn",1
      else
        ampo += 1.0,"ProjDn",j-1,"X",j,"ProjDn",j+1
      end
   end

   H = MPO(ampo, sites)


   ampo1 = AutoMPO()

   for j=1:N-1
        ampo1 += U,"ProjUp",j,"ProjUp",j+1
   end
   ampo1 += U,"ProjUp",N,"ProjUp",1

   H1 = MPO(ampo1, sites)


   ampo0 = AutoMPO()
   ampo0 += 1.0,"Id",1
   H0 = MPO(ampo0, sites)

   numlist = 1:N
   states = [isodd(numlist[n]) ? "Up" : "Dn" for n=1:N]
   psi0 = productMPS(sites,states)

   H2 = H  - (initial_energy)*H0

   sweeps = Sweeps(1)
   maxdim!(sweeps, maxD,maxD)
   cutoff!(sweeps, 1E-16)

   @show sweeps

   @printf("------------------------------------\n")
   psi = copy(psi0)
   for k=1:nk
       @show k
       psi = dmrgx(H2,psi0,psi0, sweeps)
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


       ampo1 = AutoMPO()
       ampo1 += -1.0,"ProjUp",N,"ProjUp",1
       ampo1 += 1.0,"Id",1
       H_tmp = MPO(ampo1, sites)

       temp = contract(H_tmp,psi0,maxdim = maxD, normalize=true)
       setprime!(temp,0;plev=1)

       psi0 = copy(temp)
       orthogonalize!(psi0,1)
       normalize!(psi0[1])

       psi =copy(psi0)


       Elist = inner(psi,H,psi)
       Elist1 = inner(psi,H1,psi)
       varlist = inner(H,psi,H,psi) - Elist^2
       @printf("PXP E variance=%.12f,    E=%.12f\n",varlist ,Elist)
       @printf("EXP E=%.12f\n" ,Elist1)


       varray[k]=varlist

#       if (varlist < 0.1)&&(varlist<=minvalue)
#	        H2 = H - Elist*H0
#	        minvalue = varlist
#       end


       f = h5open(string("N",string(N),"_",string(k),".h5"),"w")
       write(f,"psi",psi)
       close(f)
       minvalue = varlist
       if varlist<stop_value
           break
       end

       @printf("------------------------\n")
       for kk=1:k
              @printf("%.12f\n",varray[kk])
       end


   end

   Elist = inner(psi,H,psi)
   Elist1 = inner(psi,H1,psi)
   varlist = inner(H,psi,H,psi) - Elist^2


   @printf("------------------------\n")
   @printf("PXP E variance=%.12f,    E=%.12f\n",varlist ,Elist)
   @printf("E=%.12f\n" ,Elist1)



 end
