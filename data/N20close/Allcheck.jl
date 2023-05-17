using ITensors
using Printf
using Random
using LinearAlgebra
using ITensors.HDF5


let
  N = 20
  U =1000 
  varlist =zeros(N+1)
  Elist =zeros(N+1)
  Enn =zeros(N+1)
  Z2_overlap = zeros(N+1)
  entropy = zeros(N+1)

  for k=1:11
      f = h5open(string("N",string(N),"_",string(k),".h5"),"r")
      psi0 = read(f,"psi",MPS)
      close(f)

      sites = siteinds(psi0)

      ampo1 = AutoMPO()

      for j=1:N-1
            ampo1 += U,"ProjUp",j,"ProjUp",j+1
      end
      ampo1 += U,"ProjUp",N,"ProjUp",1

      H1 = MPO(ampo1, sites)

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

      states = [isodd(n) ? "Up" : "Dn" for n=1:N]
      Z2psi = productMPS(sites, states)

      Elist[k] = inner(psi0',H,psi0) 
      Enn[k] = inner(psi0',H1,psi0)
      varlist[k] = inner(H,psi0,H,psi0) - Elist[k]^2
      Z2_overlap[k] = inner(psi0,Z2psi)

      center = Int(N/2)
      orthogonalize!(psi0,center)
      phi = psi0[center] * psi0[center+1]
      _,S = svd(phi, (linkind(psi0, center-1), siteind(psi0, center)))
      
      SvN_initial = 0.0
      for n =1: dim(S, 1)
             p = S[n,n]^2
          if p>0.0000000001
               SvN_initial -= p * log(p)
          end
      end
     
      entropy[k] = SvN_initial
 
  end

  @printf("------------------------\n")
  for k=1:11
	  @printf("Evar=%.16f, E=%.16f, NNint=%.16f, Z2=%.16f, Entropy=%.16f, %d\n",varlist[k] ,Elist[k], Enn[k],Z2_overlap[k], entropy[k], k)
  end



end

