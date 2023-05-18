using ITensors
using Printf
using Random
using LinearAlgebra
using ITensors.HDF5


let
  N = 20
  U =1000.0
  maxD = 10
  num=1000
  num_random =1
  NN = ones(num*num_random)
  varlist =ones(num*num_random)
  Elist =zeros(num*num_random)



  mark = 0
  for k=1:num
          if isfile(string("job",string(k),"/N",string(N),"_",string(k),".h5"))
             f = h5open(string("job",string(k),"/N",string(N),"_",string(k),".h5"),"r")
             psi0 = read(f,"psi",MPS)
             close(f)
             mark = mark+1
          else
             continue
          end
          @show k

          sites = siteinds(psi0)

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
          H = MPO(ampo,sites)


          ampo1 = AutoMPO()

          for j=1:N-1
               ampo1 += U,"ProjUp",j,"ProjUp",j+1
          end

          H1 = MPO(ampo1, sites)


          Elist[mark] = inner(psi0,H,psi0)
          varlist[mark] = abs(inner(H,psi0,H,psi0) - Elist[mark]^2)
          NN[mark] = inner(psi0,H1,psi0)
  end

  varlist = varlist[1:mark]
  Elist = Elist[1:mark]

  index  = sortperm(varlist, rev=true)


  @printf("------------------------\n")
  for k=1:mark
          @printf("E variance=%.12f,    E=%.12f,  NN=%.12f,   %d,  %d\n",varlist[index[k]] ,Elist[index[k]], NN[index[k]], index[k], k)
  end



end
