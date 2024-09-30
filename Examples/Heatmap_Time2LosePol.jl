# Code to run an array of simulations and make a heatmap of how long the mechanochemical system takes to
# lose polarity presented in article "INSERT DOI". 
# One can work with pure local inhibition by changing the following mechanical parameters:
# vCTE=0   α=0   β=0   σₐ₀=0
# The system will start with flat Rho and Rac at an steady state
# Different mechanical parametes can eb changed in the main section of the code.
# Written by Andreu F Gallen working in Turlier lab and in collaboration with Orion Weiner's lab


using Gridap
using Plots
using Plots.PlotMeasures
include("Plots_RhoRacA.jl")

function plotting(ylab,po,pPNG,i)
  plot(po)
  xlabel!("ξ[μm]")
  ylabel!(ylab)
  savefig(pPNG*"t="*i*".png")
end

function conservation(sρ0,sρ,Minitial)
  return 0.1*(Minitial-(sρ0+sρ))
end
 
function threshold(x,x₀,xth)
  return  (0.5 * (tanh.(x/x₀ .- xth/x₀).+1)) 
end


function run(α,β,α₀,β₀,χ,λ⁻²,η,T,Δt,part)
  # Time discretisation parameters
  t₀  = 0.0
  t   = t₀  
  nΔt = trunc(Int,T/Δt)

  domain = (0,L)
  partition = (part)
  model = CartesianDiscreteModel(domain,partition) 
  Δx = L/partition

 # Lets copy the code in the output folder to be able to check code used for each simulation
  pPNG="./PNG/"*simulation
  mkpath(pPNG) 
  cp(@__FILE__, pPNG*split(@__FILE__, "/")[end],force=true)


  order = 1
 #Starting Boundary conditions
  x₀(x) = 0
  xₗ(x) = 0

  v₀(x) = 0
  vₗ(x) = 0

  # Creating FE space
  reffe = ReferenceFE(lagrangian,Float64,order)
  W0 = TestFESpace(model,reffe;dirichlet_tags=["tag_1"])
  WD0 = TestFESpace(model,reffe;dirichlet_tags=["tag_1","tag_2"])
  X = TrialFESpace(WD0,[x₀,xₗ])
  V = TrialFESpace(WD0,[v₀,vₗ])

  degree = 2
  #triangulate FE space
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  #Neumann and defining trial spaces and test space for ρ_b 
  Q0 = TestFESpace(model,reffe)
  RHO =  TrialFESpace(Q0)
  RHO0 =  TrialFESpace(Q0)
  Rac =  TrialFESpace(Q0)
  Rho =  TrialFESpace(Q0)

  #Building the vectors used for introducing opto influence as an increase in α and β
  α₀o(x)= 2*α₀*exp(-0.5*(x[1])^2/((L/4)^2)) 
  β₀o(x)= 2*β₀*exp(-0.5*(x[1]-L)^2/((L/4)^2))
  α₀opto = zeros(partition)
  β₀opto = zeros(partition)
  for j in 1:1:(partition)
    α₀opto[j] = α₀o(j*Δx)
    β₀opto[j] = β₀o(j*Δx)
  end

  #defining Starting conditions for some variables and dummy variables to be able to build the equations
  uh_ρ = interpolate_everywhere(kon*M0/partition/(koff+kon),RHO)
  uh_ρ_old = uh_ρ
  uh_ρ0 = interpolate_everywhere(koff*M0/partition/(koff+kon),RHO0)
  uh_ρ0_old = uh_ρ0
  uh_x = zero(X)
  uh_x_old = uh_x
  uh_v = zero(V)
  #in this code in particular we START POLARIZED
  f1(x) = 2*α₀*exp(-(x[1])^2/wrac^2)
  f2(x) = 2*β₀*exp(-(x[1]-L)^2/wrac^2)
  uh_rac = interpolate_everywhere(f1 ,Rac)
  uh_rho = interpolate_everywhere(f2 ,Rho)
  uh_rac_old = uh_rac
  uh_rho_old = uh_rho
  Minitial=sum(get_free_dof_values(uh_ρ0))+sum(get_free_dof_values(uh_ρ))

  print("Start1 sum(ρ_b) "*string(sum(get_free_dof_values(uh_ρ)))*", sum(ρ_u) "*string(sum(get_free_dof_values(uh_ρ0)))*", and sum(ρ_b+ρ_u)"*string(sum(get_free_dof_values(uh_ρ0))+sum(get_free_dof_values(uh_ρ)))*"\n")
 
 #variables to store temporal information that we would like to plot later
  tensiont = zeros(trunc(Int,T/Δt)+1,partition-2)
  vt = zeros(trunc(Int,T/Δt)+1,partition-1)
  xt = zeros(trunc(Int,T/Δt)+1,partition-1)
  αt = zeros(trunc(Int,T/Δt)+1,partition)
  βt = zeros(trunc(Int,T/Δt)+1,partition)
  Mt = zeros(trunc(Int,T/Δt))
  ract = zeros(trunc(Int,T/Δt)+1,partition+1)
  rhot = zeros(trunc(Int,T/Δt)+1,partition+1)
  #We start wriiting down the values of some parameters at time=0
  λ = 0.0
  ρt = zeros(trunc(Int,T/Δt)+1,partition+1)
  ρt[1,:] = get_free_dof_values(uh_ρ)[:]
  ρ0t = zeros(trunc(Int,T/Δt)+1,partition+1)
  ρ0t[1,:] = get_free_dof_values(uh_ρ0)[:]

  #we write down the weak form of the membrane equation
  # to do backward eurler for the time evolwe give gridap a so-called Mass Term for a   
  a(k,x,w) = ∫(k*(∇(x)⋅∇(w)))dΩ
  m(ρ_b,Δt,x,w) = ∫( (χ*ρ_b)*(x*w)/Δt )dΩ
  aₓ(ρ_b,x,w) = m(ρ_b,Δt,x,w) + a(k,x,w)
  bₓ(ρ_b,v,w) = ∫(χ*(ρ_b*v)*w)dΩ
  
  # Now for MCA bound
  # mass term for the temporal evolution de ρ_b
  mρ(Δt,ρ_b,w) = ∫( (ρ_b*w)/Δt )dΩ
  aρ(ρ_b,v,w) = mρ(Δt,ρ_b,w) + ∫(0.001*(∇(ρ_b)⋅∇(w)))dΩ  + ∫(w*(v*(∇(ρ_b)⋅VectorValue(1.0))+ρ_b*(∇(v)⋅VectorValue(1.0))))dΩ + ∫(koff*(ρ_b*w))dΩ +∫(λᵇ*((ρ_b*ρ_b*ρ_b)*w))dΩ
  bρ(w,ρ_u) = ∫(kon*(ρ_u*w)+ λ/L*(kon/(kon+koff)))dΩ + mρ(Δt,uh_ρ_old,w)  

  #Now for MCA unbound
  #mass term per la temporal evolution de ρ_b
  aρ0(ρ_u,x,x_old,w) = mρ(Δt,ρ_u,w) + ∫(kon*(ρ_u*w))dΩ +∫(D*(∇(ρ_u)⋅∇(w)))dΩ  + ∫(w*(((x-x_old)/Δt)*(∇(ρ_u)⋅VectorValue(1.0))+ρ_u*(((∇(x)-∇(x_old))/Δt)⋅VectorValue(1.0))))dΩ  
  bρ0(w,ρ_b) = ∫(koff*(ρ_b*w)+ λ/L*(koff/(kon+koff)))dΩ + mρ(Δt,uh_ρ0_old,w) 

  #Now for v
  aᵥ(ρ_b,v,w) = ∫(η*(∇(v)⋅∇(w)))dΩ + ∫( (χ*ρ_b)*(v*w) )dΩ
  bᵥ(w,uh_rho) = m(uh_ρ,Δt,uh_x,w) - m(uh_ρ,Δt,uh_x_old,w) +  ∫( w*σₐ₀*∇(uh_rho)⋅VectorValue(1.0) )dΩ #no feedback is ∫( w*∇σₐ )dΩ
  Aᵥ(v,w) = aᵥ(uh_ρ,v,w)
  Bᵥ(w) = bᵥ(w,uh_rho)
  #We can now use ρ_b and x to solve v
  op_v= AffineFEOperator(Aᵥ,Bᵥ,V,WD0)

  #SOLVING X AT t=0
  aₓ_0(x,w) = aₓ(0,x,w)
  b_0(w) =  bₓ(0,0,w)
  op_x = AffineFEOperator(aₓ_0,b_0,X,WD0)
  uh_x = solve(op_x)
  uh_x_old =  uh_x

  #SOLVE ρ_b AT t=0 vien initial velocity zero
  Aρ(ρ_b,w) = aρ(ρ_b,uh_v,w) 
  Bρ(w) = bρ(w,uh_ρ0)
  op_ρ= AffineFEOperator(Aρ,Bρ,RHO,Q0)
  uh_ρ=solve(op_ρ)
  uh_ρ_old=uh_ρ

  #SOLVE ρ_u AT t=0
  Aρ0(ρ_u,w) = aρ0(ρ_u,uh_x,uh_x_old,w)
  Bρ0(w) = bρ0(w,uh_ρ_old)
  op_ρ0 = AffineFEOperator(Aρ0,Bρ0,RHO0,Q0)
  uh_ρ0 = solve(op_ρ0)
  uh_ρ0_old = uh_ρ0

  #SOLVING V AT t=0 using x, xold and rho
  uh_v=zero(V)

  #define the production rates as vectors so that we can have heterogeneous α and β due to mechanics  
  α₀v = zeros(partition)
  α₀v[:] .= α₀
  β₀v = zeros(partition)
  β₀v[:] .= β₀

  #we define the tension for a spring
  lx=get_free_dof_values(uh_x)
  x1=circshift(lx,1)
  tension=k*(lx-x1)/Δx 
  splice!(tension,  1) 
  vt[1,:] = get_free_dof_values(uh_v)[:]
  xt[1,:] = get_free_dof_values(uh_x)[:]    
  tensiont[1,:] = tension[:]
  #ten will be used in the computation of the model, we instill the Boundary Conditions and remove negative values
  ten = zeros(partition, 1)
  ten[1:end-2] = tension[:]
  ten[end-1] = tension[end]
  ten[end] = tension[end]

  #storing information of how alpha an beta behave spatially over time
  for j in 1:1:(partition)
    αt[1,j]= α₀v[j] + α*(0.5*(1-tanh((get_free_dof_values(uh_ρ)[j])/rho0-ρth/rho0))) 
  end
  βt[1,:] = β₀v .+ β*threshold(ten,sig0,tenth) 

  #DEFINING the equations for Rac and Rho
  a_rac(rac,w,v) = (1/dᵃ)*mρ(Δt,rac,w)  + ∫(w*rac)dΩ + ∫(mTorc*w*rac*threshold(ten,ten0,10))dΩ # + ∫(1/dᵃ*w*((v)*(∇(rac)⋅VectorValue(1.0))+rac*(∇(v)⋅VectorValue(1.0))))dΩ 
  b_rac(w,rho,ρ_b) =  (1/dᵃ)*mρ(Δt,uh_rac_old,w) + ∫(w*(α₀v/(1+rho*rho) + α*(0.5*(1-tanh∘(ρ_b/rho0-ρth/rho0)))/(1+rho*rho)))dΩ  #*(tanh∘(0.9 /(ρ_b) - 0.9/0.07)+1)
 
  a_rho(rho,w,v) = (1/dᵇ)*mρ(Δt,rho,w)  + ∫(w*rho)dΩ  +∫(λʳᴬ*((rho*rho*rho)*w))dΩ  #+ ∫(1/dᵇ*w*((v)*(∇(rho)⋅VectorValue(1.0))+rho*(∇(v)⋅VectorValue(1.0))))dΩ
  b_rho(w,ten,rac) =  (1/dᵇ)*mρ(Δt,uh_rho_old,w) + ∫( w*(β₀v/(1+rac*rac) + β*threshold(ten,sig0,tenth)/(1+rac*rac)))dΩ  

  #folders to store plots to check whether the run is going well
  # mkpath(pPNG*"alpha=$α₀ beta=$β₀/rac/") 
  # mkpath(pPNG*"alpha=$α₀ beta=$β₀/rho/") 
  # mkpath(pPNG*"alpha=$α₀ beta=$β₀/0rac/") 
  # mkpath(pPNG*"alpha=$α₀ beta=$β₀/0rho/") 
 
  #SOLVE rac AT t=0
  Arac(rac,w) = a_rac(rac,w,uh_v)
  Brac(w) = b_rac(w,uh_rho,uh_ρ)
  op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
  uh_rac = solve(op_rac)
  uh_rac_old = uh_rac
  #SOLVE rho AT t=0
  Arho(rho,w) = a_rho(rho,w,uh_v)
  Brho(w) = b_rho(w,ten,uh_rac)
  op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
  uh_rho = solve(op_rho)
  uh_rho_old = uh_rho

  # Vale hem fet les condicions inicials, ara corre un for loop pel temps, nomes passa aixo pq tenim 
  # un parametre que es dot{x}(i dot rho etc) que per calcular necesitem x and x_old
  i = 0
  t=0
  dummyx0=0
 
 #give steady state as initial conditions for rac and rho 
  for ti in 1:250
    op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac
    op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
    uh_rho = solve(op_rho)
    uh_rho_old = uh_rho

    pa=get_free_dof_values(uh_rac)
    paa=get_free_dof_values(uh_rho)

    # plotting("rac",pa,pPNG*"alpha=$α₀ beta=$β₀/0rac/","$ti")
    # plotting("rho",paa,pPNG*"alpha=$α₀ beta=$β₀/0rho/","$ti")
  end
  #threshold value to start protrusion, 1.3 times the initial condition
  racth = 1.05*get_free_dof_values(uh_rac)[1]

  uh_rac_old = uh_rac
  uh_rho_old = uh_rho

  ract[1,:] = get_free_dof_values(uh_rac)[:]
  rhot[1,:] = get_free_dof_values(uh_rho)[:]
   
  for ti in t₀:Δt:(T-Δt)
    #no opto input in this code
    if t==topto
      α₀v[:] .= α₀ #.+ α₀opto[:]
      β₀v[:] .= β₀ #.+ β₀opto[:]
    end
    if t==topto
      α₀v[:] .= α₀
      β₀v[:] .= β₀ 
    end
    i1 = sum(get_free_dof_values(uh_ρ))
    i2 = sum(get_free_dof_values(uh_ρ0))
    λ = conservation(i2,i1,Minitial)
    i3 = i1+i2
    i4 = trunc(t)
    i = i + 1
    t = t + Δt  

    #boundary conditions for v and x
    rac1=get_free_dof_values(uh_rac)[1]
    vCTErac = vCTE*threshold(rac1,rac0,racth)#velocity polimerization
    v₀ = vCTErac/(1+ten[1]*ten[1]/ten0)
    #we introduce a slight relaxation for the membrane, decreases 2% x at the Boundary condition only
    dummyx0=dummyx0*0.98+Δt*vCTErac/(1+ten[1]*ten[1]/ten0)
    x₀ = dummyx0
    X = TrialFESpace(WD0,[x₀,xₗ])

    @info "Time step $i/$nΔt, time $i4, sum(ρ_b+ρ_u) $i3, and sum(ρ_b)=$i1"
    Mt[i]=i3

    ρt[i+1,:] = get_free_dof_values(uh_ρ)[:]
    ρ0t[i+1,:] = get_free_dof_values(uh_ρ0)[:]
    ract[i+1,:] = get_free_dof_values(uh_rac)[:]
    rhot[i+1,:] = get_free_dof_values(uh_rho)[:]

    # Updating v to solve ρ and x
    A(x,w) = m(uh_ρ,Δt,x,w) + a(k,x,w)
    B(w) = m(uh_ρ,Δt,uh_x,w) + bₓ(uh_ρ,uh_v,w) 
    op_x = AffineFEOperator(A,B,X,WD0)
    uh_x = solve(op_x)


    po=get_free_dof_values(uh_ρ)
    pa=get_free_dof_values(uh_rac)
    paa=get_free_dof_values(uh_rho)
    poo = get_free_dof_values(uh_ρ0)

    # plotting("rac",pa,pPNG*"alpha=$α₀ beta=$β₀/rac/","$i")
    # plotting("rho",paa,pPNG*"alpha=$α₀ beta=$β₀/rho/","$i")


    op_ρ= AffineFEOperator(Aρ,Bρ,RHO,Q0)
    uh_ρ=solve(op_ρ)
    uh_ρ_old=uh_ρ
  
    op_ρ0= AffineFEOperator(Aρ0,Bρ0,RHO0,Q0)
    uh_ρ0=solve(op_ρ0)
    uh_ρ0_old=uh_ρ0
 
    V = TrialFESpace(WD0,[v₀,vₗ])
    op_v= AffineFEOperator(Aᵥ,Bᵥ,V,WD0)
    uh_v=solve(op_v)

    op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac
    
    op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
    uh_rho = solve(op_rho)
    uh_rho_old = uh_rho

    #updating x_old per time derivarive
    uh_x_old =  uh_x 
 
    lx=get_free_dof_values(uh_x)
    x1=circshift(lx,1)
    tension=k*(lx-x1)/Δx
    splice!(tension,  1)
    ten[1:end-2] = tension[:]
    ten[end-1] = tension[end]
    ten[end] = tension[end] 
    vt[i+1,:] = get_free_dof_values(uh_v)[:]
    xt[i+1,:] = get_free_dof_values(uh_x)[:]
    tensiont[i+1,:] = tension[:]
    βt[i+1,:] = β₀v .+ β*threshold(ten,sig0,tenth) 
    for j in 1:1:(partition)
      αt[i+1,j] = α₀v[j] + α*(0.5*(1-tanh((get_free_dof_values(uh_ρ)[j])/rho0-ρth/rho0))) 
    end
    
#   IF polarity is lost, return time taken and plot data
    if t>3*topto
      auxfront = get_free_dof_values(uh_rac)[2]
      auxback = get_free_dof_values(uh_rac)[end-1]
      auxfront2 = get_free_dof_values(uh_rho)[2]
      auxback2 = get_free_dof_values(uh_rho)[end-1]
      if 1.1 > auxfront/auxback || 1.1 >  auxback2/auxfront2 
        plots_heat(i ,α₀,β₀,nΔt,Δx,vt,xt,tensiont,ρt,ρ0t,ract,rhot,Mt,λ⁻²,pPNG,α, β, dᵃ, dᵇ,αt,βt)
        return convert(Float64, t)
      end
    end

  end 
  #otherwise return and plot data when you reach time T
  plots_heat(i ,α₀,β₀,nΔt,Δx,vt,xt,tensiont,ρt,ρ0t,ract,rhot,Mt,λ⁻²,pPNG,α, β, dᵃ, dᵇ,αt,βt) 
  return convert(Float64, t)
end

###############CONSTANT MECHANICAL PARAMETERS ################
partition=100       #size of the FE system
χ = 100               # [pN s/um^3] friction coefficient
η = 10000              # [pN s/ um] 2D viscosity of the cortex
const L = 20.           # [um] System size length
k=100.0            # [pN /um]   elastic constant for the membrane
σₐ₀ =  100        # [pN /um^2] maximum active "pressure"

const koff = 1.4  # ezrin koff [1/s] toff=0.7s Fritzsche et al
const kon  = 5.0  # ezrin kon [1/s] ton =0.2s Fritzsche et al
D = 0.3   # diffusion of the membrane [um^2/s] 0.003  Fritzsche et al
M0 = 10.0 #total amount of ezrin
const wrac=L/2.5 #if we want to start polarized we use this for the Gaussian distribution
#in case one wants nonlinear deactication if the system would be unstable
λᵇ = 0   # for bound ezrin 
λʳᴬ= 0   #for Rho

#adimensional numbers that define the system
τₐ = η/σₐ₀ #for σₐ₀=100
τₑ = η/k #almost always?
Pe⁻¹ = D*η/(σₐ₀*L^2) #D*η/(σₐ₀*L^2) 0.00075 for σₐ₀=100
λ = sqrt( η*L /(χ*L*L*M0) )



rac0=0.01 # ADIMENSINAlizes rac switch for protrusion, changes the slope of the threshold function
vCTE=-0.15  #Scale for the polimerizatio velocity
ten0=10.0 #denominator for protrusion dependence on tension
tenth=5.0 #threshold for tension to activate rho
sig0 = 5.00  #adimensionalizes the tension term for rho
ρth = 0.068 #MCA/ezrin value below which protrusion starts
rho0 = 0.02 #adimensionalizes the MCA concentration term for rac

#Rac Coefficients - α and α₀ ARE DEFINED IN THE heatmap LOOP
dᵃ = 0.04  #deactivation
#α₀ = 0.50  #basal activation
#α =0 #.2*α₀ #activation through mechanics, zero for only local inhibition
Drac = 0.5 #diffusion rate
mTorc = 0.0#1
αopto = 5 # opto input for Rac

#Rho Coefficients - β and β₀ ARE DEFINED IN THE heatmap LOOP
dᵇ = 0.04 #deactivation
#β₀ = 5.50 #basal activation
#β =0.2*β₀ #activation through mechanics, zero for only local inhibition
Drho = 0.5 #diffusion rate 
βopto = 0 # opto input for Rho

#time variables
topto=50 #time to start opto signal
Δt  = 1.0 #timestep
T = 600 #time to finish simualtion
#if we wanna do multiple simulations in a row with different variables
setsv=1
setsx=2
#how many simulaons do you want the heatmap to use? total number is heatsets^2
heatsets=8

#Pe = L^2*σₐ₀/(D*η) 
len=trunc(λ,digits=2) 
const simulation= "Heatmap/sets=$heatsets T=$T ten0=$ten0 db=$dᵇ da=$dᵃ sig0=$sig0 rth=$ρth len=$len/rho0=$rho0 mTorc=$mTorc D=$D te=$τₑ ta=$τₐ M0=$M0 deltat=$Δt vCTE=$vCTE tenth=$tenth lb=$λᵇ/"

pPNG="./PNG/"*simulation
mkpath(pPNG)

xχ = zeros(setsx)
xη = zeros(setsv)
xα₀ = zeros(heatsets)
xβ₀ = zeros(heatsets)
tension = zeros(setsx,setsv,Int32(T/Δt)+1,partition-2)
v = zeros(setsx,setsv,Int32(T/Δt)+1,partition-1)
ρ_b = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
a = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
b = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
tension2 = zeros(setsx,setsv,Int32(T/Δt)+1,partition-2)
v2 = zeros(setsx,setsv,Int32(T/Δt)+1,partition-1)
ρ2 = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
αt = zeros(setsx,setsv,Int32(T/Δt)+1,partition)
βt = zeros(setsx,setsv,Int32(T/Δt)+1,partition)
time2stop = zeros(heatsets,heatsets)

#loop that will do heatsets^2 simulations to make the heatmap
for i in 1:1:heatsets
  for j in 1:1:Int8(heatsets)
    λ⁻² = 1/(λ*λ)
    # alphas and betas explored
    α₀= trunc(0.1 + 0.8*(i-1),digits=3)
    β₀= trunc(0.1 + 0.8*(j-1),digits=3)
    # this would be zero for purely local inhibition, otherwise is the scale of how mechanics afect rac and rho
    α = trunc(0.3*α₀,digits=3)
    β = trunc(0.3*β₀,digits=3)    
    @info "α₀=$α₀  β₀=$β₀    α=$α  β=$β"
    xα₀[i]=α₀
    xβ₀[j]=β₀
    lη = η
    time2stop[i,j] = run(α,β,α₀,β₀,χ,λ⁻²,lη,T,Δt,partition)
  end
end  

output = zeros(heatsets*heatsets,3)

for i in 1:1:heatsets
  for j in 1:1:heatsets
      output[i+(j-1)*heatsets,1] = xα₀[i]
      output[i+(j-1)*heatsets,2] = xβ₀[j]
      output[i+(j-1)*heatsets,3] = time2stop[i,j]
  end 
end

# Writing array to file
writedlm(pPNG*"time2stop.txt", output, ' ')

#plotting the results in a heatmap
myarray=open(readdlm,"a_betavsalpha_diagram.txt")
x = myarray[:,1]
y = myarray[:,2]
heatmap(xβ₀,xα₀,(time2stop.-50), clim=(0, T-50), margin = 5px, right_margin=17px,  size=(250, 220), framestyle = :box)
plot!(x,y, seriestype=:scatter,legend=false ,markerstrokewidth=0,  markercolor = :red,  markersize = 0.6, margin = 1px)
# ylims!(0, Ntot/L)
xlabel!("β₀")
ylabel!("α₀")
ylims!(0, maximum(xα₀))
xlims!(0, maximum(xβ₀))
# plot!(legend=:topright, legendcolumns=3)
savefig(pPNG*"alphabeta.pdf")
#plots_end()
