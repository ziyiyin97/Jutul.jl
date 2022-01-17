#=
Electro-Chemical component
A component with electric potential, concentration and temperature
The different potentials are independent (diagonal onsager matrix),
and conductivity, diffusivity is constant.
=#
using Revise
using Jutul
using MAT
using Plots

ENV["JULIA_DEBUG"] = 0;
struct SourceAtCell
    cell
    src
    function SourceAtCell(cell,src)
        new(cell,src)
    end 
end

function getHalftrans(exported_all,model1, model2, faces1,faces2,cell1,cells2)
    T_all1 = model1["operators"]["T_all"]
    T_all2 = model1["operators"]["T_all"]
    T = 1.0./((1.0./T_all1[faces1])+(1.0./T_all1[faces2]))
    return T
    #N_all = Int64.(exported["G"]["faces"]["neighbors"])
end
##
function make_system(exported,sys,bcfaces,srccells)
    T_all = exported["operators"]["T_all"]
    N_all = Int64.(exported["G"]["faces"]["neighbors"])
    isboundary = (N_all[bcfaces,1].==0) .| (N_all[bcfaces,2].==0)
    @assert all(isboundary)
    #bcfaces = exported_all["model"]["NegativeElectrode"]["CurrentCollector"]["couplingTerm"]["couplingcells"]
    bccells = N_all[bcfaces,1] + N_all[bcfaces,2]
    T_hf   = T_all[bcfaces]
    bcvaluesrc = ones(size(srccells))
    bcvaluephi = ones(size(bccells)).*0.0

    vf = []
    if haskey(exported, "volumeFraction")
        vf = exported["volumeFraction"][:, 1]
    end
    domain = exported_model_to_domain(exported, bc = bccells, b_T_hf = T_hf, vf=vf)
    G = exported["G"]    
    model = SimulationModel(domain, sys, context = DefaultContext())
    parameters = setup_parameters(model)
    parameters[:boundary_currents] = (:BCCharge, :BCMass)

    # State is dict with pressure in each cell
    phi0 = 1.0
    C0 = 1.
    T0 = 298.15
    #D = -0.7e-10 # ???
    # D = - exported_all["model"]["NegativeElectrode"]["ElectrodeActiveComponent"]["InterDiffusionCoefficient"]
    if haskey(exported,"InterDiffusionCoefficient")
        D = - exported["InterDiffusionCoefficient"]
    else
        D =  -0.0
    end
    if isa(exported["EffectiveElectricalConductivity"], Matrix)
        σ = exported["EffectiveElectricalConductivity"][1]
    else
        σ = 1.0
    end
    λ = exported["thermalConductivity"][1]

    S = model.secondary_variables
    S[:BoundaryPhi] = BoundaryPotential{Phi}()
    S[:BoundaryC] = BoundaryPotential{Phi}()
    S[:BoundaryT] = BoundaryPotential{T}()

    S[:BCCharge] = BoundaryCurrent{Charge}(srccells)
    S[:BCMass] = BoundaryCurrent{Mass}(srccells)
    S[:BCEnergy] = BoundaryCurrent{Energy}(srccells)

    init = Dict(
        :Phi                    => phi0,
        :C                      => C0,
        :T                      => T0,
        :Conductivity           => σ,
        :Diffusivity            => D,
        :ThermalConductivity    => λ,
        :BoundaryPhi            => bcvaluephi, 
        :BoundaryC              => bcvaluephi, 
        :BoundaryT              => bcvaluephi,
        :BCCharge               => bcvaluesrc.*0,#0.0227702,
        :BCMass                 => bcvaluesrc,
        :BCEnergy               => bcvaluesrc,
        )

    state0 = setup_state(model, init)
    return model, G, state0, parameters, init
end

##
function convertToIntVector(x::Float64)
    vec = Int64.(
        Vector{Float64}([x])
    )
    return vec
end
 function convertToIntVector(x::Matrix{Float64})
    vec = Int64.(   
    Vector{Float64}(x[:,1])
    )
    return vec
 end

function setup_model(exported_all)
    exported_cc = exported_all["model"]["NegativeElectrode"]["CurrentCollector"];

    sys_cc = CurrentCollector()
    
    bcfaces = convertToIntVector(exported_all["model"]["NegativeElectrode"]["CurrentCollector"]["couplingTerm"]["couplingfaces"])
    srccells = []
    (model_cc, G_cc, state0_cc, parm_cc,init_cc) = make_system(exported_cc, sys_cc, bcfaces, srccells)
    
    sys_nam = Grafite()
    exported_nam = exported_all["model"]["NegativeElectrode"]["ElectrodeActiveComponent"];
    bcfaces=[]
    srccells = []
    (model_nam, G_nam, state0_nam, parm_nam, init_nam) = 
        make_system(exported_nam, sys_nam, bcfaces, srccells)

    sys_elyte = SimpleElyte()
    exported_elyte = exported_all["model"]["Electrolyte"]
    bcfaces=[]
    srccells = []
    (model_elyte, G_elyte, state0_elyte, parm_elyte, init_elyte) = 
        make_system(exported_elyte, sys_elyte, bcfaces, srccells)


    sys_pam = NMC111()
    exported_pam = exported_all["model"]["PositiveElectrode"]["ElectrodeActiveComponent"];
    bcfaces=[]
    srccells = []
    (model_pam, G_pam, state0_pam, parm_pam, init_pam) = 
        make_system(exported_pam,sys_pam,bcfaces,srccells)   
   
    exported_pp = exported_all["model"]["PositiveElectrode"]["CurrentCollector"];
    sys_pp = CurrentCollector()
    bcfaces=[]
    srccells = []
    (model_pp, G_pp, state0_pp, parm_pp,init_pp) = 
    make_system(exported_pp,sys_pp, bcfaces, srccells)

    sys_bpp = CurrentAndVoltageSystem()
    domain_bpp = CurrentAndVoltageDomain()
    model_bpp = SimulationModel(domain_bpp, sys_bpp, context = DefaultContext())
    parm_bpp = setup_parameters(model_bpp)
    parm_bpp[:tolerances][:default] = 1e-8
    # Setup model
    

    groups = nothing
    model = MultiModel(
        (
            CC = model_cc, 
            NAM = model_nam, 
            ELYTE = model_elyte, 
            PAM = model_pam, 
            PP = model_pp,
            BPP = model_bpp
        ), 
        groups = groups)    

    state0 = exported_all["state0"]

    init_cc[:Phi] = state0["NegativeElectrode"]["CurrentCollector"]["phi"][1]           #*0
    init_pp[:Phi] = state0["PositiveElectrode"]["CurrentCollector"]["phi"][1]           #*0
    init_nam[:Phi] = state0["NegativeElectrode"]["ElectrodeActiveComponent"]["phi"][1]  #*0
    init_elyte[:Phi] = state0["Electrolyte"]["phi"][1]                                  #*0
    init_pam[:Phi] = state0["PositiveElectrode"]["ElectrodeActiveComponent"]["phi"][1]  #*0
    init_nam[:C] = state0["NegativeElectrode"]["ElectrodeActiveComponent"]["c"][1] 
    init_pam[:C] = state0["PositiveElectrode"]["ElectrodeActiveComponent"]["c"][1]
    #init_elyte[:C] = state0["Electrolyte"]["cs"][1][1]
    if haskey(state0["Electrolyte"],"cs")
        init_elyte[:C] = state0["Electrolyte"]["cs"][1][1]# for compatibility to old
    else
        init_elyte[:C] = state0["Electrolyte"]["c"][1]
    end
    init_bpp = Dict(:Phi => 1.0)

    init = Dict(
        :CC => init_cc,
        :NAM => init_nam,
        :ELYTE => init_elyte,
        :PAM => init_pam,
        :PP => init_pp,
        :BPP => init_bpp
    )

    state0 = setup_state(model, init)
    parameters = Dict(
        :CC => parm_cc,
        :NAM => parm_nam,
        :ELYTE => parm_elyte,
        :PAM => parm_pam,
        :PP => parm_pp,
        :BPP => parm_bpp
    )

    t1, t2 = exported_all["model"]["Electrolyte"]["sp"]["t"]
    z1, z2 = exported_all["model"]["Electrolyte"]["sp"]["z"]
    tDivz_eff = (t1/z1 + t2/z2)
    parameters[:ELYTE][:t] = tDivz_eff
    parameters[:ELYTE][:z] = 1
    grids = Dict(
        :CC => G_cc,
        :NAM =>G_nam,
        :ELYTE => G_elyte,
        :PAM => G_pam,
        :PP => G_pp
        )

    return model, state0, parameters, grids
end

##

function setup_coupling!(model, exported_all)
    # setup coupling CC <-> NAM charge
    target = Dict( 
        :model => :NAM,
        :equation => :charge_conservation
    )
    source = Dict( 
        :model => :CC,
        :equation => :charge_conservation
        )
    srange = Int64.(
        exported_all["model"]["NegativeElectrode"]["couplingTerm"]["couplingcells"][:,1]
        )
    trange = Int64.(
        exported_all["model"]["NegativeElectrode"]["couplingTerm"]["couplingcells"][:,2]
        )
    intersection = ( srange, trange, Cells(), Cells())
    crosstermtype = InjectiveCrossTerm
    issym = true
    coupling = MultiModelCoupling(source, target, intersection; crosstype = crosstermtype, issym = issym)
    push!(model.couplings, coupling)

    # setup coupling NAM <-> ELYTE charge
    target = Dict( 
        :model => :ELYTE,
        :equation => :charge_conservation
    )
    source = Dict( 
        :model => :NAM, 
        :equation => :charge_conservation
        )

    srange=Int64.(exported_all["model"]["couplingTerms"][1]["couplingcells"][:,1])
    trange=Int64.(exported_all["model"]["couplingTerms"][1]["couplingcells"][:,2])
    intersection = ( srange, trange, Cells(), Cells())
    crosstermtype = InjectiveCrossTerm
    issym = true
    coupling = MultiModelCoupling(source, target, intersection; crosstype = crosstermtype, issym = issym)
    push!(model.couplings, coupling)

    # setup coupling NAM <-> ELYTE mass
    target = Dict( 
        :model => :ELYTE,
        :equation => :mass_conservation
    )
    source = Dict( 
        :model => :NAM,
        :equation => :mass_conservation
        )

    srange=Int64.(exported_all["model"]["couplingTerms"][1]["couplingcells"][:,1])
    trange=Int64.(exported_all["model"]["couplingTerms"][1]["couplingcells"][:,2])
    intersection = ( srange, trange, Cells(), Cells())
    crosstermtype = InjectiveCrossTerm
    issym = true
    coupling = MultiModelCoupling(source, target, intersection; crosstype = crosstermtype, issym = issym)
    push!(model.couplings, coupling)

    # setup coupling PAM <-> ELYTE charge
    target = Dict( 
        :model => :ELYTE,
        :equation => :charge_conservation
        )
    source = Dict( 
        :model => :PAM,
        :equation => :charge_conservation
        )
    srange=Int64.(exported_all["model"]["couplingTerms"][2]["couplingcells"][:,1])
    trange=Int64.(exported_all["model"]["couplingTerms"][2]["couplingcells"][:,2])
    intersection = ( srange, trange, Cells(), Cells())
    crosstermtype = InjectiveCrossTerm
    issym = true
    coupling = MultiModelCoupling(source, target, intersection; crosstype = crosstermtype, issym = issym)
    push!(model.couplings,coupling)

    # setup coupling PAM <-> ELYTE mass
    target = Dict( 
        :model => :ELYTE,
        :equation => :mass_conservation
        )
    source = Dict( 
        :model => :PAM,
        :equation => :mass_conservation
        )
    srange=Int64.(exported_all["model"]["couplingTerms"][2]["couplingcells"][:,1])
    trange=Int64.(exported_all["model"]["couplingTerms"][2]["couplingcells"][:,2])
    intersection = ( srange, trange, Cells(), Cells())
    crosstermtype = InjectiveCrossTerm
    issym = true
    coupling = MultiModelCoupling(source, target, intersection; crosstype = crosstermtype, issym = issym)
    push!(model.couplings, coupling)

    # setup coupling PP <-> PAM charge
    target = Dict( 
        :model => :PAM,
        :equation => :charge_conservation
        )
    source = Dict( 
        :model => :PP,
        :equation => :charge_conservation
        )
    srange = Int64.(
        exported_all["model"]["PositiveElectrode"]["couplingTerm"]["couplingcells"][:,1]
        )
    trange = Int64.(
        exported_all["model"]["PositiveElectrode"]["couplingTerm"]["couplingcells"][:,2]
        )
    intersection = (srange, trange, Cells(), Cells())
    crosstermtype = InjectiveCrossTerm
    issym = true
    coupling = MultiModelCoupling(source,target, intersection; crosstype = crosstermtype, issym = issym)
    push!(model.couplings, coupling)
    
    #setup coupling PP <-> PAM charge
    target = Dict( 
        :model => :PP,
        :equation => :charge_conservation
        )
    source = Dict( 
        :model => :BPP,
        :equation => :current_equation
        )
    trange = convertToIntVector(
            exported_all["model"]["PositiveElectrode"]["CurrentCollector"]["couplingTerm"]["couplingcells"]
        )    
    srange = Int64.(ones(size(trange)))
        
    #T_all = exported["operators"]["T_all"]
    #N_all = Int64.(exported["G"]["faces"]["neighbors"])
    #isboundary = (N_all[bcfaces,1].==0) .| (N_all[bcfaces,2].==0)
    #@assert all(isboundary)
    #bcfaces = exported_all["model"]["NegativeElectrode"]["CurrentCollector"]["couplingTerm"]["couplingcells"]
    #bccells = N_all[bcfaces,1] + N_all[bcfaces,2]
    #T_hf   = T_all[bcfaces]
   
    intersection = (srange, trange, Cells(), Cells())
    crosstermtype = InjectiveCrossTerm
    issym = true
    coupling = MultiModelCoupling(source,target, intersection; crosstype = crosstermtype, issym = issym)
    
    push!(model.couplings, coupling)

    # source = Dict( 
    #     :model => :PP,
    #     :equation => :charge_conservation
    #     )
    # target = Dict( 
    #     :model => :BPP,
    #     :equation => :current_equation
    #     )
    # srange = Int64.(
    #     Vector{Float64}([10.0])
    #     )
    # trange = Int64.(
    #     Vector{Float64}([1])
    #     )
    # intersection = (srange, trange, Cells(), Cells())
    # crosstermtype = InjectiveCrossTerm
    # issym = false
    # coupling = MultiModelCoupling(source,target, intersection; crosstype = crosstermtype, issym = issym)
    
    
    #push!(model.couplings, coupling)
end


##
function currentFun(t,inputI)
    #inputI = 9.4575
    tup = 0.1
    if ( t<= tup)
        val = sineup(0, inputI, 0, tup, t) 
    else
        val = inputI;
    end
end

function test_battery(name)
##
    #name="model1d_notemp"
    fn = string(dirname(pathof(Jutul)), "/../data/models/", name, ".mat")
    exported_all = MAT.matread(fn)

    model, state0, parameters, grids = setup_model(exported_all)    
    setup_coupling!(model, exported_all)
    #inputI = 9.4575
    inputI = exported_all["states"][10]["PositiveElectrode"]["CurrentCollector"]["I"]
    cFun(time) = currentFun(time, inputI)
    forces_pp = nothing 
    #forces_pp = (src = SourceAtCell(10,9.4575*0.0),)
    currents = Dict( :current => cFun)
    forces = Dict(
        :CC => nothing,
        :NAM => nothing,
        :ELYTE => nothing,
        :PAM => nothing,
        :PP => forces_pp,
        :BPP => currents
    )
    for (k, p) in parameters
        p[:tolerances][:default] = 1e-3
    end
    #parameters[:PP][:charge_conservation]=1e-3
    sim = Simulator(model, state0 = state0, parameters = parameters, copy_state = true)
    #timesteps = exported_all["schedule"]["step"]["val"][1:27]
    timesteps = exported_all["schedule"]["step"]["val"][1:29]
    cfg = simulator_config(sim)
    cfg = simulator_config(sim)
    cfg[:linear_solver] = nothing
    cfg[:info_level] = 2
    cfg[:debug_level] = 0
    cfg[:max_residual] = 1e20
    cfg[:min_nonlinear_iterations] = 3
    cfg[:extra_timing] = true
##

    states, report = simulate(sim, timesteps, forces = forces, config = cfg)
    stateref = exported_all["states"]

    return states, grids, state0, stateref, parameters, exported_all, model, timesteps, cfg, report
end

##
name="model1d_notemp"
name="model1D_50"
#name="model1D_500"
#name="model1D_5000"
#name="model2D_1100"
#name="sector_1656"
#name="model3D_492"#.mat
#name="sector_1656"
#name="model3D_3936"
#name="sector_55200" #Tobig for direct linear_solver
#name="spiral_16560"
states, grids, state0, stateref, parameters, exported_all, model, timesteps, cfg, report = test_battery(name);
steps = size(states, 1)
E = Matrix{Float64}(undef,steps,2)
for step in 1:steps
    phi = states[step][:PP][:Phi][10]
    E[step,1] = phi
    phi_ref = stateref[step]["PositiveElectrode"]["ElectrodeActiveComponent"]["phi"][10]
    E[step,2] = phi_ref
end
plot1 = Plots.plot([], []; title = "E", size=(1000, 800))
plot!(plot1,cumsum(timesteps),E)
closeall()
display(plot1)
##
error()
using Plots

function plot_phi()
    plot1 = Plots.plot([], []; title = "Phi", size=(1000, 800))

    p = plot!(plot1, legend = false)
    submodels = (:CC, :NAM, :ELYTE, :PAM, :PP)
    # submodels = (:NAM, :ELYTE, :PAM)
    # submodels = (:CC, :NAM)
    # submodels = (:ELYTE,)
    # submodels = (:PP, :PAM)

    var = :Phi
    steps = size(states, 1)
    for i in 1:steps
        for mod in submodels
            x = grids[mod]["cells"]["centroids"]
            plot!(plot1, x, states[i][mod][var], lw=2, color=RGBA(0.5, 0.5, 0.5, 0.5))
        end
    end
    return plot1
end
plot1 = plot_phi()
closeall()
display(plot1)

##

#for (i, dt) in enumerate(timesteps)
for i in 1:15   
    refstep=i
    sim_step=i

    p1 = Plots.plot(title="Phi", size=(1000, 800))
    p2 = Plots.plot(title="Flux", size=(1000, 800))
    p3 = Plots.plot(title="C", size=(1000, 800))

    fields = ["CurrentCollector","ElectrodeActiveComponent"]
    components = ["NegativeElectrode","PositiveElectrode"]
    #components = ["NegativeElectrode"]
    #components = ["PositiveElectrode"]
    #components = []
    for component = components
        for field in fields
            G = exported_all["model"][component][field]["G"]
            x = G["cells"]["centroids"]
            xf= G["faces"]["centroids"][end]
            xfi= G["faces"]["centroids"][2:end-1]

            state = stateref[refstep][component]
            phi_ref = state[field]["phi"]
            j_ref = state[field]["j"]

            Plots.plot!(p1,x,phi_ref;linecolor="red")
            Plots.plot!(p2,xfi,j_ref;linecolor="red")
            if haskey(state[field],"c")
                c = state[field]["c"]
                Plots.plot!(p3,x,c;linecolor="red")
            end
        end
    end

    fields = [] 
    fields = ["Electrolyte"]

    for field in fields
        G = exported_all["model"][field]["G"]
        x = G["cells"]["centroids"]
        xf= G["faces"]["centroids"][end]
        xfi= G["faces"]["centroids"][2:end-1]

        state = stateref[refstep]
        phi_ref = state[field]["phi"]
        j_ref = state[field]["j"]

        Plots.plot!(p1,x,phi_ref;linecolor="red")
        Plots.plot!(p2,xfi,j_ref;linecolor="red")
        if haskey(state[field],"cs")
            c = state[field]["cs"][1]
            Plots.plot!(p3,x,c;linecolor="red")
        end
    end

    ##


    mykeys = [:CC, :NAM] # :ELYTE]
    mykeys = [:PP, :PAM]
    #mykeys = [:ELYTE]
    mykeys =  keys(grids)
    for key in mykeys
        G = grids[key]
        x = G["cells"]["centroids"]
        xf= G["faces"]["centroids"][end]
        xfi= G["faces"]["centroids"][2:end-1]     
        p = plot(p1, p2, layout = (1, 2), legend = false)
        phi = states[sim_step][key][:Phi]
        Plots.plot!(
            p1, x, phi; markershape=:circle, linestyle=:dot, seriestype = :scatter
            )
        
        if haskey(states[sim_step][key], :TotalCurrent)
            j = states[sim_step][key][:TotalCurrent][1:2:end-1]
        else
            j = -states[sim_step][key][:TPkGrad_Phi][1:2:end-1]
        end
        
        Plots.plot!(p2, xfi, j; markershape=:circle,linestyle=:dot, seriestype = :scatter)
        if(haskey(states[sim_step][key], :C))
            cc = states[sim_step][key][:C]
            Plots.plot!(p3, x, cc; markershape=:circle, linestyle=:dot, seriestype = :scatter)
        end
    end

    display(plot!(p1, p2, p3,layout = (3, 1), legend = false))
end
error()
##
E = Matrix{Float64}(undef,27,2)
for step in 1:27
    phi = states[step][:PP][:Phi][10]
    E[step,1] = phi
    phi_ref = stateref[step]["PositiveElectrode"]["ElectrodeActiveComponent"]["phi"][10]
    E[step,2] = phi_ref
end

##
function print_diff_j(s, sref, n)
    k = :TotalCurrent
    if haskey(s[n], k)
        Δ = abs.(1 .+ (sref[n][k]) ./ s[n][k][2:2:end])
    else
        Δ = abs.(1 .+ (-sref[n][k]) ./ s[n][:TPkGrad_Phi][2:2:end])
    end
    println("k = $k, n = $n")
    println("rel.diff = $(maximum(Δ))")
end
##
EAC = "ElectrodeActiveComponent"
PE = "PositiveElectrode"
NE = "NegativeElectrode"

# Transelate bw states from matlab and julia
j2m = Dict{Symbol, String}(
    # :C                  => "cs",
    :T                  => "T",
    :Phi                => "phi", 
    :Conductivity       => "conductivity",
    :Diffusivity        => "D",
    :TotalCurrent       => "j",
    :ChargeCarrierFlux  => "LiFlux", # This is not correct - CC is more than Li
    :ELYTE              => "Electrolyte"
)
m2j = Dict(value => key for (key, value) in j2m)

rs = stateref[:, 1]
rs_elyte = [s[j2m[:ELYTE]] for s in rs];
rs_pam = [s[PE][EAC] for s in rs];

states_pam = [s[:PAM] for s in states];
states_elyte = [s[:ELYTE] for s in states];
##

states_comp = states_elyte
ref_states = get_ref_states(j2m, rs_elyte);
for (n, state) in enumerate(states_comp)
    print_diff_j(states_comp, ref_states, n)
end

plot(E)
##
i=4;plot([states[i][:ELYTE][:C], stateref[i]["Electrolyte"]["cs"][:,1]])
##
i=10;plot([states[i][:PAM][:C], stateref[i]["Electrolyte"]["cs"][:,1]])
i=10;plot([states[i][:NAM][:C], stateref[i]["NegativeElectrode"]["ElectrodeActiveComponent"]["c"][:,1]])
i=10;plot([states[i][:PAM][:C], stateref[i]["PositiveElectrode"]["ElectrodeActiveComponent"]["c"][:,1]])

i=10;plot([states[i][:CC][:Phi],stateref[i]["NegativeElectrode"]["CurrentCollector"]["phi"]])
i=10;plot([states[i][:NAM][:Phi],stateref[i]["NegativeElectrode"]["ElectrodeActiveComponent"]["phi"]])
i=10;plot([states[i][:ELYTE][:Phi], stateref[i]["Electrolyte"]["phi"]])