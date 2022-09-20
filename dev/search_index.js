var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Jutul","category":"page"},{"location":"#Jutul","page":"Home","title":"Jutul","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Jutul.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Jutul]","category":"page"},{"location":"#Jutul.AMGPreconditioner","page":"Home","title":"Jutul.AMGPreconditioner","text":"AMG on CPU (Julia native)\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.BlockMajorLayout","page":"Home","title":"Jutul.BlockMajorLayout","text":"Same as EntityMajorLayout, but the nzval is a matrix\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.CartesianMesh","page":"Home","title":"Jutul.CartesianMesh","text":"CartesianMesh(dims, [Δ, [origin]])\n\nCreate a Cartesian mesh with dimensions specified by the Tuple dims.\n\nArguments\n\ndims::Tuple: Number of grid cells in each direction. For example, (nx, ny) will give a 2D grids with nx cells in the x-direction.\nΔ::Tuple=ones(length(dims)): Equal length to dims. First option: A Tuple of scalars where each entry is the length of each cell in that direction. For\n\nexample, specifying (Δx, Δy) for a uniform grid with each grid cell having area ofΔx*Δy. Second option: ATuple` of vectors where each entry contains the cell sizes in the direction.\n\norigin=zeros(length(dims)): The origin of the first corner in the grid.\n\nExamples\n\nGenerate a uniform 3D mesh that discretizes a domain of 2 by 3 by 5 units with 3 by 5 by 2 cells:\n\njulia> CartesianMesh((3, 5, 2), (2.0, 3.0, 5.0))\nCartesianMesh (3D) with 3x5x2=30 cells\n\nGenerate a non-uniform 2D mesh:\n\njulia> CartesianMesh((2, 3), ([1.0, 2.0], [0.1, 3.0, 2.5]))\nCartesianMesh (3D) with 3x5x2=30 cells\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.CompactAutoDiffCache","page":"Home","title":"Jutul.CompactAutoDiffCache","text":"Cache that holds an AD vector/matrix together with their positions.\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.ConstantVariables","page":"Home","title":"Jutul.ConstantVariables","text":"A set of constants, repeated over the entire set of Cells or some other entity\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.DampedJacobiPreconditioner","page":"Home","title":"Jutul.DampedJacobiPreconditioner","text":"Damped Jacobi preconditioner on CPU\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.DefaultContext","page":"Home","title":"Jutul.DefaultContext","text":"Default context\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.EntityMajorLayout","page":"Home","title":"Jutul.EntityMajorLayout","text":"Domain entities sequentially in rows:\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.EquationMajorLayout","page":"Home","title":"Jutul.EquationMajorLayout","text":"Equations are stored sequentially in rows, derivatives of same type in columns:\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.FlowDiscretization","page":"Home","title":"Jutul.FlowDiscretization","text":"Discretization of kgradp + upwind\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.GroupWisePreconditioner","page":"Home","title":"Jutul.GroupWisePreconditioner","text":"Multi-model preconditioners\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.ILUZeroPreconditioner","page":"Home","title":"Jutul.ILUZeroPreconditioner","text":"ILU(0) preconditioner on CPU\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.JutulAutoDiffCache","page":"Home","title":"Jutul.JutulAutoDiffCache","text":"An AutoDiffCache is a type that holds both a set of AD values and a map into some global Jacobian.\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.JutulDiscretization","page":"Home","title":"Jutul.JutulDiscretization","text":"Ask discretization for entry i for specific entity\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.LUPreconditioner","page":"Home","title":"Jutul.LUPreconditioner","text":"Full LU factorization as preconditioner (intended for smaller subsystems)\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.ParallelCSRContext","page":"Home","title":"Jutul.ParallelCSRContext","text":"Context that uses threads etc to accelerate loops\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.SPU","page":"Home","title":"Jutul.SPU","text":"Single-point upwinding.\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.TPFA","page":"Home","title":"Jutul.TPFA","text":"Two-point flux approximation.\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.TwoPointFiniteVolumeGeometry","page":"Home","title":"Jutul.TwoPointFiniteVolumeGeometry","text":"TwoPointFiniteVolumeGeometry(neighbors, areas, volumes, normals, cell_centers, face_centers)\n\nStore two-point geometry information for a given list of neighbors specified as a 2 by n matrix where n is the number of faces such that face i connectes cells N[1, i] and N[2, i].\n\nThe two-point finite-volume geometry contains the minimal set of geometry information required to compute standard finite-volume discretizations.\n\n\n\n\n\n","category":"type"},{"location":"#Jutul.UnaryTabulatedVariable","page":"Home","title":"Jutul.UnaryTabulatedVariable","text":"\n\n\n\n","category":"type"},{"location":"#Jutul.absolute_increment_limit-Tuple{JutulVariables}","page":"Home","title":"Jutul.absolute_increment_limit","text":"Absolute allowable change for variable during a nonlinear update.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.align_to_jacobian!-NTuple{4, Any}","page":"Home","title":"Jutul.align_to_jacobian!","text":"Update an equation so that it knows where to store its derivatives in the Jacobian representation.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.align_to_jacobian!-Tuple{Jutul.ConservationLawTPFAStorage, ConservationLaw, Any, Any, Cells}","page":"Home","title":"Jutul.align_to_jacobian!","text":"Update positions of law's derivatives in global Jacobian\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.allocate_array_ad-Tuple{AbstractMatrix}","page":"Home","title":"Jutul.allocate_array_ad","text":"allocate_array_ad(v::AbstractMatrix, ...)\n\nConvert matrix to AD matrix.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.allocate_array_ad-Tuple{AbstractVector}","page":"Home","title":"Jutul.allocate_array_ad","text":"allocate_array_ad(v::AbstractVector, ...)\n\nConvert vector to AD vector.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.allocate_array_ad-Union{Tuple{Vararg{R}}, Tuple{R}} where R<:Integer","page":"Home","title":"Jutul.allocate_array_ad","text":"allocate_array_ad(n[, m]; <keyword arguments>)\n\nAllocate vector or matrix as AD with optionally provided context and a specified non-zero on the diagonal.\n\nArguments\n\nn::Integer: number of entries in vector, or number of rows if m is given.\nm::Integer: number of rows (optional)\n\nKeyword arguments\n\nnpartials = 1: Number of partials derivatives to allocate for each element\ndiag_pos = nothing: Indices of where to put entities on the diagonal (if any)\n\nOther keyword arguments are passed onto get_ad_entity_scalar.\n\nExamples:\n\nAllocate a vector with a single partial:\n\njulia> allocate_array_ad(2)\n2-element Vector{ForwardDiff.Dual{nothing, Float64, 1}}:\n Dual{nothing}(0.0,0.0)\n Dual{nothing}(0.0,0.0)\n\nAllocate a vector with two partials, and set the first to one:\n\njulia> allocate_array_ad(2, diag_pos = 1, npartials = 2)\n2-element Vector{ForwardDiff.Dual{nothing, Float64, 2}}:\n Dual{nothing}(0.0,1.0,0.0)\n Dual{nothing}(0.0,1.0,0.0)\n\nSet up a matrix with two partials, where the first column has partials [1, 0] and the second [0, 1]:\n\njulia> allocate_array_ad(2, 2, diag_pos = [1, 2], npartials = 2)\n2×2 Matrix{ForwardDiff.Dual{nothing, Float64, 2}}:\n Dual{nothing}(0.0,1.0,0.0)  Dual{nothing}(0.0,1.0,0.0)\n Dual{nothing}(0.0,0.0,1.0)  Dual{nothing}(0.0,0.0,1.0)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.apply_forces!-NTuple{4, Any}","page":"Home","title":"Jutul.apply_forces!","text":"Apply a set of forces to all equations. Equations that don't support a given force will just ignore them, thanks to the power of multiple dispatch.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.apply_forces_to_equation!-NTuple{7, Any}","page":"Home","title":"Jutul.apply_forces_to_equation!","text":"Update an equation with the effect of a force. The default behavior for any force we do not know about is to assume that the force does not impact this particular equation.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.as_value-Union{Tuple{AbstractArray{D}}, Tuple{D}} where D<:ForwardDiff.Dual","page":"Home","title":"Jutul.as_value","text":"Create a mapped array that produces only the values when indexed.\n\nOnly useful for AD arrays, otherwise it does nothing.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.associated_entity-Tuple{JutulEquation}","page":"Home","title":"Jutul.associated_entity","text":"Return the domain entity the equation is associated with\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.associated_entity-Tuple{JutulVariables}","page":"Home","title":"Jutul.associated_entity","text":"The entity a variable is associated with, and can hold partial derivatives with respect to.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.cell_dims-Tuple{Any, Any}","page":"Home","title":"Jutul.cell_dims","text":"cell_dims(g, pos)::Tuple\n\nGet physical cell dimensions of cell with index pos for grid g.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.cell_index-Tuple{Any, Tuple}","page":"Home","title":"Jutul.cell_index","text":"cell_index(g, pos)\n\nGet linear (scalar) index of mesh cell from provided IJK tuple pos.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.convert_state_ad","page":"Home","title":"Jutul.convert_state_ad","text":"Convert a state containing variables as arrays of doubles to a state where those arrays contain the same value as Dual types. The dual type is currently taken from ForwardDiff.\n\n\n\n\n\n","category":"function"},{"location":"#Jutul.coord_offset-Tuple{Any, AbstractFloat}","page":"Home","title":"Jutul.coord_offset","text":"Lower corner for one dimension, without any transforms applied\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.declare_pattern-Tuple{Any, Any, Any, Any, Vararg{Any}}","page":"Home","title":"Jutul.declare_pattern","text":"Give out source, target arrays of equal length for a given equation attached to the given model.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.declare_sparsity","page":"Home","title":"Jutul.declare_sparsity","text":"Give out I, J arrays of equal length for a given equation attached to the given model.\n\n\n\n\n\n","category":"function"},{"location":"#Jutul.degrees_of_freedom_per_entity-Tuple{Any, ConstantVariables}","page":"Home","title":"Jutul.degrees_of_freedom_per_entity","text":"Constant variables hold no degrees of freedom.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.degrees_of_freedom_per_entity-Tuple{Any, ScalarVariable}","page":"Home","title":"Jutul.degrees_of_freedom_per_entity","text":"Number of independent primary variables / degrees of freedom per computational entity.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.dim-Tuple{AbstractJutulMesh}","page":"Home","title":"Jutul.dim","text":"dim(g)::Integer\n\nGet the dimension of a mesh.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.full_cell-Tuple{Any, Jutul.TrivialGlobalMap}","page":"Home","title":"Jutul.full_cell","text":"Inner cell to local cell (full set)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.get_ad_entity_scalar-Union{Tuple{T}, Tuple{T, Any}, Tuple{T, Any, Any}} where T<:Real","page":"Home","title":"Jutul.get_ad_entity_scalar","text":"get_ad_entity_scalar(v::Real, npartials, diag_pos = nothing; <keyword_arguments>)\n\nGet scalar with partial derivatives as AD instance.\n\nArguments\n\nv::Real: Value of AD variable.\nnpartials: Number of partial derivatives each AD instance holds.\ndiag_pos = nothing: Position(s) of where to set 1 as the partial derivative instead of zero.\n\nKeyword arguments\n\ntag = nothing: Tag for AD instance. Two AD values of the different tag cannot interoperate to avoid perturbation confusion (see ForwardDiff documentation).\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.get_dependencies-Tuple{Any, Any}","page":"Home","title":"Jutul.get_dependencies","text":"Get dependencies of variable when viewed as a secondary variable. Normally autogenerated with @jutul_secondary\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.get_entity_tag-Tuple{Any, Any}","page":"Home","title":"Jutul.get_entity_tag","text":"Combine a base tag (which can be nothing) with a entity to get a tag that captures base tag + entity tag for use with AD initialization.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.get_entries-Tuple{CompactAutoDiffCache}","page":"Home","title":"Jutul.get_entries","text":"Get entries of autodiff cache. Entries are AD vectors that hold values and derivatives.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.get_entries-Tuple{JutulEquation}","page":"Home","title":"Jutul.get_entries","text":"Get the entries of the main autodiff cache for an equation.\n\nNote: This only gets the .equation field's entries.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.get_primary_variable_ordered_entities-Tuple{SimulationModel}","page":"Home","title":"Jutul.get_primary_variable_ordered_entities","text":"Get only the entities where primary variables are present, sorted by their order in the primary variables.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.global_cell-Tuple{Any, Jutul.TrivialGlobalMap}","page":"Home","title":"Jutul.global_cell","text":"Local cell -> global cell (full set)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.global_face-Tuple{Any, Jutul.TrivialGlobalMap}","page":"Home","title":"Jutul.global_face","text":"Local face -> global face (full set)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.initialize_context!-NTuple{4, Any}","page":"Home","title":"Jutul.initialize_context!","text":"Initialize context when setting up a model\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.initialize_extra_state_fields!-Tuple{Any, JutulModel}","page":"Home","title":"Jutul.initialize_extra_state_fields!","text":"Add variables that need to be in state, but are never AD variables (e.g. phase status flag)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.initialize_storage!-Tuple{Any, JutulModel}","page":"Home","title":"Jutul.initialize_storage!","text":"Initialize the already allocated storage at the beginning of a simulation. Use this to e.g. set up extra stuff in state0 needed for initializing the simulation loop.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.initialize_variable_value!-Tuple{Any, Any, GroupedVariables, Symbol, AbstractVector}","page":"Home","title":"Jutul.initialize_variable_value!","text":"Initializer for the value of non-scalar primary variables\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.interior_cell-Tuple{Any, Jutul.TrivialGlobalMap}","page":"Home","title":"Jutul.interior_cell","text":"Local cell in full set -> inner cell (or zero)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.local_ad-Tuple{Any, Any, Any}","page":"Home","title":"Jutul.local_ad","text":"local_ad(state::T, index::I, ad_tag::∂T) where {T, I<:Integer, ∂T}\n\nCreate localad for state for index I of AD tag of type adtag\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.local_cell-Tuple{Any, Jutul.TrivialGlobalMap}","page":"Home","title":"Jutul.local_cell","text":"Global cell -> local cell (full set)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.local_face-Tuple{Any, Jutul.TrivialGlobalMap}","page":"Home","title":"Jutul.local_face","text":"Global face -> local face (full set)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.maximum_value-Tuple{JutulVariables}","page":"Home","title":"Jutul.maximum_value","text":"Upper (inclusive) limit for variable.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.minimum_value-Tuple{JutulVariables}","page":"Home","title":"Jutul.minimum_value","text":"Lower (inclusive) limit for variable.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_cells-Tuple{AbstractJutulMesh}","page":"Home","title":"Jutul.number_of_cells","text":"number_of_cells(g)::Integer\n\nGet the number of cells in a mesh.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_degrees_of_freedom-Tuple{JutulModel}","page":"Home","title":"Jutul.number_of_degrees_of_freedom","text":"Total number of degrees of freedom for a model, over all primary variables and all entities.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_entities-Tuple{Any, JutulEquation}","page":"Home","title":"Jutul.number_of_entities","text":"Get the number of entities (e.g. the number of cells) that the equation is defined on.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_entities-Tuple{Any, JutulVariables}","page":"Home","title":"Jutul.number_of_entities","text":"Number of entities (e.g. Cells, Faces) a variable is defined on. By default, each primary variable exists on all cells of a discretized domain\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_entities-Tuple{JutulAutoDiffCache}","page":"Home","title":"Jutul.number_of_entities","text":"Get number of entities a cache is defined on.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_entities-Tuple{T} where T<:(AbstractVector)","page":"Home","title":"Jutul.number_of_entities","text":"Number of entities for vector stored in state (just the number of elements)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_entities-Tuple{T} where T<:AbstractArray","page":"Home","title":"Jutul.number_of_entities","text":"Number of entities for matrix stored in state (convention is number of columns)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_equations-Tuple{Any, JutulEquation}","page":"Home","title":"Jutul.number_of_equations","text":"Get the total number of equations on the domain of model.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_equations_per_entity-Tuple{JutulModel, JutulEquation}","page":"Home","title":"Jutul.number_of_equations_per_entity","text":"Get the number of equations per entity. For example, mass balance of two components will have two equations per grid cell (= entity)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_faces-Tuple{Any}","page":"Home","title":"Jutul.number_of_faces","text":"number_of_faces(g)::Integer\n\nGet the number of faces in a mesh.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_partials_per_entity-Tuple{JutulEquation}","page":"Home","title":"Jutul.number_of_partials_per_entity","text":"Get the number of partials\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.number_of_values","page":"Home","title":"Jutul.number_of_values","text":"Total number of values for a model, for a given type of variables over all entities\n\n\n\n\n\n","category":"function"},{"location":"#Jutul.read_results-Tuple{Any}","page":"Home","title":"Jutul.read_results","text":"Read results from a given outputpath provded to simulate or simulatorconfig\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.relative_increment_limit-Tuple{JutulVariables}","page":"Home","title":"Jutul.relative_increment_limit","text":"Relative allowable change for variable during a nonlinear update. A variable with value |x| and relative limit 0.2 cannot change more than |x|*0.2.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.replace_variables!-Tuple{Any}","page":"Home","title":"Jutul.replace_variables!","text":"Replace a variable that already exists (either primary or secondary)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.set_parameters!-Tuple{Any}","page":"Home","title":"Jutul.set_parameters!","text":"Set a parameter (adding if it does not exist)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.set_primary_variables!-Tuple{Any}","page":"Home","title":"Jutul.set_primary_variables!","text":"Set a primary variable (adding if it does not exist)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.set_secondary_variables!-Tuple{Any}","page":"Home","title":"Jutul.set_secondary_variables!","text":"Set a secondary variable (adding if it does not exist)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.setup_adjoint_storage-Tuple{Any}","page":"Home","title":"Jutul.setup_adjoint_storage","text":"setup_adjoint_storage(model; state0 = setup_state(model), parameters = setup_parameters(model))\n\nSet up storage for use with solve_adjoint_sensitivities!.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.setup_state!","page":"Home","title":"Jutul.setup_state!","text":"Initialize primary variables and other state fields, given initial values as a Dict\n\n\n\n\n\n","category":"function"},{"location":"#Jutul.setup_state-Tuple{JutulModel, Vararg{Any}}","page":"Home","title":"Jutul.setup_state","text":"Set up a state. You likely want to overload setup_state! instead of this one.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.setup_storage!-Tuple{Any, JutulModel}","page":"Home","title":"Jutul.setup_storage!","text":"Allocate storage for a given model. The storage consists of all dynamic quantities used in the simulation. The default implementation allocates properties, equations and linearized system.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.setup_storage-Tuple{JutulModel}","page":"Home","title":"Jutul.setup_storage","text":"Allocate storage for the model. You should overload setup_storage! if you have a custom definition.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.solve_adjoint_sensitivities!-NTuple{6, Any}","page":"Home","title":"Jutul.solve_adjoint_sensitivities!","text":"solve_adjoint_sensitivities!(∇G, storage, states, state0, timesteps, G; forces = setup_forces(model))\n\nNon-allocating version of solve_adjoint_sensitivities.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.solve_adjoint_sensitivities-NTuple{4, Any}","page":"Home","title":"Jutul.solve_adjoint_sensitivities","text":"solve_adjoint_sensitivities(model, states, reports, G; extra_timing = false, state0 = setup_state(model), forces = setup_forces(model), raw_output = false, kwarg...)\n\nCompute sensitivities of model parameter with name target for objective function G.\n\nSolves the adjoint equations: For model equations F the gradient with respect to parameters p is     ∇ₚG = Σₙ (∂Fₙ / ∂p)ᵀ λₙ where n ∈ [1, N]. Given Lagrange multipliers λₙ from the adjoint equations     (∂Fₙ / ∂xₙ)ᵀ λₙ = - (∂J / ∂xₙ)ᵀ - (∂Fₙ₊₁ / ∂xₙ)ᵀ λₙ₊₁ where the last term is omitted for step n = N and G is the objective function.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.solve_numerical_sensitivities-NTuple{5, Any}","page":"Home","title":"Jutul.solve_numerical_sensitivities","text":"solve_numerical_sensitivities(model, states, reports, G, target;\n                                            forces = setup_forces(model),\n                                            state0 = setup_state(model),\n                                            parameters = setup_parameters(model),\n                                            epsilon = 1e-8)\n\nCompute sensitivities of model parameter with name target for objective function G.\n\nThis method uses numerical perturbation and is primarily intended for testing of solve_adjoint_sensitivities.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.spu_upwind-Union{Tuple{I}, Tuple{R}, Tuple{I, I, R, AbstractArray{R}}} where {R<:Real, I<:Integer}","page":"Home","title":"Jutul.spu_upwind","text":"Perform single-point upwinding based on signed potential.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.spu_upwind_mult-NTuple{4, Any}","page":"Home","title":"Jutul.spu_upwind_mult","text":"Perform single-point upwinding based on signed potential, then multiply the result with that potential\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.spu_upwind_mult_index-NTuple{4, Any}","page":"Home","title":"Jutul.spu_upwind_mult_index","text":"Perform single-point upwinding based on signed potential, then multiply the result with that potential\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.synchronize-Tuple{JutulContext}","page":"Home","title":"Jutul.synchronize","text":"Synchronize backend after allocations.\n\nSome backends may require notification that storage has been allocated.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.tpfv_geometry","page":"Home","title":"Jutul.tpfv_geometry","text":"tpfv_geometry(g)\n\nGenerate two-point finite-volume geometry for a given grid, if supported.\n\nSee also TwoPointFiniteVolumeGeometry.\n\n\n\n\n\n","category":"function"},{"location":"#Jutul.transfer-Tuple{Any, Any}","page":"Home","title":"Jutul.transfer","text":"Transfer v to the representation expected by a given context.\n\nFor the defalt context, the transfer function does nothing. For other context such as the CUDA version, it may convert integers and floats to other types (e.g. Float32) and Arrays to CuArrays.\n\nYou will likely have to implement some transfer operators for your own types if you want to simulate with a non-default context.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.two_point_potential_drop-NTuple{5, Real}","page":"Home","title":"Jutul.two_point_potential_drop","text":"Two-point potential drop with gravity (generic)\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.two_point_potential_drop_half_face-Tuple{Any, Any, AbstractVector, Any, Any}","page":"Home","title":"Jutul.two_point_potential_drop_half_face","text":"Two-point potential drop (with derivatives only respect to \"c_self\")\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.update_before_step!-NTuple{4, Any}","page":"Home","title":"Jutul.update_before_step!","text":"\n\n\n\n","category":"method"},{"location":"#Jutul.update_equation!-Tuple{Any, JutulEquation, Any, Any, Any}","page":"Home","title":"Jutul.update_equation!","text":"Update equation based on currently stored properties\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.update_linearized_system_equation!-Tuple{Any, Any, Any, JutulEquation, CompactAutoDiffCache}","page":"Home","title":"Jutul.update_linearized_system_equation!","text":"Update a linearized system based on the values and derivatives in the equation.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.update_secondary_variable!-NTuple{5, Any}","page":"Home","title":"Jutul.update_secondary_variable!","text":"Update a secondary variable. Normally autogenerated with @jutul_secondary\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.update_state_dependents!-Tuple{Any, JutulModel, Any, Any}","page":"Home","title":"Jutul.update_state_dependents!","text":"Perform updates of everything that depends on the state.\n\nThis includes properties, governing equations and the linearized system\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.update_values!-Tuple{AbstractArray, AbstractArray}","page":"Home","title":"Jutul.update_values!","text":"update_values!(x, dx)\n\nReplace values of x in-place by y, leaving x with the avlues of y and the partials of x.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.value-Tuple{AbstractDict}","page":"Home","title":"Jutul.value","text":"value(d::Dict)\n\nCall value on all elements of some Dict.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.value-Tuple{Any}","page":"Home","title":"Jutul.value","text":"Take value of AD.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.values_per_entity-Tuple{Any, JutulVariables}","page":"Home","title":"Jutul.values_per_entity","text":"Number of values held by a primary variable. Normally this is equal to the number of degrees of freedom, but some special primary variables are most conveniently defined by having N values and N-1 independent variables.\n\n\n\n\n\n","category":"method"},{"location":"#Jutul.@jutul_secondary-Tuple{Any}","page":"Home","title":"Jutul.@jutul_secondary","text":"Designate the function as updating a secondary variable.\n\nThe function is then declared, in addition to helpers that allows checking what the dependencies are and unpacking the dependencies from state.\n\nIf we define the following function annotated with the macro: @jutulsecondary function updateas_secondary!(target, var::MyVarType, model, a, b, c)     @. target = a + b / c end\n\nThe macro also defines:  function get_dependencies(var::MyVarType, model)    return [:a, :b, :c] end\n\nfunction updatesecondaryvariable!(arraytarget, var::MyVarType, model, state)     updateassecondary!(arraytarget, var, model, state.a, state.b, state.c) end\n\nNote that the names input arguments beyond the parameter dict matter, as these will be fetched from state.\n\n\n\n\n\n","category":"macro"}]
}
