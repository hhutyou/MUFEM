# This is a benchmark code for 2D elastic deformation, which can be used to verify the correctness of the implementation of the finite element method. The code is organized as follows:
# 1. Include necessary packages and user-defined functions.
# 2. Define parameters for the 2D elastic problem.
# 3. Consider heterogeneity if needed, and define material properties for each element.
# 4. Define parameters for loading, computation, and output.
# 5. Solve the problem using the defined solvers.
# 6. Post-process the results if needed, including outputting results to files and plotting.
# Note: The code is designed to be flexible and can be modified to include different constitutive models, physical fields, heterogeneity, and other features as needed. The user can also add additional functions for post-processing or other purposes.
# 
using Revise # for development
using LinearAlgebra, SparseArrays, DelimitedFiles, Meshes # for computations

#----------- Include necessary functions for 2D elastic problem -----------#
include("../../src/Types.jl")
include("../../src/BaseFunctions.jl")
include("../../src/interpolationFromIntegrationPointToNode.jl")

#----------- Generate 1m×1m mesh and define boundary conditions -----------#
const celltype = Quadrilateral
const nx_els, ny_els = 10, 10          # number of quad elements in x and y directions
const Lx, Ly = 1.0, 1.0               # domain dimensions [m]
grid = CartesianGrid((0.0, 0.0), (Lx, Ly); dims = (nx_els, ny_els))

# node: (nnode, 2) matrix of [x, y] coordinates
node_grid = Matrix{Float64}(undef, nvertices(grid), 2)
for (i, v) ∈ enumerate(vertices(grid))
	c = Meshes.coordinates(v)
	node_grid[i, 1] = c[1]
	node_grid[i, 2] = c[2]
end
const node = node_grid

# element: (nel, 4) matrix of node indices for each quadrilateral
topo = topology(grid)
elem_conn = Matrix{Int64}(undef, nelements(topo), 4)
for (i, conn) ∈ enumerate(elements(topo))
	elem_conn[i, :] .= indices(conn)
end
const element = elem_conn
const nnode = size(node, 1)
const nel = size(element, 1)
@info "Mesh generated: $nel quad elements, $nnode nodes"

# generate shape functions and doftables
include("../../src/ShapeFunctions.jl")
include("../../src/DofTables.jl")
@info "Formulating the shape functions and doftables takes"
@time begin
	const Ns, Nu, Bu, detjacob = shapeFunctions(Quadrilateral())
	const iK, jK, edofMat = dofTables(Quadrilateral())
	const freq = zeros(Int64, nnode)
	for i ∈ 1:nnode
		freq[i] = size(findall(element[:, :] .== i), 1)
	end
end

# find the boundary nodes
const tol_mesh = 1e-12 * max(Lx, Ly)   # floating-point tolerance for coordinate comparisons
const nodes_left = findall(x -> x < tol_mesh, node[:, 1])          # x ≈ 0
const nodes_right = findall(x -> x > Lx - tol_mesh, node[:, 1])          # x ≈ Lx
const nodes_bottom = findall(x -> x < tol_mesh, node[:, 2])          # y ≈ 0
const nodes_top = findall(x -> x > Ly - tol_mesh, node[:, 2])          # y ≈ Ly
@info "Boundary nodes — left: $(length(nodes_left)), right: $(length(nodes_right)), bottom: $(length(nodes_bottom)), top: $(length(nodes_top))"
const fix_x = nodes_bottom
const fix_y = nodes_bottom
const load_x = Int64[]
const load_y = nodes_top

# generate boundary doftables
const fixeddofs = union(ydirect(fix_y), xdirect(fix_x))
const loaddofs = union(xdirect(load_x), ydirect(load_y))
const freedofs = setdiff(1:2size(node, 1), fixeddofs, loaddofs)

#----------- Material model: von Mises with isotropic linear hardening -----------#
# E = 30e9 Pa, ν = 0.2, σ_y = 5e8 Pa (initial yield), H = 3e8 Pa (hardening modulus)
const mat = VonMisesHardening(30e9, 0.2, 5e8, 3e8)

#----------- Parameters for loading, computation and output -----------#
const tol_u, maxit = 1e-12, 100 # convergence tolerance for N-R solver of displacement
const u_increment, u_increment_total = 0.01, 0.1
const traction = 0.0
const output = joinpath(@__DIR__, "output") * "/"
isdir(output) || mkpath(output)
const output_num = 1 # Output results per output_num steps
################################
#----------- Sovle -----------#
include("../../src/ProcessLib/Mechanics/Mechanics.jl") # turn on the main mechanics process
Fload = []
Uload = []
Fload, Uload = solvers(Quadrilateral(), mat)
