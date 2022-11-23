# This function defines the degree of freedom (DOF) table for quadrilateral elements in a finite element mesh. The DOF table is used to map the local element DOFs to the global DOFs in the system. In this case, each node has two DOFs (displacement in x and y directions), and the function constructs the edofMat matrix that contains the global DOF indices for each element. The iK and jK arrays are also generated for assembling the global stiffness matrix efficiently.
function dofTables(::Quadrilateral)
	edofMat = zeros(Int64, nel, 2 * size(element, 2))
	for iel ∈ 1:nel  #loop within element
		edofMat[iel, 1] = 2 * element[iel, 1] - 1
		edofMat[iel, 2] = 2 * element[iel, 1]
		edofMat[iel, 3] = 2 * element[iel, 2] - 1
		edofMat[iel, 4] = 2 * element[iel, 2]
		edofMat[iel, 5] = 2 * element[iel, 3] - 1
		edofMat[iel, 6] = 2 * element[iel, 3]
		edofMat[iel, 7] = 2 * element[iel, 4] - 1
		edofMat[iel, 8] = 2 * element[iel, 4]
	end
	iK::Array{Int64} = reshape(kron(edofMat, ones(Int64, 8, 1))', 64 * nel)
	jK::Array{Int64} = reshape(kron(edofMat, ones(Int64, 1, 8))', 64 * nel)
	return iK, jK, edofMat
end

