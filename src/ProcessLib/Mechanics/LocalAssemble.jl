#----------- Assemble of the displacement equation -----------#
function localAssemble(::Quadrilateral, mat::LinearElastic, u::Array{T}, fu_ext::Array{T}) where {T <: Float64}
	dim = DisplacementDim(Quadrilateral())
	if dim != 2
		error("Currently only 2D problems are supported for Mechanics process. Please check the mesh and cell type.")
	end
	sigma = zeros(Float64, 3, 4, nel)
	strain = zeros(Float64, 3, 4, nel)
	residual = zeros(Float64, 2 * nnode)
	# temporary variables
	Jac_ele = zeros(Float64, 64, nel)
	rhs_ele = zeros(Float64, 8, nel)
	u_ele_total = u[edofMat]
	#----------- update stress and other quantities at integration points -----------#
	for iel ∈ 1:nel
		#------------ initiate variable on single element---------#
		u_ele = view(u_ele_total, iel, :)
		epsilon_ele = reshape(view(Bu,:,:,iel) * u_ele, 3, 4)
		# element-wise parameters
		De_ele = view(De,:,:,iel)
		Ke_ele = zeros(Float64, 3, 3, 4) # element stiffness matrix at each integration point
		for ip ∈ 1:4
			# update stress, Jacobian matrix on single integration point
			epsilon_ip = epsilon_ele[:, ip]
			################################
			sigma[:, ip, iel] .= De_ele * epsilon_ip
			strain[:, ip, iel] .= epsilon_ip
			Ke_ele[:, :, ip] .= De_ele
		end
		rhs_ele[:, iel] .= kron(view(detjacob, :, iel)', ones(8, 3)) .* (view(Bu,:,:,iel)') * (vec(view(sigma,:,:,iel)))
		Jac_ele[:, iel] = vec(kron(detjacob[:, iel]', ones(8, 3)) .* Bu[:, :, iel]' * blockdiag(sparse(Ke_ele[:, :, 1]), sparse(Ke_ele[:, :, 2]),
								  sparse(Ke_ele[:, :, 3]), sparse(Ke_ele[:, :, 4])) * Bu[:, :, iel])
	end
	# compute residual of the displacement equation
	fu_int = sparse(vec(edofMat'), ones(Int64, 8 * nel), vec(rhs_ele))
	Ku = sparse(iK, jK, vec(Jac_ele))
	residual .= fu_ext .- fu_int

	return residual, (Ku .+ Ku') ./ 2, Array(sigma), Array(strain)
end
