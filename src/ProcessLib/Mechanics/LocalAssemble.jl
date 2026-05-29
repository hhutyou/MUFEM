#----------- Assemble: linear elastic -----------#
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

#----------- Assemble: von Mises isotropic linear hardening -----------#
# eps_p_step : plastic strain at start of load step  (4, 4, nel)
# alpha_step : equiv. plastic strain at start of step (4, nel)
# Returns the updated eps_p_new and alpha_new alongside the usual outputs.
function localAssemble(
	::Quadrilateral,
	mat     :: VonMisesHardening,
	u       :: Array{T},
	fu_ext  :: Array{T},
	eps_p_step :: Array{Float64, 3},
	alpha_step :: Array{Float64, 2},
) where {T <: Float64}

	λ      = mat.E * mat.ν / ((1.0 + mat.ν) * (1.0 - 2.0mat.ν))
	μ      = mat.E / (2.0 * (1.0 + mat.ν))
	K_bulk = λ + 2.0μ / 3.0
	e4     = [1.0, 1.0, 1.0, 0.0]
	De_4   = [λ+2μ λ λ 0.0; λ λ+2μ λ 0.0; λ λ λ+2μ 0.0; 0.0 0.0 0.0 μ]

	sigma     = zeros(Float64, 3, 4, nel)   # [σ11, σ22, σ12] for output
	strain    = zeros(Float64, 3, 4, nel)   # [ε11, ε22, γ12] for output
	eps_p_new = zeros(Float64, 4, 4, nel)
	alpha_new = zeros(Float64, 4, nel)
	residual  = zeros(Float64, 2 * nnode)
	Jac_ele   = zeros(Float64, 64, nel)
	rhs_ele   = zeros(Float64, 8, nel)
	u_ele_total = u[edofMat]

	for iel ∈ 1:nel
		u_ele    = view(u_ele_total, iel, :)
		# 3-component strain at each IP: rows 3ip-2:3ip → reshape to (3, 4)
		eps_3_ele = reshape(view(Bu, :, :, iel) * u_ele, 3, 4)

		Ke_ele = zeros(Float64, 3, 3, 4)   # condensed 3×3 tangent per IP

		for ip ∈ 1:4
			eps_3 = eps_3_ele[:, ip]
			# augment to 4-component strain: insert εzz = 0 at index 3
			eps_4 = [eps_3[1], eps_3[2], 0.0, eps_3[3]]

			sigma_4, ep_ip, alp_ip, D3 = returnMappingVonMises(
				eps_4,
				eps_p_step[:, ip, iel],
				alpha_step[ip, iel],
				De_4, μ, K_bulk,
				mat.σ_y, mat.H,
			)

			sigma[:, ip, iel]    .= [sigma_4[1], sigma_4[2], sigma_4[4]]
			strain[:, ip, iel]   .= eps_3
			eps_p_new[:, ip, iel] = ep_ip
			alpha_new[ip, iel]    = alp_ip
			Ke_ele[:, :, ip]      = D3
		end

		rhs_ele[:, iel] .= kron(view(detjacob, :, iel)', ones(8, 3)) .* (view(Bu, :, :, iel)') * vec(view(sigma, :, :, iel))
		Jac_ele[:, iel]  = vec(
			kron(detjacob[:, iel]', ones(8, 3)) .* Bu[:, :, iel]' *
			blockdiag(sparse(Ke_ele[:, :, 1]), sparse(Ke_ele[:, :, 2]),
			          sparse(Ke_ele[:, :, 3]), sparse(Ke_ele[:, :, 4])) *
			Bu[:, :, iel]
		)
	end

	fu_int   = sparse(vec(edofMat'), ones(Int64, 8 * nel), vec(rhs_ele))
	Ku       = sparse(iK, jK, vec(Jac_ele))
	residual .= fu_ext .- fu_int

	return residual, (Ku .+ Ku') ./ 2, Array(sigma), Array(strain), eps_p_new, alpha_new
end
