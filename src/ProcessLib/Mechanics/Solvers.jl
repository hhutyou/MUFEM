function solvers(::Quadrilateral, mat::LinearElastic)
	#----------- Initiation -----------#
	dim = DisplacementDim(Quadrilateral())
	if dim != 2
		error("Currently only 2D problems are supported for Mechanics process. Please check the mesh and cell type.")
	end
	u = zeros(Float64, 2 * nnode)
	fu_ext = zeros(Float64, 2 * nnode)
	fu_int = zeros(Float64, 2 * nnode)
	sigma = zeros(Float64, 3, 4, nel)
	strain = zeros(Float64, 3, 4, nel)
	#----------- Start loading -----------#
	fu_ext .= fu_ext!(fu_ext, traction)
	Δu = 0.0
	for stp ∈ 1:typemax(Int64)
		Δu = Δu + u_increment
		abs(Δu) > abs(u_increment_total) ? (print("Computation finished."); break) : (@info "Step = $stp, Δu = $Δu.")
		u[loaddofs] .+= u_increment
		# temporary container
		u_old = zeros(Float64, 2 * nnode)
		#
		nit = 0
		while (nit < maxit)
			nit += 1
			@info "-iteration: $nit"
			# solve displacement (N-R)
			residual, Jac, sigma, strain = localAssemble(Quadrilateral(), mat, u, fu_ext)
			u[freedofs] .= u[freedofs] .+ Jac[freedofs, freedofs] \ residual[freedofs]
			norm(u) > eps() ? err_u = norm(u .- u_old) / norm(u) : err_u = 0.0
			@info "--err_u = $err_u."
			err_u < tol_u ? break : nothing
			u_old .= u
		end
		# Post-processing
		@info "Outputing: displacement, stress and strain."
		begin
			if mod(stp, output_num) == 0
				#output u ux uy uz
				open(output * "u_step_$stp.dat", "w") do io
					write(
						io,
						"TITLE=\"Displacement\" VARIABLES=\"X\",\"Y\",\"u\",\"ux\",\"uy\" ZONE N=$nnode,E=$nel,F=FEPOINT,ET=QUADRILATERAL, \n",
					)
					writedlm(
						io,
						[node[1:nnode, :] sqrt.(
							u[1:2:end] .^ 2 .+ u[2:2:end] .^ 2,
						) u[1:2:end] u[2:2:end]],
					)
					writedlm(io, element)
				end
				#output σx σy σxy
				open(output * "sigma_step_$stp.dat", "w") do io
					write(
						io,
						"TITLE=\"Stress\" VARIABLES=\"X\",\"Y\",\"stress_x\",\"stress_y\",\"stress_xy\" ZONE N=$nnode,E=$nel,F=FEPOINT,ET=QUADRILATERAL, \n",
					)
					writedlm(io, [node interpolationFromIntegrationPointToNode(sigma)'])
					writedlm(io, element)
				end
				#output epsx epsy epsxy
				open(output * "strain_step_$stp.dat", "w") do io
					write(
						io,
						"TITLE=\"Strain\" VARIABLES=\"X\",\"Y\",\"strain_x\",\"strain_y\",\"strain_xy\" ZONE N=$nnode,E=$nel,F=FEPOINT,ET=QUADRILATERAL, \n",
					)
					writedlm(io, [node interpolationFromIntegrationPointToNode(strain)'])
					writedlm(io, element)
				end
			end
			#
			push!(Fload, sum(fu_int[loaddofs]))
			push!(Uload, sum(u[loaddofs]) / size(u[loaddofs], 1))
		end
	end
	return Fload, Uload
end
