function interpolationFromIntegrationPointToNode(x::Array{Float64, 2})
	## matrix used in caculating node stress
	ixf = Int.(element'[:])
	jxf = Int.(ones(4 * nel))
	###
	for iel ∈ 1:nel
		Tmatrix = inv(Ns[:, :, iel])
		x[:, iel] = x[:, iel]' * Tmatrix'  ## extrapolate value to nodes
	end
	return sparse(ixf, jxf, x[:]) ./ freq
end
##
function interpolationFromIntegrationPointToNode(x::Array{Float64, 3})
	## matrix used in caculating node stress
	ig = Int.(kron(ones(4 * nel), [1, 2, 3]))
	jg = Int.(kron(element'[:], ones(3)))
	###
	for iel ∈ 1:nel
		Tmatrix = inv(Ns[:, :, iel])
		x[:, :, iel] = x[:, :, iel] * Tmatrix'  ## extrapolate value to nodes
	end
	return sparse(ig, jg, x[:]) ./ kron(freq', ones(3))
end
