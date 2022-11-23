# This file defines the shape functions and their derivatives for quadrilateral elements.
function shapeFunctions(::Quadrilateral)
	##1.feisoq
	function feisoq(s::Float64, t::Float64)
		#derivation
		dNds::Array{Float64} = [-(1 - t) / 4, (1 - t) / 4, (1 + t) / 4, -(1 + t) / 4]
		dNdt::Array{Float64} = [-(1 - s) / 4, -(1 + s) / 4, (1 + s) / 4, (1 - s) / 4]
		return dNds, dNdt
	end
	##2.fejacob
	function fejacob(nnel::Int64, dNds::Array{Float64}, dNdt::Array{Float64}, Xcoor::Array{Float64}, Ycoor::Array{Float64})
		jacob = zeros(Float64, 2, 2)
		for g ∈ 1:nnel
			jacob[1, 1] = jacob[1, 1] + dNds[g] * Xcoor[g]
			jacob[1, 2] = jacob[1, 2] + dNds[g] * Ycoor[g]

			jacob[2, 1] = jacob[2, 1] + dNdt[g] * Xcoor[g]
			jacob[2, 2] = jacob[2, 2] + dNdt[g] * Ycoor[g]
		end
		return jacob
	end
	##3.federiv
	function federiv(nnel::Int64, dNds::Array{Float64}, dNdt::Array{Float64}, invjacob::Array{Float64})
		dNdx = zeros(Float64, nnel)
		dNdy = zeros(Float64, nnel)
		for s ∈ 1:nnel
			dNdx[s] = invjacob[1, 1] * dNds[s] + invjacob[1, 2] * dNdt[s]
			dNdy[s] = invjacob[2, 1] * dNds[s] + invjacob[2, 2] * dNdt[s]
		end
		return dNdx, dNdy
	end
	##Assemble
	nnel::Int64 = size(element, 2)
	Ns = zeros(Float64, 4, 4, nel)
	Nu = zeros(Float64, 8, 8, nel)
	Bu = zeros(Float64, 12, 8, nel)
	detjacob = zeros(Float64, 4, nel)
	for iel ∈ 1:nel
		#节点坐标
		nd = zeros(Int64, nnel)
		Xcoor = zeros(Float64, nnel)
		Ycoor = zeros(Float64, nnel)
		for i ∈ 1:nnel
			nd[i] = element[iel, i]
			Xcoor[i] = node[nd[i], 1]
			Ycoor[i] = node[nd[i], 2]
		end
		gauss = [-0.5774 -0.5774; -0.5774 0.5774; 0.5774 -0.5774; 0.5774 0.5774]
		for i ∈ 1:4
			s::Float64 = gauss[i, 1]
			t::Float64 = gauss[i, 2]
			N1::Float64 = 0.25 * (1 - s) * (1 - t)
			N2::Float64 = 0.25 * (1 + s) * (1 - t)
			N3::Float64 = 0.25 * (1 + s) * (1 + t)
			N4::Float64 = 0.25 * (1 - s) * (1 + t)
			Ns[i, :, iel] = [N1 N2 N3 N4]
			dNds, dNdt = feisoq(s, t)
			jacob = fejacob(nnel, dNds, dNdt, Xcoor, Ycoor) #雅可比矩阵
			detjacob[i, iel] = det(jacob)
			invjacob = inv(jacob)
			dNdx, dNdy = federiv(nnel, dNds, dNdt, invjacob)
			Nu[(2*i-1):(2*i), :, iel] = [[N1 0; 0 N1] [N2 0; 0 N2] [N3 0; 0 N3] [N4 0; 0 N4]]
			Bu[(3*i-2):(3*i), :, iel] = [[dNdx[1] 0.0; 0.0 dNdy[1]; dNdy[1] dNdx[1]] [dNdx[2] 0; 0 dNdy[2]; dNdy[2] dNdx[2]] [
				dNdx[3] 0.0; 0.0 dNdy[3]; dNdy[3] dNdx[3]] [dNdx[4] 0; 0 dNdy[4]; dNdy[4] dNdx[4]]]
		end
	end
	return Ns, Nu, Bu, detjacob
end
