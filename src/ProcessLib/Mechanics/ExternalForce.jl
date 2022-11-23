# This function assemble the external force to the global force vector
function fu_ext!(fu_ext::Array{T2}, traction::T2) where {T2 <: Float64}
	if traction == 0.0
		return fu_ext
	end
	Up_conf = -traction
	##
	elenum_top::Int64 = size(ymax, 1) .- 1
	topNodes = sort(map(tuple, ymax, node[ymax, 1]), by = x -> x[2])
	topNodes = [topNodes[i][1] for i ∈ 1:size(ymax, 1)]
	##
	for i ∈ 1:elenum_top
		tsctr = 2 * [topNodes[i], topNodes[i+1]]
		f_ext[tsctr] = f_ext[tsctr] + [0.5, 0.5] * Up_conf * abs(node[topNodes[i+1], 1] - node[topNodes[i], 1])
	end
	return f_ext
end