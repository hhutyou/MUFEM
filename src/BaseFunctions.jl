## dofs given the node number
function xdirect(x::T)::T where {T <: Array{Int}}
	x = 2 .* x .- 1
end
function ydirect(x::T) where {T <: Array{Int}}
	x = 2 .* x
end