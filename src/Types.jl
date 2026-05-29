#
# Dispatch axis: coupled physical process
abstract type ProcessType end
struct Mechanics <: ProcessType end
struct Hydromechanics <: ProcessType end

# Dispatch axis: element topology
abstract type CellType end
struct Triangle <: CellType end
struct Quadrilateral <: CellType end
struct Tetrahedron <: CellType end
struct Hexahedron <: CellType end

# spatial dimension is a type-level property, not per-instance data
DisplacementDim(::Union{Triangle, Quadrilateral}) = 2
DisplacementDim(::Union{Tetrahedron, Hexahedron}) = 3

# Dispatch axis + data: constitutive model
abstract type ConstitutiveModel end

struct LinearElastic <: ConstitutiveModel
	E::Float64   # Young's modulus [Pa]
	ν::Float64   # Poisson's ratio [-]
end

struct ElastoPlastic <: ConstitutiveModel
	E::Float64   # Young's modulus [Pa]
	ν::Float64   # Poisson's ratio [-]
	φ::Float64   # friction angle [rad]
	c::Float64   # cohesion [Pa]
end

struct VonMisesHardening <: ConstitutiveModel
	E::Float64    # Young's modulus [Pa]
	ν::Float64    # Poisson's ratio [-]
	σ_y::Float64  # initial uniaxial yield stress [Pa]
	H::Float64    # isotropic linear hardening modulus [Pa]
end
