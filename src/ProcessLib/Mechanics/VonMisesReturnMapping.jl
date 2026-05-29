# Return-mapping algorithm for J2 plasticity with isotropic linear hardening
# (von Mises yield criterion, plane-strain, small-strain).
#
# Voigt convention used throughout:
#   stress/strain vectors: [11, 22, 33, 12]
#   strain[4] = γ12 = 2·ε12  (engineering shear strain)
#   De_4 (4×4): [λ+2μ λ λ 0; λ λ+2μ λ 0; λ λ λ+2μ 0; 0 0 0 μ]
#
# State variables per integration point:
#   eps_p  – plastic strain 4-vector [ε11_p, ε22_p, ε33_p, γ12_p]
#   alpha  – equivalent (accumulated) plastic strain scalar
#
# Inputs
#   eps_4     : total strain 4-vector at current iteration  [ε11, ε22, 0, γ12]
#   eps_p_old : plastic strain at start of load step
#   alpha_old : equiv. plastic strain at start of load step
#   De_4      : 4×4 plane-strain elastic stiffness
#   μ         : shear modulus
#   K_bulk    : bulk modulus  K = λ + 2μ/3
#   σ_y       : initial uniaxial yield stress
#   H         : isotropic linear hardening modulus
#
# Returns
#   sigma_4_new  : updated 4-component stress
#   eps_p_new    : updated 4-component plastic strain
#   alpha_new    : updated equivalent plastic strain
#   D_ep_3x3     : 3×3 consistent algorithmic tangent  (rows/cols {1,2,4} of D_ep_4x4,
#                  i.e. relating [dσ11,dσ22,dσ12] to [dε11,dε22,dγ12] at εzz=0)

function returnMappingVonMises(
    eps_4     :: Vector{Float64},
    eps_p_old :: Vector{Float64},
    alpha_old :: Float64,
    De_4      :: Matrix{Float64},
    μ         :: Float64,
    K_bulk    :: Float64,
    σ_y       :: Float64,
    H         :: Float64,
)
    idx3 = [1, 2, 4]   # indices for condensed 3-component Voigt
    e4   = [1.0, 1.0, 1.0, 0.0]

    # --- trial elastic stress ---
    sigma_trial = De_4 * (eps_4 - eps_p_old)

    # --- deviatoric trial stress ---
    p = (sigma_trial[1] + sigma_trial[2] + sigma_trial[3]) / 3.0
    s = [sigma_trial[1]-p, sigma_trial[2]-p, sigma_trial[3]-p, sigma_trial[4]]

    # --- von Mises equivalent stress  q = √(3/2 · sᵢⱼsᵢⱼ) in tensor notation ---
    # s:s (tensor double contraction) = s11²+s22²+s33²+2·s12²
    q = sqrt(max(1.5 * (s[1]^2 + s[2]^2 + s[3]^2 + 2.0*s[4]^2), 0.0))

    # --- current yield stress ---
    q_y = σ_y + H * alpha_old

    # --- condensed elastic tangent (3×3) ---
    De_3x3 = De_4[idx3, idx3]

    # --- elastic step: no correction needed ---
    if q <= q_y || q < 1e-14
        return sigma_trial, copy(eps_p_old), alpha_old, De_3x3
    end

    # --- plastic step: closed-form radial return ---
    Δγ = (q - q_y) / (3.0μ + H)
    r  = 3.0μ * Δγ / q          # = 3μΔγ/q  (deviatoric scaling factor)

    # stress update:  σ_new = σ_trial − (3μΔγ/q)·s
    sigma_new = sigma_trial - r .* s

    # plastic strain update (engineering shear for index 4):
    #   Δεᵢᵢ_p = Δγ·(3sᵢᵢ)/(2q)    (normal)
    #   Δγ12_p = Δγ·(3s12)/q        (engineering shear = 2·tensor shear)
    Δeps_p = Δγ * (1.5 / q) * [s[1], s[2], s[3], 2.0*s[4]]
    eps_p_new = eps_p_old + Δeps_p

    alpha_new = alpha_old + Δγ

    # --- consistent algorithmic tangent (4×4) ---
    # Derived by implicit differentiation of the return-mapping equations:
    #   D_ep = (1−r)·De + r·K·(e⊗eᵀ) − β·(s⊗sᵀ)
    # where β = 9μ²·q_y / (q³·(3μ+H))
    β_ep    = 9.0μ^2 * q_y / (q^3 * (3.0μ + H))
    D_ep_4x4 = (1.0-r) .* De_4 .+ r*K_bulk .* (e4 * e4') .- β_ep .* (s * s')

    # condense to 3×3 by fixing εzz=0 (plane-strain constraint)
    D_ep_3x3 = D_ep_4x4[idx3, idx3]

    return sigma_new, eps_p_new, alpha_new, D_ep_3x3
end
