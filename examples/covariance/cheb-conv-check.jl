using LinearAlgebra, Plots

include("utils.jl")

M = load_sparse_matrix(
    reinterpret(Float64, read("M_data.bin")), 
    reinterpret(Int64, read("M_colind.bin")), 
    reinterpret(Int64, read("M_rowptr.bin"))
    )
L = load_sparse_matrix(
    reinterpret(Float64, read("L_data.bin")), 
    reinterpret(Int64, read("L_colind.bin")), 
    reinterpret(Int64, read("L_rowptr.bin"))
    )

##

R = cholesky(Matrix(M + M')/2).U
A = R' \ (Matrix(L) / R)
F = eigen(A)
Phi = R \ F.vectors
Lam = F.values

##

M_l = Diagonal(sum(eachrow(M)))
A_l = sqrt.(M_l) \ (Matrix(L) / sqrt.(M_l))
F_l = eigen(A_l)
Phi_l = sqrt.(M_l) \ F_l.vectors
Lam_l = F_l.values

##

mesh_file = "../../../../butterfly-LBO-models/bunny0.obj"
meshname  = split(split(mesh_file, "/")[end], ".")[1]
_, mesh, vertices, n = load_mesh(mesh_file)

k  = 1e-4
nu = 0.0
g(l) = matern_sdf(k, nu, l)

ps   = Int64.(readdlm("ps.txt")[1:end-1])[1:4]
tols = Float64.(readdlm("tols.txt")[1:end-1])

cheb_matvecs = Vector{Matrix{Float64}}(undef, length(ps))
lbo_matvecs  = Vector{Matrix{Float64}}(undef, length(tols))

filestring = @sprintf("cheb_p%i_kappa%.1e_nu%.1e", ps[1], k, nu)
randvecs = reshape(reinterpret(Float64, read("randvecs_$filestring.bin")), :, n)

## read covariance matvecs

for (i, p) in enumerate(ps)
    filestring = @sprintf("cheb_p%i_kappa%.1e_nu%.1e", p, k, nu)
    cheb_matvecs[i] = reshape(reinterpret(Float64, read("matvecs_$filestring.bin")), :, n)
end

for (i, tol) in enumerate(tols)
    filestring = @sprintf("lbo_tol%.0e_kappa%.1e_nu%.1e", tol, k, nu)
    lbo_matvecs[i] = reshape(reinterpret(Float64, read("matvecs_$filestring.bin")), :, n)
end

##

C = Phi * Diagonal(g.(Lam).^2) * Phi'
eigs = eigvals(C)

C_l = Phi_l * Diagonal(g.(Lam_l).^2) * Phi_l'
eigs_l = eigvals(C_l)

C_cheb    = [cheb_matvecs[i]/randvecs for i in eachindex(cheb_matvecs)]
eigs_cheb = eigvals.(C_cheb)

C_lbo    = [lbo_matvecs[i]/randvecs  for i in eachindex(lbo_matvecs)]
eigs_lbo = eigvals.(C_lbo)

##

gr(size=(500,500))
plot(
    [Lam Lam_l], 
    labels=["true" "lumped"], 
    yscale=:log10, ylims=[1, maximum(Lam)],
    xlabel="ℓ", ylabel="λℓ",
    linewidth=2
    )

pl = plot(
    yscale=:log10,
    xlabel="ℓ", ylabel="eigenvalue ℓ of C",
    legend=:topright
)
for (i, p) in enumerate(ps)
    v = reverse(real.(eigs_cheb[i]))
    v[v .< 1e-16] .= NaN
    plot!(pl,
        v, label="p=$p", linestyle=:dash, color=palette(:default)[i]
    )
end
for (i, tol) in enumerate(tols[[1,3,4]])
    v = reverse(real.(eigs_lbo[[1,3,4]][i]))
    v[v .< 1e-16] .= NaN
    plot!(pl,
        v, label="tol=$tol", color=palette(:default)[i]
    )
end

v = reverse(real.(eigs_l))
v[v .< 1e-16] .= NaN
plot!(pl,
    v, label="lumped", linestyle=:dash, color=:black #, linewidth=2
)

v = reverse(real.(eigs))
v[v .< 1e-16] .= NaN
plot!(pl,
    v, label="true", color=:black
)

display(pl)
