using FileIO, LinearAlgebra, SparseArrays, Printf, DelimitedFiles, Meshes, GeometryBasics, MeshIO, Plots, LaTeXStrings
import GLMakie as Mke

include("utils.jl")

## load mesh
println("Loading mesh...")

mesh_file = ARGS[1]
meshname  = split(split(mesh_file, "/")[end], ".")[1]
_, mesh, vertices, n = load_mesh(mesh_file)

## read other parameters

k  = parse(Float64, ARGS[2])
nu = parse(Float64, ARGS[3])

ps   = Int64.(readdlm("ps.txt")[1:end-1])
tols = Float64.(readdlm("tols.txt")[1:end-1])

num_ps   = length(ps)
num_tols = length(tols)-1 # last is reference

cheb_matvecs = Vector{Matrix{Float64}}(undef, num_ps)
lbo_matvecs  = Vector{Matrix{Float64}}(undef, num_tols)

## read covariance matvecs

for (i, p) in enumerate(ps)
    filestring = @sprintf("cheb_p%i_kappa%.1e_nu%.1e", p, k, nu)
    cheb_matvecs[i] = reshape(reinterpret(Float64, read("matvecs_$filestring.bin")), :, n)
end

for (i, tol) in enumerate(tols)
    filestring = @sprintf("lbo_tol%.0e_kappa%.1e_nu%.1e", tol, k, nu)
    if i==length(tols)
        global ref_matvecs = reshape(reinterpret(Float64, read("matvecs_$filestring.bin")), :, n)
    else
        lbo_matvecs[i] = reshape(reinterpret(Float64, read("matvecs_$filestring.bin")), :, n)
    end
end

## read timings

perf = readdlm(@sprintf("performance_kappa%.1e_nu%.1e.txt", k, nu))

lbo_timings  = Float64.(perf[2:2+num_tols-1, 4])
cheb_timings = Float64.(perf[num_tols+4:end, 2])

## estimate errors using Hutchinson

function hutchinson_frobenius(matvecs; corr_matvecs=nothing)
    nmv = size(matvecs, 1) 
    zs  = norm.(eachrow(matvecs)).^2
    mean_z = sum(zs) / nmv
    var_z  = sum((zs .- mean_z).^2) / (nmv-1)

    if isnothing(corr_matvecs)
        return mean_z, var_z
    else
        nmv = size(corr_matvecs, 1) 
        ws  = norm.(eachrow(corr_matvecs)).^2
        mean_w = sum(ws) / nmv
        var_w  = sum((ws .- mean_w).^2) / (nmv-1)

        rho = sum((zs .- mean_z).*(ws .- mean_w)) / sqrt(var_z*var_w) / (nmv-1)

        return mean_z, var_z, rho
    end
end

function ratio_variance(m_y, v_y, m_x, v_x, rho, nmv) 
    return (v_y/m_x^2 - 2m_y*rho/m_x^3 + m_y^2*v_x/m_x^4) / nmv
end

nmv = size(lbo_matvecs[1], 1)

cheb_frobsq = Vector{Vector{Float64}}(undef, num_ps)
lbo_frobsq  = Vector{Vector{Float64}}(undef, num_tols)

m_x, v_x = hutchinson_frobenius(ref_matvecs)

for i in eachindex(cheb_matvecs)
    m_y, v_y, rho = hutchinson_frobenius(
        cheb_matvecs[i] - ref_matvecs[1:nmv, :],
        corr_matvecs=ref_matvecs
        )
    r   = m_y / m_x
    v_r = ratio_variance(m_y, v_y, m_x, v_x, rho, nmv) 
    cheb_frobsq[i] = [r, v_r, m_y, v_y, rho]
end

for i in eachindex(lbo_matvecs)
    m_y, v_y, rho = hutchinson_frobenius(
        lbo_matvecs[i] - ref_matvecs[1:nmv, :],
        corr_matvecs=ref_matvecs
        )
    r   = m_y / m_x
    v_r = ratio_variance(m_y, v_y, m_x, v_x, rho, nmv) 
    lbo_frobsq[i] = [r, v_r, m_y, v_y, rho]
end

# norm squared
cheb_normsq = getindex.(cheb_frobsq, 1)
lbo_normsq  = getindex.(lbo_frobsq, 1)

# norms
cheb_norm = sqrt.(cheb_normsq)
@show lbo_norm  = sqrt.(lbo_normsq)

# standard deviation of norm squared
cheb_stdsq = sqrt.(getindex.(cheb_frobsq, 2))
lbo_stdsq  = sqrt.(getindex.(lbo_frobsq, 2))

# upper and lower differences from norm
cheb_uppers = sqrt.(cheb_normsq + 2*cheb_stdsq) - cheb_norm
cheb_lowers = cheb_norm - sqrt.(cheb_normsq - 2*cheb_stdsq)
lbo_uppers  = sqrt.(lbo_normsq + 2*lbo_stdsq) - lbo_norm
lbo_lowers  = lbo_norm - sqrt.(lbo_normsq - 2*lbo_stdsq)

##

gr(size=(400,300))
pl = plot(
    # title=@sprintf(
    #     "%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
    #     split(mesh_file, "/")[end], n, k, nu
    #     ),
    xscale=:log10,
    yscale=:log10,
    xlabel="time per sample (s)",
    ylabel=L"$|| \Sigma - \tilde{\Sigma} \ || / || \Sigma ||$",
    xlims=[1e-3, 3],
    ylims=[6e-4, 1.1]
    )
plot!(pl,
    cheb_timings, cheb_norm, 
    yerror=(cheb_lowers, cheb_uppers),
    label="Chebyshev", line=(:red, 2, :dash), 
    marker=4, markerstrokewidth=0, markercolor=:red
    )
plot!(pl,
    lbo_timings, lbo_norm, 
    yerror=(lbo_lowers, lbo_uppers),
    label="Butterfly", line=(:blue, 1, :dash, 0.5), 
    marker=4, msw=1, markercolor=:blue
    )
savefig(pl, @sprintf("output/%s_est_error_kappa%.1e_nu%.1e.png", meshname, k, nu))
savefig(pl, @sprintf("output/%s_est_error_kappa%.1e_nu%.1e.pdf", meshname, k, nu))