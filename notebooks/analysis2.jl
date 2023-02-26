### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ d4fe6472-5c91-11ed-15a1-6916a5210bde
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")

    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 1e60da4f-ddd8-4f57-8efb-302a195bf69e
TableOfContents()

# ╔═╡ 4e9129c9-4ac7-4d00-965b-d11d4f1ad21c
path = joinpath(dirname(pwd()), "results", "masses.csv")

# ╔═╡ d75a0d1a-5c4c-48bf-a5f8-5a6aa4d93491
df = DataFrame(CSV.File(path))

# ╔═╡ 3e404ad0-bb46-4de2-be54-b01381acc73f
df_set1, df_set2 = df[1:27, :], df[28:54, :]

# ╔═╡ 2db18926-4636-43f0-aa16-527cac6973bb
begin
	gt_mass_hd = copy(df_set1[!, :ground_truth_mass_hd]) .* 3
	gt_mass_md = copy(df_set1[!, :ground_truth_mass_md]) .* 3
	gt_mass_ld = copy(df_set1[!, :ground_truth_mass_ld]) .* 3
end

# ╔═╡ 14d68ff4-5627-41cf-aa58-008754550af9
md"""
# Plot
"""

# ╔═╡ 34c904ee-deb3-4da8-a00e-be6882c080ad
medphys_theme = Theme(;
    Axis=(
        backgroundcolor=:white,
        xgridcolor=:gray,
        xgridwidth=0.5,
        xlabelfont=:Helvetica,
        xticklabelfont=:Helvetica,
        xlabelsize=20,
        xticklabelsize=20,
        # xminorticksvisible = true,
        ygridcolor=:gray,
        ygridwidth=0.5,
        ylabelfont=:Helvetica,
        yticklabelfont=:Helvetica,
        ylabelsize=20,
        yticklabelsize=20,
        # yminortickvisible = true,
        bottomsplinecolor=:black,
        leftspinecolor=:black,
        titlefont=:Helvetica,
        titlesize=30,
    ),
)

# ╔═╡ 2fad3bf4-323c-463e-853a-27ffe5239c1f
function lin_reg()
	f = Figure()

	## A ##
	ax1 = Axis(f[1,1])
	scatter!(ax1, df_set1[!, :ground_truth_mass_hd], df_set1[!, :predicted_mass_hd])
	scatter!(ax1, df_set1[!, :ground_truth_mass_md], df_set1[!, :predicted_mass_md])
	scatter!(ax1, df_set1[!, :ground_truth_mass_ld], df_set1[!, :predicted_mass_ld])
	lines!(ax1, [-1000, 1000], [-1000, 1000])
	xlims!(ax1; low=-10, high=200)
    ylims!(ax1; low=-10, high=200)
    ax1.xticks = [0, 50, 100, 150, 200]
    ax1.yticks = [0, 50, 100, 150, 200]

	ax1.xlabel = "Known Mass (mg)"
    ax1.ylabel = "Calculated Mass (mg)"
    ax1.title = "Set 1 (52 - 733 mg/cc)"

	## B ##
	ax2 = Axis(f[2, 1])
	sc1 = scatter!(ax2, df_set2[!, :ground_truth_mass_hd], df_set2[!, :predicted_mass_hd])
	sc2 = scatter!(ax2, df_set2[!, :ground_truth_mass_md], df_set2[!, :predicted_mass_md])
	sc3 = scatter!(ax2, df_set2[!, :ground_truth_mass_ld], df_set2[!, :predicted_mass_ld])
	ln1 = lines!(ax2, [-1000, 1000], [-1000, 1000])
	xlims!(ax2; low=-10, high=75)
    ylims!(ax2; low=-10, high=75)
    ax2.xticks = [0, 25, 50, 75]
    ax2.yticks = [0, 25, 50, 75]

	ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Set 2 (27 - 797 mg/cc)"

	##-- LABELS --##
    f[1:2, 2] = Legend(
        f,
        [sc1, sc2, sc3, ln1],
        ["High Density", "Medium Density", "Low Density", "Unity"];
        framevisible=false,
    )

    for (label, layout) in zip(["A", "B"], [f[1, 1], f[2, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            textsize=25,
            padding=(0, 60, 25, 0),
            halign=:right,
        )
    end
	
	f
end

# ╔═╡ e5fe65f3-6b6e-4d20-86b0-64d8032ef6a9
function lin_reg2()
	f = Figure()

	## A ##
	ax1 = Axis(f[1,1])
	scatter!(ax1, gt_mass_hd, df_set1[!, :predicted_mass_hd])
	scatter!(ax1, gt_mass_md, df_set1[!, :predicted_mass_md])
	scatter!(ax1, gt_mass_ld, df_set1[!, :predicted_mass_ld])
	lines!(ax1, [-1000, 1000], [-1000, 1000])
	xlims!(ax1; low=-10, high=200)
    ylims!(ax1; low=-10, high=200)
    ax1.xticks = [0, 50, 100, 150, 200]
    ax1.yticks = [0, 50, 100, 150, 200]

	ax1.xlabel = "Known Mass (mg)"
    ax1.ylabel = "Calculated Mass (mg)"
    ax1.title = "Set 1 Corrected (52 - 733 mg/cc)"

	## B ##
	ax2 = Axis(f[2, 1])
	sc1 = scatter!(ax2, df_set2[!, :ground_truth_mass_hd], df_set2[!, :predicted_mass_hd])
	sc2 = scatter!(ax2, df_set2[!, :ground_truth_mass_md], df_set2[!, :predicted_mass_md])
	sc3 = scatter!(ax2, df_set2[!, :ground_truth_mass_ld], df_set2[!, :predicted_mass_ld])
	ln1 = lines!(ax2, [-1000, 1000], [-1000, 1000])
	xlims!(ax2; low=-10, high=75)
    ylims!(ax2; low=-10, high=75)
    ax2.xticks = [0, 25, 50, 75]
    ax2.yticks = [0, 25, 50, 75]

	ax2.xlabel = "Known Mass (mg)"
    ax2.ylabel = "Calculated Mass (mg)"
    ax2.title = "Set 2 (27 - 797 mg/cc)"

	##-- LABELS --##
    f[1:2, 2] = Legend(
        f,
        [sc1, sc2, sc3, ln1],
        ["High Density", "Medium Density", "Low Density", "Unity"];
        framevisible=false,
    )

    for (label, layout) in zip(["A", "B"], [f[1, 1], f[2, 1]])
        Label(
            layout[1, 1, TopLeft()],
            label;
            textsize=25,
            padding=(0, 60, 25, 0),
            halign=:right,
        )
    end
	
	f
end

# ╔═╡ c21d19f1-cd3d-4962-a662-edecd0e5d5c6
with_theme(medphys_theme) do
    lin_reg()
end

# ╔═╡ 3a264f06-00bf-4920-a25b-69e558ae761e
with_theme(medphys_theme) do
    lin_reg2()
end

# ╔═╡ Cell order:
# ╠═d4fe6472-5c91-11ed-15a1-6916a5210bde
# ╠═1e60da4f-ddd8-4f57-8efb-302a195bf69e
# ╠═4e9129c9-4ac7-4d00-965b-d11d4f1ad21c
# ╠═d75a0d1a-5c4c-48bf-a5f8-5a6aa4d93491
# ╠═3e404ad0-bb46-4de2-be54-b01381acc73f
# ╠═2db18926-4636-43f0-aa16-527cac6973bb
# ╟─14d68ff4-5627-41cf-aa58-008754550af9
# ╟─34c904ee-deb3-4da8-a00e-be6882c080ad
# ╟─2fad3bf4-323c-463e-853a-27ffe5239c1f
# ╟─e5fe65f3-6b6e-4d20-86b0-64d8032ef6a9
# ╟─c21d19f1-cd3d-4962-a662-edecd0e5d5c6
# ╟─3a264f06-00bf-4920-a25b-69e558ae761e
