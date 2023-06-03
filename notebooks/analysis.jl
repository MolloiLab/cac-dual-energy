### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ a746ff38-1cde-4853-b001-7f2dfd6a693c
using DrWatson

# ╔═╡ 7cad0e61-2cde-44ef-903a-531f10b322f5
# ╠═╡ show_logs = false
@quickactivate "cac-dual-energy"

# ╔═╡ b169b3f6-819c-4f29-b32d-d5b28b9ae85f
using PlutoUI, Statistics, CSV, DataFrames, CairoMakie, Colors, GLM, MLJBase, Printf

# ╔═╡ 4ed99e12-126a-4728-a700-6942c678b403
using StatsBase: quantile!, rmsd

# ╔═╡ debe4f43-1698-461d-95c6-17f6052eee36
include(srcdir("plot_utils.jl")); include(srcdir("helper_functions.jl")); 

# ╔═╡ 7f4c4f06-9875-4fc5-bf0a-77fe243067ec
TableOfContents()

# ╔═╡ 7a4b45d8-2fbd-41d1-9df6-c31d15361327
root_path = datadir("output")

# ╔═╡ 58f66f99-006d-4a14-a2a4-3342e764c01d
medphys_theme = Theme(
    Axis = (
        backgroundcolor = :white,
		xgridcolor = :gray,
		xgridwidth = 0.1,
		xlabelsize = 15,
		xticklabelsize = 15,
		ygridcolor = :gray,
		ygridwidth = 0.1,
		ylabelsize = 15,
		yticklabelsize = 15,
		bottomsplinecolor = :black,
		leftspinecolor = :black,
		titlesize = 25
	)
);

# ╔═╡ b4bc009b-c9dc-40f4-9444-f12cae95289c
md"""
# Accuracy
"""

# ╔═╡ d9071ce3-0109-4979-9474-a9256acc1b58
df_m = CSV.read(datadir("results", "material_decomposition.csv"), DataFrame);

# ╔═╡ 8212d99a-a5c8-4e44-a2e4-8a6e598b9367
df_a = CSV.read(datadir("results", "agatston.csv"), DataFrame);

# ╔═╡ 20080b7f-bc37-453a-a426-c478ec805076
df_v = CSV.read(datadir("results", "volume_fraction.csv"), DataFrame);

# ╔═╡ 28cfcdf4-abd7-4699-bcf2-7d4af05340e1
co_1, r_squared_1, rms_values_1, pred_1 = calculate_coefficients(df_m);

# ╔═╡ 684af87a-a9f2-42fd-be26-40a12719bf22
co_2, r_squared_2, rms_values_2, pred_2 = calculate_coefficients(df_a);

# ╔═╡ 3d5c68d6-9e5c-430d-b17a-552bbd4d254d
co_3, r_squared_3, rms_values_3, pred_3 = calculate_coefficients(df_v);

# ╔═╡ ebd72800-16d2-42fc-9b7d-816efcbb675a
function accuracy()
	f = Figure()

	##-- A --##
	ax = Axis(
		f[1, 1],
		xticks = [0, 25, 50, 75, 100, 125],
		yticks = [0, 25, 50, 75, 100, 125],
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Material Decomposition",
	)
	
	df = df_m
	sc1=scatter!(df[!, :ground_truth_mass_hd], df[!, :predicted_mass_hd])
	sc2=scatter!(df[!, :ground_truth_mass_md], df[!, :predicted_mass_md])
	sc3=scatter!(df[!, :ground_truth_mass_ld], df[!, :predicted_mass_ld], color=:red)
	ln1=lines!([-1000, 1000], [-1000, 1000])
	ln2=lines!(collect(1:1000), pred_1, linestyle=:dashdot)
	create_textbox(f[1, 1], co_1, r_squared_1, rms_values_1)
	
	xlims!(ax, low=0, high=125)
	ylims!(ax, low=0, high=125)

	##-- B --##
	ax = Axis(
		f[2, 1],
		xticks = [0, 25, 50, 75, 100, 125],
		yticks = [0, 25, 50, 75, 100, 125],
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Agatson Scoring",
	)
	
	df = df_a
	scatter!(df[!, :ground_truth_mass_hd], df[!, :predicted_mass_hd])
	scatter!(df[!, :ground_truth_mass_md], df[!, :predicted_mass_md])
	scatter!(df[!, :ground_truth_mass_ld], df[!, :predicted_mass_ld], color=:red)
	lines!([-1000, 1000], [-1000, 1000],)
	lines!(collect(1:1000), pred_2, linestyle=:dashdot)
	create_textbox(f[2, 1], co_2, r_squared_2, rms_values_2)
	xlims!(ax, low=0, high=125)
	ylims!(ax, low=0, high=125)

	##-- C --##
	ax = Axis(
		f[3, 1],
		xticks = [0, 25, 50, 75, 100, 125],
		yticks = [0, 25, 50, 75, 100, 125],
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Volume Fraction",
	)
	df = df_v
	scatter!(df[!, :ground_truth_mass_hd], df[!, :predicted_mass_hd])
	scatter!(df[!, :ground_truth_mass_md], df[!, :predicted_mass_md])
	scatter!(df[!, :ground_truth_mass_ld], df[!, :predicted_mass_ld], color=:red)
	lines!([-1000, 1000], [-1000, 1000],)
	lines!(collect(1:1000), pred_3, linestyle=:dashdot)
	create_textbox(f[3, 1], co_3, r_squared_3, rms_values_3)
	xlims!(ax, low=0, high=125)
	ylims!(ax, low=0, high=125)

	#-- LABELS --##
	f[2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip(["A", "B","C"], [f[1,1], f[2,1], f[3,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end
	
	save(plotsdir("linear_reg.png"), f)
	f
end

# ╔═╡ 190d1636-41ab-4e98-ae6c-5d872a742c0f
with_theme(accuracy, medphys_theme)

# ╔═╡ 7ffa22cc-8214-4d39-9ff0-20b0b63f06cb
md"""
# Sensitivity and Specificity
"""

# ╔═╡ db529cfe-fd35-4a35-b9fd-8367bd7fb6fe
md"""
False-positive expectations
- 0.0 std => ~50% false-positive
- 0.5 std => ~30% false-positive
- 1.0 std => ~16% false-positive
- 1.5 std => ~7% false-positive
- 2.0 std => ~2% false-positive
"""

# ╔═╡ 5eb31f2b-bfaa-41d0-9375-9a076ba789d1
std_level = 1.5

# ╔═╡ 1f39b6b5-b839-4dff-bbf0-10b3aa7b84ef
md"""
## False Negative
"""

# ╔═╡ 2cff3bc9-1f37-4e56-a02f-b61c4cceb6f1
md"""
#### Agatston
"""

# ╔═╡ d2b062b3-30dc-454a-b2fa-56c41312cc38
array_a = hcat(df_a[!, :calculated_agat_large], df_a[!, :calculated_agat_medium], df_a[!, :calculated_agat_small]);

# ╔═╡ 9b52009a-4ef7-4f7f-aadb-89320a32f906
total_cac = length(array_a)

# ╔═╡ 47aa7f6c-067d-40c1-94b5-5bfdd30bb607
num_zero_a = length(findall(x -> x <= 0, array_a))

# ╔═╡ 17dad848-6491-45e4-8eec-3ed3e08c5318
length(findall(x -> x <= 0, df_a[!, :calculated_agat_large])), length(findall(x -> x <= 0, df_a[!, :calculated_agat_medium])), length(findall(x -> x <= 0, df_a[!, :calculated_agat_small]))

# ╔═╡ 12edad23-47d8-4c30-b595-0e48380868dd
40/360, 38/360, 25/360

# ╔═╡ fb9da5d1-a8cc-4485-9cf6-ab9c98c11816
begin
	df_a_ld, df_a_md, df_a_hd = groupby(df_a, :insert_sizes)
	
	length(findall(x -> x <= 0, hcat(df_a_ld[!, :calculated_agat_large], df_a_ld[!, :calculated_agat_medium], df_a_ld[!, :calculated_agat_small]))), length(findall(x -> x <= 0, hcat(df_a_md[!, :calculated_agat_large], df_a_md[!, :calculated_agat_medium], df_a_md[!, :calculated_agat_small]))), length(findall(x -> x <= 0, hcat(df_a_hd[!, :calculated_agat_large], df_a_hd[!, :calculated_agat_medium], df_a_hd[!, :calculated_agat_small])))
end

# ╔═╡ b8da1153-ae0a-4baa-9576-dffd76d2267b
md"""
#### Material Decomposition
"""

# ╔═╡ fdd43925-fc6a-4e26-8aea-8beb6c2c21fe
begin
	false_negative_m = []
	for i in 1:3:nrow(df_m)-2
		mean_m, std_m = mean(df_m[i:i+2, :bkg_mass]), std(df_m[i:i+2, :bkg_mass]) * std_level 
		array_m = hcat(df_m[i:i+2, :predicted_mass_hd], df_m[i:i+2, :predicted_mass_md], df_m[i:i+2, :predicted_mass_ld]);
		neg = length(findall(x -> x <= mean_m + std_m, array_m))
		push!(false_negative_m, neg)
	end
end

# ╔═╡ 0691ad9a-10c4-489b-a6bf-03dd46ef7f09
total_zero_m = sum(false_negative_m)

# ╔═╡ 0e9db9e2-5a41-48a6-9d11-7c38f172c131
md"""
#### Volume Fraction
"""

# ╔═╡ 4234bb6c-d94a-418f-a6cd-2bf1df8b7aff
begin
	false_negative_v = []
	for i in 1:3:nrow(df_v)-2
		mean_v, std_v = mean(df_v[i:i+2, :bkg_mass]), std(df_v[i:i+2, :bkg_mass]) * std_level 
		array_v = hcat(df_v[i:i+2, :predicted_mass_hd], df_v[i:i+2, :predicted_mass_md], df_v[i:i+2, :predicted_mass_ld]);
		neg = length(findall(x -> x <= mean_v + std_v, array_v))
		push!(false_negative_v, neg)
	end
end

# ╔═╡ 99ea3200-1d7d-4f6f-9ec0-8dd9427aadfc
total_zero_v = sum(false_negative_v)

# ╔═╡ 364ade03-b649-43be-9580-85f0e2fae7bd
total_zero_m, total_zero_v, num_zero_a

# ╔═╡ 881488ad-3eaf-46ed-8a86-7207df316ca4
md"""
## False Positive
"""

# ╔═╡ 5302a029-1bd0-4582-99bb-d60832cc696c
md"""
#### Agatston
"""

# ╔═╡ f127ec74-b7d3-4c83-bacd-40fcda4e14c7
array_a_pos = df_a[!, :bkg_mass]

# ╔═╡ 934bb2d3-9c17-4549-9dfb-e16ebb5963bc
total_cac_pos = length(array_a_pos)

# ╔═╡ bb1c2f6f-0d5a-49c8-9133-2aa65485a405
total_zero_a_pos = length(findall(x -> x > 0, array_a_pos))

# ╔═╡ 78a2bdd7-54c6-46a0-a9b4-a28f3c97448a
md"""
#### Material Decomposition
"""

# ╔═╡ 7be46fdc-0554-48b2-aa8e-55770c68f16b
begin
	false_positive_m = []
	for i in 1:3:nrow(df_m)-2
		mean_m, std_m = mean(df_m[i:i+2, :bkg_mass]), std(df_m[i:i+2, :bkg_mass]) * std_level
		array_m_pos = df_m[i:i+2, :bkg_mass]
		pos = length(findall(x -> x > mean_m + std_m, array_m_pos))
		push!(false_positive_m, pos)
	end
end

# ╔═╡ 0e428d01-39de-4274-a83d-1e2369ae35f5
total_zero_m_pos = sum(false_positive_m)

# ╔═╡ af083168-9989-407d-980b-9fbde3c8397c
md"""
#### Volume Fraction
"""

# ╔═╡ 89445e61-9398-477f-afc0-5c69ea6e3212
begin
	false_positive_v = []
	for i in 1:3:nrow(df_v)-2
		mean_v, std_v = mean(df_v[i:i+2, :bkg_mass]), std(df_v[i:i+2, :bkg_mass]) * std_level
		array_v_pos = df_v[i:i+2, :bkg_mass]
		pos = length(findall(x -> x > mean_v + std_v, array_v_pos))
		push!(false_positive_v, pos)
	end
end

# ╔═╡ 7ec8ae80-bdb9-46ea-adb5-0f76e861e50d
total_zero_v_pos = sum(false_positive_v)

# ╔═╡ 3560a592-71fa-4817-82b1-3946842d6f34
function sensitivity_specificity()
    f = Figure()
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:3, ["Material Decomposition", "Volume Fraction", "Agatston"]),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_zero_m / total_cac) * 100
	h2 = (total_zero_v / total_cac) * 100
	h3 = (num_zero_a / total_cac) * 100 
	heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[1:3], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)

	##-- B --##
    ax = Axis(
		f[2, 1]; 
		xticks = (1:3, ["Material Decomposition", "Volume Fraction", "Agatston"]),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_zero_m_pos / total_cac_pos) * 100
	h2 = (total_zero_v_pos / total_cac_pos) * 100
	h3 = (total_zero_a_pos / total_cac_pos) * 100
    heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[1:3], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	save(plotsdir("sensitivity_specificity.eps"), f)

    return f
end

# ╔═╡ d8a3e7f0-65d8-47d7-885a-f3fa6ba63739
with_theme(sensitivity_specificity, medphys_theme)

# ╔═╡ 07bab099-72e6-4315-aeac-ee7acc73f3dc
total_zero_m_pos, total_zero_a_pos

# ╔═╡ Cell order:
# ╠═a746ff38-1cde-4853-b001-7f2dfd6a693c
# ╠═7cad0e61-2cde-44ef-903a-531f10b322f5
# ╠═b169b3f6-819c-4f29-b32d-d5b28b9ae85f
# ╠═4ed99e12-126a-4728-a700-6942c678b403
# ╠═debe4f43-1698-461d-95c6-17f6052eee36
# ╠═7f4c4f06-9875-4fc5-bf0a-77fe243067ec
# ╠═7a4b45d8-2fbd-41d1-9df6-c31d15361327
# ╠═58f66f99-006d-4a14-a2a4-3342e764c01d
# ╟─b4bc009b-c9dc-40f4-9444-f12cae95289c
# ╠═d9071ce3-0109-4979-9474-a9256acc1b58
# ╠═8212d99a-a5c8-4e44-a2e4-8a6e598b9367
# ╠═20080b7f-bc37-453a-a426-c478ec805076
# ╠═28cfcdf4-abd7-4699-bcf2-7d4af05340e1
# ╠═684af87a-a9f2-42fd-be26-40a12719bf22
# ╠═3d5c68d6-9e5c-430d-b17a-552bbd4d254d
# ╟─190d1636-41ab-4e98-ae6c-5d872a742c0f
# ╟─ebd72800-16d2-42fc-9b7d-816efcbb675a
# ╟─7ffa22cc-8214-4d39-9ff0-20b0b63f06cb
# ╟─db529cfe-fd35-4a35-b9fd-8367bd7fb6fe
# ╠═5eb31f2b-bfaa-41d0-9375-9a076ba789d1
# ╟─1f39b6b5-b839-4dff-bbf0-10b3aa7b84ef
# ╟─2cff3bc9-1f37-4e56-a02f-b61c4cceb6f1
# ╠═d2b062b3-30dc-454a-b2fa-56c41312cc38
# ╠═9b52009a-4ef7-4f7f-aadb-89320a32f906
# ╠═47aa7f6c-067d-40c1-94b5-5bfdd30bb607
# ╠═17dad848-6491-45e4-8eec-3ed3e08c5318
# ╠═12edad23-47d8-4c30-b595-0e48380868dd
# ╠═fb9da5d1-a8cc-4485-9cf6-ab9c98c11816
# ╟─b8da1153-ae0a-4baa-9576-dffd76d2267b
# ╠═fdd43925-fc6a-4e26-8aea-8beb6c2c21fe
# ╠═0691ad9a-10c4-489b-a6bf-03dd46ef7f09
# ╟─0e9db9e2-5a41-48a6-9d11-7c38f172c131
# ╠═4234bb6c-d94a-418f-a6cd-2bf1df8b7aff
# ╠═99ea3200-1d7d-4f6f-9ec0-8dd9427aadfc
# ╠═364ade03-b649-43be-9580-85f0e2fae7bd
# ╟─881488ad-3eaf-46ed-8a86-7207df316ca4
# ╟─5302a029-1bd0-4582-99bb-d60832cc696c
# ╠═f127ec74-b7d3-4c83-bacd-40fcda4e14c7
# ╠═934bb2d3-9c17-4549-9dfb-e16ebb5963bc
# ╠═bb1c2f6f-0d5a-49c8-9133-2aa65485a405
# ╟─78a2bdd7-54c6-46a0-a9b4-a28f3c97448a
# ╠═7be46fdc-0554-48b2-aa8e-55770c68f16b
# ╠═0e428d01-39de-4274-a83d-1e2369ae35f5
# ╟─af083168-9989-407d-980b-9fbde3c8397c
# ╠═89445e61-9398-477f-afc0-5c69ea6e3212
# ╠═7ec8ae80-bdb9-46ea-adb5-0f76e861e50d
# ╟─3560a592-71fa-4817-82b1-3946842d6f34
# ╟─d8a3e7f0-65d8-47d7-885a-f3fa6ba63739
# ╠═07bab099-72e6-4315-aeac-ee7acc73f3dc
