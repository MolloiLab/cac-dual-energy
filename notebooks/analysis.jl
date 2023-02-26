### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 321db384-5cbe-49b5-bf1e-80cc1a603699
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ a35241f6-289d-433c-9e1e-432d26fcc34b
using Markdown; using InteractiveUtils

# ╔═╡ 68b40210-7a3c-403e-a206-bcea384a7d27
begin
	using PlutoUI, Statistics, CSV, DataFrames, CairoMakie, Colors, GLM, MLJBase, Printf
	using StatsBase: quantile!, rmsd
end

# ╔═╡ debe4f43-1698-461d-95c6-17f6052eee36
include(srcdir("plot_utils.jl"));

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
		xlabelsize = 20,
		xticklabelsize = 20,
		ygridcolor = :gray,
		ygridwidth = 0.1,
		ylabelsize = 20,
		yticklabelsize = 20,
		bottomsplinecolor = :black,
		leftspinecolor = :black,
		titlesize = 30
	)
);

# ╔═╡ b4bc009b-c9dc-40f4-9444-f12cae95289c
md"""
# Accuracy
"""

# ╔═╡ d9071ce3-0109-4979-9474-a9256acc1b58
df_m = CSV.read(string(datadir("results"),"/","masses_matdecomp.csv"), DataFrame);

# ╔═╡ 92efa7e8-ad6e-4e84-a4ab-4773f8fa2be7
let
	df = df_m
	gt_array = vec(hcat(df[!, :ground_truth_mass_hd], df[!, :ground_truth_mass_md], df[!, :ground_truth_mass_ld]))
	calc_array = vec(hcat(df[!, :predicted_mass_hd], df[!, :predicted_mass_md], df[!, :predicted_mass_ld]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_i1
	model_i1 = lm(@formula(Y ~ X), data)
	global r2_1
	r2_1 = GLM.r2(model_i1)
	global rms_values1
	rms_values1 = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_i1))
	]
end

# ╔═╡ d837159a-b6f2-40b8-822c-ad6349ea081a
begin
	newX1 = DataFrame(X=collect(1:1000));
	pred_i1 = GLM.predict(model_i1, newX1)
end

# ╔═╡ 004faee6-9239-450d-995a-98037baf6536
co1 = coef(model_i1)

# ╔═╡ 8212d99a-a5c8-4e44-a2e4-8a6e598b9367
df_a = CSV.read(string(datadir("results"),"/","masses_agat.csv"), DataFrame);

# ╔═╡ cc32f947-3c21-490a-9fd6-625349e68843
let
	df = df_a
	gt_array = vec(hcat(df[!, :ground_truth_mass_hd], df[!, :ground_truth_mass_md], df[!, :ground_truth_mass_ld]))
	calc_array = vec(hcat(df[!, :calculated_agat_large], df[!, :calculated_agat_medium], df[!, :calculated_agat_small]))
	data = DataFrame(
		X = gt_array,
		Y= calc_array
	)
	global model_a
	model_a = lm(@formula(Y ~ X), data)
	global r2a
	r2a = GLM.r2(model_a)
	global rms_valuesa
	rms_valuesa = [
		rms(data[!, :X], data[!, :Y]),
		rmsd(data[!, :Y], GLM.predict(model_a))
	]
end

# ╔═╡ 8f93a665-5293-4c4e-9a7e-72c4c9094ba8
begin
	newX3 = DataFrame(X=collect(1:1000));
	pred_a = GLM.predict(model_a, newX3)
end

# ╔═╡ f2210674-17da-40b1-be01-e0dfc5c35d1b
co3 = coef(model_a)

# ╔═╡ fea99e7b-3e03-4b68-a8dd-7532a275fe27
function accuracy()
	f = Figure()

	##-- A --##
	ax = Axis(
		f[1, 1],
		xticks = [0, 25, 50, 75, 100, 125],
		yticks = [0, 25, 50, 75, 100, 125],
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Integrated Calcium Mass",
	)
	
	df = df_m
	sc1=scatter!(df[!, :ground_truth_mass_hd], df[!, :predicted_mass_hd])
	errorbars!(df[!, :ground_truth_mass_hd], df[!, :predicted_mass_hd], rms(df[!, :ground_truth_mass_hd], df[!, :predicted_mass_hd]))
	sc2=scatter!(df[!, :ground_truth_mass_md], df[!, :predicted_mass_md])
	errorbars!(df[!, :ground_truth_mass_md], df[!, :predicted_mass_md], rms(df[!, :ground_truth_mass_md], df[!, :predicted_mass_md]))
	sc3=scatter!(df[!, :ground_truth_mass_ld], df[!, :predicted_mass_ld], color=:red)
	errorbars!(df[!, :ground_truth_mass_ld], df[!, :predicted_mass_ld], rms(df[!, :ground_truth_mass_ld], df[!, :predicted_mass_ld]))
	ln1=lines!([-1000, 1000], [-1000, 1000])
	ln2=lines!(collect(1:1000), pred_i1, linestyle=:dashdot)
	create_textbox(f[1, 1], co1, r2_1, rms_values1)
	
	xlims!(ax, low=0, high=125)
	ylims!(ax, low=0, high=125)

	ax = Axis(
		f[2, 1],
		xticks = [0, 25, 50, 75, 100, 125],
		yticks = [0, 25, 50, 75, 100, 125],
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Agatson Scoring",
	)

	##-- B --##
	df = df_a
	scatter!(df[!, :ground_truth_mass_hd], df[!, :calculated_agat_large])
	errorbars!(df[!, :ground_truth_mass_hd], df[!, :calculated_agat_large], rms(df[!, :ground_truth_mass_hd], df[!, :calculated_agat_large]))
	scatter!(df[!, :ground_truth_mass_md], df[!, :calculated_agat_medium])
	errorbars!(df[!, :ground_truth_mass_md], df[!, :calculated_agat_medium], rms(df[!, :ground_truth_mass_md], df[!, :calculated_agat_medium]))
	scatter!(df[!, :ground_truth_mass_ld], df[!, :calculated_agat_small], color=:red)
	errorbars!(df[!, :ground_truth_mass_ld], df[!, :calculated_agat_small], rms(df[!, :ground_truth_mass_ld], df[!, :calculated_agat_small]))
	lines!([-1000, 1000], [-1000, 1000],)
	lines!(collect(1:1000), pred_a, linestyle=:dashdot)
	create_textbox(f[2, 1], co3, r2a, rms_valuesa)
	xlims!(ax, low=0, high=125)
	ylims!(ax, low=0, high=125)

	#-- LABELS --##
	f[1:2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end
	
	save(plotsdir("linear_reg.eps"), f)
	f
end

# ╔═╡ 190d1636-41ab-4e98-ae6c-5d872a742c0f
with_theme(medphys_theme) do
    accuracy()
end

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

# ╔═╡ 78e7c2b7-c7e7-4fd5-8a2a-e8d67e39d829
df_a

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
	df_a_ld, df_a_md, df_a_hd = groupby(df_a, :inserts)
	
	length(findall(x -> x <= 0, hcat(df_a_ld[!, :calculated_agat_large], df_a_ld[!, :calculated_agat_medium], df_a_ld[!, :calculated_agat_small]))), length(findall(x -> x <= 0, hcat(df_a_md[!, :calculated_agat_large], df_a_md[!, :calculated_agat_medium], df_a_md[!, :calculated_agat_small]))), length(findall(x -> x <= 0, hcat(df_a_hd[!, :calculated_agat_large], df_a_hd[!, :calculated_agat_medium], df_a_hd[!, :calculated_agat_small])))
end

# ╔═╡ b8da1153-ae0a-4baa-9576-dffd76d2267b
md"""
#### Integrated
"""

# ╔═╡ fdd43925-fc6a-4e26-8aea-8beb6c2c21fe
begin
	false_negative_i = []
	for i in 1:3:nrow(df_m)-2
		mean_i, std_i = mean(df_m[i:i+2, :mass_bkg]), std(df_m[i:i+2, :mass_bkg])*std_level 
		array_i = hcat(df_m[i:i+2, :predicted_mass_hd], df_i[i:i+2, :predicted_mass_md], df_i[i:i+2, :predicted_mass_ld]);
		neg = length(findall(x -> x <= mean_i + std_i, array_i))
		push!(false_negative_i, neg)
	end
end

# ╔═╡ 0691ad9a-10c4-489b-a6bf-03dd46ef7f09
total_zero_i = sum(false_negative_i)

# ╔═╡ 364ade03-b649-43be-9580-85f0e2fae7bd
total_zero_i, num_zero_a

# ╔═╡ 881488ad-3eaf-46ed-8a86-7207df316ca4
md"""
## False Positive
"""

# ╔═╡ 5302a029-1bd0-4582-99bb-d60832cc696c
md"""
#### Agatston
"""

# ╔═╡ f127ec74-b7d3-4c83-bacd-40fcda4e14c7
array_a_pos = df_a[!, :mass_bkg]

# ╔═╡ 934bb2d3-9c17-4549-9dfb-e16ebb5963bc
total_cac_pos = length(array_a_pos)

# ╔═╡ bb1c2f6f-0d5a-49c8-9133-2aa65485a405
total_zero_a_pos = length(findall(x -> x > 0, array_a_pos))

# ╔═╡ 78a2bdd7-54c6-46a0-a9b4-a28f3c97448a
md"""
#### Integrated
"""

# ╔═╡ 7be46fdc-0554-48b2-aa8e-55770c68f16b
begin
	false_positive_i = []
	for i in 1:3:nrow(df_i)-2
		mean_i, std_i = mean(df_i[i:i+2, :mass_bkg]), std(df_i[i:i+2, :mass_bkg])*std_level
		array_i_pos = df_i[i:i+2, :mass_bkg]
		pos = length(findall(x -> x > (mean_i + std_level), array_i_pos))
		push!(false_positive_i, pos)
	end
end

# ╔═╡ 0e428d01-39de-4274-a83d-1e2369ae35f5
total_zero_i_pos = sum(false_positive_i)

# ╔═╡ 3560a592-71fa-4817-82b1-3946842d6f34
function sensitivity_specificity()
    f = Figure()
    colors = Makie.wong_colors()

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (1:2, ["Integrated", "Agatston"]),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2]
	h1 = (total_zero_i / total_cac) * 100
	h2 = (num_zero_a / total_cac) * 100 
	heights1 = [h1, h2]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
    barplot!(table, heights1; color=colors[1:2], bar_labels=[l1, l2])

    ylims!(ax; low=0, high=100)

	##-- B --##
	ax = Axis(
		f[2, 1]; 
		xticks = (1:2, ["Integrated", "Agatston"]),
		title = "False-Positive (CAC>0)",
		ylabel = "False-Positive (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2]
	h1 = (total_zero_i_pos / total_cac_pos) * 100
	h2 = (total_zero_a_pos / total_cac_pos) * 100
    heights1 = [h1, h2]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
    barplot!(table, heights1; color=colors[1:2], bar_labels=[l1, l2])

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
with_theme(medphys_theme) do
    sensitivity_specificity()
end

# ╔═╡ 07bab099-72e6-4315-aeac-ee7acc73f3dc
total_zero_i_pos, total_zero_a_pos

# Cell order:
# ╠═f36e510f-ac3a-4986-bf56-d3b80a068820
# ╠═d750f7ce-a801-11ed-04b2-1575dfb3a26e
# ╠═f1ae7024-bc5e-43c0-ad67-d58787f67004
# ╠═6bd3b773-d4cd-430b-b63a-91cc04ff1c95
# ╠═c3ee18d6-847b-4eb8-ac80-888a2e69c12f
# ╠═4f43a4f5-4dd2-4366-ab5b-b8d0089b1dd9
# ╟─9cf7599f-667f-45cf-a111-95149a7b2b0a
# ╟─2d2bdc53-e225-4369-a82b-f56ac153263e
# ╟─aa9185b1-e0d6-46db-9867-d6e44bcf70f8
# ╠═f4579ec7-9013-4ed6-aa32-5fb6a7de82b9
# ╠═07b01b11-a479-4a1d-88e2-6bb3ce1952c4
# ╠═22e2507a-3153-4cde-b95c-ca50108dd78a
# ╠═38d3637a-3407-46d9-a91c-5526198f6ce9
# ╠═12c47081-bb76-4ba0-ba5d-af15b21991bd
# ╠═22994bc2-8b11-48a2-bf56-68e592b0d702
# ╠═6e9e50e7-66b8-460a-8ebf-d5cf83a2443d
# ╠═5fe1ae64-d676-479f-8d2d-463740c48ae5
# ╟─45d93a55-8f32-4a44-bdd7-abac89bb41cc
# ╟─452d14b2-cc8f-4814-9b0b-679492f8426a
# ╟─2915b2a1-416e-4365-b0b4-fd2b07c8c1a7
# ╠═2abe98b8-91a9-42ed-b18a-192e7922efd2
# ╠═124809eb-61aa-4151-8608-3907d235560f
# ╟─9dba29d0-0402-48da-94ff-34ef7b64f6e7
# ╠═781a4de9-6b5b-4791-87f1-25207398fefc
# ╟─0db59538-fc2f-438c-9d2e-858714ace83b
# ╟─dc47bcc1-66a1-4137-908a-8301021a2f72
# ╠═b3f2e899-36dc-45d7-a61d-de14f5a8836c
# ╠═9fbd3b35-6266-438e-97bd-579592f65727
# ╠═3ce554d6-28f6-42b0-8926-ad2f7543cf4f
# ╠═4d55c88e-9506-4f23-8d50-6251db07841d
# ╠═070879f7-2b29-44a2-9644-49e02f48c51e
# ╠═a43bdcbb-8abc-4e8c-83b4-29547f41ad39
# ╠═3424c91c-c592-4d1d-b226-b4d1ae4b3563
# ╟─ce8a2deb-76f5-4b58-a51e-7e02c77782b2
# ╠═20533ef0-7222-41af-a8ad-6744cf159b78
# ╠═273414ff-5e2a-4338-8527-58f0952d52ba
# ╟─ec1d2254-040f-40de-83b9-6306bc869a89
# ╟─5a3a85b4-545f-4514-88d2-5df807a25aff
# ╠═02d80c8b-0ccf-4e79-8cb0-27d97a8aa390
# ╠═a483459c-da31-4ae5-bf92-d585824c4c6a
# ╠═0544b3db-432e-406c-bc88-c3ebcc4c34a3
# ╟─2b5c1c6d-46cc-4474-8927-8309d769d617
# ╠═a497ba6b-61ee-48c6-b7af-e8454dbb613e
# ╠═233f365d-6ed2-4a4d-9533-ffde89972401

# ╔═╡ Cell order:
# ╠═a35241f6-289d-433c-9e1e-432d26fcc34b
# ╠═321db384-5cbe-49b5-bf1e-80cc1a603699
# ╠═68b40210-7a3c-403e-a206-bcea384a7d27
# ╠═debe4f43-1698-461d-95c6-17f6052eee36
# ╠═7f4c4f06-9875-4fc5-bf0a-77fe243067ec
# ╠═7a4b45d8-2fbd-41d1-9df6-c31d15361327
# ╠═58f66f99-006d-4a14-a2a4-3342e764c01d
# ╟─b4bc009b-c9dc-40f4-9444-f12cae95289c
# ╠═d9071ce3-0109-4979-9474-a9256acc1b58
# ╠═92efa7e8-ad6e-4e84-a4ab-4773f8fa2be7
# ╠═d837159a-b6f2-40b8-822c-ad6349ea081a
# ╠═004faee6-9239-450d-995a-98037baf6536
# ╠═8212d99a-a5c8-4e44-a2e4-8a6e598b9367
# ╠═cc32f947-3c21-490a-9fd6-625349e68843
# ╠═8f93a665-5293-4c4e-9a7e-72c4c9094ba8
# ╠═f2210674-17da-40b1-be01-e0dfc5c35d1b
# ╠═fea99e7b-3e03-4b68-a8dd-7532a275fe27
# ╠═190d1636-41ab-4e98-ae6c-5d872a742c0f
# ╠═7ffa22cc-8214-4d39-9ff0-20b0b63f06cb
# ╠═db529cfe-fd35-4a35-b9fd-8367bd7fb6fe
# ╠═5eb31f2b-bfaa-41d0-9375-9a076ba789d1
# ╠═1f39b6b5-b839-4dff-bbf0-10b3aa7b84ef
# ╠═2cff3bc9-1f37-4e56-a02f-b61c4cceb6f1
# ╠═78e7c2b7-c7e7-4fd5-8a2a-e8d67e39d829
# ╠═d2b062b3-30dc-454a-b2fa-56c41312cc38
# ╠═9b52009a-4ef7-4f7f-aadb-89320a32f906
# ╠═47aa7f6c-067d-40c1-94b5-5bfdd30bb607
# ╠═17dad848-6491-45e4-8eec-3ed3e08c5318
# ╠═12edad23-47d8-4c30-b595-0e48380868dd
# ╠═fb9da5d1-a8cc-4485-9cf6-ab9c98c11816
# ╠═b8da1153-ae0a-4baa-9576-dffd76d2267b
# ╠═fdd43925-fc6a-4e26-8aea-8beb6c2c21fe
# ╠═0691ad9a-10c4-489b-a6bf-03dd46ef7f09
# ╠═364ade03-b649-43be-9580-85f0e2fae7bd
# ╠═881488ad-3eaf-46ed-8a86-7207df316ca4
# ╠═5302a029-1bd0-4582-99bb-d60832cc696c
# ╠═f127ec74-b7d3-4c83-bacd-40fcda4e14c7
# ╠═934bb2d3-9c17-4549-9dfb-e16ebb5963bc
# ╠═bb1c2f6f-0d5a-49c8-9133-2aa65485a405
# ╠═78a2bdd7-54c6-46a0-a9b4-a28f3c97448a
# ╠═7be46fdc-0554-48b2-aa8e-55770c68f16b
# ╠═0e428d01-39de-4274-a83d-1e2369ae35f5
# ╠═3560a592-71fa-4817-82b1-3946842d6f34
# ╠═d8a3e7f0-65d8-47d7-885a-f3fa6ba63739
# ╠═07bab099-72e6-4315-aeac-ee7acc73f3dc
