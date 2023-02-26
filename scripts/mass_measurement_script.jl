### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 110dbe7c-a5ac-433a-bea9-87cc9130d359
begin
    let
        using Pkg
        Pkg.activate(mktempdir())
        Pkg.Registry.update()
        Pkg.add("PlutoUI")
        Pkg.add("CairoMakie")
        Pkg.add("Statistics")
        Pkg.add("StatsBase")
        Pkg.add("ImageMorphology")
        Pkg.add("ImageFiltering")
		Pkg.add("TestImages")
        Pkg.add("CSV")
        Pkg.add("DataFrames")
        Pkg.add("GLM")
		Pkg.add("HypothesisTests")
        Pkg.add(; url="https://github.com/JuliaHealth/DICOM.jl")
        Pkg.add(; url="https://github.com/Dale-Black/DICOMUtils.jl")
        Pkg.add("ImageComponentAnalysis")
        Pkg.add("Tables")
    end

    using PlutoUI
    using CairoMakie
    using Statistics
    using StatsBase: quantile!
    using ImageMorphology
    using ImageFiltering
	using TestImages
    using CSV
    using DataFrames
    using GLM
	using DICOM
    using DICOMUtils
    using ImageComponentAnalysis
    using Tables
	using HypothesisTests
end

# ╔═╡ 9e2582eb-4500-45a7-bcee-b2ae25b15a6b
TableOfContents()

# ╔═╡ d1efb736-1148-437b-b13e-f771ee92c3cd
md"""
## Script
"""

# ╔═╡ f2b7e4a2-894c-4727-bd46-8ca3729a4622
function predict_concentration(x, y, p)
	A = p[1] + (p[2] * x) + (p[3] * y) + (p[4] * x^2) + (p[5] * x * y) + (p[6] * y^2)
	B = 1 + (p[7] * x) + (p[8] * y)
	F = A / B
end;

# ╔═╡ a25d3dbc-a715-41ac-8c0d-0d3d328fc859
begin
	r = "/Users/xings/Google Drive/Research/dual energy/Segmentation/small/mask_L_HD.csv"
	mask1 = Array(CSV.read(r, DataFrame; header=false))
end;

# ╔═╡ 636a8595-8e77-4471-9dd5-7af56be2a39a
begin
	SIZES = ["small","medium","large"]
	base_root = "/Users/xings/Google Drive/Research/dual energy/"
	KVs = [80, 135]
	#calc = [733,411,151,669,370,90,552,222,52]
	calc = [797,101,37,403,48,32,199,41,27]
	dfs = []
	for SIZE in SIZES
		intensityS135 = Vector{Float64}()  # high energy 
		intensityS80 = Vector{Float64}()   # low energy
		intensityM135 = Vector{Float64}()
		intensityM80 = Vector{Float64}()
		intensityL135 = Vector{Float64}()
		intensityL80 = Vector{Float64}()
		
		vol_S = Vector{Float64}()   # total volume of VOI
		vol_M = Vector{Float64}()
		vol_L = Vector{Float64}()
		
		# Load calibration parameters
		param_root = string(base_root,"Calibration/",SIZE, "Calibration.csv")
		params_df = DataFrame(CSV.File(param_root))

		# Load masks
		mask_root = string(base_root,"dilate_segmentation/",SIZE,"/")
		mask_L_HD = Array(CSV.read(string(mask_root, "mask_L_HD.csv"), DataFrame; header=false))
		mask_M_HD = Array(CSV.read(string(mask_root, "mask_M_HD.csv"), DataFrame; header=false))
	    mask_S_HD = Array(CSV.read(string(mask_root, "mask_S_HD.csv"), DataFrame; header=false))
		mask_L_MD = Array(CSV.read(string(mask_root, "mask_L_MD.csv"), DataFrame; header=false))
		mask_M_MD = Array(CSV.read(string(mask_root, "mask_M_MD.csv"), DataFrame; header=false))
		mask_S_MD = Array(CSV.read(string(mask_root, "mask_S_MD.csv"), DataFrame; header=false))
	    mask_L_LD = Array(CSV.read(string(mask_root, "mask_L_LD.csv"), DataFrame; header=false))
	    mask_M_LD = Array(CSV.read(string(mask_root, "mask_M_LD.csv"), DataFrame; header=false))
	    mask_S_LD = Array(CSV.read(string(mask_root, "mask_S_LD.csv"), DataFrame; header=false))


		
		# Draw VOI and Count number of voxels of VOI

		for ENERGY in KVs
			pth = string(base_root, "validation_originalsize_dcms/", SIZE,"1/", ENERGY, "/")
			header, dcm_array, _ = dcm_reader(pth)
			pixel_size = DICOMUtils.get_pixel_size(header)

			for i in 1:size(dcm_array,3)
				d = dcm_array[:,:,i]
				if ENERGY == 80
					push!(intensityS80,mean(d[mask_S_HD]))
					push!(intensityM80,mean(d[mask_M_HD]))
					push!(intensityL80,mean(d[mask_L_HD]))
						
					push!(intensityS80,mean(d[mask_S_MD]))
					push!(intensityM80,mean(d[mask_M_MD]))
					push!(intensityL80,mean(d[mask_L_MD]))
						
					push!(intensityS80,mean(d[mask_S_LD]))
					push!(intensityM80,mean(d[mask_M_LD]))
					push!(intensityL80,mean(d[mask_L_LD]))

					push!(vol_S,count(mask_S_HD)*pixel_size[1]*pixel_size[2]*pixel_size[3])
					push!(vol_S,count(mask_S_MD)*pixel_size[1]*pixel_size[2]*pixel_size[3])
					push!(vol_S,count(mask_S_LD)*pixel_size[1]*pixel_size[2]*pixel_size[3])
					
					push!(vol_M,count(mask_M_HD)*pixel_size[1]*pixel_size[2]*pixel_size[3])
					push!(vol_M,count(mask_M_MD)*pixel_size[1]*pixel_size[2]*pixel_size[3])
					push!(vol_M,count(mask_M_LD)*pixel_size[1]*pixel_size[2]*pixel_size[3])

					push!(vol_L,count(mask_L_HD)*pixel_size[1]*pixel_size[2]*pixel_size[3])
					push!(vol_L,count(mask_L_MD)*pixel_size[1]*pixel_size[2]*pixel_size[3])
					push!(vol_L,count(mask_L_LD)*pixel_size[1]*pixel_size[2]*pixel_size[3])
					
				else
					push!(intensityS135,mean(d[mask_S_HD]))
					push!(intensityM135,mean(d[mask_M_HD]))
					push!(intensityL135,mean(d[mask_L_HD]))
						
					push!(intensityS135,mean(d[mask_S_MD]))
					push!(intensityM135,mean(d[mask_M_MD]))
					push!(intensityL135,mean(d[mask_L_MD]))
						
					push!(intensityS135,mean(d[mask_S_LD]))
					push!(intensityM135,mean(d[mask_M_LD]))
					push!(intensityL135,mean(d[mask_L_LD]))
				end
			end
		end

		# Predict Density
		pred_ρ_S = Vector{Float64}()
		pred_ρ_M = Vector{Float64}()
		pred_ρ_L = Vector{Float64}()

		for j in 1:size(intensityS135,1)
			ρ_S = predict_concentration(intensityS80[j],intensityS135[j],params_df[!,1])
			ρ_M = predict_concentration(intensityM80[j],intensityM135[j],params_df[!,1])
			ρ_L = predict_concentration(intensityL80[j],intensityL135[j],params_df[!,1])
	
			push!(pred_ρ_S,ρ_S)
			push!(pred_ρ_M, ρ_M)
			push!(pred_ρ_L,ρ_L)
		end


		header, _, _ = dcm_reader(string(base_root, "validation_originalsize_dcms/", SIZE,"1/", 80, "/"))
		pixel_size = DICOMUtils.get_pixel_size(header)
		
		# Calculate ground truth mass
		known_mass_L = pixel_size[3]*2.5^2*π*calc*1e-03
		known_mass_M = pixel_size[3]*1.5^2*π*calc*1e-03
		known_mass_S = pixel_size[3]*0.5^2*π*calc*1e-03
		
		# Predict Mass
		pred_mass_S = Vector{Float64}()
		pred_mass_M = Vector{Float64}()
		pred_mass_L = Vector{Float64}()

		for k in 1:size(pred_ρ_S,1)
			push!(pred_mass_S,pred_ρ_S[k]*vol_S[k]*1e-03)
			push!(pred_mass_M,pred_ρ_M[k]*vol_M[k]*1e-03)
			push!(pred_mass_L,pred_ρ_L[k]*vol_L[k]*1e-03)
		end

		df = DataFrame(:known_density=>calc,:known_mass_S=>known_mass_S,:predicted_mass_S=>pred_mass_S,
		:known_mass_M=>known_mass_M,:predicted_mass_M=>pred_mass_M,
		:known_mass_L=>known_mass_L,:predicted_mass_L=>pred_mass_L)

		push!(dfs,df)

	end
end

# ╔═╡ 04e3639d-14e4-4db1-ab3e-3321b3690ee5
md"""
## Save dfs
"""

# ╔═╡ 85b9f16e-ffd5-474a-8e30-2f49fe731449
begin
	for i = 1:3
		output_path = string(base_root,"validation results/",SIZES[i],"_mass_originalSize.csv")
		#CSV.write(output_path, dfs[i])
	end
end

# ╔═╡ eaa17af2-209b-4240-98d3-72ecd70db329
begin
	small_df = dfs[1]
	med_df = dfs[2]
	large_df = dfs[3]
end;

# ╔═╡ 5a70a000-b714-4730-ac94-4d869aac5074
md"""
## Plots
"""

# ╔═╡ 203035d0-9f8b-4e93-9a4e-ae1bb5346c12
md"""
### Small Phantom
"""

# ╔═╡ 1109c410-18f7-40f0-a09f-f3c77d15c1a5
begin
	fig01 = Figure()
	ax01 = Axis(fig01[1,1])

	scatter!(small_df[!,:known_mass_S],small_df[!,:predicted_mass_S],label = "small insert")
	scatter!(small_df[!,:known_mass_M],small_df[!,:predicted_mass_M],label = "medium insert")
	scatter!(small_df[!,:known_mass_L],small_df[!,:predicted_mass_L],label = "large insert")
	lines!(0:50,0:50,color=:red,label = "unity")

	ax01.title = "Measurement vs Truth (All Inserts)"
	ax01.xlabel = "Ground Truth Mass (mg)"
	ax01.ylabel = "Predicted Mass (mg)"


	fig01[1, 2] = Legend(fig01, ax01, framevisible = false)
	
	fig01
end

# ╔═╡ a8433d82-48dc-483c-8acb-957012db0f73
begin
	linearRegressor11 = lm(@formula(predicted_mass_S ~ known_mass_S), small_df)
	linearFit11 = predict(linearRegressor11)
	m11 = linearRegressor11.model.pp.beta0[2]
	b11 = linearRegressor11.model.pp.beta0[1]
end;

# ╔═╡ 36e3e91a-1960-479d-bf16-48a20c782b43
begin
	fig11 = Figure()
	ax11 = Axis(fig11[1,1])

	scatter!(small_df[!,:known_mass_S],small_df[!,:predicted_mass_S])
	lines!(small_df[!,:known_mass_S],small_df[!,:known_mass_S],color=:red,label = "unity")
	lines!(small_df[!,:known_mass_S],linearFit11, color =:green, linestyle=:dashdot,label = "fitted line")

	ax11.title = "Measurement vs Truth (small insert)"
	ax11.xlabel = "Ground Truth Mass (mg)"
	ax11.ylabel = "Predicted Mass (mg)"

	fig11[1, 2] = Legend(fig11, ax11, framevisible = false)
	
	fig11
end

# ╔═╡ 9c574d58-188d-4695-964c-304588c803cb
begin
	linearRegressor21 = lm(@formula(predicted_mass_M ~ known_mass_M), small_df)
	linearFit21 = predict(linearRegressor21)
	m21 = linearRegressor21.model.pp.beta0[2]
	b21 = linearRegressor21.model.pp.beta0[1]
end;

# ╔═╡ f923b5d6-dde8-4298-ae92-41bc67bf397c
begin
	fig21 = Figure()
	ax21 = Axis(fig21[1,1])

	scatter!(small_df[!,:known_mass_M],small_df[!,:predicted_mass_M])
	lines!(small_df[!,:known_mass_M],small_df[!,:known_mass_M],color=:red,label = "unity")
	lines!(small_df[!,:known_mass_M],linearFit21, color =:green, linestyle=:dashdot,label = "fitted line")

	ax21.title = "Measurement vs Truth (medium insert)"
	ax21.xlabel = "Ground Truth Mass (mg)"
	ax21.ylabel = "Predicted Mass (mg)"

	fig21[1, 2] = Legend(fig21, ax21, framevisible = false)
	
	fig21
end

# ╔═╡ 80edf668-a19f-4130-bc47-9703f8392aee
begin
	linearRegressor31 = lm(@formula(predicted_mass_L ~ known_mass_L), small_df)
	linearFit31 = predict(linearRegressor31)
	m31 = linearRegressor31.model.pp.beta0[2]
	b31 = linearRegressor31.model.pp.beta0[1]
end;

# ╔═╡ 1a4b6690-b923-4db6-ab9d-87ac0e033b89
begin
	fig31 = Figure()
	ax31 = Axis(fig31[1,1])

	scatter!(small_df[!,:known_mass_L],small_df[!,:predicted_mass_L])
	lines!(small_df[!,:known_mass_L],small_df[!,:known_mass_L],color=:red,label = "unity")
	lines!(small_df[!,:known_mass_L],linearFit31, color =:green, linestyle=:dashdot,label = "fitted line")

	ax31.title = "Measurement vs Truth (large insert)"
	ax31.xlabel = "Ground Truth Mass (mg)"
	ax31.ylabel = "Predicted Mass (mg)"

	fig31[1, 2] = Legend(fig31, ax31, framevisible = false)
	
	fig31
end

# ╔═╡ a4293798-afad-4e82-a727-928e661132eb
md"""
### Medium Phantom
"""

# ╔═╡ 9e4ab6e4-f34c-4fda-a98b-b75f60955833
begin
	fig02 = Figure()
	ax02 = Axis(fig02[1,1])

	scatter!(med_df[!,:known_mass_S],med_df[!,:predicted_mass_S],label = "small insert")
	scatter!(med_df[!,:known_mass_M],med_df[!,:predicted_mass_M],label = "medium insert")
	scatter!(med_df[!,:known_mass_L],med_df[!,:predicted_mass_L],label = "large insert")
	lines!(0:50,0:50,color=:red,label = "unity")

	ax02.title = "Measurement vs Truth (All Inserts)"
	ax02.xlabel = "Ground Truth Mass (mg)"
	ax02.ylabel = "Predicted Mass (mg)"


	fig02[1, 2] = Legend(fig02, ax02, framevisible = false)
	
	fig02
end

# ╔═╡ 1680428c-45bb-4e12-9e50-b11304c2fbe7
begin
	linearRegressor12 = lm(@formula(predicted_mass_S ~ known_mass_S), med_df)
	linearFit12 = predict(linearRegressor12)
	m12 = linearRegressor12.model.pp.beta0[2]
	b12 = linearRegressor12.model.pp.beta0[1]
end;

# ╔═╡ 8f72ad01-919d-4452-8238-4256bde4c1f1
begin
	fig12 = Figure()
	ax12 = Axis(fig12[1,1])

	scatter!(med_df[!,:known_mass_S],med_df[!,:predicted_mass_S])
	lines!(med_df[!,:known_mass_S],med_df[!,:known_mass_S],color=:red,label = "unity")
	lines!(med_df[!,:known_mass_S],linearFit12, color =:green, linestyle=:dashdot,label = "fitted line")

	ax12.title = "Measurement vs Truth (small insert)"
	ax12.xlabel = "Ground Truth Mass (mg)"
	ax12.ylabel = "Predicted Mass (mg)"

	fig12[1, 2] = Legend(fig12, ax12, framevisible = false)
	
	fig12
end

# ╔═╡ ff8f655f-4b76-400e-84bb-d72a57e32e3e
begin
	linearRegressor22 = lm(@formula(predicted_mass_M ~ known_mass_M), med_df)
	linearFit22 = predict(linearRegressor22)
	m22 = linearRegressor22.model.pp.beta0[2]
	b22 = linearRegressor22.model.pp.beta0[1]
end;

# ╔═╡ 03cdecd0-fd78-4c8f-a991-78c616b1e805
begin
	fig22 = Figure()
	ax22 = Axis(fig22[1,1])

	scatter!(med_df[!,:known_mass_M],med_df[!,:predicted_mass_M])
	lines!(med_df[!,:known_mass_M],med_df[!,:known_mass_M],color=:red,label = "unity")
	lines!(med_df[!,:known_mass_M],linearFit22, color =:green, linestyle=:dashdot,label = "fitted line")

	ax22.title = "Measurement vs Truth (medium insert)"
	ax22.xlabel = "Ground Truth Mass (mg)"
	ax22.ylabel = "Predicted Mass (mg)"

	fig22[1, 2] = Legend(fig22, ax22, framevisible = false)
	
	fig22
end

# ╔═╡ 1de48ee3-2a45-4519-8976-a127250484a3
begin
	linearRegressor32 = lm(@formula(predicted_mass_L ~ known_mass_L), med_df)
	linearFit32 = predict(linearRegressor32)
	m32 = linearRegressor32.model.pp.beta0[2]
	b32 = linearRegressor32.model.pp.beta0[1]
end;

# ╔═╡ b1009439-ceaa-4b03-a9b6-612650c43b6b
begin
	fig32 = Figure()
	ax32 = Axis(fig32[1,1])

	scatter!(med_df[!,:known_mass_L],med_df[!,:predicted_mass_L])
	lines!(med_df[!,:known_mass_L],med_df[!,:known_mass_L],color=:red,label = "unity")
	lines!(med_df[!,:known_mass_L],linearFit32, color =:green, linestyle=:dashdot,label = "fitted line")

	ax32.title = "Measurement vs Truth (large insert)"
	ax32.xlabel = "Ground Truth Mass (mg)"
	ax32.ylabel = "Predicted Mass (mg)"

	fig32[1, 2] = Legend(fig32, ax32, framevisible = false)
	
	fig32
end

# ╔═╡ eaaf03b7-8a77-411c-8462-8b12123feffe
md"""
### Large Phantom
"""

# ╔═╡ 3ebb0112-5ac4-462e-b290-ea3f5e885d39
begin
	fig03 = Figure()
	ax03 = Axis(fig03[1,1])

	scatter!(large_df[!,:known_mass_S],large_df[!,:predicted_mass_S],label = "small insert")
	scatter!(large_df[!,:known_mass_M],large_df[!,:predicted_mass_M],label = "medium insert")
	scatter!(large_df[!,:known_mass_L],large_df[!,:predicted_mass_L],label = "large insert")
	lines!(0:50,0:50,color=:red,label = "unity")

	ax03.title = "Measurement vs Truth (All Inserts)"
	ax03.xlabel = "Ground Truth Mass (mg)"
	ax03.ylabel = "Predicted Mass (mg)"


	fig03[1, 2] = Legend(fig03, ax03, framevisible = false)
	
	fig03
end

# ╔═╡ 798de15d-92f9-43aa-ad38-4e2c363776ce
begin
	linearRegressor13 = lm(@formula(predicted_mass_S ~ known_mass_S), large_df)
	linearFit13 = predict(linearRegressor13)
	m13 = linearRegressor13.model.pp.beta0[2]
	b13 = linearRegressor13.model.pp.beta0[1]
end;

# ╔═╡ f6dd68e0-2733-440d-b764-f8219511f267
begin
	fig13 = Figure()
	ax13 = Axis(fig13[1,1])

	scatter!(large_df[!,:known_mass_S],large_df[!,:predicted_mass_S])
	lines!(large_df[!,:known_mass_S],large_df[!,:known_mass_S],color=:red,label = "unity")
	lines!(large_df[!,:known_mass_S],linearFit13, color =:green, linestyle=:dashdot,label = "fitted line")

	ax13.title = "Measurement vs Truth (small insert)"
	ax13.xlabel = "Ground Truth Mass (mg)"
	ax13.ylabel = "Predicted Mass (mg)"

	fig13[1, 2] = Legend(fig13, ax13, framevisible = false)
	
	fig13
end

# ╔═╡ 85cdf99e-3df9-44bb-b8b5-c60414883a55
begin
	linearRegressor23 = lm(@formula(predicted_mass_M ~ known_mass_M), large_df)
	linearFit23 = predict(linearRegressor23)
	m23 = linearRegressor23.model.pp.beta0[2]
	b23 = linearRegressor23.model.pp.beta0[1]
end;

# ╔═╡ 56d16745-c163-42cc-9276-6b8e342431ee
begin
	fig23 = Figure()
	ax23 = Axis(fig23[1,1])

	scatter!(large_df[!,:known_mass_M],large_df[!,:predicted_mass_M])
	lines!(large_df[!,:known_mass_M],large_df[!,:known_mass_M],color=:red,label = "unity")
	lines!(large_df[!,:known_mass_M],linearFit23, color =:green, linestyle=:dashdot,label = "fitted line")

	ax23.title = "Measurement vs Truth (medium insert)"
	ax23.xlabel = "Ground Truth Mass (mg)"
	ax23.ylabel = "Predicted Mass (mg)"

	fig23[1, 2] = Legend(fig23, ax23, framevisible = false)
	
	fig23
end

# ╔═╡ 9be749be-9c81-4b9c-8d75-4410657f3587
begin
	linearRegressor33 = lm(@formula(predicted_mass_L ~ known_mass_L), large_df)
	linearFit33 = predict(linearRegressor33)
	m33 = linearRegressor33.model.pp.beta0[2]
	b33 = linearRegressor33.model.pp.beta0[1]
end;

# ╔═╡ 7b5a5038-5730-4f2f-883c-9e4295b1b231
begin
	fig33 = Figure()
	ax33 = Axis(fig33[1,1])

	scatter!(large_df[!,:known_mass_L],large_df[!,:predicted_mass_L])
	lines!(large_df[!,:known_mass_L],large_df[!,:known_mass_L],color=:red,label = "unity")
	lines!(large_df[!,:known_mass_L],linearFit33, color =:green, linestyle=:dashdot,label = "fitted line")

	ax33.title = "Measurement vs Truth (large insert)"
	ax33.xlabel = "Ground Truth Mass (mg)"
	ax33.ylabel = "Predicted Mass (mg)"

	fig33[1, 2] = Legend(fig33, ax33, framevisible = false)
	
	fig33
end

# ╔═╡ 0a3d783b-ca45-46ba-80f4-4874edf911e3
md"""
## RMSE & RMSD
"""

# ╔═╡ 5dd818a6-c3f9-4277-b9f8-500833491fbc
md"""
### Small Phantom
"""

# ╔═╡ 08e17c5a-4fce-4db5-9fd5-be7484017f71
begin
	RMSE_L1 = 0
	RMSE_M1 = 0
	RMSE_S1 = 0
	RMSD_L1 = 0
	RMSD_M1 = 0
	RMSD_S1 = 0
	RMSD1 = 0
	for i in 1:size(small_df[!,:known_mass_L],1)
		RMSE_L1 = RMSE_L1+(small_df[!,:known_mass_L][i]-small_df[!,:predicted_mass_L][i])^2
		RMSE_M1 = RMSE_M1+(small_df[!,:known_mass_M][i]-small_df[!,:predicted_mass_M][i])^2
		RMSE_S1 = RMSE_S1+(small_df[!,:known_mass_S][i]-small_df[!,:predicted_mass_S][i])^2

		RMSD_L1 = RMSD_L1+(small_df[!,:predicted_mass_L][i]-m31*small_df[!,:known_mass_L][i]+b31)^2
		RMSD_M1 = RMSD_M1+(small_df[!,:predicted_mass_M][i]-m21*small_df[!,:known_mass_M][i]+b21)^2
		RMSD_S1 = RMSD_S1+(small_df[!,:predicted_mass_S][i]-m11*small_df[!,:known_mass_S][i]+b11)^2

		RMSD1 = RMSD1+(small_df[!,:predicted_mass_L][i]-1.0002*small_df[!,:known_mass_L][i]-0.1642)^2+(small_df[!,:predicted_mass_L][i]-1.0002*small_df[!,:known_mass_L][i]-0.1642)^2+(small_df[!,:predicted_mass_L][i]-1.0002*small_df[!,:known_mass_L][i]-0.1642)^2
	end
	RMSE1 = sqrt((RMSE_L1+RMSE_M1+RMSE_S1)/(size(small_df[!,:known_mass_L],1)*3))
	RMSE_L1 = sqrt(RMSE_L1/size(small_df[!,:known_mass_L],1))
	RMSE_M1 = sqrt(RMSE_M1/size(small_df[!,:known_mass_M],1))
	RMSE_S1 = sqrt(RMSE_S1/size(small_df[!,:known_mass_S],1))

	RMSD1 = sqrt(RMSD1/(size(small_df[!,:known_mass_L],1)*3))
	RMSD_L1 = sqrt(RMSD_L1/size(small_df[!,:known_mass_L],1))
	RMSD_M1 = sqrt(RMSD_M1/size(small_df[!,:known_mass_M],1))
	RMSD_S1 = sqrt(RMSD_S1/size(small_df[!,:known_mass_S],1))

end;

# ╔═╡ 2660f2a1-11b4-4a0c-b13a-96191ef44e8f
md"""
### Medium Phantom
"""

# ╔═╡ 44edd386-ed6c-4616-a3fe-a8dd22166a60
begin
	RMSE_L2 = 0
	RMSE_M2 = 0
	RMSE_S2 = 0
	RMSD_L2 = 0
	RMSD_M2 = 0
	RMSD_S2 = 0
	RMSD2 = 0
	for i in 1:size(med_df[!,:known_mass_L],1)
		RMSE_L2 = RMSE_L2+(med_df[!,:known_mass_L][i]-med_df[!,:predicted_mass_L][i])^2
		RMSE_M2 = RMSE_M2+(med_df[!,:known_mass_M][i]-med_df[!,:predicted_mass_M][i])^2
		RMSE_S2 = RMSE_S2+(med_df[!,:known_mass_S][i]-med_df[!,:predicted_mass_S][i])^2

		RMSD_L2 = RMSD_L2+(med_df[!,:predicted_mass_L][i]-m32*med_df[!,:known_mass_L][i]+b32)^2
		RMSD_M2 = RMSD_M2+(med_df[!,:predicted_mass_M][i]-m22*med_df[!,:known_mass_M][i]+b22)^2
		RMSD_S2 = RMSD_S2+(med_df[!,:predicted_mass_S][i]-m12*med_df[!,:known_mass_S][i]+b12)^2

		RMSD2 = RMSD2+(med_df[!,:predicted_mass_L][i]-1.0032*med_df[!,:known_mass_L][i]+0.0579)^2+(med_df[!,:predicted_mass_L][i]-1.0032*med_df[!,:known_mass_L][i]+0.0579)^2+(med_df[!,:predicted_mass_L][i]-1.0032*med_df[!,:known_mass_L][i]+0.0579)^2
	end
	RMSE2 = sqrt((RMSE_L2+RMSE_M2+RMSE_S2)/(size(med_df[!,:known_mass_L],1)*3))
	RMSE_L2 = sqrt(RMSE_L2/size(med_df[!,:known_mass_L],1))
	RMSE_M2 = sqrt(RMSE_M2/size(med_df[!,:known_mass_M],1))
	RMSE_S2 = sqrt(RMSE_S2/size(med_df[!,:known_mass_S],1))

	RMSD2 = sqrt(RMSD2/(size(med_df[!,:known_mass_L],1)*3))
	RMSD_L2 = sqrt(RMSD_L2/size(med_df[!,:known_mass_L],1))
	RMSD_M2 = sqrt(RMSD_M2/size(med_df[!,:known_mass_M],1))
	RMSD_S2 = sqrt(RMSD_S2/size(med_df[!,:known_mass_S],1))

end;

# ╔═╡ eeca22ad-10ed-4410-a910-8694d478b612
md"""
### Large Phantom
"""

# ╔═╡ 19b193eb-4eb1-4fc9-8a72-35537a157736
begin
	RMSE_L3 = 0
	RMSE_M3 = 0
	RMSE_S3 = 0
	RMSD_L3 = 0
	RMSD_M3 = 0
	RMSD_S3 = 0
	RMSD3 = 0
	for i in 1:size(large_df[!,:known_mass_L],1)
		RMSE_L3 = RMSE_L3+(large_df[!,:known_mass_L][i]-large_df[!,:predicted_mass_L][i])^2
		RMSE_M3 = RMSE_M3+(large_df[!,:known_mass_M][i]-large_df[!,:predicted_mass_M][i])^2
		RMSE_S3 = RMSE_S3+(large_df[!,:known_mass_S][i]-large_df[!,:predicted_mass_S][i])^2

		RMSD_L3 = RMSD_L3+(large_df[!,:predicted_mass_L][i]-m33*large_df[!,:known_mass_L][i]+b33)^2
		RMSD_M3 = RMSD_M3+(large_df[!,:predicted_mass_M][i]-m23*large_df[!,:known_mass_M][i]+b23)^2
		RMSD_S3 = RMSD_S3+(large_df[!,:predicted_mass_S][i]-m13*large_df[!,:known_mass_S][i]+b13)^2

		RMSD3 = RMSD3+(large_df[!,:predicted_mass_L][i]-0.9993*large_df[!,:known_mass_L][i]-0.5624)^2+(large_df[!,:predicted_mass_L][i]-0.9993*large_df[!,:known_mass_L][i]-0.5624)^2+(large_df[!,:predicted_mass_L][i]-0.9993*large_df[!,:known_mass_L][i]-0.5624)^2
	end
	RMSE3 = sqrt((RMSE_L3+RMSE_M3+RMSE_S3)/(size(large_df[!,:known_mass_L],1)*3))
	RMSE_L3 = sqrt(RMSE_L3/size(large_df[!,:known_mass_L],1))
	RMSE_M3 = sqrt(RMSE_M3/size(large_df[!,:known_mass_M],1))
	RMSE_S3 = sqrt(RMSE_S3/size(large_df[!,:known_mass_S],1))

	RMSD3 = sqrt(RMSD3/(size(large_df[!,:known_mass_L],1)*3))
	RMSD_L3 = sqrt(RMSD_L3/size(large_df[!,:known_mass_L],1))
	RMSD_M3 = sqrt(RMSD_M3/size(large_df[!,:known_mass_M],1))
	RMSD_S3 = sqrt(RMSD_S3/size(large_df[!,:known_mass_S],1))

end;

# ╔═╡ Cell order:
# ╠═110dbe7c-a5ac-433a-bea9-87cc9130d359
# ╠═9e2582eb-4500-45a7-bcee-b2ae25b15a6b
# ╟─d1efb736-1148-437b-b13e-f771ee92c3cd
# ╠═f2b7e4a2-894c-4727-bd46-8ca3729a4622
# ╠═a25d3dbc-a715-41ac-8c0d-0d3d328fc859
# ╠═636a8595-8e77-4471-9dd5-7af56be2a39a
# ╟─04e3639d-14e4-4db1-ab3e-3321b3690ee5
# ╠═85b9f16e-ffd5-474a-8e30-2f49fe731449
# ╠═eaa17af2-209b-4240-98d3-72ecd70db329
# ╟─5a70a000-b714-4730-ac94-4d869aac5074
# ╟─203035d0-9f8b-4e93-9a4e-ae1bb5346c12
# ╟─1109c410-18f7-40f0-a09f-f3c77d15c1a5
# ╠═a8433d82-48dc-483c-8acb-957012db0f73
# ╟─36e3e91a-1960-479d-bf16-48a20c782b43
# ╠═9c574d58-188d-4695-964c-304588c803cb
# ╟─f923b5d6-dde8-4298-ae92-41bc67bf397c
# ╠═80edf668-a19f-4130-bc47-9703f8392aee
# ╟─1a4b6690-b923-4db6-ab9d-87ac0e033b89
# ╟─a4293798-afad-4e82-a727-928e661132eb
# ╟─9e4ab6e4-f34c-4fda-a98b-b75f60955833
# ╠═1680428c-45bb-4e12-9e50-b11304c2fbe7
# ╟─8f72ad01-919d-4452-8238-4256bde4c1f1
# ╠═ff8f655f-4b76-400e-84bb-d72a57e32e3e
# ╟─03cdecd0-fd78-4c8f-a991-78c616b1e805
# ╠═1de48ee3-2a45-4519-8976-a127250484a3
# ╟─b1009439-ceaa-4b03-a9b6-612650c43b6b
# ╟─eaaf03b7-8a77-411c-8462-8b12123feffe
# ╟─3ebb0112-5ac4-462e-b290-ea3f5e885d39
# ╠═798de15d-92f9-43aa-ad38-4e2c363776ce
# ╟─f6dd68e0-2733-440d-b764-f8219511f267
# ╠═85cdf99e-3df9-44bb-b8b5-c60414883a55
# ╟─56d16745-c163-42cc-9276-6b8e342431ee
# ╠═9be749be-9c81-4b9c-8d75-4410657f3587
# ╟─7b5a5038-5730-4f2f-883c-9e4295b1b231
# ╟─0a3d783b-ca45-46ba-80f4-4874edf911e3
# ╟─5dd818a6-c3f9-4277-b9f8-500833491fbc
# ╠═08e17c5a-4fce-4db5-9fd5-be7484017f71
# ╟─2660f2a1-11b4-4a0c-b13a-96191ef44e8f
# ╠═44edd386-ed6c-4616-a3fe-a8dd22166a60
# ╟─eeca22ad-10ed-4410-a910-8694d478b612
# ╠═19b193eb-4eb1-4fc9-8a72-35537a157736
