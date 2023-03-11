### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 89014bd0-3d2f-4ed5-a4ec-bfa36a1b722d
# ╠═╡ show_logs = false
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ 71621e0e-0b10-4398-95c1-c73cab737e5b
# ╠═╡ show_logs = false
begin
    using PlutoUI, CairoMakie, MAT, DICOM, DICOMUtils
end

# ╔═╡ 3982c500-52e3-4dad-ac99-e6423ecc3282
TableOfContents()

# ╔═╡ ee37a9bf-545e-4c68-b22d-0c3fce8df0ba
energies = [80, 135]

# ╔═╡ 03049b0b-9ac2-4d0e-b105-e36055cc6f42
densities = ["Density1", "Density2", "Density3"]

# ╔═╡ 753bf413-2094-4642-8c13-26d70c8f99f8
sizes_folders = ["Small", "Medium", "Large"]

# ╔═╡ ef8cd6e6-2ec2-4f57-b85d-54b489139812
sizes_folders1 = ["Small1","Medium1","Large1"]

# ╔═╡ 98f46c6e-da8c-4b1d-817e-b3285452f210
sizes = ["small", "medium", "large"]

# ╔═╡ 4eb56f0e-fbcf-40f4-befe-a5ae4a1c4a49
file_nums = [1,2,3]

# ╔═╡ 6e8c4d1a-08ad-4cf7-976c-e4b718951b53
begin
	## For files in SIZE ONLY
	for size_folder in sizes_folders
		for density in densities
    		for energy in energies
				for file_num in file_nums
					
					@info size_folder, density, energy, file_num
					local file_size
					if(size_folder == "Small")
						file_size = "small"
					elseif(size_folder == "Medium")
						file_size = "medium"
					else
						file_size = "large"
					end
					
					## Path to the original dataset of .mat files
					home_root = datadir(
						"mat_measurement_bone_marrow",
						"SIZE",
						size_folder
					)
					
					path_root = string(home_root,
						"/",
						density,
						"energy",
						string(energy),
						file_size,
						".mat"
					)
					
					vars1 = matread(path_root)
					array1 = vars1[string("I")]
					array1 = Int16.(round.(array1))
					
					## Path to known DICOM file				
					dcm_path = datadir("sample.dcm")
					
					dcm = dcm_parse(dcm_path)
					dcm[tag"Pixel Data"] = array1
					dcm[tag"Instance Number"] = file_num
					dcm[tag"Rows"] = size(array1, 1)
					dcm[tag"Columns"] = size(array1, 2)
	
					## Path to output the newly creted DICOM files
					output_root = datadir(
						"dcms_measurement_new",
						size_folder,
						density,
						string(energy)
					)
					
					if !isdir(output_root)
						mkpath(output_root)
					end
					output_path = datadir(
						"dcms_measurement_new",
						size_folder,
						density,
						string(energy),
						string(file_num)*".dcm"
					)
					dcm_write(output_path, dcm)
				end
			end
        end
    end

	## For files in SIZE1 ONLY
	for size_folder1 in sizes_folders1
		for density in densities
    		for energy in energies
				for file_num in file_nums

					@info size_folder1, density, energy, file_num
					local file_size
					if(size_folder1 == "Small1")
						file_size = "small"
					elseif(size_folder1 == "Medium1")
						file_size = "medium"
					else
						file_size = "large"
					end
					
					## Path to the original dataset of .mat files
					path_root = datadir(
						"mat_measurement_bone_marrow",
						"SIZE1",
						size_folder1,
						density*"energy"*string(energy)*file_size*".mat"	
					)
					
					vars1 = matread(path_root)
					array1 = vars1[string("I")]
					array1 = Int16.(round.(array1))
					
					## Path to known DICOM file
					dcm_path = datadir("sample.dcm")
					
					dcm = dcm_parse(dcm_path)
					dcm[tag"Pixel Data"] = array1
					dcm[tag"Instance Number"] = file_num
					dcm[tag"Rows"] = size(array1, 1)
					dcm[tag"Columns"] = size(array1, 2)
	
					## Path to output the newly creted DICOM files
					output_root = datadir(
						"dcms_measurement_new",
						size_folder1,
						density,
						string(energy)
					)
					output_path = joinpath(output_root, string(file_num) * ".dcm")
					dcm_write(output_path, dcm)
				end
			end
        end
    end
end

# ╔═╡ 40095bc7-a7a5-4455-9198-52a383e5434e
md"""
## Check DICOM image(s)
"""

# ╔═╡ 9162cae4-c86b-4aa4-aee7-5fce03a71747
output_path = datadir("dcms_measurement_new")

# ╔═╡ 36946509-a59e-4218-b501-cb9fc401fcc9
dcmdir_combined = dcmdir_parse(dirname(output_path));

# ╔═╡ 05a9b1a9-26c7-4472-9f3a-44e131cc84d1
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ c442c24c-2ab7-4dab-9bc1-95ec68d69996
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ f1a6dd2e-e9d2-4bda-aee6-0219b25a4098
heatmap(transpose(vol_combined[:, :, c]); colormap=:grays)

# ╔═╡ Cell order:
# ╠═89014bd0-3d2f-4ed5-a4ec-bfa36a1b722d
# ╠═71621e0e-0b10-4398-95c1-c73cab737e5b
# ╠═3982c500-52e3-4dad-ac99-e6423ecc3282
# ╠═ee37a9bf-545e-4c68-b22d-0c3fce8df0ba
# ╠═03049b0b-9ac2-4d0e-b105-e36055cc6f42
# ╠═753bf413-2094-4642-8c13-26d70c8f99f8
# ╠═ef8cd6e6-2ec2-4f57-b85d-54b489139812
# ╠═98f46c6e-da8c-4b1d-817e-b3285452f210
# ╠═4eb56f0e-fbcf-40f4-befe-a5ae4a1c4a49
# ╠═6e8c4d1a-08ad-4cf7-976c-e4b718951b53
# ╟─40095bc7-a7a5-4455-9198-52a383e5434e
# ╠═9162cae4-c86b-4aa4-aee7-5fce03a71747
# ╠═36946509-a59e-4218-b501-cb9fc401fcc9
# ╠═05a9b1a9-26c7-4472-9f3a-44e131cc84d1
# ╟─c442c24c-2ab7-4dab-9bc1-95ec68d69996
# ╠═f1a6dd2e-e9d2-4bda-aee6-0219b25a4098
