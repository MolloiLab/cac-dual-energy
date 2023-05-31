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

# ╔═╡ fe172646-1144-4313-a385-989f6a5e5784
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ 30a78917-86fe-4e7b-92f1-12746e0abbe3
begin
    using PlutoUI, CairoMakie, MAT, DICOM, DICOMUtils
end

# ╔═╡ c828ad36-6fa9-4a01-b7c7-0ffb9ef6a5cb
TableOfContents()

# ╔═╡ 542fbaa9-ab5f-4ef5-b48b-7177f26883f4
energies = [80, 135]

# ╔═╡ 3433372b-7bf6-4cc3-9203-f33b0b1485d2
densities = ["Density1", "Density2", "Density3"]

# ╔═╡ 4930fca9-77a7-4d8f-bb62-f3d6f21e4654
sizes_folders = ["Small", "Medium", "Large"]

# ╔═╡ 61010e80-bdbd-4a70-82c8-c02bc359cd74
sizes_folders1 = ["Small1","Medium1","Large1"]

# ╔═╡ 9f2e2eb5-0fa0-4c1d-a351-8951fd0beef3
sizes = ["small", "medium", "large"]

# ╔═╡ 439f6a90-5b03-42eb-9e34-2a2d29b17a8e
file_nums = 1:3

# ╔═╡ dd2f247c-900e-4691-b5c4-249fd2aa46ec
chop(lowercase(sizes_folders[1]))

# ╔═╡ 96d399f6-253d-4bbc-abae-3ba635d88b6f
begin
	## For files in SIZE ONLY
	for size_folder in sizes_folders
		for density in densities
    		for energy in energies
				for file_num in file_nums
					
					@info size_folder, density, energy, file_num
					file_size = lowercase(size_folder)

					## ------------- Inserts ------------- ##
					## Path to the original dataset of .mat files
					path = datadir(
						"mat_measurement_bone_marrow",
						"SIZE",
						size_folder,
						string(density, "energy", string(energy), file_size, ".mat")
					)
					vars1 = matread(path)
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
					output_path = datadir(
						"dcms_measurement_combined",
						size_folder,
						density,
						string(energy)
					)
					
					if !isdir(output_path)
						mkpath(output_path)
					end
					dcm_write(string(output_path, "/", file_num, ".dcm"), dcm)

					## ------------- Calibration Rod ------------- ##
					local pixel
					if size_folder == "Small"
						pixel = 30
					elseif size_folder == "Medium"
						pixel = 35
					elseif size_folder == "Large"
						pixel = 40
					end
					file_num = file_num + 3
					path = datadir(
						"mat_calibration_bone_marrow",
						size_folder,
						string("200rod", string(energy), "kV", pixel, ".mat")
					)
					vars1 = matread(path)
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
					output_path = datadir(
						"dcms_measurement_combined",
						size_folder,
						density,
						string(energy)
					)
					
					if !isdir(output_path)
						mkpath(output_path)
					end
					dcm_write(string(output_path, "/", file_num, ".dcm"), dcm)
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
					file_size = chop(lowercase(size_folder1))
					
					## ------------- Inserts ------------- ##
					## Path to the original dataset of .mat files
					path = datadir(
						"mat_measurement_bone_marrow",
						"SIZE1",
						size_folder1,
						string(density, "energy", string(energy), file_size, ".mat")
					)
					vars1 = matread(path)
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
					output_path = datadir(
						"dcms_measurement_combined",
						size_folder1,
						density,
						string(energy)
					)
					
					if !isdir(output_path)
						mkpath(output_path)
					end
					dcm_write(string(output_path, "/", file_num, ".dcm"), dcm)

					## ------------- Calibration Rod ------------- ##
					local pixel
					if size_folder1 == "Small1"
						pixel = 30
					elseif size_folder1 == "Medium1"
						pixel = 35
					elseif size_folder1 == "Large1"
						pixel = 40
					end
					file_num = file_num + 3
					path = datadir(
						"mat_calibration_bone_marrow",
						chop(size_folder1),
						string("200rod", string(energy), "kV", pixel, ".mat")
					)
					vars1 = matread(path)
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
					output_path = datadir(
						"dcms_measurement_combined",
						size_folder1,
						density,
						string(energy)
					)
					
					if !isdir(output_path)
						mkpath(output_path)
					end
					dcm_write(string(output_path, "/", file_num, ".dcm"), dcm)
				end
			end
        end
    end
end

# ╔═╡ 1e6b7cb2-1e3d-4ae8-b998-db1d0ecdd11c
md"""
## Check DICOM image(s)
"""

# ╔═╡ 9ab0c2dc-3cdd-4cb4-a12b-dc293b80e392
begin
	SIZES = ["Small", "Medium", "Large"]
	SIZE = SIZES[3]
end

# ╔═╡ cf4bd188-ff24-4576-90df-5c9c061f7ae4
output_path = datadir("dcms_measurement_combined", SIZE, "Density1", "80")

# ╔═╡ e5f772e8-cd4e-4bcc-b439-0e1955f6d214
dcmdir_combined = dcmdir_parse(output_path);

# ╔═╡ 2916512c-8a9c-463d-98b3-843c208ea234
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ 4715a5e1-5bdd-45af-a7ec-339582ed2d31
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=4, show_value=true)

# ╔═╡ b5a4ebe4-da2f-40a0-9e52-8ee3369b7cd8
if SIZE == "Small"
	centers = [187, 318]
elseif SIZE == "Medium"
	centers = [238, 365]
elseif SIZE == "Large"
	centers = [285, 415]
end

# ╔═╡ 619f9ba7-0e03-4fb3-82c6-94664346fa55
let
	f = Figure()
	ax = Axis(f[1, 1])
	heatmap!(vol_combined[:, :, c]; colormap=:grays)
	scatter!(centers[1]:centers[1], centers[2]:centers[2], ;markersize=5)
	f
end

# ╔═╡ Cell order:
# ╠═fe172646-1144-4313-a385-989f6a5e5784
# ╠═30a78917-86fe-4e7b-92f1-12746e0abbe3
# ╠═c828ad36-6fa9-4a01-b7c7-0ffb9ef6a5cb
# ╠═542fbaa9-ab5f-4ef5-b48b-7177f26883f4
# ╠═3433372b-7bf6-4cc3-9203-f33b0b1485d2
# ╠═4930fca9-77a7-4d8f-bb62-f3d6f21e4654
# ╠═61010e80-bdbd-4a70-82c8-c02bc359cd74
# ╠═9f2e2eb5-0fa0-4c1d-a351-8951fd0beef3
# ╠═439f6a90-5b03-42eb-9e34-2a2d29b17a8e
# ╠═dd2f247c-900e-4691-b5c4-249fd2aa46ec
# ╠═96d399f6-253d-4bbc-abae-3ba635d88b6f
# ╟─1e6b7cb2-1e3d-4ae8-b998-db1d0ecdd11c
# ╠═9ab0c2dc-3cdd-4cb4-a12b-dc293b80e392
# ╠═cf4bd188-ff24-4576-90df-5c9c061f7ae4
# ╠═e5f772e8-cd4e-4bcc-b439-0e1955f6d214
# ╠═2916512c-8a9c-463d-98b3-843c208ea234
# ╟─4715a5e1-5bdd-45af-a7ec-339582ed2d31
# ╠═b5a4ebe4-da2f-40a0-9e52-8ee3369b7cd8
# ╠═619f9ba7-0e03-4fb3-82c6-94664346fa55
