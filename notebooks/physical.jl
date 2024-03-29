### A Pluto.jl notebook ###
# v0.19.36

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

# ╔═╡ 27f94bff-b27c-4ff2-aecb-ce1ff005320e
# ╠═╡ show_logs = false
using Pkg; Pkg.instantiate()

# ╔═╡ 1e230606-179f-4875-aa86-ddb63f62b5ac
using DrWatson

# ╔═╡ aa37c780-d6c7-4b42-ab26-0ce4b8b2d82d
# ╠═╡ show_logs = false
@quickactivate "cac-dual-energy"

# ╔═╡ 6648546a-35b1-40df-8493-79284a46f609
using PlutoUI: bind, TableOfContents, Slider

# ╔═╡ cd5ede11-8d2b-43a9-bb6f-5768dd3e7f33
# ╠═╡ show_logs = false
using CairoMakie: Figure, Axis, heatmap!, scatter!, hist!, heatmap

# ╔═╡ 33307f19-7302-4c79-9c69-573996c678d6
using ImageMorphology: label_components, component_centroids, component_boxes

# ╔═╡ f5fe79b0-eede-49b4-bec8-6ce8d61026b5
using ImageMorphology: dilate

# ╔═╡ 9d453f1c-dc43-4176-a7fa-4ea11c7ae3ee
using Statistics: mean, std

# ╔═╡ 00f3a2e5-9a32-4358-8611-b06c5b0b566b
using StatsBase: countmap, quantile

# ╔═╡ 76908202-063f-4a4c-b4aa-c1349da6adb7
using Clustering: kmeans, assignments

# ╔═╡ c41f7045-cf79-44e8-9bc5-6c7c01af0196
using LinearAlgebra: norm, normalize, dot, atan

# ╔═╡ 5f0d1e1a-d7ba-4cf0-94f4-c2c160545e6b
using StaticArrays: SVector

# ╔═╡ 1ca1fc19-9eb4-4f3c-a8e5-9b059bec02bf
using CalciumScoring: VolumeFraction, Agatston, score

# ╔═╡ d9d39c0c-a1fc-4294-b5bf-ba2bfdcc3068
using DICOM: dcmdir_parse

# ╔═╡ f0ef6fa8-330e-43f7-82b7-9c8e3d95f3f9
using DICOM # this is for tag".." functionality

# ╔═╡ 775b9988-5929-4386-a1fd-e83c792f5d70
using MaterialDecomposition: fit_calibration

# ╔═╡ 7d2e63c6-b81f-4e59-b52e-984fd621f643
using DataFrames: DataFrame

# ╔═╡ baecaf71-d7e0-4421-a182-daa7b6516fb7
using MaterialDecomposition: quantify

# ╔═╡ de0183b1-a794-4260-b3dc-a84ac4d889f5
# ╠═╡ show_logs = false
include(srcdir("active_contours.jl"));

# ╔═╡ 64de6299-c69c-409d-bcb9-dd5ce09ae0b5
TableOfContents()

# ╔═╡ d7008c7d-1d26-4f62-8d01-2a32242b7e54
md"""
# Calibration
"""

# ╔═╡ cbf70317-1ec3-45e0-bee9-831df955519b
md"""
## Load
"""

# ╔═╡ 9ab1fb7e-f7b0-4d55-84a1-41178bf92e34
begin
	cal_path_le = "/Volumes/USB DISK/calibration/Spectral 0.5  50 Mono/4"
	cal_path_he = "/Volumes/USB DISK/calibration/Spectral 0.5  70 Mono/7"
end

# ╔═╡ 37982cc4-d1c8-4485-b1fe-23704ffed31b
begin
	dcms_cal_le = dcmdir_parse(cal_path_le)
	dcms_cal_he = dcmdir_parse(cal_path_he)
end

# ╔═╡ cd00f90b-93ac-4b8d-8695-f70c1143418e
begin
	header_cal_le = dcms_cal_le[1].meta
	header_cal_he = dcms_cal_he[1].meta
end;

# ╔═╡ 011d3dcc-32f4-45f4-a15c-1589fad4d9ce
begin
	dcm_cal_le = header_cal_le[tag"PixelData"]
	dcm_cal_he = header_cal_he[tag"PixelData"]
end;

# ╔═╡ 943402ca-853a-49cb-a3a7-e19a414584ab
md"""
## Visualize
"""

# ╔═╡ 87c38e85-d9dd-4de2-af20-855689c56a7d
@bind z1_cal Slider(axes(dcm_cal_le, 3), default=div(size(dcm_cal_le, 3), 2), show_value=true)

# ╔═╡ 91305d5b-d176-4d3e-8645-70073e7301d4
let
	f = Figure(size = (700, 500))

	ax = Axis(
		f[1, 1],
		title = "Low Energy"
	)
	heatmap!(dcm_cal_le[:, :, z1_cal], colormap=:grays)

	ax = Axis(
		f[1, 2],
		title = "High Energy"
	)
	heatmap!(dcm_cal_he[:, :, z1_cal], colormap=:grays)
	# heatmap!(calibration_cv[:, :, z1_cal], colormap=(:jet, 0.5))

	f
end

# ╔═╡ 3540a1ab-d70b-4b0c-91dc-ff5ff7d101c5
md"""
## Segment
"""

# ╔═╡ 6a84ea41-9548-4595-88b2-482fd2120d95
hw = div(size(dcm_cal_le, 3), 2)

# ╔═╡ 5a063a55-2fe4-4a71-b672-365e1ff43464
offset = 10

# ╔═╡ 99b110a4-0698-4fd1-b4ef-45749453d991
begin
	overlaid_le = dcm_cal_le[:, :, (hw - offset):(hw + offset)]
	overlaid_slice_le = dcm_cal_le[:, :, hw]

	overlaid_he = dcm_cal_he[:, :, (hw - offset):(hw + offset)]
	overlaid_slice_he = dcm_cal_he[:, :, hw]
end;

# ╔═╡ 07c64082-b299-4dcd-b931-ed193357dc64
function find_octagon_center(A, B, C)
    mid_AB = (A .+ B) ./ 2
    mid_BC = (B .+ C) ./ 2
    
    # Handle vertical slopes and potential division by zero.
    slope_AB = (B[1] != A[1]) ? (B[2] - A[2]) / (B[1] - A[1]) : Inf
    slope_BC = (C[1] != B[1]) ? (C[2] - B[2]) / (C[1] - B[1]) : Inf
    
    perp_slope_AB = -1 / slope_AB
    perp_slope_BC = -1 / slope_BC
    
    c_AB = mid_AB[2] - perp_slope_AB * mid_AB[1]
    c_BC = mid_BC[2] - perp_slope_BC * mid_BC[1]
    
    A_matrix = [1 -perp_slope_AB; 1 -perp_slope_BC]
    b_vector = [c_AB; c_BC]
    
    center = A_matrix \ b_vector
    return center
end

# ╔═╡ 0d273cb4-58e5-4356-beb7-e7898f8e88a3
rad = 7

# ╔═╡ 93e6bb48-b6bd-4f0f-982b-0148b2dbddc5
function create_sphere_mask(img::AbstractArray{T, 3}, centroid, radius) where T
    # Initialize mask with all zeros
    mask = zeros(Bool, size(img))

    # Define the center of the sphere (x0, y0, z0)
    x0, y0, z0 = centroid..., div(size(img, 3), 2)

    # Set all voxels inside the sphere to 1 and all voxels outside the sphere to 0
    for x in axes(img, 1), y in axes(img, 2), z in axes(img, 3)
        if ((x - x0)^2 + (y - y0)^2 + (z - z0)^2) <= radius^2
            mask[x, y, z] = true
        end
    end
    return mask
end

# ╔═╡ d9eba74f-9709-4725-b0fa-f54e618b4cea
md"""
## Fit Calibration
"""

# ╔═╡ d08522ac-cf9c-41c2-9dcf-d2875e9987c2
densities_cal = [0.0, 0.050, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600]

# ╔═╡ 5494e7c6-4f57-4e20-a2ec-964580b69740
md"""
# Load DICOM(s)
"""

# ╔═╡ 0b3c923a-2440-446f-8169-1c3813012d5d
begin
	output_path_le = "/Volumes/USB DISK/e/Spectral 0.5  50 Mono/7"
	output_path_he = "/Volumes/USB DISK/e/Spectral 0.5  70 Mono/12"
end

# ╔═╡ f68dd63d-fec7-48cb-b49d-536fa1eb0a14
begin
	dcms_le = dcmdir_parse(output_path_le)
	dcms_he = dcmdir_parse(output_path_he)
end

# ╔═╡ 118433ed-6b47-4906-b8fe-1b4602aa50be
begin
	header_le = dcms_le[1].meta
	header_he = dcms_he[1].meta
end;

# ╔═╡ 59acf8dd-0157-44c2-a7e8-b4c45abe00b2
begin
	dcm_arr_le = header_le[tag"PixelData"]
	dcm_arr_he = header_he[tag"PixelData"]
end;

# ╔═╡ d9030aa1-1ce1-4ed4-8863-7c681784ca2e
pixel_size = [0.5, 0.5, 0.5]

# ╔═╡ 4511769f-07af-4582-8f3d-f43ee8de48b4
md"""
## Visualize
"""

# ╔═╡ a4ef8928-7c3f-447f-b2d0-540d7e2d0fab
@bind a Slider(axes(dcm_arr_le, 3), default=160, show_value=true)

# ╔═╡ 7acaae9a-d31a-4e7f-9084-f97b71e9275a
let
	f = Figure(size = (800, 500))

	ax = Axis(
		f[1, 1],
		title = "Low Energy"
	)
	heatmap!(dcm_arr_le[:, :, a], colormap=:grays)

	ax = Axis(
		f[1, 2],
		title = "High Energy"
	)
	heatmap!(dcm_arr_le[:, :, a], colormap=:grays)


	f
end

# ╔═╡ 2ff1c49a-7df0-499b-8f41-9c519dd31fb1
md"""
# Mask Heart Function
"""

# ╔═╡ 7ec1add5-7c6e-4bf6-8125-6efa5fe7c48e
function erode_mask(img, num_erosions)
	new_img = img
	i = 0
	while i < num_erosions
		new_img = erode(new_img)
		i += 1
	end
	return new_img
end

# ╔═╡ 77cdf367-1fe5-4458-a584-a21e26c3de44
function create_circle_mask(img::AbstractMatrix, centroids, radius)
    # initialize mask with all zeros
    mask = zeros(size(img))

    # define the center of the circle
    x0, y0 = centroids[1], centroids[2]

    # set all pixels inside the circle to 1 and all pixels outside the circle to 0
    for x in axes(img, 1), y in axes(img, 2)
        if ((x - x0)^2 + (y - y0)^2) <= radius^2
            mask[x, y] = 1
        else
            mask[x, y] = 0
        end
    end
    return Bool.(mask)
end

# ╔═╡ 300e3a7b-dc1a-43b3-a50a-c955cebdb09b
function centroids_from_mask(mask)
	cc_labels = label_components(mask)
	largest_connected_component, _ = sort(collect(pairs(countmap(cc_labels[cc_labels .!= 0]))), by=x->x[2], rev=true)
	largest_connected_indices = findall(cc_labels .== largest_connected_component[1])

	new_mask = zeros(size(mask))
	for i in largest_connected_indices
		new_mask[i] = 1
	end
	centroids = Int.(round.(component_centroids(label_components(new_mask))[end]))
end

# ╔═╡ 5d1e36c6-6941-486f-a5ab-242c843a3f06
function calculate_octagon_centroids(
	image; 
	thresh1 = 1900, thresh2 = 1500, thresh3 = 1100
)
    cal_mask1 = image .> thresh1
    cal_components1 = label_components(cal_mask1)
    centroid1 = centroids_from_mask(cal_components1)

    cal_mask2 = (image .> thresh2) .& (image .< thresh1 + 100)
    cal_components2 = label_components(cal_mask2)
    centroid2 = centroids_from_mask(cal_components2)

    cal_mask3 = (image .> thresh3) .& (image .< thresh2 + 200)
    cal_components3 = label_components(cal_mask3)
    centroid3 = centroids_from_mask(cal_components3)

    octagon_center = find_octagon_center(centroid1, centroid2, centroid3)
    radius = norm(centroid1 .- octagon_center)
    centroids_cal = [centroid1, centroid2, centroid3]

	last_angle = atan(centroid3[2] - octagon_center[2], centroid3[1] - octagon_center[1])
    for i in 1:5
		new_angle = last_angle - (π / 4) * (i)
        new_centroid_x = octagon_center[1] + radius * cos(new_angle)
        new_centroid_y = octagon_center[2] + radius * sin(new_angle)
        rounded_centroid_x = round(Int, new_centroid_x)
        rounded_centroid_y = round(Int, new_centroid_y)
        push!(centroids_cal, (rounded_centroid_x, rounded_centroid_y))
    end
    
    return centroids_cal
end

# ╔═╡ 2b8d8c6d-0aee-4691-b9bc-e4ffe041dd1a
centroids_cal = calculate_octagon_centroids(overlaid_slice_le);

# ╔═╡ 0b5565c2-4978-4484-8cdd-a39977beb623
masks_cal = [create_sphere_mask(overlaid_le, centroid, rad) for centroid in centroids_cal];

# ╔═╡ c9fbcb1f-5316-4b07-a166-67a6545d453c
begin
	intensities_cal_le = reverse([mean(overlaid_le[i]) for i in masks_cal])
	intensities_cal_he = reverse([mean(overlaid_he[i]) for i in masks_cal])

	# intensities_cal_le[1] = 60
	# intensities_cal_he[1] = 20
end

# ╔═╡ 9889a63b-6b31-45b4-9953-e9f112910851
df_cal = DataFrame(
    "Density Calibration (mg/cm^3)" => densities_cal,
    "Intensity Calibration (Low Energy)" => intensities_cal_le,
    "Intensity Calibration (High Energy)" => intensities_cal_he
)

# ╔═╡ 16c2b182-29cc-4d7d-9c83-ba41e3b3c5b3
begin
	p0 = zeros(8)
	intensities_total_cal = hcat(intensities_cal_le, intensities_cal_he)
	fit_params = fit_calibration(intensities_total_cal, densities_cal, p0)
	# fit_params = fit_calibration(intensities_total_cal[2:end, :], densities_cal[2:end], p0)
end

# ╔═╡ c7bbdbd9-6dac-4dd7-8615-72b701f75b86
md"""
Slice: $(@bind z3_cal Slider(axes(overlaid_le, 3), default=div(size(overlaid_le, 3), 2), show_value=true))

Mask Number: $(@bind i_cal Slider(eachindex(masks_cal), show_value=true))
"""

# ╔═╡ 5b7ac3e2-fda1-44ed-bec4-338ee6428a6f
let
	f = Figure()

	ax = Axis(
		f[1, 1],
		title = "Calibration & Masks"
	)
	heatmap!(overlaid_le[:, :, z3_cal], colormap=:grays)
	heatmap!(masks_cal[i_cal][:, :, z3_cal], colormap = (:jet, 0.50))
	# for i in eachindex(centroids_cal)
	# 	scatter!(centroids_cal[i])
	# end

	# scatter!(centroids_cal[i_cal])
	
	f
end

# ╔═╡ 5de49946-026c-48d0-a7d0-70c293032524
begin
	half_x, half_y = size(dcm_arr_le, 1) ÷ 2, size(dcm_arr_le, 2) ÷ 2
	init_circle = create_circle_mask(dcm_arr_le[:, :, 3], (half_x, half_y), 140)

	init_mask = BitArray(undef, size(dcm_arr_le))
	for z in axes(dcm_arr_le, 3)
		init_mask[:, :, z] = init_circle
	end

	init_mask = init_mask .* initial_level_set(size(init_mask))
end;

# ╔═╡ 3b0c49e1-b162-471f-ae2e-355c39356812
heart_cv = chan_vese(dcm_arr_le; init_level_set = init_mask);

# ╔═╡ 67d8627b-e486-4d99-805a-91c871e45a73
centroids = centroids_from_mask(heart_cv)

# ╔═╡ 3518dc88-328b-48a9-a4dd-a5d1ed919382
heart_mask = create_circle_mask(dcm_arr_le[:, :, 3], centroids, 100);

# ╔═╡ 7cea2363-611b-47b6-b0d0-a68a84127d81
# idxs = getindex.(findall(isone, heart_mask), [1 2]);

# ╔═╡ b4b62a45-e7c3-4f8b-9962-ee7584469d83
md"""
## Visualize
"""

# ╔═╡ a3e7c703-f71d-4f23-bc16-6c399c210e0d
@bind z2 Slider(axes(dcm_arr_le, 3), default=130, show_value=true)

# ╔═╡ 50e5afe8-e276-4e56-9034-f566869b5021
let
	f = Figure()

	ax = Axis(f[1, 1])
	heatmap!(dcm_arr_le[:, :, z2], colormap=:grays)
	heatmap!(heart_mask, colormap = (:jet, 0.3))

	f
end

# ╔═╡ 1eb04f5a-593e-466c-bc32-ea3cffe260d4
md"""
# Inserts
"""

# ╔═╡ 9777cf9f-a34f-4db4-9a85-c04d7f0ec34a
begin
	dcm_heart_le = dcm_arr_le .* heart_mask
	dcm_heart_he = dcm_arr_he .* heart_mask
end;

# ╔═╡ c2450d27-f296-4bbb-b958-34a3778b3c2f
@bind d Slider(axes(dcm_heart_le, 3), default=130, show_value=true)

# ╔═╡ 27cb6246-d6f7-4c21-9844-31d259a5568c
let
	f = Figure(size = (800, 500))

	Axis(
		f[1, 1],
		title = "Low Energy"
	)
	heatmap!(dcm_heart_le[:, :, d], colormap=:grays)

	Axis(
		f[1, 2],
		title = "High Energy"
	)
	heatmap!(dcm_heart_he[:, :, d], colormap=:grays)
	
	f
end

# ╔═╡ 1f91f759-a6b3-4095-9eb4-93944a28741f
md"""
## Endpoints
"""

# ╔═╡ eac3a8f2-7f1b-4e7c-b85e-a90151309602
function find_heart_endpoints(dcm_heart, air_threshold=-2000)
	# Find the indices of the elements in the array that are between 30 and 60
	selected_inds = findall(dcm_heart .<= air_threshold)
	
	# Create boolean array from cartesian indices
	local bool_arr4 = zeros(size(dcm_heart))
	for i in selected_inds
		bool_arr4[i] = 1
	end

	local cc_labels 
	while true
		# Use connected component labeling to identify and label all connected components
		cc_labels = label_components(bool_arr4)
		if length(unique(cc_labels)) == 3
			break
		end
		
		bool_arr4 = erode(bool_arr4)
		if (length(unique(cc_labels)) <= 2)
			@warn "Erosion failed; only background is left"
			break
		end
	end
	
	begining_slice = findlast(isone, cc_labels)[3] + 1
	end_slice = findfirst(x->x==2, cc_labels)[3] - 1
	return begining_slice, end_slice
end

# ╔═╡ 2d6c800c-10e6-47fb-a5da-600d4c294283
begining_slice, end_slice = find_heart_endpoints(dcm_heart_le)

# ╔═╡ 6f9baa30-8314-4f6c-a21f-1cb24ae924da
md"""
## Plane Fitting
"""

# ╔═╡ 194a806d-74c9-45a7-9078-f73400df421b
function kmeans_cleaning(points)
	points_t = points'
	kmeans_result = kmeans(points_t, 2)
	labels = assignments(kmeans_result)
	mean_distances = [mean([norm(points_t[:, i] - kmeans_result.centers[:, j]) for i in findall(x -> x == j, labels)]) for j in 1:2]
	main_cluster = argmin(mean_distances)
	cleaned_points = points_t[:, labels .== main_cluster]
	cleaned_points = cleaned_points'
end

# ╔═╡ b8ba8056-a03c-4e96-ac51-e206ead16257
function find_heart_plane(dcm_heart, endpoints; air_threshold = -300, std_threshold = 0.25)
	# Find the indices of the elements in the array that are less than
	# air_threshold and filter to exclude beginning and end slices
	remove = [collect(1:endpoints[1])..., collect(endpoints[2]:320)...]
	selected_inds = findall(dcm_heart .<= air_threshold)
	for r in remove
		selected_inds = [i for i in selected_inds if i.I[3] != r]
	end

	selected_inds = getindex.(selected_inds, [1 2 3])

	# Clean points
	filtered_points = kmeans_cleaning(selected_inds)

	# Create boolean array from cartesian indices
	bool_arr = zeros(size(dcm_heart))
	for i in eachrow(filtered_points)
		bool_arr[i...] = 1
	end

	# Use connected component labeling to identify and label all connected components
	cc_labels = label_components(bool_arr)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	return filtered_points, cc_labels
end

# ╔═╡ 6217bc62-64a3-404b-8b2a-5ebab33e6e4d
idx_plane, cc_labels = find_heart_plane(dcm_heart_le, (begining_slice, end_slice));

# ╔═╡ b8c0bc97-bf63-4a31-b500-a39ca8e77b3d
function extract_plane_points(cc_labels)
	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	# # Find the value with the most occurrences
	most_common_value, _ = sort(collect(pairs(counts)), by=x->x[2], rev=true)

	# # Find the indices of the most common value in the original array
	most_common_indices = findall(cc_labels .== most_common_value[1])
	most_common_indices = getindex.(most_common_indices, [1 2 3])

	pts = most_common_indices[1+100, :], most_common_indices[Int(round(end/2)), :], most_common_indices[end-100, :];
end;

# ╔═╡ c495ba88-fc99-4a78-ab85-c324cbca495b
centroids

# ╔═╡ 02f19951-5077-4744-9d18-4b99d6ec33e3
pts = extract_plane_points(cc_labels)

# ╔═╡ cea55170-b870-42dc-93c5-ce796039cc10
md"""
## Segment Calcium Inserts
"""

# ╔═╡ dece0067-9ddd-47c5-85f8-a016e8db620f
function get_insert_centers(mpr, threshold)
	# z = div(size(mpr, 3), 2)
	# mpr_slice = mpr[:, :, z]
	# mpr_slice_thresh = mpr_slice .> threshold
	mpr_slice_thresh = mpr .> threshold
	
	# Use connected component labeling to identify and label all connected components
	cc_labels = label_components(mpr_slice_thresh)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	# Find the value with the most occurrences
	most_common_value_a, most_common_value_b = sort(collect(pairs(counts)), by=x->x[2], rev=true)

	# Find the indices of the most common value in the original array
	most_common_indices_a = findall(cc_labels .== most_common_value_a[1])

	# Create boolean array from new cartesian indices
	bool_arr_a = zeros(size(mpr_slice_thresh))
	for i in most_common_indices_a
		bool_arr_a[i] = 1
	end
	centroids_a = Int.(round.(component_centroids(label_components(bool_arr_a))[end]))
	box_a = component_boxes(label_components(bool_arr_a))

	# Find the indices of the most common value in the original array
	most_common_indices_b = findall(cc_labels .== most_common_value_b[1])

	# Create boolean array from new cartesian indices
	bool_arr_b = zeros(size(mpr_slice_thresh))
	for i in most_common_indices_b
		bool_arr_b[i] = 1
	end
	centroids_b = Int.(round.(component_centroids(label_components(bool_arr_b))[end]))

	# centers_a, centers_b = [centroids_a..., z], [centroids_b..., z]
	centers_a, centers_b = centroids_a, centroids_b
	return centers_a, centers_b
	
end

# ╔═╡ 8bb9db58-d77d-4ed3-94a4-78894a0baa13
centers_a, centers_b = get_insert_centers(dcm_heart_le, 200);

# ╔═╡ 11328d2a-cec1-4c81-8558-686950717939
function _in_cylinder(
	pt::SVector{3, Int}, pt1::SVector{3, Float64}, pt2::SVector{3, Float64}, radius
)
    v = pt2 - pt1
    w = pt - pt1

    # Compute the dot product
    c1 = dot(w, v)
    if c1 <= 0
        return norm(w) <= radius
    end

    c2 = dot(v, v)
    if c2 <= c1
        return norm(pt - pt2) <= radius
    end

    # Compute the perpendicular distance
    b = c1 / c2
    pb = pt1 + b * v
    return norm(pt - pb) <= radius
end

# ╔═╡ 3af82fb8-a8b1-4f43-95ba-9260bc7beb5c
function create_cylinder(array, pt1, pt2, radius, offset)
    # Convert the points to static arrays
    pt1 = SVector{3, Float64}(pt1)
    pt2 = SVector{3, Float64}(pt2)

    # Compute the unit vector in the direction from pt1 to pt2
    direction = normalize(pt2 - pt1)

    # Adjust the endpoints of the cylinder by the offset
    pt1 = pt1 - offset * direction
    pt2 = pt2 + offset * direction

    # Initialize the 3D array
    cylinder = zeros(Int, size(array)...)
    # Iterate over the 3D array
    for k in axes(cylinder, 3)
        for j in axes(cylinder, 2)
            for i in axes(cylinder, 1)
                # Create a static vector for the current point
                pt = SVector{3, Int}(i, j, k)

                # Check if the current point is inside the cylinder
                if _in_cylinder(pt, pt1, pt2, radius)
                    cylinder[i, j, k] = 1
                end
            end
        end
    end
    return Bool.(cylinder)
end

# ╔═╡ 32c5fe98-c82a-44c8-ab88-2dce146fa28a
cylinder = create_cylinder(dcm_heart_le, centers_a, centers_b, 8, -25);

# ╔═╡ 31e877f0-b96d-467d-8c47-9733d69f2346
begin
	br = create_cylinder(dcm_heart_le, centers_a, centers_b, 12, -25);
	background_ring = Bool.(br .- cylinder)
end;

# ╔═╡ 08b9e104-57fa-4807-843e-5ef958079375
@bind z Slider(axes(dcm_heart_le, 3), default=div(size(dcm_heart_le, 3), 2), show_value=true)

# ╔═╡ 313a9618-f56d-49e3-aa07-4741421beb65
let
	f = Figure(size = (800, 800))

	ax = Axis(
		f[1, 1],
		title = "Insert Mask (Low Energy)"
	)
	heatmap!(dcm_heart_le[:, :, z]; colormap = :grays)
	heatmap!(cylinder[:, :, z]; colormap = (:jet, 0.25))

	ax = Axis(
		f[1, 2],
		title = "Background Mask (Low Energy)"
	)
	heatmap!(dcm_heart_le[:, :, z]; colormap = :grays)
	heatmap!(background_ring[:, :, z]; colormap = (:viridis, 0.25))

	ax = Axis(
		f[2, 1],
		title = "Insert Mask (High Energy)"
	)
	heatmap!(dcm_heart_he[:, :, z]; colormap = :grays)
	heatmap!(cylinder[:, :, z]; colormap = (:jet, 0.25))

	ax = Axis(
		f[2, 2],
		title = "Background Mask (High Energy)"
	)
	heatmap!(dcm_heart_he[:, :, z]; colormap = :grays)
	heatmap!(background_ring[:, :, z]; colormap = (:viridis, 0.25))


	f
end

# ╔═╡ 8ff22fa9-5891-4705-8f17-d6d95a36e4c2
md"""
## Segment Calibration Insert
"""

# ╔═╡ c8ed7211-5b67-44ce-bf8f-f6b65806fc57
begin
	binary_calibration = falses(size(dcm_heart_le))
	binary_calibration[centers_a...] = true
	binary_calibration = dilate(binary_calibration)
end;

# ╔═╡ 01dbc828-4fa7-48d4-8a0a-edb275e8f715
md"""
## Remove Outliers (Air)
"""

# ╔═╡ ce8494e6-43ab-40a5-82dd-ea64710184f9
function remove_outliers(vector)
    Q1 = quantile(vector, 0.25)
    Q3 = quantile(vector, 0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    return [x for x in vector if x > lower_bound]
end

# ╔═╡ 02fb9d18-4686-4c73-a2eb-eaf8b0f5c012
begin
	dcm_heart_clean_le = remove_outliers(dcm_heart_le[cylinder])
	dcm_heart_clean_he = remove_outliers(dcm_heart_he[cylinder])
end;

# ╔═╡ db4672ce-933a-4aa8-9220-1e263ea540dd
let
	f = Figure(size = (700, 1000))
	
	ax = Axis(
		f[1, 1],
		title = "Original (Low Energy)"
	)
	hist!(dcm_heart_le[cylinder])

	ax = Axis(
		f[1, 2],
		title = "Clean (Low Energy)"
	)
	hist!(dcm_heart_clean_le)

	ax = Axis(
		f[2, 1],
		title = "Original (High Energy)"
	)
	hist!(dcm_heart_he[cylinder])

	ax = Axis(
		f[2, 2],
		title = "Clean (High Energy)"
	)
	hist!(dcm_heart_clean_he)

	f
end

# ╔═╡ f2bd880f-e79d-4831-9f8d-65b5f4e1ba68
md"""
# Score
"""

# ╔═╡ 6eadda32-c20e-4429-ab4f-6480facd988e
md"""
## Ground Truth
"""

# ╔═╡ 83b3c3c3-4c91-48ed-909b-aafdce1427f9
begin
	gt_density = 0.10 # mg/mm^3

	# π * (diameter/2)^2 * length
	gt_volume = π * (5/2)^2 * 7 # mm3
	gt_mass = gt_density * gt_volume
end

# ╔═╡ 1f674c41-4c28-45f7-a9df-2fbb784ee115
md"""
## Material Decomposition
"""

# ╔═╡ 54fcf8a4-b397-4d38-a953-fdf575c44f71
voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]

# ╔═╡ 1a6b0fb7-6b00-45a4-b034-9d7767b6b386
begin
	intensity_le = mean(dcm_heart_le[cylinder])
	intensity_he = mean(dcm_heart_he[cylinder])
end

# ╔═╡ fff6da7e-6eb3-4e87-a48a-2309d001dc57
vol_dil = sum(cylinder[:, :, :, 1]) * voxel_size

# ╔═╡ 1ef192cb-4349-4ce3-8bd8-0590b7b8923f
md_mass = quantify(intensity_le, intensity_he, fit_params, vol_dil)

# ╔═╡ 32ef19be-5d24-495f-b1be-aa1f3b15531c
md_mass

# ╔═╡ be7a5329-1147-4889-a2d7-d07005369e4b
md"""
## Volume Fraction
"""

# ╔═╡ 4c130872-3a57-40a3-8c43-51e46b963b69
hu_calcium_400 = mean(dcm_heart_le[binary_calibration])

# ╔═╡ a52ce3c2-3a82-4124-a260-0e1497038cd2
std(dcm_heart_le[binary_calibration])

# ╔═╡ 22fcd352-feaa-450d-9aa4-0f8863a81b50
ρ_calcium_400 = 0.400 # mg/mm^3

# ╔═╡ f869c063-e18f-4936-92e7-7cadc22768de
begin
	hu_heart_tissue_bkg_le = mean(dcm_heart_le[background_ring])
	hu_heart_tissue_bkg_he = mean(dcm_heart_he[background_ring])
end;

# ╔═╡ 3e54e379-9513-4741-9ade-657741be2615
begin
	vf_mass_le = score(dcm_heart_clean_le, hu_calcium_400, hu_heart_tissue_bkg_le, voxel_size, ρ_calcium_400, VolumeFraction())
	vf_mass_he = score(dcm_heart_clean_he, hu_calcium_400, hu_heart_tissue_bkg_he, voxel_size, ρ_calcium_400, VolumeFraction())
end

# ╔═╡ 913d3972-e9f4-486d-8a82-973de4271d8c
vf_mass_le, vf_mass_he

# ╔═╡ 728cf85d-7200-41b7-850b-cc68f17f68cb
md"""
## Agatston
"""

# ╔═╡ 94919901-f29b-4b30-ad9a-5f5b0507aad3
begin
	overlayed_le = dcm_arr_le .* cylinder
	overlayed_he = dcm_arr_he .* cylinder
end;

# ╔═╡ fdcb9783-0652-4651-8900-84bc81a686a7
begin
	kV_le = 80 # 50 monoenergetic => 80 kV
	kV_he = 135 # 70 monoenergetic => 135 kV
end

# ╔═╡ 5ad5b33e-4f1a-40ad-851e-0fd170651006
mass_cal_factor = ρ_calcium_400 / hu_calcium_400

# ╔═╡ 8d7d7e3f-9e23-4c01-9b8b-41b1e10f2bf8
begin
	a_agatston_le, a_volume_le, a_mass_le = score(overlayed_le, pixel_size, mass_cal_factor, Agatston(); kV=kV_le)
	a_agatston_he, a_volume_he, a_mass_he = score(overlayed_he, pixel_size, mass_cal_factor, Agatston(); kV=kV_he)
end

# ╔═╡ a6d47dee-b7a0-4d84-af95-200ba3eeab28
a_mass_le, a_mass_he

# ╔═╡ f0c86451-a6be-47c1-95a4-967e5a90b4ce
md"""
# Results
"""

# ╔═╡ Cell order:
# ╠═1e230606-179f-4875-aa86-ddb63f62b5ac
# ╠═aa37c780-d6c7-4b42-ab26-0ce4b8b2d82d
# ╠═27f94bff-b27c-4ff2-aecb-ce1ff005320e
# ╠═6648546a-35b1-40df-8493-79284a46f609
# ╠═cd5ede11-8d2b-43a9-bb6f-5768dd3e7f33
# ╠═33307f19-7302-4c79-9c69-573996c678d6
# ╠═f5fe79b0-eede-49b4-bec8-6ce8d61026b5
# ╠═9d453f1c-dc43-4176-a7fa-4ea11c7ae3ee
# ╠═00f3a2e5-9a32-4358-8611-b06c5b0b566b
# ╠═76908202-063f-4a4c-b4aa-c1349da6adb7
# ╠═c41f7045-cf79-44e8-9bc5-6c7c01af0196
# ╠═5f0d1e1a-d7ba-4cf0-94f4-c2c160545e6b
# ╠═1ca1fc19-9eb4-4f3c-a8e5-9b059bec02bf
# ╠═d9d39c0c-a1fc-4294-b5bf-ba2bfdcc3068
# ╠═f0ef6fa8-330e-43f7-82b7-9c8e3d95f3f9
# ╠═de0183b1-a794-4260-b3dc-a84ac4d889f5
# ╠═64de6299-c69c-409d-bcb9-dd5ce09ae0b5
# ╟─d7008c7d-1d26-4f62-8d01-2a32242b7e54
# ╟─cbf70317-1ec3-45e0-bee9-831df955519b
# ╠═9ab1fb7e-f7b0-4d55-84a1-41178bf92e34
# ╠═37982cc4-d1c8-4485-b1fe-23704ffed31b
# ╠═cd00f90b-93ac-4b8d-8695-f70c1143418e
# ╠═011d3dcc-32f4-45f4-a15c-1589fad4d9ce
# ╟─943402ca-853a-49cb-a3a7-e19a414584ab
# ╟─87c38e85-d9dd-4de2-af20-855689c56a7d
# ╟─91305d5b-d176-4d3e-8645-70073e7301d4
# ╟─3540a1ab-d70b-4b0c-91dc-ff5ff7d101c5
# ╠═6a84ea41-9548-4595-88b2-482fd2120d95
# ╠═5a063a55-2fe4-4a71-b672-365e1ff43464
# ╠═99b110a4-0698-4fd1-b4ef-45749453d991
# ╠═07c64082-b299-4dcd-b931-ed193357dc64
# ╠═5d1e36c6-6941-486f-a5ab-242c843a3f06
# ╠═2b8d8c6d-0aee-4691-b9bc-e4ffe041dd1a
# ╠═0d273cb4-58e5-4356-beb7-e7898f8e88a3
# ╠═93e6bb48-b6bd-4f0f-982b-0148b2dbddc5
# ╠═0b5565c2-4978-4484-8cdd-a39977beb623
# ╠═c9fbcb1f-5316-4b07-a166-67a6545d453c
# ╠═32ef19be-5d24-495f-b1be-aa1f3b15531c
# ╟─c7bbdbd9-6dac-4dd7-8615-72b701f75b86
# ╟─5b7ac3e2-fda1-44ed-bec4-338ee6428a6f
# ╟─d9eba74f-9709-4725-b0fa-f54e618b4cea
# ╠═775b9988-5929-4386-a1fd-e83c792f5d70
# ╠═7d2e63c6-b81f-4e59-b52e-984fd621f643
# ╠═d08522ac-cf9c-41c2-9dcf-d2875e9987c2
# ╠═9889a63b-6b31-45b4-9953-e9f112910851
# ╠═16c2b182-29cc-4d7d-9c83-ba41e3b3c5b3
# ╟─5494e7c6-4f57-4e20-a2ec-964580b69740
# ╠═0b3c923a-2440-446f-8169-1c3813012d5d
# ╠═f68dd63d-fec7-48cb-b49d-536fa1eb0a14
# ╠═118433ed-6b47-4906-b8fe-1b4602aa50be
# ╠═59acf8dd-0157-44c2-a7e8-b4c45abe00b2
# ╠═d9030aa1-1ce1-4ed4-8863-7c681784ca2e
# ╟─4511769f-07af-4582-8f3d-f43ee8de48b4
# ╟─a4ef8928-7c3f-447f-b2d0-540d7e2d0fab
# ╟─7acaae9a-d31a-4e7f-9084-f97b71e9275a
# ╟─2ff1c49a-7df0-499b-8f41-9c519dd31fb1
# ╠═7ec1add5-7c6e-4bf6-8125-6efa5fe7c48e
# ╠═77cdf367-1fe5-4458-a584-a21e26c3de44
# ╠═300e3a7b-dc1a-43b3-a50a-c955cebdb09b
# ╠═5de49946-026c-48d0-a7d0-70c293032524
# ╠═3b0c49e1-b162-471f-ae2e-355c39356812
# ╠═67d8627b-e486-4d99-805a-91c871e45a73
# ╠═3518dc88-328b-48a9-a4dd-a5d1ed919382
# ╠═7cea2363-611b-47b6-b0d0-a68a84127d81
# ╟─b4b62a45-e7c3-4f8b-9962-ee7584469d83
# ╟─a3e7c703-f71d-4f23-bc16-6c399c210e0d
# ╟─50e5afe8-e276-4e56-9034-f566869b5021
# ╟─1eb04f5a-593e-466c-bc32-ea3cffe260d4
# ╠═9777cf9f-a34f-4db4-9a85-c04d7f0ec34a
# ╟─c2450d27-f296-4bbb-b958-34a3778b3c2f
# ╟─27cb6246-d6f7-4c21-9844-31d259a5568c
# ╟─1f91f759-a6b3-4095-9eb4-93944a28741f
# ╠═eac3a8f2-7f1b-4e7c-b85e-a90151309602
# ╠═2d6c800c-10e6-47fb-a5da-600d4c294283
# ╟─6f9baa30-8314-4f6c-a21f-1cb24ae924da
# ╠═194a806d-74c9-45a7-9078-f73400df421b
# ╠═b8ba8056-a03c-4e96-ac51-e206ead16257
# ╠═6217bc62-64a3-404b-8b2a-5ebab33e6e4d
# ╠═b8c0bc97-bf63-4a31-b500-a39ca8e77b3d
# ╠═c495ba88-fc99-4a78-ab85-c324cbca495b
# ╠═02f19951-5077-4744-9d18-4b99d6ec33e3
# ╟─cea55170-b870-42dc-93c5-ce796039cc10
# ╠═dece0067-9ddd-47c5-85f8-a016e8db620f
# ╠═8bb9db58-d77d-4ed3-94a4-78894a0baa13
# ╠═11328d2a-cec1-4c81-8558-686950717939
# ╠═3af82fb8-a8b1-4f43-95ba-9260bc7beb5c
# ╠═32c5fe98-c82a-44c8-ab88-2dce146fa28a
# ╠═31e877f0-b96d-467d-8c47-9733d69f2346
# ╟─08b9e104-57fa-4807-843e-5ef958079375
# ╟─313a9618-f56d-49e3-aa07-4741421beb65
# ╟─8ff22fa9-5891-4705-8f17-d6d95a36e4c2
# ╠═c8ed7211-5b67-44ce-bf8f-f6b65806fc57
# ╟─01dbc828-4fa7-48d4-8a0a-edb275e8f715
# ╠═ce8494e6-43ab-40a5-82dd-ea64710184f9
# ╠═02fb9d18-4686-4c73-a2eb-eaf8b0f5c012
# ╟─db4672ce-933a-4aa8-9220-1e263ea540dd
# ╟─f2bd880f-e79d-4831-9f8d-65b5f4e1ba68
# ╟─6eadda32-c20e-4429-ab4f-6480facd988e
# ╠═83b3c3c3-4c91-48ed-909b-aafdce1427f9
# ╟─1f674c41-4c28-45f7-a9df-2fbb784ee115
# ╠═baecaf71-d7e0-4421-a182-daa7b6516fb7
# ╠═54fcf8a4-b397-4d38-a953-fdf575c44f71
# ╠═1a6b0fb7-6b00-45a4-b034-9d7767b6b386
# ╠═fff6da7e-6eb3-4e87-a48a-2309d001dc57
# ╠═1ef192cb-4349-4ce3-8bd8-0590b7b8923f
# ╟─be7a5329-1147-4889-a2d7-d07005369e4b
# ╠═4c130872-3a57-40a3-8c43-51e46b963b69
# ╠═a52ce3c2-3a82-4124-a260-0e1497038cd2
# ╠═22fcd352-feaa-450d-9aa4-0f8863a81b50
# ╠═f869c063-e18f-4936-92e7-7cadc22768de
# ╠═3e54e379-9513-4741-9ade-657741be2615
# ╠═913d3972-e9f4-486d-8a82-973de4271d8c
# ╟─728cf85d-7200-41b7-850b-cc68f17f68cb
# ╠═94919901-f29b-4b30-ad9a-5f5b0507aad3
# ╠═fdcb9783-0652-4651-8900-84bc81a686a7
# ╠═5ad5b33e-4f1a-40ad-851e-0fd170651006
# ╠═8d7d7e3f-9e23-4c01-9b8b-41b1e10f2bf8
# ╠═a6d47dee-b7a0-4d84-af95-200ba3eeab28
# ╟─f0c86451-a6be-47c1-95a4-967e5a90b4ce
