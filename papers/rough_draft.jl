### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ cae3f654-00ed-11ee-3ce3-01366403ae53
using DrWatson

# ╔═╡ 1d2eaa0b-b92d-4fc3-a3cc-75dbc906893e
# ╠═╡ show_logs = false
@quickactivate "cac-dual-energy"

# ╔═╡ 063c2f86-93d6-4e68-987b-b8dc2a992489
using PlutoUI, DataFrames, Images

# ╔═╡ 3c7267f4-f8d2-4a3c-9692-74cb56927652
TableOfContents()

# ╔═╡ f7db83da-2959-4ec8-a142-33cb3a59b72c
md"""
# Abstract
"""

# ╔═╡ 459c113c-f773-462d-bbb9-d760cfaf7abc
md"""
Agatston scoring is limited by factors such as arbitrary thresholding and does not include all the calcium information in computed tomography scans. Agatston scoring’s limitations lead to problems with reproducibility and accuracy and, therefore, may not correlate well with known calcium mass. A scoring technique that removes the need for thresholding and quantifies calcium mass more accurately and reproducibly is needed.


This study adapted dual energy material decomposition to create a novel calcium scoring approach that include all the calcium information within an image. Dual energy material decomposition was used to calculate calcium mass and compared to known calcium mass, along with Agatston scoring on simulated CT scans across different patient sizes, energies, and motion levels. 

The simulation was created, based on previously validated software, to match the scanning parameters of a 320-slice CT scanner. To simulate different patient sizes, additional fat rings were added, which resulted in a small-sized phantom of 30x20cm, a medium-sized phantom of 35x25cm, and a large-sized phantom of 40x30cm. Three calcification inserts of different diameters (1, 3, and 5 mm) and different hydroxyapatite (HA) densities were placed within the phantom. All calcium scoring measurements were repeated across different kVs (80-135 kV), patient sizes (small, medium, large), calcium insert sizes (1-5 mm), and calcium insert densities []. Physical phantom scans, acquired by a previous group, were then used to validate the results of the simulation study.

[RESULTS SUMMARY HERE]
"""

# ╔═╡ 5ce128bc-7ea7-401a-af5f-50430fdf3622
md"""
# 1. Introduction
"""

# ╔═╡ ae603790-83b6-44e0-9e26-de9c87c791f5
md"""
Coronary artery calcification (CAC) is a standard atherosclerotic marker and an essential predictor of coronary heart disease, which is the leading cause of death in the United States [CITE]. 

Agatston scoring is the most common CAC scoring technique and is a good predictor of major adverse cardiac events (MACE)[CITE]. However, Agatston scoring has known limitations, such as inaccuracy [CITE], lack of reproducibility [CITE] and low sensitivity for microcalcifcations [CITE]. These limitations are likely due, in part, to the arbitrary thresholding requirement of Agatston scoring. Agastson scoring excludes all voxels below (typically) 130 Hounsfield units. Hence, any potential calcifications in these voxels are discarded. Therefore low density calcifications and microcalcifications are missed [CITE]. Similarly, different vendors and scanners can result in varying calcium scores because of the non-specific threshold [CITE]. 

Dual-energy computed tomography (DECT) is a promising technology that has been shown to improve the visualization of calcified plaque components and enhance plaque assessment [CITE]. DECT exploits the fact that different tissues have different mass attenuation coefficients when interacting with X-rays of different energies [CITE]. DECT, therefore, allows for material specifc analysis of coronary CT imaging, which can be potentially be used to directly quantify the calcium mass in DECT scans.

In this study, we investigate a new calcium mass quantification technique based on dual energy material decompositoion. We compare this method and traditional Agatston scoring against the ground truth calcium mass, in simulated calcification phantoms.
"""

# ╔═╡ 612a13c7-f6ff-4a73-9da5-a7f980c99911
md"""
# 2. Method
"""

# ╔═╡ 6fc47b71-d107-423a-a60a-474e809077a8
md"""
## 2.1 - Simulation
"""

# ╔═╡ 49c6753c-04fd-4f1b-bb2c-815e7415499c
md"""
The simulation study was set to match the scanning parameters of the 320-slice CT scanner (Canon Aquilion One, Canon America Medical Systems, Tustin, CA), as previously reported [CITE]. The X-ray spectrum was created with an interpolating polynomial model [CITE]. The linear attenuation coefficients were made from their chemical composition [CITE]. Poisson noise was added to simulate quantum noise. The simulation did not include Compton scatter, the dominant attenuation mechanism in CT imaging due to the interaction of free electrons with the incoming X-ray, but beam hardening was included. A3200x2200 pixel digital phantom was designed based on an anthropomorphic thorax phantom with a size of 30x20 cm2 (QRM-Thorax, QRM, Mӧhrendorf, Germany). To simulate different patient sizes, additional fat rings emulated by a mixture of 20% water and 80% lipid were added, which resulted in a medium-sized phantom of 35x25 cm2 and a large-sized phantom of 40x30 cm2. There were nine calcification inserts within the thorax with different densities and sizes. Three calcification inserts of different diameters (1, 3, and 5 mm), each with a length of 1.5 mm, and different hydroxyapatite (HA) densities were placed within each phantom. A combination of HA and myocardium was used to vary the calcification densities. For the normal-density study, the HA densities in the inserts were [INSERT HERE] mgHAcm-3. For the low-density study, the densities were changed to [INSERT HERE] mgHAcm-3. Each phantom also contained a 10 mm diameter calibration rod. All phantom sizes and density levels were scanned using 80/135 (dual energy) and 120 kV tube voltages. Dual energy material decomposition analysis was performed on the 80/135 kV images, whereas Agatston scoring was performed on the 120 kV images. For small, medium, and large patient sizes, the exposure value was adjusted to 0.9, 2.0, and 5.4 mR, respectively, resulting in similar noise levels for different-sized phantoms.

Simulation materials and geometries are shown in Figure 1. Acquisition and reconstruction parameters for the simulated and physical phantoms are shown in Table 1. The calibration rods were all 10 mm in diameter. All calcium scoring measurements were repeated across each kV, patient size, calcium insert size, and calcium insert density.

Segmenting regions of interest (ROIs) is important in calcium measurement. For this study, segmentations were done automatically based on previous work by Praagh et al. [CITE] and adapted for simulated phantoms. The automatic segmentation approach effectively segments calcium inserts based on the known geometry of the simulated phantom without requiring manual intervention. 
"""

# ╔═╡ 338d1ce2-bdf9-4f03-acb3-0d89e9760575
load(plotsdir("phantom materials.png"))

# ╔═╡ 2e6a0d86-d6ee-4371-a70d-f9d68c39201a
md"""
Fig. 1 Shows a sketch of the simulated phantom with the colors highlighting the different materials in the simulated phantoms.
"""

# ╔═╡ e81bfb10-c574-4811-8db6-fbd3cb557532
table1 = DataFrame(
	"Parameter" => [
		"Manufacturer"
		"CT System"
		"Reconstruction"
		"Tube Voltage (kV)"
		"Exposure Small (mR)"
		"Exposure Medium (mR)"
		"Exposure Large (mR)"
		"Slice Thickness (mm)"
		"Matrix Size Small (pixels)"
		"Matrix Size Medium (pixels)"
		"Matrix Size Large (pixels)"
		"Detector Element Width (mm)"
		"Detector Thickness (mm)"
	],
	"Simulation" => [
		"Canon"
		"Aquilion One Vision"
		"FBP"
		"80, 100, 120, 135"
		"0.9"
		"2"
		"5.4"
		"0.5"
		"640 x 440"
		"740 x 540"
		"840 x 640"
		"0.5"
		"0.5"
	]
)

# ╔═╡ 1efb1bca-b94a-4c04-980f-321a8adb1b21
md"""
Table 1. Acquisition and reconstruction parameters for the simulated phantoms.
"""

# ╔═╡ 5779b09b-28a1-4bd7-a887-db0407219a25
md"""
## 2.2 - Agatston Scoring
"""

# ╔═╡ 4c3bc9ee-0f51-4cb7-b211-e08910c784bd
md"""
Agatston scoring is defined at a tube voltage of 120 kV [CITE], but recent papers have shown how Agatston scoring can be adjusted for use in low-dose scans (70, 80, and 100 kV) [CITE]. For this study, we assumed an exponentially decreasing trendline and extrapolated beyond to account for a higher tube voltage of 135 kV. This results in an extrapolation formula shown in Equation 1, where ``y_{thresh}`` corresponds to the extrapolated threshold (HU) and ``TV`` corresponds to the tube voltage (kV). All kV-specific thresholds used in this study are shown in Table 2.
"""

# ╔═╡ dbb02eca-9938-44ac-8333-8b1abf53591c
md"""

```math
\begin{equation}
y_{thresh}=(378) e^{-0.009(TV)}
\end{equation}
\tag{1}
```

"""

# ╔═╡ a3ddda21-a9bd-4e8d-8b07-634d0ca6158d
table2 = DataFrame(
	"Tube Voltage (kV)" => [80, 120, 135],
	"Threshold (HU)" => [177, 130, 112]
)

# ╔═╡ 3824a49a-50ac-4653-a17b-85210c34d3d1
md"""
Table 2. Tube voltage adapted thresholds for Agatston scoring.
"""

# ╔═╡ 51d4b9dc-ff3a-426f-8ebc-83fb75a20308
md"""
## 2.3 - Material Decomposition
"""

# ╔═╡ dcc77efa-583e-4019-bd6b-207626ee1e75
md"""
## 2.4 - Statistical Analysis
"""

# ╔═╡ 2aa250f9-85c4-415c-aae5-048f055fe53e
md"""
All calcium scoring calculations and analyses were performed using the programming language Julia [CITE]. Root-mean-square error (RMSE) and root-mean-square deviation (RMSD) were calculated for all linear regression measurements to test for accuracy (RMSE) and precision (RMSD). Equation 7 shows how to calculate RMSE and RMSD. ``N`` is the total number of data points, ``\hat{y}`` is the calculated calcium masses, and ``y`` is either the ground truth calcium masses (RMSE) or the linear regression-based calcium masses (RMSD), which is computed based on the calculated calcium masses.

"""

# ╔═╡ d8789bd9-968f-4b4b-8b66-f7031a8ea537
md"""
```math
\begin{equation}
RMS = \sqrt{\frac{\sum{|y-{\hat{y}|}^2}}{N}}
\end{equation}
\tag{7}
```
"""

# ╔═╡ a6f51821-4fb0-4476-8973-543308342e03
md"""
# 3. Results
"""

# ╔═╡ 6ab70428-8020-4103-9496-5cdac8432070
md"""
## 3.1 - Accuracy
"""

# ╔═╡ 53ef1987-367e-497b-892f-e43f1221558c
md"""
### 3.1.1 - Low-Density
"""

# ╔═╡ 9a55d632-800f-4235-8587-2300b99aca0a
md"""
### 3.1.2 - Normal-Density
"""

# ╔═╡ 540ce011-98c3-4a85-97a7-9ef52272d7fb
md"""
## 3.2 - Reproducibility
"""

# ╔═╡ adb94eff-1c57-4991-8794-54d5d65b67f1
md"""
### 3.2.1 - Low-Density
"""

# ╔═╡ cd9e6eb8-84c8-4fbb-9e57-c08cc0429082
md"""
### 3.2.2 - Normal-Density
"""

# ╔═╡ 3389a24e-fdf7-402b-882b-076fb350c04d
md"""
## 3.3 - Sensitivity and Specificity
"""

# ╔═╡ a2595481-79da-478b-a1ab-dfcc1c9d4c8d
md"""
### 3.3.1 - Low-Density
"""

# ╔═╡ 8ba9f9ef-ce24-45b3-9bca-b877a56eba8a
md"""
### 3.3.2 - Normal-Density
"""

# ╔═╡ b8aaea3d-fe7c-41d1-96ad-f53a8c97f2da
md"""
# 4. Discussion
"""

# ╔═╡ db408b4c-027a-46eb-b0ff-89f485f112fa
md"""
Calcium mass was measured for different patient sizes, calcium sizes, and calcium densities on simulated phantoms. The calculated calcium mass was compared against the known mass of the calcium inserts.

The results indicate that material decomposition calcium mass is more accurate, reproducible, sensitive, and specific than Agatston scoring.

Agatston scoring has commonly been used in the past for predicting patient outcomes. However, a limitation of Agatston scoring is that it's only defined at 120 kVp, and a threshold of 130 HU is commonly used for calcium detection. However, the calcium attenuation coefficient is energy dependent, which makes scoring challenging when images are acquired at lower kVps, to reduce patient radiation dose. Recent reports have introduced correction factors for calcium measurements at lower kVps [CITE]. Another limitation of a thresholding approach for calcium measurement is that it is affected by partial volume effect and motion. We have introduced a new method for calcium mass quantifications based on dual energy material decomposition that addresses the above limitations.

Previous studies have shown that up to 5% of patients with suspected zero CAC (CAC=0) will experience MACE, despite no detectible calcium by Agatston scoring [CITE]. One question arises as to whether these patients had calcium that is not detectible by traditional Agatston scoring or simply no calcium. Dual energy material decomposition attempts to address this concern by removing the intensity thresholding requirements of standard Agatston scoring. This study shows that dual energy material decomposition is more sensitive to low-density calcifications than Agatston scoring. The results showed that the percentage of false-negative (CAC=0) scores on the stationary simulated phantom were [] and [] for dual energy material decomposition and Agatston scoring, respectively. The substantial reduction in false-negative zero calcium scores for dual energy material decomposition compared to the existing techniques will help address the current limitation for patients with false-negative (CAC=0) scores.

Dual energy material decomposition can detect calcifications that are currently indetectable by the Agatston scoring approach due to its thresholding requirement. Furthermore, a previous study has shown that calcium volume was positively and independently associated with major adverse cardiac event risk, and calcium density was inversely associated with major adverse cardiac event risk [CITE]. Another study has shown that calcium density score was the strongest positive independent predictor of major adverse cardiac events, compared to Agatston score, mass score, and volume score [CITE]. Disagreements between these studies are possibly related to the thresholding approach of Agatston scoring and poor reproducibility of Agatston scoring, which is also a limitation of all the traditional calcium (mass, volume, and density) scoring approaches based on the Agatston technique. Dual energy material decomposition provides a more accurate, reproducible, and quantitative approach to calcium measurement. Future studies on patient data comparing Agatston scoring with dual energy material decomposition might help explain these seemingly contradictory results better

When acquiring the mass of calcium in the QRM phantom using Agatston scoring, it has been shown the mass scores can be significantly different from the known physical mass [CITE]. Our results show that dual energy material decomposition measures calcium mass more accurately than Agatston scoring. We also see an improvement in reproducibility when comparing dual energy material decomposition against Agatston scoring. 

Although many parameters were considered in the simulation, only one type of reconstruction was utilized, which prevented us from addressing the effect that reconstruction type has on these calcium scoring techniques. Previous studies have shown good agreement between Agatston scores based on filtered back-projection, hybrid iterative reconstruction, and deep learning-based reconstruction [CITE]. 

Slice thickness plays an important role in calcium scoring, and traditional Agatston scoring is only defined at a slice thickness of 3 mm. Recent studies have shown that the accuracy and sensitivity of Agatston scoring are improved when slice thickness is decreased [CITE]. Our simulation was limited to 0.5 mm slice thickness which is expected to provide more accurate and sensitive comparisons for Agatston scoring. Nonetheless, future studies might provide insights by varying the slice thickness. This study was also limited by the lack of realistic cardiac hardware, such as stents, that are known to cause blooming or motion artifacts [CITE]. Future studies should account for this by including physical phantoms with realistic cardiac hardware or a robust patient data set. 

Previous studies show that Agatston scoring consistently underestimates calcium density and volume, with even further underestimation for low-density and motion-affected plaques [CITE]. Werf et el. indicates that low-density calcifications might fall below the 130 HU threshold because of blurring from motion which artificially reduces the Agatston score [CITE]. Our study is consistent with these results; we showed that Agatston scoring produced the most false-negative (CAC=0) classifications for the simulated data and the physical phantom scans. Future studies are warranted in physical phantoms with lower-density calcification inserts (< 200 mg/cm3) to understand how dual energy material decomposition compares to Agatston mass within the low-density regime on physical data. Very high coronary artery calcium density (> 1000 mg/cm3) is quite rare, and Agatston scoring has already been shown to be a good predictor of cardiovascular disease within this subset of patients [CITE]. 

Tzolos et al. showed that Agatston scoring struggles to detect small calcifications in the coronary arteries of patients due to the threshold requirement of 130 HU and the minimum connected component requirement of 1 mm [33]. Our results are consistent with this study and indicate that the size of the calcium insert is a critical variable in accounting for false-negative (CAC=0) Agatston scores. The small inserts resulted in [] false-negative (CAC=0) scores out of [] total scores, whereas the medium and large inserts accounted for only [] and [] false-negative (CAC=0) scores, respectively. Density was also an important factor. Low-density calcifications resulted in [] false-negative (CAC=0) scores out of [] total scores, whereas the medium-density and high-density calcifications resulted in [] and [] false-negative (CAC=0) scores, respectively. Based on our results, integrated calcium mass and volume fraction calcium mass improved upon many of the issues associated with Agatston scoring and resulted in fewer false-negative (CAC=0) scores.
"""

# ╔═╡ 39f702fb-44c4-4af9-8086-e525f3939a49
md"""
# 5. Conclusion
"""

# ╔═╡ 4d28c1b8-642b-4e05-9868-9907e6d12718
md"""
This study shows dual energy material decomposition improves sensitivity, accuracy, and reproducibility of calcium mass measurement as compared to traditional Agatston scoring. This improvement in calcium scoring can potentially improve risk stratification for patients undergoing calcium scoring and further improve outcomes compared to Agatston scoring.
"""

# ╔═╡ 4dabc904-0b14-426b-8691-2c77a12a8306
md"""
# Acknowledgements
"""

# ╔═╡ 20baeecb-dfc7-41c3-b5a9-c49fa63b4095
md"""
# References
"""

# ╔═╡ Cell order:
# ╠═cae3f654-00ed-11ee-3ce3-01366403ae53
# ╠═1d2eaa0b-b92d-4fc3-a3cc-75dbc906893e
# ╠═063c2f86-93d6-4e68-987b-b8dc2a992489
# ╠═3c7267f4-f8d2-4a3c-9692-74cb56927652
# ╟─f7db83da-2959-4ec8-a142-33cb3a59b72c
# ╟─459c113c-f773-462d-bbb9-d760cfaf7abc
# ╟─5ce128bc-7ea7-401a-af5f-50430fdf3622
# ╟─ae603790-83b6-44e0-9e26-de9c87c791f5
# ╟─612a13c7-f6ff-4a73-9da5-a7f980c99911
# ╟─6fc47b71-d107-423a-a60a-474e809077a8
# ╟─49c6753c-04fd-4f1b-bb2c-815e7415499c
# ╟─338d1ce2-bdf9-4f03-acb3-0d89e9760575
# ╟─2e6a0d86-d6ee-4371-a70d-f9d68c39201a
# ╟─e81bfb10-c574-4811-8db6-fbd3cb557532
# ╟─1efb1bca-b94a-4c04-980f-321a8adb1b21
# ╟─5779b09b-28a1-4bd7-a887-db0407219a25
# ╟─4c3bc9ee-0f51-4cb7-b211-e08910c784bd
# ╟─dbb02eca-9938-44ac-8333-8b1abf53591c
# ╟─a3ddda21-a9bd-4e8d-8b07-634d0ca6158d
# ╟─3824a49a-50ac-4653-a17b-85210c34d3d1
# ╟─51d4b9dc-ff3a-426f-8ebc-83fb75a20308
# ╟─dcc77efa-583e-4019-bd6b-207626ee1e75
# ╟─2aa250f9-85c4-415c-aae5-048f055fe53e
# ╟─d8789bd9-968f-4b4b-8b66-f7031a8ea537
# ╟─a6f51821-4fb0-4476-8973-543308342e03
# ╟─6ab70428-8020-4103-9496-5cdac8432070
# ╟─53ef1987-367e-497b-892f-e43f1221558c
# ╟─9a55d632-800f-4235-8587-2300b99aca0a
# ╟─540ce011-98c3-4a85-97a7-9ef52272d7fb
# ╟─adb94eff-1c57-4991-8794-54d5d65b67f1
# ╟─cd9e6eb8-84c8-4fbb-9e57-c08cc0429082
# ╟─3389a24e-fdf7-402b-882b-076fb350c04d
# ╟─a2595481-79da-478b-a1ab-dfcc1c9d4c8d
# ╟─8ba9f9ef-ce24-45b3-9bca-b877a56eba8a
# ╟─b8aaea3d-fe7c-41d1-96ad-f53a8c97f2da
# ╟─db408b4c-027a-46eb-b0ff-89f485f112fa
# ╟─39f702fb-44c4-4af9-8086-e525f3939a49
# ╟─4d28c1b8-642b-4e05-9868-9907e6d12718
# ╟─4dabc904-0b14-426b-8691-2c77a12a8306
# ╟─20baeecb-dfc7-41c3-b5a9-c49fa63b4095
