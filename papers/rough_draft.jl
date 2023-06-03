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


This study adapted two techniques to create two novel calcium scoring approaches that include all the calcium information within an image. The integrated Hounsfield technique was adjusted for calcium mass quantification and the volume fraction technique. Integrated calcium mass and volume fraction calcium mass were calculated and compared to known calcium mass, along with Agatston scoring and spatially weighted calcium scoring on simulated CT scans across different patient sizes, energies, and motion levels. 

The simulation was created, based on previously validated software, to match the scanning parameters of a 320-slice CT scanner. To simulate different patient sizes, additional fat rings were added, which resulted in a small-sized phantom of 30x20cm, a medium-sized phantom of 35x25cm, and a large-sized phantom of 40x30cm. Three calcification inserts of different diameters (1, 3, and 5 mm) and different hydroxyapatite (HA) densities were placed within the phantom. All calcium scoring measurements were repeated across different kVs (80-135 kV), patient sizes (small, medium, large), calcium insert sizes (1-5 mm), and calcium insert densities (25-800 mg/cm3). Physical phantom scans, acquired by a previous group, were then used to validate the results of the simulation study.

Both integrated and volume fraction calcium mass yielded lower root mean squared error (RMSE) and deviation (RMSD) values than Agatston scoring in all accuracy comparisons. Integrated and volume fraction calcium mass also yielded lower RMSE and RMSD values than Agatston or spatially weighted calcium scoring in every reproducibility comparison. Specifically, the RMSE and RMSD of the low-density calcifications for the integrated calcium mass technique vs. known calcium mass were 0.495 mg and 0.494 mg, respectively. The RMSE and RMSD values for the volume fraction calcium mass technique were 0.585 mg and 0.575 mg, respectively. Compared to Agatston scoring, which produced RMSE and RMSD values of 3.509 mg and 2.240 mg, respectively. Reproducibility across all patient sizes, calcium densities, and tube voltages, integrated calcium mass produced RMSE and RMSD values of 0.708 mg and 0.692 mg, respectively. Volume fraction calcium mass produced RMSE and RMSD values of 0.783 mg and 0.776 mg, respectively. Spatially weighted calcium scoring produced RMSE and RMSD values of 4.722 mg and 4.714 mg, respectively. Agatston scoring produced RMSE and RMSD values of 0.982 mg and 0.980 mg, respectively. The percentage of false-negative (CAC=0) scores produced by integrated calcium mass, volume fraction calcium mass, spatially weighted calcium scoring, and Agatston scoring were 11.111%, 11.111%, 42.593%, and 38.889%, respectively. The percentage of false-positive (CAC>0) scores produced by integrated calcium mass, volume fraction calcium mass, spatially weighted calcium scoring, and Agatston scoring were 11.111%, 11.111%, 16.667, and 0.0%, respectively. 
"""

# ╔═╡ 5ce128bc-7ea7-401a-af5f-50430fdf3622
md"""
# 1. Introduction
"""

# ╔═╡ ae603790-83b6-44e0-9e26-de9c87c791f5
md"""
The coronary artery calcification (CAC) score is a standard atherosclerotic marker [CITE]. CAC scoring is a test performed using computed tomography (CT) that measures the amount of calcium buildup within the walls of the heart's arteries and is an essential predictor of coronary heart disease [CITE]. The leading cause of death in the United States is coronary heart disease, killing 659,000 people annually, and it is the third leading cause of mortality worldwide [CITE].

Agatston scoring is the most common CAC scoring technique [CITE] and is a good predictor of major adverse cardiac events (MACE) [CITE]. Although, studies have shown that a significant number of patients have been characterized as having no calcium (CAC=0) while still developing MACE [CITE]. This is possibly due in part to the intensity thresholding requirements associated with the Agatston scoring technique, so other approaches like spatially weighted calcium scoring have been put forth as an alternative to Agatston scoring. Spatially weighted calcium scoring improves upon traditional calcium scoring by avoiding thresholding and has been shown to predict MACE more accurately [CITE], [CITE]. Spatially weighted calcium scoring is still limited in distinguishing non-zero CAC from noise. It also lacks quantitative insight as it is an arbitrary score without any direct physical association.

The Agatston technique can be used to estimate calcium volume and calcium mass. The calcium mass score acquired via the Agaston technique is fundamentally similar to traditional Agatston calcium scoring and suffers from many of the same limitations inherent to the Agatston scoring approach [CITE], [CITE]. A previous study also showed that calcium mass quantification via the Agatston approach produces up to 50% underestimation of calcium mass for large patients [CITE]. 

A CAC measurement technique that improves the accuracy, reproducibility, sensitivity, and specificity of the measured calcium mass would help to improve coronary heart disease diagnosis and stratify risk more accurately, especially for undetectable and nearly undetectable calcium levels. A more quantitative calcium mass calculation would likewise help improve risk stratification and diagnosis of coronary heart disease. The integrated intensity or Hounsfield unit (HU) technique recently demonstrated the ability to accurately assess coronary artery cross-sectional area, even past the visible threshold [CITE]. Similarly, volume fraction mass quantification has demonstrated the ability to accurately measure the mass of lung tissue in pulmonary computed tomography (CT) scans. These techniques have not been adapted for calcium mass quantification and CAC scoring. 


This study adapted the integrated HU and volume fraction mass quantification techniques to quantify coronary artery calcium in CT scans by accounting for the partial volume effect. These new calcium mass quantification techniques, integrated calcium mass, and volume fraction calcium mass were evaluated on simulated CT scans where the total calcium mass was known. These techniques were compared to the known calcium mass to determine which technique was the most robust and accurate compared to the gold standard Agatston scoring technique and the recently introduced spatially weighted calcium scoring technique. 

Throughout this study, Agatston scoring refers to the calcium mass score derived from the Agatston approach. Spatially weighted calcium scoring refers only to the calcium score without associated physical units. Integrated calcium mass and volume fraction calcium mass refer only to the calcium mass calculated via the integrated HU approach and volume fraction mass quantification technique, respectively. This study's calcium scoring algorithms are publicly available at https://github.com/Dale-Black/CalciumScoring.jl [13].

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
The simulation study was set to match the scanning parameters of the 320-slice CT scanner (Canon Aquilion One, Canon America Medical Systems, Tustin, CA), as previously reported 12. The X-ray spectrum was created with an interpolating polynomial model 13. The linear attenuation coefficients were made from their chemical composition 14. Poisson noise was added to simulate quantum noise. The simulation did not include Compton scatter, the dominant attenuation mechanism in CT imaging due to the interaction of free electrons with the incoming X-ray, but beam hardening was included. A3200x2200 pixel digital phantom was designed based on an anthropomorphic thorax phantom with a size of 30x20 cm2 (QRM-Thorax, QRM, Mӧhrendorf, Germany). To simulate different patient sizes, additional fat rings emulated by a mixture of 20% water and 80% lipid were added, which resulted in a medium-sized phantom of 35x25 cm2 and a large-sized phantom of 40x30 cm2. There were nine calcification inserts within the thorax with different densities and sizes. Three calcification inserts of different diameters (1, 3, and 5 mm), each with a length of 1.5 mm, and different hydroxyapatite (HA) densities were placed within each phantom. A combination of HA and myocardium was used to vary the calcification densities. For the normal-density study, the HA densities in the inserts were 200, 400, and 800 mgHAcm-3. For the low-density study, the densities were changed to 25, 50, and 100 mgHAcm-3. Each phantom also contained a 10 mm diameter calibration rod. All phantom sizes and density levels were scanned using 80, 100, 120, and 135 kV tube voltages. For small, medium, and large patient sizes, the exposure value was adjusted to 0.9, 2.0, and 5.4 mR, respectively, resulting in similar noise levels for different-sized phantoms.

Simulation materials and geometries are shown in Figure 1. Acquisition and reconstruction parameters for the simulated and physical phantoms are shown in Table 1. The calibration rods were all 10 mm in diameter. All calcium scoring measurements were repeated across each kV, patient size, calcium insert size, and calcium insert density.

Segmenting regions of interest (ROIs) is important in calcium measurement. For this study, segmentations were done automatically based on previous work by Praagh et al. 15 and adapted for simulated phantoms. The automatic segmentation approach effectively segments calcium inserts based on the known geometry of the simulated phantom without requiring manual intervention. 
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
	"Tube Voltage (kV)" => [80, 100, 120, 135],
	"Threshold (HU)" => [177, 145, 130, 112]
)

# ╔═╡ 3824a49a-50ac-4653-a17b-85210c34d3d1
md"""
Table 2. Tube voltage adapted thresholds for Agatston scoring.
"""

# ╔═╡ 51d4b9dc-ff3a-426f-8ebc-83fb75a20308
md"""
## 2.3 - Material Decomposition
"""

# ╔═╡ 9fca9de4-04cb-4c86-9b3a-58429fc1922b
md"""
## 2.4 - Volume Fraction Mass
"""

# ╔═╡ 34e7d7a9-6205-4e57-ad62-5410e0649c1f
md"""
The volume fraction calcium mass technique is similar to integrated calcium mass, but instead of calculating the calcium within an entire ROI, the percent of calcium contained within one voxel is calculated. The percent calcium contained within each voxel is then summed up within a ROI to obtain the total percentage of calcium for that ROI (Eq. 5). Given the known size of the ROI and density of the calibration rod, the volume and mass of calcium can be calculated (Eq. 6). 

Equation 5 shows how to calculate the percentage of calcium contained within each voxel. ``i`` is the intensity of one voxel of interest (HU), ``S_{Obj}`` is the intensity of pure background, which can be obtained from a ring-like object as seen in Fig. 3B, and ``S_{Obj}`` is the intensity of pure calcium which can be obtained from a calibration rod of any density. The result, ``k_{i}``, is then the percentage of calcium within one voxel ``i``. The entire region of voxels is then summed to give ``K`` which is the total percentage of calcium contained within the ROI.

Equation 6 shows how to calculate the volume of calcium ``V_{Obj}``, given the total percentage of calcium within an ROI ``K`` and the known volume of that ROI ``V_{ROI}``. This can then be converted into a calcium mass ``M_{Obj}`` given the density of the calibration rod. Fig. 3 shows the parameters needed to calculate mass via the integrated calcium mass technique (Fig. 3A) and the volume fraction calcium mass technique (Fig. 3B).
"""

# ╔═╡ df7ea29b-d887-4f78-8c83-9011a6b8bf29
md"""
```math
\begin{aligned}
k_{i} &= \frac{i - S_{Bkg}}{S_{Obj} - S_{Bkg}} \\
K &= \sum_{i} k_{i}
\end{aligned}
\tag{5}
```

\

```math
\begin{aligned}
V_{Obj} &= K \times V_{ROI} \\
M_{Obj} &= V_{Obj} \times \rho_{S_{Obj}}
\end{aligned}
\tag{6}
```
"""

# ╔═╡ 644a5004-9822-4598-8ad2-dc6a5e24f739
load(plotsdir("volume fraction.png"))

# ╔═╡ 04996a71-4c48-4f50-a95a-75e5fabff034
md"""
Fig 3. Shows two identical simulated vessel lumens with ROIs for calcium measurement.(A) Shows the ROIs needed for the integrated calcium mass technique. The central ROI that yields ``S_{Obj}`` and the background ROI that yields ``S_{Bkg}`` are unaffected by the partial volume effect, while the object ROI used to calculate ``I`` is affected by the partial volume effect. (B) Shows the ROIs needed for the volume fraction technique and a zoomed-in portion of the simulated vessel lumen, where ``i`` is the individual voxel intensity. ``S_{Bkg}`` is a ring-like background ROI unaffected by the partial volume effect. ``S_{Obj}`` is the mean intensity of known calcium density obtained from a calibration rod unaffected by the partial volume effect.
"""

# ╔═╡ dcc77efa-583e-4019-bd6b-207626ee1e75
md"""
## 2.5 - Statistical Analysis
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
Calcium mass was measured for different patient sizes, kVps, calcium sizes, and calcium densities, with and without motion on simulated phantoms, and validated in a subset of physical phantom scans acquired by Praagh et al. and shared with us for analysis. The calculated calcium mass was compared against the known mass of the calcium inserts.

The results indicate that integrated calcium mass and volume fraction calcium mass are more accurate, reproducible, sensitive, and specific than Agatston scoring (Table 4).

Agatston scoring has commonly been used in the past for predicting patient outcomes. However, a limitation of Agatston scoring is that it's only defined at 120 kVp, and a threshold of 130 HU is commonly used for calcium detection. However, the calcium attenuation coefficient is energy dependent, which makes scoring challenging when images are acquired at lower kVps, to reduce patient radiation dose. Recent reports have introduced correction factors for calcium measurements at lower kVps [CITE]. Another limitation of a thresholding approach for calcium measurement is that it is affected by partial volume effect and motion. We have introduced two new methods for calcium mass quantifications based on integrated intensity (HU) and volume fraction that address the above limitations.

Previous studies have shown that up to 5% of patients with suspected zero CAC (CAC=0) will experience MACE, despite no detectible calcium by Agatston scoring [CITE]. One question arises as to whether these patients had calcium that is not detectible by traditional Agatston scoring or simply no calcium. Integrated calcium mass, volume fraction calcium mass, and spatially weighted calcium scoring attempt to address this concern by removing the intensity thresholding requirements of standard Agatston scoring. This study shows that integrated calcium mass and volume fraction calcium mass are more sensitive to low-density calcifications than spatially weighted calcium scoring and Agatston scoring. The results showed that the percentage of false-negative (CAC=0) scores on the stationary simulated phantom were 11.111, 11.111, 42.593, and 38.889 for integrated calcium mass, volume fraction calcium mass, spatially weighted calcium scoring, and Agatston scoring, respectively (Figure 6). The substantial reduction in false-negative zero calcium scores for the integrated calcium mass and volume fraction calcium mass techniques compared to the existing techniques will help address the current limitation for patients with false-negative (CAC=0) scores.

The integrated calcium mass and volume fraction calcium mass techniques can detect calcifications that are currently indetectable by the Agatston scoring approach due to its thresholding requirement. Furthermore, a previous study has shown that calcium volume was positively and independently associated with major adverse cardiac event risk, and calcium density was inversely associated with major adverse cardiac event risk [CITE]. Another study has shown that calcium density score was the strongest positive independent predictor of major adverse cardiac events, compared to Agatston score, mass score, and volume score [CITE]. Disagreements between these studies are possibly related to the thresholding approach of Agatston scoring and poor reproducibility of Agatston scoring, which is also a limitation of all the traditional calcium (mass, volume, and density) scoring approaches based on the Agatston technique. Integrated calcium mass and volume fraction calcium mass provide more accurate, reproducible, and quantitative approaches to calcium measurement. Future studies on patient data comparing Agatston scoring with integrated calcium mass and volume fraction calcium mass might help explain these seemingly contradictory results better

When acquiring the mass of calcium in the QRM phantom using Agatston scoring, it has been shown the mass scores can be significantly different from the known physical mass [CITE]. Our results (Figures 4 and 7) show that integrated calcium mass and volume fraction calcium mass measure calcium mass more accurately than Agatston scoring, on both simulated and physical phantoms. We also see an improvement in reproducibility when comparing integrated calcium mass and volume fraction calcium mass against Agatston scoring and spatially weighted calcium scoring (Fig. 5). 

Although many parameters were considered in the simulation, only one type of reconstruction was utilized, which prevented us from addressing the effect that reconstruction type has on these calcium scoring techniques. Previous studies have shown good agreement between Agatston scores based on filtered back-projection, hybrid iterative reconstruction, and deep learning-based reconstruction [CITE]. 

Slice thickness plays an important role in calcium scoring, and traditional Agatston scoring is only defined at a slice thickness of 3 mm. Recent studies have shown that the accuracy and sensitivity of Agatston scoring are improved when slice thickness is decreased [CITE]. Our simulation was limited to 0.5 mm slice thickness which is expected to provide more accurate and sensitive comparisons for Agatston scoring. Nonetheless, future studies might provide insights by varying the slice thickness. This study was also limited by the lack of realistic cardiac hardware, such as stents, that are known to cause blooming or motion artifacts [CITE]. Future studies should account for this by including physical phantoms with realistic cardiac hardware or a robust patient data set. 

Fewer total scans were included in the physical phantom analysis compared to the simulated analysis, resulting in fewer calcium scores. This is a limitation and becomes more pronounced in the sensitivity and specificity analysis, where false-negative (CAC=0) and false-positive (CAC>0) results are rare, to begin with. In addition, motion was not included in the physical phantom analysis. Likewise, the simulated motion is based on motion in magnetic resonance imaging. More physical phantom scans, including scans affected by motion, need to be included in future studies. This would provide more robust accuracy, reproducibility, sensitivity, and specificity measurements and would result in a more reliable expected percentage of false-negative (CAC=0) and false-positive (CAC>0) scores for each calcium quantification technique.

Previous studies show that Agatston scoring consistently underestimates calcium density and volume, with even further underestimation for low-density and motion-affected plaques [CITE]. Werf et el. indicates that low-density calcifications might fall below the 130 HU threshold because of blurring from motion which artificially reduces the Agatston score [CITE]. Our study is consistent with these results (Fig. 4C and 7C); we showed that Agatston scoring produced the most false-negative (CAC=0) classifications for the simulated data and the physical phantom scans. Future studies are warranted in physical phantoms with lower-density calcification inserts (< 200 mg/cm3) to understand how integrated calcium mass and volume fraction calcium mass compares to Agatston mass within the low-density regime on physical data. Very high coronary artery calcium density (> 1000 mg/cm3) is quite rare, and Agatston scoring has already been shown to be a good predictor of cardiovascular disease within this subset of patients [CITE]. 

Tzolos et al. showed that Agatston scoring struggles to detect small calcifications in the coronary arteries of patients due to the threshold requirement of 130 HU and the minimum connected component requirement of 1 mm [33]. Our results are consistent with this study and indicate that the size of the calcium insert is a critical variable in accounting for false-negative (CAC=0) Agatston scores. The small inserts resulted in 40 false-negative (CAC=0) scores out of 216 total scores, whereas the medium and large inserts accounted for only 24 and 20 false-negative (CAC=0) scores, respectively. Density was also an important factor. Low-density calcifications resulted in 40 false-negative (CAC=0) scores out of 216 total scores, whereas the medium-density and high-density calcifications resulted in 32 and 12 false-negative (CAC=0) scores, respectively. Based on our results, integrated calcium mass and volume fraction calcium mass improved upon many of the issues associated with Agatston scoring and resulted in fewer false-negative (CAC=0) scores.
"""

# ╔═╡ 39f702fb-44c4-4af9-8086-e525f3939a49
md"""
# 5. Conclusion
"""

# ╔═╡ 4d28c1b8-642b-4e05-9868-9907e6d12718
md"""
This study shows that both the integrated calcium mass and volume fraction calcium mass techniques improve sensitivity, accuracy, and reproducibility of calcium mass measurement as compared to traditional Agatston scoring. This improvement in calcium scoring can potentially improve risk stratification for patients undergoing calcium scoring and further improve outcomes compared to Agatston scoring.
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
# ╟─9fca9de4-04cb-4c86-9b3a-58429fc1922b
# ╟─34e7d7a9-6205-4e57-ad62-5410e0649c1f
# ╟─df7ea29b-d887-4f78-8c83-9011a6b8bf29
# ╠═644a5004-9822-4598-8ad2-dc6a5e24f739
# ╟─04996a71-4c48-4f50-a95a-75e5fabff034
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
# ╠═db408b4c-027a-46eb-b0ff-89f485f112fa
# ╠═39f702fb-44c4-4af9-8086-e525f3939a49
# ╠═4d28c1b8-642b-4e05-9868-9907e6d12718
# ╠═4dabc904-0b14-426b-8691-2c77a12a8306
# ╠═20baeecb-dfc7-41c3-b5a9-c49fa63b4095
