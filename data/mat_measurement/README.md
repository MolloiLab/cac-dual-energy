# Simulated Measurement Images

## Folder Structure
`.mat` files of mass measurement phantom arranged by phantom sizes, densities and KV \
`PHANTOMSIZE/IMAGES` \
[SIZE] are first set of images (good result), [SIZE1] are second set of images with more focus on lower density (sysmetic error)

## Image Geometry
QRM phantom
* 9 inserts in one image with 3 different densities of 3 different sizes
* Diameter of inserts:
  * **Small**: 1 mm
  * **Medium**: 3 mm
  * **Large**: 5 mm

## Calcium Density(ies)

| Folder Name | Density Ref    | High Density (mg/cc) | Medium Density (mg/cc) | Low Density (mg/cc) |
| ----------- | -------------- | -------------------- | ---------------------- | ------------------- |
| Large       | Density1       | 733                  | 411                    | 151                 |
| Large       | Density2       | 669                  | 370                    | 90                  |
| Large       | Density3       | 552                  | 222                    | 52                  |
| Large1      | Density1       | 797                  | 101                    | 37                  |
| Large1      | Density2       | 403                  | 48                     | 32                  |
| Large1      | Density3       | 199                  | 41                     | 27                  |
| Medium      | Density1       | 733                  | 411                    | 151                 |
| Medium      | Density2       | 669                  | 370                    | 90                  |
| Medium      | Density3       | 552                  | 222                    | 52                  |
| Medium1     | Density1       | 797                  | 101                    | 37                  |
| Medium1     | Density2       | 403                  | 48                     | 32                  |
| Medium1     | Density3       | 199                  | 41                     | 27                  |
| Small       | Density1       | 733                  | 411                    | 151                 |
| Small       | Density2       | 669                  | 370                    | 90                  |
| Small       | Density3       | 552                  | 222                    | 52                  |
| Small1      | Density1       | 797                  | 101                    | 37                  |
| Small1      | Density2       | 403                  | 48                     | 32                  |
| Small1      | Density3       | 199                  | 41                     | 27                  |


## Simulated Computed Tomography Parameter
**Exposure** 
Exposure values are same as the one used in integrated technique. Increasing the exposure can help reduce the noise and possibly contribute to a better measurement.
| Phantom Size | Exposure (mR) |
| ------------ | ------------- |
| Large        | 5.4           |
| Medium       | 2.0           |
| Small        | 0.9           |
