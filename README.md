# **dluskora**
*for metamorphic petrologists*

**dluskola** is a viable tool for simulating trace-elements zoning patterns in a
prograde garnet porphyroblast in low-*T* eclogite.
This is based on the theory given by **Skora et al. (2006)**.
Using this tool, we can discuss the kinetics of REE diffusion and garnet growth in subduction zones.

## Features
This package encompasses two functions: **dluskora1** and **dluskora2**. By choosing more appropriate one and setting some parameters, you can easily draw theoretical core-to-rim REE profiles.
**dluskora1** is suitable for a situation under which a garnet radius can be approximated by a linear function of time,
while **dluskora2** (beta version) is useful when a radius–time trajectory is not expressed as a straight line.
Moreover, these functions are also compatible with the idea of parameter reduction by **Fukushima et al. (*Island Arc*, in prep)**.
## Requirement
R 3.6.0
## Installation
Please type the following command in your R console.

`remotes::install_github("Ryo-fkushima/dluskora")`
## Usage
`dluskora1(fac1=(1.78 *10^(21)), Q=300000, syR=0.60, c_ave=30,
D_0=(4 * 10^(13)), R=8.3, T_1=450, T_2=600, K_d=15, Mr=100,
Mt=45, fg=1, ft=2, garsize=0.27)`

Output is generated with the *plot* function. If you have REE profiles in natural garnet samples, you can compare them to the calculated curve by using the *lines* or *points* functions.

## Author
Ryo Fukushima (Tohoku University, Sendai, Japan)

## References
Fukushima R, Tsujimori T, Aoki S, Aoki K (in prep) A new scheme of diffusion-limited REE-uptake model for prograde-zoned garnets in low-temperature eclogite: Principle and application. Isl Arc

Skora S, Baumgartner LP, Mahlen NJ, Johnson CM, Pilet S, Hellebrand E (2006) Diffusion-limited REE uptake by eclogite garnets
and its consequences for Lu–Hf and Sm–Nd geochronology. Contrib Mineral Petrol 152:703-720



