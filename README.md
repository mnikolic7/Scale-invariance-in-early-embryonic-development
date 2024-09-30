# Scale-invariance-in-early-embryonic-development
Data and analysis accompanying the paper:  "Scale invariance in early embryonic development"   Milos Nikolic, Victoria Antonetti, Feng Liu, Gentian Muhaxheri, Mariela D. Petkova, Martin Scheeler, Eric M. Smith, William Bialek and Thomas Gregor. 

[arxiv.org/abs/2312.17684](https://arxiv.org/abs/2312.17684)

------------------------------------------------------------------------

Here we provide the data we used in our analysis, and our MATLAB
analysis code. For more details, see the reference.

We analyze gene expression patterns from three data sets: 
1. the pair-rule genes [Petkova et al. 2019], 
2. the gap genes [Petkova et al. 2019],
3. maternal morphogen, Bicoid [Liu et al. 2013]. 
	note: Cephalic furrow position (CF) is also from [Liu et al. 2013]

We are only concerned with the WT expression data in 1. and 2.
For 3. we use data only from the 2XA fly line which used as the 
reference line because it has almost identical dosage of Bcd compared 
to WT. 

Sources:
- Petkova, Mariela D., Gašper Tkačik, William Bialek, Eric F. Wieschaus,
and Thomas Gregor. "Optimal decoding of cellular identities in a genetic
network." Cell 176, no. 4 (2019): 844-855.

- Liu, Feng, Alexander H. Morrison, and Thomas Gregor. "Dynamic 
interpretation of maternal inputs by the Drosophila segmentation gene 
network." Proceedings of the National Academy of Sciences 110, no. 17 
(2013): 6724-6729.

Additionally, we use some analysis protocols that are described in 
more detail in:
- Dubuis, Julien O., Reba Samanta, and Thomas Gregor. "Accurate 
measurements of dynamics and reproducibility in small genetic networks."
Molecular systems biology 9.1 (2013): 639.

------------------------------------------------------------------------
Pair rule gene data [Petkova et al. 2019]
------------------------------------------------------------------------
DATA DESCRIPTION:

Every rawProfiles_"gene".mat file contains a data structure with the
following fields: 

.index = embryo image # (for internal purposes) 
.dist = depth of the cellularization membrane, in microns. 
.age =  embryo age in minutes after the start of nc14
."gene_name"  =  raw fluorescence gene expression profiles taken along 
	the dorsal side of the embryo in scaled coordinates (x_s=x/L). 
	Each profile vector is of length 1000 pixels, 
	where 0 corresponds to the anterior (A) and 1000 to the 
	posterior (P) of the embryo.  
.L = embryo length in microns, extracted from the raw images. 

------------------------------------------------------------------------
Gap gene data [Petkova et al. 2019]
------------------------------------------------------------------------
DATA DESCRIPTION:

rawProfiles_gapGenes_Hb_Gt_Kni_Kr.mat contains the gap gene data.

.index = embryo image # (for internal purposes) 
.dist = depth of the cellularization membrane, in microns. 
.age =  embryo age in minutes after the start of nc14
."gene_name"  =  raw fluorescence gene expression profiles taken along 
	the dorsal side of the embryo in scaled coordinates (x_s=x/L). 
	Each profile vector is of length 900 pixels, 
	where first entry corresponds to x/L=0.051 and the last entry
	corresponds to the x/L=0.950.
.xs = the vector of relative positions xs=x/L from 0.05 to 0.95. 
.L = embryo length in microns, extracted from the raw images. 


The dataset gapGenes_timeCorrected_N301.mat contains the subset of the  
gap gene expression profiles (without the imaging/staining sessions that
had higher background, see SFig.3 in reference). The gap genes are also 
time corrected as in Dubuis et al. 2013. It uses the same format. 

------------------------------------------------------------------------
Bicoid data [Liu et al. 2013]
------------------------------------------------------------------------
DATA DESCRIPTION:

Fig5_rawProfiles_Bcd.mat contains three variables:
L - embryo lengths
xs - relative positions
profiles_bcd - 90 by 582 matrix, each column is a single Bcd expresssion
	profile. Profiles are normalized to the maximum of the mean of all 
	profiles, in the same way as the gap genes. 

------------------------------------------------------------------------
Cephalic furrow data [Liu et al. 2013]
------------------------------------------------------------------------
Fig1_cephalicFurrow.mat contains two variables:

CF - position of CF in relative coordinates
L_CF - length of the corresponding embryos. 

