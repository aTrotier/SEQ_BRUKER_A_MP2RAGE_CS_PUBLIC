# Introduction

This part of the repository contains code in order to find an optimal Compressed-Sensing reconstruction from a fully dataset with undersampled poisson disc distribution.



This code :

1. Read a bruker fully k-space 
2. Generate a poisson disk mask and also create and undersample k-space with it
3. Reconstruct an image in a subfunction with various CS parameters. Then compare the results with the fully image. This function keeps the image of the lowest RMSE and display the parameters used to obtain it.


**Warning** Depending of the range of reconstruction you want to test, this function can be very long.



# Run the code

First you need to compile the mex file (see next section : poisson disc mask)

The main script to use is **script_opti_reco**. Parameters can be modify in section :

* **Create Mask** to change the undersampling of the kspace
* **Find best reco for W in spatial direction and TV** : change the CS parameters list tested when you performe the Wavelet regularization along the spatial dimension  and the Total Variation regularization along the Time of Inversion (TI)
* **Find best reco for W in spatial direction only** : change the CS parameters list tested when you performe the Wavelet regularization along the spatial dimension



The last section : **Do your own reco** can be used to try different reconstruction one by one. You can launch it after section **Create fully and undersample dataset**



# Poisson disc mask

## Informations

Undersampled poisson mask is created with a C function : poissonDiscMex.c and poissonDiscMex.h. 

The original code is available in the espirit toolbox developped by Lustig M. and available at this adress : **https://people.eecs.berkeley.edu/~mlustig/Software.html**. 

I have done some modifications   :

- To make the code usable for Bruker (using global variable)
- To start the seed always at the same value in order to obtain reproducible mask
- Integrate a fully encoded part at the center of the kspace in order to generate the sensitivity map



This function should generate exactly the same distribution than the one included in the a_MP2RAGE_CS sequence. If some modifications is done in the function poissonDisc.h/.c in the sequence it should also been include in that file.

## How to use

In order to use the code under matlab it should be compiled with the following command  under matlab :

`mex poissondDiscMex.c`

Then it can be called under matlab with the function :

`[mask] = poissonDiscMex(fovy,fovz,sky,skz,ry,rz,ncal,cutcorners, pp);`

where 

- **mask** is a binary mask along ky-kz direction
- **fovy/fovz** are the field of view (both should be in the same units)
- **sky/skz** is the size of the k-space
- **ry/rz** is the acceleration rate along ky and kz direction
- **ncal** is the size of the fully central part
- **cutcorners** if a flag equal to 0 or 1 if you want to obtain an elliptic kspace (it cut the corners)
- **pp** is a double value corresponding to the rate of sparsity along the radial dimension (try it before using higher values in your sequence)