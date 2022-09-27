# Bruker Sequence : MP2RAGE + Compressed Sensing acceleration

This repository regroups the compiled sequence and matlab reconstruction of the MP2RAGE+CS sequence. This repository is still in **beta** don't hesitate to open a issue and share your comments.



Sequence principle is described in this MRM publication :  [![DOI](https://zenodo.org/badge/DOI/10.1002/mrm.27438.svg)](https://doi.org/10.1002/mrm.27438)

**Source code is available as a private [submodule](https://github.com/aTrotier/SEQ_BRUKER_A_MP2RAGE_CS_PUBLIC)** if you want the source code, send a request to : <aurelien.trotier@rmsb.u-bordeaux.fr>



## Plan

- [Bruker Sequence : MP2RAGE + Compressed Sensing acceleration](#bruker-sequence--mp2rage--compressed-sensing-acceleration)
  - [Plan](#plan)
  - [Folder structure](#folder-structure)
  - [Sequence installation](#sequence-installation)
  - [Acquisition](#acquisition)

## Folder structure

* **SEQ** : Binary the sequence on Bruker scanner (PV6.0.1).
* **RECO_a_MP2RAGE_CS** : Reconstruction available in matlab and Julia
  * **RECO_JULIA_A_MP2RAGE_CS** : Offline and online reconstruction through macro
  * **RECO_MATLAB_A_MP2RAGE_CS**
    * **script_T1MAP_MP2RAGE.m** is the main script to reconstruct the T1 maps
    * **simu folder** : regroup script to find optimal parameters for acquisition and reconstruction
      * **opti_reco** folder : script to find the optimal parameters for reconstruction
      * **opti_acq** folder : various script to analyze the acquisition parameters effect


## Sequence installation

Sequence has been developped for Paravision **PV6.0.1**. Minor modification are required for PV6.0 compatibility (feel free to contact us)

**Installation step :**

* Copy the binary sequence under the folder
  `/opt/PV6.0.1/share/`

* To install : `File -> Import -> Binary Method ` and select the sequence in the share folder.


Sequence is now install and available under the **Palette** tab/Explorer tab/Scan Programs & Protocols :

```
Object : AnyObject
Region : AnyRegion
Application : UserMethods
```

To use it drag and drop to an exam card.

## Acquisition

We will only described the MP2RAGE specific parameters. The other ones are standard (Sequence is based on the Bruker FLASH sequence). The parameter name in the bruker card is in **bold**. In brackets are the corresponding names in the publication.

* ROUTINE tab

  * **Repetition Time** (TR) : is the time between 2 excitations within GRE blocks

  * **Recovery Time** (MP2RAGETR) : is the time between 2 inversions pulses

  * **Flip Angle TI1** and **Flip Angle TI2** : correspond of the excitation flip angle of the first and second GRE Blocks

  * **Train Echo**  : is the number of echos per GRE readout

  * **Effective TI** (TI1/TI2) : correspond to the time between inversion and the middle of each GRE blocks.

  * **Acceleration Type** :

    * **Fully** : no acceleration, all k-space is encoded (if the **Echo Train is equal to the ky matrix size. The image will be directly reconstructed under Bruker)
    * **Grappa** (WIP -> supposed to work, just check the mask to be sure)
    * **Compressed sensing** :
      * **Mask_With_File** : *deprecated* (if you need to use your own mask contact us)
      * **Mask_With_algorithm** : if you change any parameters don't forget to click on **CalcMask**

    > Both Grappa and Compressed use **Acceleration CS Y/Z** so choose your acceleration and click on **CalcMask**. It indicates the total acceleration **Acceleration Totale** (French :D).
    >
    > **Size of center mask** indicates the size of the squared fully sampled central part of the k-space (Used later for coil sensitivity calibration)
    >
    > **GenpointMask** : the mask algorithm generates a number of ky/kz lines which can not be a multiple of the **Train Echo**. Thus, we do not acquire some lines at the bottom of the k-space, in order to remain a multiple. The number of non-acquired lines is called the **GenpointMask**

    **AGAIN : Don't forget to click on CalcMask**

* Contrast tab / Sel IR : You can change the inversion pulse parameters

  * Usually we use a sech of 10 ms with a slab Thickness of 200%

* Sequence tab / Main : You can change the excitation pulse parameters

  * Usually we use sinc10H of 1ms
