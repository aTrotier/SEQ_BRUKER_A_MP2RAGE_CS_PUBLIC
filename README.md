# Bruker Sequence : MP2RAGE + Compressed Sensing acceleration

This repository regroups the compiled sequence and matlab reconstruction of the MP2RAGE+CS sequence. This repository is still in **beta** don't hesitate to open a issue and share your comments.

Sequence principle is described in this MRM publication :  [![DOI](https://zenodo.org/badge/DOI/10.1002/mrm.27438.svg)](https://doi.org/10.1002/mrm.27438)

If you want the source code, send a request to : <aurelien.trotier@rmsb.u-bordeaux.fr>



## Plan

* [Folder structure](#Folder-structure)
* [Sequence installation](#Sequence-installation)
* [Acquisition](#Acquisition)
* [Reconstruction](#Reconstruction)
* [Simulation](#Simulation)
* [NNotes](#Sequence installation)

## Folder structure

* **a_MP2RAGE_CS** : Includes the file to install the sequence on Bruker scanner (PV6.0.1).
* **reco-a_MP2RAGE_cs** : Matlab script for reconstruction
  * **script_T1MAP_MP2RAGE.m** is the main script to reconstruct the T1 maps
  * **simu folder** : regroup script to find optimal parameters for acquisition and reconstruction
    * **opti_reco** folder : script to find the optimal parameters for reconstruction
    * **opti_acq** folder : various script to analyze the acquisition parameters effect


## Sequence installation

Sequence has been developped for Paravision **PV6.0.1**. Minor modification are required for PV6.0 compatibility (feel free to contact us)

**Installation step :**

* Copy the folder **a_MP2RAGE_CS** to the location (Change the {USER} by your user name) : 
  `/opt/PV6.0.1/prog/curdir/{USER}/Paravision/methods/src/`

* In the Paravision Workspace Explorer go to : `Method Development/User Methods ({USER})`

  * right-click on the sequence a_MP2RAGE_CS
  * Build/install
  * Let all the build options selected -> OK
  * You should see these last lines in popup window **output** :

  ```
  Link a_MP2RAGE_CS.so
  Install /opt/PV6.0.1/prog/curdir/{USER}/Paravision/methods/a_MP2RAGE_CS.so
  Install /opt/PV6.0.1/prog/curdir/{USER}/Paravision/methods/a_MP2RAGE_CS.xml
  ```



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



## Recontruction

**Requirements :** 

* Matlab (tested on version > 2016b
* [BART toolbox](https://mrirecon.github.io/bart/) (Tested for BART 0.5).

Before launching any script add bart to your matlab path (launch the startup.m file in the bart folder).

To perform the offline recontruction : 

* Download the bruker dataset located in `/opt/PV6.0.1/data/{USER}/`

* Add to matlab path the folder and subfolder  **reco-a_MP2RAGE_cs**

* Launch the script : **script_T1MAP_MP2RAGE.m**

  * A popup window ask for the bruker dataset : select the scan folder you want to reconstruct (it is a number)
  * A matlab object : **OBJ_MP2RAGE_RECO** is created (here called **param_in**) which regroup the parameter that will be used for the Compressed-sensing reconstruction.
  * You can change the reconstruction parameter, for example :

  ```matlab
  param_in.ITER = 30;
  param_in.OPTION_RECO = "NR_RECO";
  param_in.W_LAMBDA = 0.02;
  ```

  * To run the reconstruction pass the **OBJ_MP2RAGE_RECO** to the function **a_MP2RAGE_CS_bart**

  ```matlab
  s_out = a_MP2RAGE_CS_bart(param_in);
  ```

  * If all goes right, multiples figures will popup :
    * ky/kz mask
    * LookupTable
  * The reconstruction data are stored in matlab structure (here called **s_out**) which include
    * **T1map**
    * **MP2RAGE_mask** : MP2RAGE image modified with a filter that nulls the salt and-pepper background 
    * **imTI** : Ti1 and Ti2 images obtained after the CS reconstruction
    * **LUT** : structure that stores the parameters and image used to create the T1map
    * **param_in** are also store

## SIMULATION

* opti_reco : **documentation todo**
* opti_acq : 
  * **simu_lut_MP2RAGE.m** : plot the lookuptable (important to check the range of T1 value measurable)
  * **opti_CNR** : set 2 T1 tissues at the beginning of the script. It will determined what is the best parameters to choose to get the highest CNR. It takes into account the acquisition time (divide by sqrt(**MP2RAGE_TR**) and/or sqrt(**Echo train**))
