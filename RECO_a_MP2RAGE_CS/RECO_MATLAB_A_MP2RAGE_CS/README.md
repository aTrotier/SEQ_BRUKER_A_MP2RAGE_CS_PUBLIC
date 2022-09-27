
# Matlab reconstruction

## Folder structure

* **RECO_MATLAB_A_MP2RAGE_CS**
  * **script_T1MAP_MP2RAGE.m** is the main script to reconstruct the T1 maps
  * **simu folder** : regroup script to find optimal parameters for acquisition and reconstruction
    * **opti_reco** folder : script to find the optimal parameters for reconstruction
    * **opti_acq** folder : various script to analyze the acquisition parameters effect


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