# RadPathFusionLung
Repository contains code that allows the registration of histology slices and CT in the context of lung cancer. 

The code is provided here as suplementary material for the manuscript by Rusu et. al. "Co-registration of pre-operative CT with ex vivo surgically excised ground glass nodules to define spatial extent of invasive adenocarcinoma on in vivo imaging: a proof-of-concept study" which was published in European Radiology 2017 [10.1007/s00330-017-4813-0](https://www.doi.org/10.1007/s00330-017-4813-0)

Author: Mirabela Rusu (mirabela.rusu@gmail.com)

## Dependences: 

1. Matlab: This code was developed under version R2014a
2. Elastix: 
 * download from: [http://elastix.isi.uu.nl/](http://elastix.isi.uu.nl/)
 * This code was tested under version 4.8
 * The code assumes Elastix can be found in the path: "C:\Program Files\elastix_v4.8\\". If needed the path needs to be edited on line 5 of the file [matlab/run_all_steps.m](matlab/run_all_steps.m)
 
 
## How to run

1. Go to the folder matlab
2. Open run_all_steps.m in matlab
3. Run script 
4. When promted: 'File "...\run_all_steps.m" is not found in the current folder or the matbal path.' Choose, "Add to path"
  The scipt will run the 5 individual steps: 
    * step1_reconstrucHistology
    * step2_exhaustiveSearch
    * step3_refineHistology
    * step4_refineCTToHistology
    * step5_mapHistologyOntoCT

## Input: 
The input data has to availble in a folder structure as showed by the example:
```
1. CT
   |- CT.mha
   |- CT_crop.mha
   |- CT_crop_blood_label.mha
   | ....
   
2. histology 
   |- imgs
     |- Hist1.tif
	 |- Hist2.tif
	 |- Hist3.tif
   |- masks   
     |- lesion 
       |- Hist1.tif
	   |- Hist2.tif
	   |- Hist3.tif
     |- invasion 
       |- Hist1.tif	   
	   |- Hist2.tif
	   |- Hist3.tif
 ```
The file [matlab/run_all_steps.m] needs to be edited at line 11 to update the name of the masks folders.  
 
## Output: 

1. The code will create two folders example/output_3 and example/output_4. The former contains the results when assuming a distance between slices of 3mm, while the later conains the results of the reconstruction when assuming a distance between slices of 4mm. 
2. Each folder will contain 5 subfolder with the results of 5 individual steps
 
