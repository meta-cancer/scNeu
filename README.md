# scNeu: A pipeline for identify neutrophils from single-cell data

Neutrophils are the most abundant leukocytes in human blood and are essential components of innate immunity. Until recently, neutrophils were considered homogeneous and transcriptionally inactive cells, but both concepts are being challenged. Single-cell RNA sequencing (scRNA-seq) offers an unbiased view of cells along a continuum of transcriptional states. However, the use of scRNA-seq to characterize neutrophils has proven technically difficult, explaining in part the paucity of published single-cell data on neutrophils. We have found that modifications to the data analysis pipeline, rather than to the existing scRNA-seq chemistries, can significantly increase the detection of human neutrophils in scRNA-seq.

scNeu provides a comprehensive tool for identify neutrophils from 10X Genomics single-cell data.

## Pipeline Overview

scNeu pipeline includes three steps:

### Step1: Prefilter from 10X genomics single-cell raw matrix (in R)

```shell
singularity exec --nv \
    --bind ${software}:${software} \
    --bind ${inPth}:${inPth} \
    --bind ${outPth}:${outPth} \
    ${sif_path} \
    Rscript /scNeu/step1_raw_prefilter.R ${i} ${inPth} ${outPth}${i}/step1_prefilter ${software}
```

### Step2: Process and predict with SVM model (in Python)

In step2, we can choose whether to use GPU to accelerate the calculation of the model according to the actual situation.

```shell
singularity exec --nv \
    --bind ${software}:${software} \
    --bind ${inPth}:${inPth} \
    --bind ${outPth}:${outPth} \
    ${sif_path} \
    python /scNeu/step2_auto_predict.py -s ${i} \
        -i ${outPth}${i}/step1_prefilter \
        -o ${outPth}${i}/step2_predict \
        --model /scNeu/model.joblib \
        --scaler /scNeu/scaler.joblib \
        --geneList /scNeu/geneList.csv \
        --usegpu ${useGPU}
```

### Step3: Perform final neutrophil annotations (in R)

```shell
singularity exec --nv \
    --bind ${software}:${software} \
    --bind ${inPth}:${inPth} \
    --bind ${outPth}:${outPth} \
    ${sif_path} \
    Rscript /scNeu/step3_final_annotation.R ${i} ${outPth}${i}/step1_prefilter ${outPth}${i}/step3_annotation ${outPth}${i}/step2_predict ${outPth}${i} ${software}
```

Please cite:

Here are two usages of this pipelines: singularity and conda/mamba environment.

## Using singularity

1. First, download the ".sif" files in release pages (build with singularity 4.0.1).

2. Execute pipeline:

   ```shell
   #!/bin/bash
   
   #################### args ####################
   export sif_path="/0.db/Pan/scNeu/scNeu.sif"
   export software="/0.db/Pan/scNeu/"
   export sampleList="/0.db/Pan/scNeu/example_sampleList.txt"
   export inPth="/0.db/Pan/scNeu/example_rawdata/"
   export outPth="/0.db/Pan/scNeu/example_output/"
   
   ################### config ###################
   useGPU=True # Boolean: Choose True or False
   
   ##################### main #####################
   mkdir -p ${outPth}
   cd ${outPth}
   for i in `cat $sampleList`
   do
   	echo "!!!!!!!!!!!!!!Start:"${i}"!!!!!!!!!!!!!!!!!!!"
   	mkdir -p ${outPth}${i}
   	cd ${outPth}${i}
   	
   	#### step1：prefilter from raw matrix ####
   	echo "!!!!!!!!!!!!!! "${i}": Step1__Prefilter from raw matrix!!!!!!!!!!!!!!!!!!!"
   	mkdir -p ./step1_prefilter
   	singularity exec --nv \
           --bind ${software}:${software} \
           --bind ${inPth}:${inPth} \
           --bind ${outPth}:${outPth} \
           ${sif_path} \
           Rscript /scNeu/step1_raw_prefilter.R ${i} ${inPth} ${outPth}${i}/step1_prefilter ${software}
   	
   	### step2：process and predict ####
   	echo "!!!!!!!!!!!!!! "${i}": Step2__Process and predict !!!!!!!!!!!!!!!!!!!"
   	mkdir -p ./step2_predict
       singularity exec --nv \
           --bind ${software}:${software} \
           --bind ${inPth}:${inPth} \
           --bind ${outPth}:${outPth} \
           ${sif_path} \
           python /scNeu/step2_auto_predict.py -s ${i} \
           -i ${outPth}${i}/step1_prefilter \
           -o ${outPth}${i}/step2_predict \
   		--model /scNeu/model.joblib \
   		--scaler /scNeu/scaler.joblib \
   		--geneList /scNeu/geneList.csv \
   		--usegpu ${useGPU}
   
   	### step3：perform final annotations ####
   	echo "!!!!!!!!!!!!!! "${i}": Step3__Perform final annotations !!!!!!!!!!!!!!!!!!!"
   	mkdir -p ./step3_annotation
       singularity exec --nv \
           --bind ${software}:${software} \
           --bind ${inPth}:${inPth} \
           --bind ${outPth}:${outPth} \
           ${sif_path} \
           Rscript /scNeu/step3_final_annotation.R ${i} ${outPth}${i}/step1_prefilter ${outPth}${i}/step3_annotation ${outPth}${i}/step2_predict ${outPth}${i} ${software}
   	
   	cd ${outPth}
   	echo "!!!!!!!!!!!!!!Done:"${i}"!!!!!!!!!!!!!!!!!!!"
   done
   ```
   
   and then execute:
   
   ```
   bash ./sig_scNeu.sh
   ```
   
   

## Using conda/mamba

1. First, create the conda/mamba environment.

   ```shell
   mamba create -n scNeu -c rapidsai -c nvidia -c conda-forge -c bioconda -c defaults python=3.8 scanpy=1.9.1 cuml=22.12.00 cupy=9.6.0 rmm=22.12.00 dask-cudf=22.12.01 r-seurat=4.1.1 r-dplyr=1.0.9 r-qs=0.26.3 bioconductor-singlecellexperiment=1.20.0 numpy=1.22.4 matplotlib=3.6.3 numba=0.56.2 rpy2=3.5.11 r-base=4.2.3 r-stringi=1.7.8 lcms2=2.12 icu=70.1 r-pheatmap=1.0.12 bioconductor-biocparallel=1.32.5 ipykernel=6.28.0 anndata2ri=1.3.1 backports.zoneinfo=0.2.1 parallel importlib-resources=6.4.0
   ```

   or use conda "environment.yml":

   ```shell
   mamba env create -f environment.yml
   ```
   
   and then:
   
   ```shell
   sed -i '17s/.*/'"$(sed -n '19p' "${CONDA_PREFIX}/lib/python3.8/site-packages/setuptools/_vendor/jaraco/context.py")"'/' "${CONDA_PREFIX}/lib/python3.8/site-packages/setuptools/_vendor/jaraco/context.py"
   ```
   
   
   
2. Execute pipeline (replace "args" to actual path in line 7~10):

   ```shell
   #!/bin/bash
   # mamba activate scNeu
   export R_HOME=${CONDA_PREFIX}/lib/R
   R_LIB_PATH=${CONDA_PREFIX}/lib/R/library
   
   #################### args ####################
   export software="/0.db/Pan/scNeu/"
   export sampleList="/0.db/Pan/scNeu/example_sampleList.txt"
   export inPth="/0.db/Pan/scNeu/example_rawdata/"
   export outPth="/0.db/Pan/scNeu/example_output/"
   
   ################### config ###################
   useGPU=True # Boolean: Choose True or False
   
   ##################### main #####################
   model=${software}"model.joblib"
   scaler=${software}"scaler.joblib"
   geneList=${software}"geneList.csv"
   
   mkdir -p ${outPth}
   cd ${outPth}
   for i in `cat $sampleList`
   do
   	echo "!!!!!!!!!!!!!!Start:"${i}"!!!!!!!!!!!!!!!!!!!"
   	mkdir -p ${outPth}${i}
   	cd ${outPth}${i}
   	
   	#### step1：prefilter from raw matrix ####
   	echo "!!!!!!!!!!!!!! "${i}": Step1__Prefilter from raw matrix!!!!!!!!!!!!!!!!!!!"
   	mkdir -p ./step1_prefilter
   	Rscript -e ".libPaths('$R_LIB_PATH'); source('${software}step1_raw_prefilter.R')" ${i} ${inPth} ${outPth}${i}/step1_prefilter ${software}
   	
   	### step2：process and predict ####
   	echo "!!!!!!!!!!!!!! "${i}": Step2__Process and predict !!!!!!!!!!!!!!!!!!!"
   	mkdir -p ./step2_predict
   	python ${software}step2_auto_predict.py -s ${i} \
   	-i ${outPth}${i}/step1_prefilter \
   	-o ${outPth}${i}/step2_predict \
   	--model ${model} \
   	--scaler ${scaler} \
   	--geneList ${geneList} \
   	--usegpu ${useGPU}
   
   	### step3：perform final annotations ####
   	echo "!!!!!!!!!!!!!! "${i}": Step3__Perform final annotations !!!!!!!!!!!!!!!!!!!"
   	mkdir -p ./step3_annotation
   	Rscript -e ".libPaths('$R_LIB_PATH'); source('${software}step3_final_annotation.R')" ${i} ${outPth}${i}/step1_prefilter ${outPth}${i}/step3_annotation ${outPth}${i}/step2_predict ${outPth}${i}
   	
   	cd ${outPth}
   	echo "!!!!!!!!!!!!!!Done:"${i}"!!!!!!!!!!!!!!!!!!!"
   done
   ```
   
   or choose parallel execution (alter"-j 3" parameter to change the number of parallel tasks) script: "para_scNeu.sh",and then execute scripts:
   
   ```shell
   bash ./run_scNeu.sh
   ```
   
   

