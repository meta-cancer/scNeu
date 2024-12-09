#!/bin/bash

#################### args ####################
export sif_path="/0.db/Pan/scNeu/scNeu.sif"
export software="/0.db/Pan/scNeu/"
export sampleList="/0.db/Pan/scNeu/example_sampleList.txt"
export inPth="/0.db/Pan/scNeu/example_rawdata/"
export outPth="/0.db/Pan/scNeu/example_output/"

################### config ###################
useGPU=False # Boolean: Choose True or False

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


