#!/bin/bash
# mamba activate scNeu
export R_HOME=${CONDA_PREFIX}/lib/R
export R_LIB_PATH=${CONDA_PREFIX}/lib/R/library

#################### args ####################
export software="/0.db/Pan/scNeu/"
export sampleList="/0.db/Pan/scNeu/example_sampleList.txt"
export inPth="/0.db/Pan/scNeu/example_rawdata/"
export outPth="/0.db/Pan/scNeu/example_output/"

################### config ###################
export useGPU=True  # Boolean: Choose True or False

##################### main #####################
mkdir -p ${outPth}
cd ${outPth}

run_sample() {
    i=$1
    echo "!!!!!!!!!!!!!!Start:"${i}"!!!!!!!!!!!!!!!!!!!"
    mkdir -p ${outPth}${i}
    cd ${outPth}${i}
    
	echo "Software path: ${software}"
	echo "inPth path: ${inPth}"
	echo "outPth path: ${outPth}"
	echo "R_LIB_PATH path: ${R_LIB_PATH}"
	echo "R_HOME path: ${R_HOME}"
	model=${software}"model.joblib"
	scaler=${software}"scaler.joblib"
	geneList=${software}"geneList.csv"
    
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
    Rscript -e ".libPaths('$R_LIB_PATH'); source('${software}step3_final_annotation.R')" ${i} ${outPth}${i}/step1_prefilter ${outPth}${i}/step3_annotation ${outPth}${i}/step2_predict ${outPth}${i} ${software}
    
    cd ${outPth}
    echo "!!!!!!!!!!!!!!Done:"${i}"!!!!!!!!!!!!!!!!!!!"
}

export -f run_sample

# Using GNU Parallel to run in parallel
cat $sampleList | parallel --env software,R_LIB_PATH,inPth,outPth,R_HOME,useGPU -j 3 run_sample
