# required: mamba activate scNeu or using singularity
import argparse
import configparser
import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
from rpy2.robjects import r
import anndata2ri
# import rpy2.robjects as ro
# from rpy2.robjects import pandas2ri

# 激活转换器
anndata2ri.activate()
# pandas2ri.activate()
r('.libPaths("$CONDA_PREFIX/lib/R/library")')
r('library(qs)')
r('library(dplyr)')
r('library(Seurat)')
r('library(SingleCellExperiment)')

def process_sample(file_path,geneList_path):
    """
    处理单个样本
    返回numpy格式的特征矩阵X和标签y
    """
    r(f'seurat_object <- qread("{file_path}")')
    r('sce_object <- as.SingleCellExperiment(seurat_object)')
    adata = r('sce_object')

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 归一化 nCount_RNA 和 nFeature_RNA 除以总和,再乘以目标总计数(这里是10000) 最后取log(1+x)
    total_counts = adata.obs['nCount_RNA'].sum()
    adata.obs['nCount'] = np.log1p((adata.obs['nCount_RNA'] / total_counts) * 1e4)
    total_feature = adata.obs['nFeature_RNA'].sum()
    adata.obs['nFeature'] = np.log1p((adata.obs['nCount_RNA'] / total_feature) * 1e4)

    # 筛选基因/细胞子集
    gene_list_df = pd.read_csv(geneList_path)
    gene_list = gene_list_df.iloc[:, 0].tolist()
    adata = adata[:, adata.var_names.isin(gene_list)]
    
    X = pd.DataFrame(adata.X.toarray(),columns=adata.var_names)
    X['nCount'] = adata.obs['nCount'].values
    X['nFeature'] = adata.obs['nFeature'].values
    y = adata.obs['celltypes']
    y_binary = np.where(y == 'Neu', 1, 0)
    y = y_binary
    unique, counts = np.unique(y, return_counts=True)
    category_counts = dict(zip(unique, counts))
    print("类别及其计数：", category_counts)

    return X.values, y
    

def fitSVM(file_path,model_path,scaler_path): 
    import rmm
    import cupy as cp
    from cuml.preprocessing import StandardScaler
    from cuml.svm import SVC
    from cuml.metrics.accuracy import accuracy_score
    from joblib import load

    # rmm 初始化 统一管理内存显存
    rmm.reinitialize(managed_memory=True)
    cp.cuda.set_allocator(rmm.rmm_cupy_allocator)

    # 加载模型和标准化器
    model = load(model_path)
    scaler = load(scaler_path)

    # 加载特征数据
    X = np.load(os.path.join(file_path, 'features.npy'))
    # 检查是否存在标签数据
    labels_path = os.path.join(file_path, 'labels.npy')
    y = np.load(labels_path) if os.path.exists(labels_path) else None

    # 使用加载的标准化器对外部数据进行标准化
    X_scaled = scaler.transform(X)

    # 使用模型进行预测
    y_pred = model.predict(X_scaled)

    if y is not None:
        # 评估模型在外部数据集上的性能
        accuracy = accuracy_score(y, y_pred)
        print("Accuracy:", accuracy)
        return y_pred,accuracy
    else:
        print("No labels provided; returning predictions.")
        return y_pred,None

def fitSVM_CPU(file_path,model_path,scaler_path): 
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    from sklearn.metrics import accuracy_score
    from joblib import load

    # 加载模型和标准化器
    model = load(model_path)
    scaler = load(scaler_path)

    # 加载特征数据
    X = np.load(os.path.join(file_path, 'features.npy'))
    # 检查是否存在标签数据
    labels_path = os.path.join(file_path, 'labels.npy')
    y = np.load(labels_path) if os.path.exists(labels_path) else None

    # 使用加载的标准化器对外部数据进行标准化
    X_scaled = scaler.transform(X)

    # 使用模型进行预测
    y_pred = model.predict(X_scaled)

    if y is not None:
        # 评估模型在外部数据集上的性能
        accuracy = accuracy_score(y, y_pred)
        print("Accuracy:", accuracy)
        return y_pred,accuracy
    else:
        print("No labels provided; returning predictions.")
        return y_pred,None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="<<<scNeu: A pipeline for identify neutrophils from single-cell data>>>")
    # 定义参数
    parser.add_argument("-s", "--sample", required=True, help="Sample name or ID")
    parser.add_argument("-i", "--input_path", required=True, help="Input data dir path ")
    parser.add_argument("-o", "--output_path", required=True, help="Output dir path")
    parser.add_argument("-m", "--model", required=True, help="Model")
    parser.add_argument("-c", "--scaler", required=True, help="Scaler")
    parser.add_argument("-g", "--geneList", required=True, help="GeneList")
    parser.add_argument("-u", "--usegpu", type=bool, required=True, help="Use NVIDIA GPU or not (boolean)")

    args = parser.parse_args()
    sample_name = args.sample
    input_path = args.input_path
    output_path = args.output_path
    model_path = args.model
    scaler_path = args.scaler
    geneList_path = args.geneList
    use_gpu = args.usegpu
    
    # 显示提示信息
    print(f"!!!!!  <INFO> Sample Name: {sample_name}  !!!!!")
    print(f"!!!!!  <INFO> Input  Path: {input_path}   !!!!!")
    print(f"!!!!!  <INFO> Output Path: {output_path}  !!!!!")
    print(f"!!!!!  <INFO> Model Path: {model_path}  !!!!!")
    print(f"!!!!!  <INFO> Scaler Path: {scaler_path}  !!!!!")
    print(f"!!!!!  <INFO> GeneList Path: {geneList_path}  !!!!!")
    print(f"!!!!!  <INFO> GPU Used: {use_gpu}  !!!!!")
    
    # 检查输入目录
    if not os.path.exists(args.input_path):
        print(f"Error: Input path doesn't '{args.input_path}' exist!", file=sys.stderr)
        sys.exit(1)
    
    # 检查输出目录
    if not os.path.exists(args.output_path):
        print(f"Output path '{args.output_path}' doesn't exist, creating...")
        os.makedirs(args.output_path)
        print(f"Creating Complete!")
    
    # Step1：根据.prefilter.qs文件生成features和labels的numpy数组
    # qs_path = next(os.path.join(dp, f) for dp, dn, filenames in os.walk(input_path) for f in filenames if f.endswith('.prefilter.qs'))
    qs_path = os.path.join(input_path, f"{sample_name}.prefilter.qs")
    X, y = process_sample(qs_path,geneList_path)
    
    # Step2: 保存features和labels
    np.save(os.path.join(output_path, 'features.npy'), X)
    np.save(os.path.join(output_path, 'labels.npy'), y)
    
    # Step2: 根据features矩阵生成预测的predict数组
    if use_gpu == True:
        y_pred, acc = fitSVM(output_path, model_path,scaler_path)
    else:
        y_pred, acc = fitSVM_CPU(output_path, model_path,scaler_path)
    
    # Step3: 统计与保存
    pred_path = os.path.join(output_path, 'predictions.csv')
    pd.DataFrame(y_pred).to_csv(pred_path, index=False, header=False)
    
    if acc is not None:
        count_filename = f"{acc}"
        count_filepath = os.path.join(output_path, count_filename)
        with open(count_filepath, 'w') as f:
            f.write(f"Count of '0's: {acc}\n")
    print(f"Data and prediction file saved to {output_path}")
    
