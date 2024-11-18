# TESLA Cell Type Annotation
import os,csv,re,time
import pickle
import random
import warnings
warnings.filterwarnings('ignore')
import cv2
import numpy as np
import pandas as pd
import scanpy as sc
import TESLA as tesla
import json
import slideio
import matplotlib.colors as clr
import matplotlib.pyplot as plt
from scipy import stats
from scanpy import read_10x_h5
from scipy.sparse import issparse

cnt_color = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
os.chdir("/Users/jinghuang/Library/CloudStorage/Dropbox/Jian_Jing/Jindan_collaboration")

#================================== cell type markers ==================================
d_g={'T_CD8': ['PTPRC','CD3E','CD8A','CD8B1'], 
	 'T_CD4': ['PTPRC','CD3E','CD4'], 
	 'Macrophages': ['PTPRC','ITGAM','ADGRE1','MAFB','LYZ2','C1QB'], 
	 'Monocytes':['PTPRC','ITGAM','CD14','S100A8','S100A9'],
	 'Tumor': ['EPCAM']
}

#Macrophage gene 'C1QB' is highly in PF18.

d_num={'T_CD8': 3, 
	   'T_CD4': 3, 
	   'Macrophages': 6, 
	   'Monocytes': 5,
	   'Tumor': 1
}

d_s_g={
	'T_CD8': ['HAVCR2','PDCD1','GZMB'],
	'T_CD4': ['HAVCR2','PDCD1','FOXP3','RORC','TBX21','CCR7'],
	'Macrophages': ['NOS2','ARG1','MRC1','CD86'],
	'Monocytes': ['TGFB1','ARG1','CD274']
}

#================================== sample information ==================================
sample_list=["P12","PF12","P12_upd","P18","PF18","XP_upd","XPF"]
sample="P18"
data_dir="./data/"+sample+"/"

#read in data
res=50
cnt_name="cnt_scanx" #"cnt_scanx" for P12, PF12, P18, PF18, and XPF; while "cnt_cv2" for XP_upd.
enhanced_exp_adata=sc.read(data_dir+sample+"_sudo_log_s=1_res="+str(res)+"_nbr=10_k=2_"+cnt_name+".h5ad")
img=cv2.imread(data_dir+"he.tif")
resize_factor=1000/np.min(img.shape[0:2])
resize_width=int(img.shape[1]*resize_factor)
resize_height=int(img.shape[0]*resize_factor)
binary_scanx=np.load(data_dir+sample+"_scanx_cnt_binary.npy")
binary_cv2=np.load(data_dir+sample+"_cv2_cnt_binary.npy")
binary_resized_scanx=cv2.resize(binary_scanx, ((resize_width, resize_height)))
binary_resized_cv2=cv2.resize(binary_cv2, ((resize_width, resize_height)))

overall_ratio={}
cell_type_ratio={}

#----------------------------------------------
# P12 upd sample
sample="P12_upd"
data_dir="./data/"+sample+"/"

#read in data
res=50
cnt_name="cnt_scany" #"cnt_scanx" for P12, PF12, P18, PF18, and XPF; "cnt_scany" for P12_upd; and "cnt_cv2" for XP_upd.
enhanced_exp_adata=sc.read(data_dir+sample+"_sudo_log_s=1_res="+str(res)+"_nbr=10_k=2_"+cnt_name+".h5ad")

# read in image file 
P12_upd_slide = slideio.open_slide(data_dir+"he.tif", "GDAL")
scene=P12_upd_slide.get_scene(0)
scene.size
img=scene.read_block()
img=cv2.cvtColor(img,cv2.COLOR_RGB2BGR) #RGB -> BGR
img.shape

resize_factor=1000/np.min(img.shape[0:2])
resize_width=int(img.shape[1]*resize_factor)
resize_height=int(img.shape[0]*resize_factor)

# contour 
binary_scany=np.load(data_dir+sample+"_cnt_scany_binary.npy") # naming is not consistent as the above
binary_resized_scany=cv2.resize(binary_scany, ((resize_width, resize_height)))

overall_ratio={}
cell_type_ratio={}

#========================================================================================
#---------------------------------------- T_CD8+ ----------------------------------------
#========================================================================================
#---------- cell type ----------
cell_type="T_CD8"
genes=d_g[cell_type]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
num=d_num[cell_type]
print(genes)
print(num)
min_UMI=0

d_r={}
#target_size can be set to "small" or "large"
binary_detect=binary_scany
binary_visual=binary_cv2
pred_refined1, target_clusters1, c_m1=tesla.annotation(img=img, 
													#binary=binary_scanx,
													binary=binary_detect,
													sudo_adata=enhanced_exp_adata, 
													genes=genes, 
													resize_factor=resize_factor,
													num_required=num, 
													target_size="small",
													min_UMI=min_UMI)

#check target clusters
for cluster in target_clusters1:
	ret_img_cluster=tesla.visualize_annotation(img=img, 
										binary=binary_visual, #correct the boundary issue
										resize_factor=resize_factor,
										pred_refined=pred_refined1, 
										target_clusters=[cluster], 
										c_m=c_m1)
	cv2.imwrite("./results/"+sample+"/"+cell_type+"/annotation/"+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_cluster="+str(cluster)+".jpg", ret_img_cluster)

#update the detected clusters
target_clusters1=[]

#plot clusters
ret_img1=tesla.visualize_annotation(img=img, 
							#binary=binary_cv2, #correct the boundary issue
							binary=binary_visual,
							resize_factor=resize_factor,
							pred_refined=pred_refined1, 
							target_clusters=target_clusters1, 
							c_m=c_m1)

cv2.imwrite('./results/'+sample+'/'+cell_type+'/annotation/'+sample+'_'+cell_type+'_res='+str(res)+'_num='+str(num)+'_umi='+str(min_UMI)+'.jpg', ret_img1)

#cell type ratio
cell_type_ind=pd.Series(pred_refined1).astype(int).isin(target_clusters1).values
ratio=np.sum(cell_type_ind)/np.sum(binary_resized_scany)
cell_type_ratio[cell_type]=np.round(ratio*100,1)
print(cell_type_ratio)

#save the results
d_r["pred"]=pred_refined1
d_r["target_clusters"]=target_clusters1
d_r["c_m"]=c_m1
res_path=data_dir+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_pred_results.pkl"
with open(res_path,'wb') as f: pickle.dump(d_r,f)


#---------- subtype ----------
subtype_ratio={}
genes=d_s_g[cell_type]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
print(genes)
min_UMI=0

for g in genes:
	pred_refined2, target_clusters2, c_m2=tesla.annotation(img=img,
														#binary=binary_scanx,
														binary=binary_detect,
														sudo_adata=enhanced_exp_adata,
														genes=[g],
														resize_factor=resize_factor,
														num_required=1,
														target_size="small",
														min_UMI=min_UMI)
	pred_refined2[cell_type_ind==False]=-1
	ratio=np.sum(pd.Series(pred_refined2).astype(int).isin(target_clusters2))/np.sum(cell_type_ind)
	subtype_ratio[g]=np.round(ratio*100,1)
	ret_img2=tesla.visualize_annotation(img=ret_img1,
									binary=binary_resized_cv2,
									resize_factor=1,
									pred_refined=pred_refined2,
									target_clusters=target_clusters2,
									c_m=c_m2,
									cnt_color=clr.LinearSegmentedColormap.from_list('blue',["#EAE7CC", '#0000BA'], N=256))
	text=sample+": "+cell_type+" "+g+" ("+str(subtype_ratio[g])+"%)"
	ret_img2_upd=cv2.putText(ret_img2,text,(50,50),cv2.FONT_ITALIC,1.2,(0,0,0),2,cv2.LINE_AA)
	cv2.imwrite("./results/"+sample+"/"+cell_type+"/annotation/"+sample+"_"+cell_type+"_subtype_"+g+".jpg", ret_img2_upd)

print(subtype_ratio)
overall_ratio[cell_type]=subtype_ratio


#========================================================================================
#---------------------------------------- T_CD4+ ----------------------------------------
#========================================================================================
#---------- cell type ----------
cell_type="T_CD4"
genes=d_g[cell_type]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
num=d_num[cell_type]
print(genes)
print(num)
min_UMI=0

d_r={}
#target_size can be set to "small" or "large"
binary_detect=binary_scanx
binary_visual=binary_cv2
pred_refined1, target_clusters1, c_m1=tesla.annotation(img=img, 
													#binary=binary_scanx,
													binary=binary_detect,
													sudo_adata=enhanced_exp_adata, 
													genes=genes, 
													resize_factor=resize_factor,
													num_required=num, 
													target_size="small",
													min_UMI=min_UMI)

#check target clusters
for cluster in target_clusters1:
	ret_img_cluster=tesla.visualize_annotation(img=img, 
										binary=binary_visual, #correct the boundary issue
										resize_factor=resize_factor,
										pred_refined=pred_refined1, 
										target_clusters=[cluster], 
										c_m=c_m1)
	cv2.imwrite("./results/"+sample+"/"+cell_type+"/annotation/"+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_cluster="+str(cluster)+".jpg", ret_img_cluster)

#update the detected clusters
target_clusters1=[]

#plot clusters
ret_img1=tesla.visualize_annotation(img=img, 
							#binary=binary_cv2, #correct the boundary issue
							binary=binary_visual,
							resize_factor=resize_factor,
							pred_refined=pred_refined1, 
							target_clusters=target_clusters1, 
							c_m=c_m1)

cv2.imwrite('./results/'+sample+'/'+cell_type+'/annotation/'+sample+'_'+cell_type+'_res='+str(res)+'_num='+str(num)+'_umi='+str(min_UMI)+'.jpg', ret_img1)

#cell type ratio
cell_type_ind=pd.Series(pred_refined1).astype(int).isin(target_clusters1).values
ratio=np.sum(cell_type_ind)/np.sum(binary_resized_scanx)
cell_type_ratio[cell_type]=np.round(ratio*100,1)
print(cell_type_ratio)

#save the results
d_r["pred"]=pred_refined1
d_r["target_clusters"]=target_clusters1
d_r["c_m"]=c_m1
res_path=data_dir+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_pred_results.pkl"
with open(res_path,'wb') as f: pickle.dump(d_r,f)


#---------- subtype ----------
#add in one subtype gene of "CCR7"
#cell_type="T_CD4"
#genes=["CCR7"]
#genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
#num=d_num[cell_type]
#print(genes)
#print(num)#

#min_UMI=0#

#binary_detect=binary_scanx #binary_scanx for P12, PF12, P18, PF18, and XPF; while binary_cv2 for XP_upd.
#binary_visual=binary_cv2
#res_path=data_dir+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_pred_results.pkl"
#with open(res_path, 'rb') as f: d_r=pickle.load(f)#

#pred_refined1=d_r["pred"]
#target_clusters1=d_r["target_clusters"]
#ret_img1=cv2.imread('./results/'+sample+'/'+cell_type+'/annotation/'+sample+'_'+cell_type+'_res='+str(res)+'_num='+str(num)+'_umi='+str(min_UMI)+'.jpg')
#cell_type_ind=pd.Series(pred_refined1).astype(int).isin(target_clusters1).values#

#subtype_ratio={}


subtype_ratio={}
genes=d_s_g[cell_type]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
print(genes)
min_UMI=0

for g in genes:
	pred_refined2, target_clusters2, c_m2=tesla.annotation(img=img,
														#binary=binary_scanx,
														binary=binary_detect,
														sudo_adata=enhanced_exp_adata,
														genes=[g],
														resize_factor=resize_factor,
														num_required=1,
														target_size="small",
														min_UMI=min_UMI)
	pred_refined2[cell_type_ind==False]=-1
	ratio=np.sum(pd.Series(pred_refined2).astype(int).isin(target_clusters2))/np.sum(cell_type_ind)
	subtype_ratio[g]=np.round(ratio*100,1)
	ret_img2=tesla.visualize_annotation(img=ret_img1,
									binary=binary_resized_cv2,
									resize_factor=1,
									pred_refined=pred_refined2,
									target_clusters=target_clusters2,
									c_m=c_m2,
									cnt_color=clr.LinearSegmentedColormap.from_list('blue',["#EAE7CC", '#0000BA'], N=256))
	text=sample+": "+cell_type+" "+g+" ("+str(subtype_ratio[g])+"%)"
	ret_img2_upd=cv2.putText(ret_img2,text,(50,50),cv2.FONT_ITALIC,1.2,(0,0,0),2,cv2.LINE_AA)
	cv2.imwrite("./results/"+sample+"/"+cell_type+"/annotation/"+sample+"_"+cell_type+"_subtype_"+g+".jpg", ret_img2_upd)

print(subtype_ratio)
overall_ratio[cell_type]=subtype_ratio


#========================================================================================
#---------------------------------------- Macrophages ----------------------------------------
#========================================================================================
#---------- cell type ----------
cell_type="Macrophages"
genes=d_g[cell_type]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
num=d_num[cell_type]
print(genes)
print(num)
min_UMI=0

d_r={}
#target_size can be set to "small" or "large"
binary_detect=binary_scanx
binary_visual=binary_cv2
pred_refined1, target_clusters1, c_m1=tesla.annotation(img=img, 
													#binary=binary_scanx,
													binary=binary_detect,
													sudo_adata=enhanced_exp_adata, 
													genes=genes, 
													resize_factor=resize_factor,
													num_required=num, 
													target_size="small",
													min_UMI=min_UMI)

#check target clusters
for cluster in target_clusters1:
	ret_img_cluster=tesla.visualize_annotation(img=img, 
										binary=binary_visual, #correct the boundary issue
										resize_factor=resize_factor,
										pred_refined=pred_refined1, 
										target_clusters=[cluster], 
										c_m=c_m1)
	cv2.imwrite("./results/"+sample+"/"+cell_type+"/annotation/"+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_cluster="+str(cluster)+".jpg", ret_img_cluster)

#update the detected clusters
target_clusters1=[]

#plot clusters
ret_img1=tesla.visualize_annotation(img=img, 
							#binary=binary_cv2, #correct the boundary issue
							binary=binary_visual,
							resize_factor=resize_factor,
							pred_refined=pred_refined1, 
							target_clusters=target_clusters1, 
							c_m=c_m1)

cv2.imwrite('./results/'+sample+'/'+cell_type+'/annotation/'+sample+'_'+cell_type+'_res='+str(res)+'_num='+str(num)+'_umi='+str(min_UMI)+'.jpg', ret_img1)

#cell type ratio
cell_type_ind=pd.Series(pred_refined1).astype(int).isin(target_clusters1).values
ratio=np.sum(cell_type_ind)/np.sum(binary_resized_scanx)
cell_type_ratio[cell_type]=np.round(ratio*100,1)
print(cell_type_ratio)

#save the results
d_r["pred"]=pred_refined1
d_r["target_clusters"]=target_clusters1
d_r["c_m"]=c_m1
res_path=data_dir+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_pred_results.pkl"
with open(res_path,'wb') as f: pickle.dump(d_r,f)


#---------- subtype ----------
#add in two subtype genes of "MRC1" and "CD86"
cell_type="Macrophages"
genes=["MRC1","CD86"]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
num=d_num[cell_type]
print(genes)
print(num)

min_UMI=0
binary_detect=binary_scanx #binary_scanx for P12, PF12, P18, PF18, and XPF; while binary_cv2 for XP_upd.
binary_visual=binary_cv2
res_path=data_dir+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_pred_results.pkl"
with open(res_path, 'rb') as f: d_r=pickle.load(f)

pred_refined1=d_r["pred"]
target_clusters1=d_r["target_clusters"]
ret_img1=cv2.imread('./results/'+sample+'/'+cell_type+'/annotation/'+sample+'_'+cell_type+'_res='+str(res)+'_num='+str(num)+'_umi='+str(min_UMI)+'.jpg')
cell_type_ind=pd.Series(pred_refined1).astype(int).isin(target_clusters1).values

subtype_ratio={}


subtype_ratio={}
genes=d_s_g[cell_type]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
print(genes)
min_UMI=0

for g in genes:
	pred_refined2, target_clusters2, c_m2=tesla.annotation(img=img,
														#binary=binary_scanx,
														binary=binary_detect,
														sudo_adata=enhanced_exp_adata,
														genes=[g],
														resize_factor=resize_factor,
														num_required=1,
														target_size="small",
														min_UMI=min_UMI)
	pred_refined2[cell_type_ind==False]=-1
	ratio=np.sum(pd.Series(pred_refined2).astype(int).isin(target_clusters2))/np.sum(cell_type_ind)
	subtype_ratio[g]=np.round(ratio*100,1)
	ret_img2=tesla.visualize_annotation(img=ret_img1,
									binary=binary_resized_cv2,
									resize_factor=1,
									pred_refined=pred_refined2,
									target_clusters=target_clusters2,
									c_m=c_m2,
									cnt_color=clr.LinearSegmentedColormap.from_list('blue',["#EAE7CC", '#0000BA'], N=256))
	text=sample+": "+cell_type+" "+g+" ("+str(subtype_ratio[g])+"%)"
	ret_img2_upd=cv2.putText(ret_img2,text,(50,50),cv2.FONT_ITALIC,1.2,(0,0,0),2,cv2.LINE_AA)
	cv2.imwrite("./results/"+sample+"/"+cell_type+"/annotation/"+sample+"_"+cell_type+"_subtype_"+g+".jpg", ret_img2_upd)

print(subtype_ratio)
overall_ratio[cell_type]=subtype_ratio


#========================================================================================
#---------------------------------------- Monocytes ----------------------------------------
#========================================================================================
#---------- cell type ----------
cell_type="Monocytes"
genes=d_g[cell_type]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
num=d_num[cell_type]
print(genes)
print(num)
min_UMI=0

d_r={}
#target_size can be set to "small" or "large"
binary_detect=binary_scanx
binary_visual=binary_cv2
pred_refined1, target_clusters1, c_m1=tesla.annotation(img=img, 
													#binary=binary_scanx,
													binary=binary_detect,
													sudo_adata=enhanced_exp_adata, 
													genes=genes, 
													resize_factor=resize_factor,
													num_required=num, 
													target_size="small",
													min_UMI=min_UMI)

#check target clusters
for cluster in target_clusters1:
	ret_img_cluster=tesla.visualize_annotation(img=img, 
										binary=binary_visual, #correct the boundary issue
										resize_factor=resize_factor,
										pred_refined=pred_refined1, 
										target_clusters=[cluster], 
										c_m=c_m1)
	cv2.imwrite("./results/"+sample+"/"+cell_type+"/annotation/"+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_cluster="+str(cluster)+".jpg", ret_img_cluster)

#update the detected clusters
target_clusters1=[]

#plot clusters
ret_img1=tesla.visualize_annotation(img=img, 
							#binary=binary_cv2, #correct the boundary issue
							binary=binary_visual,
							resize_factor=resize_factor,
							pred_refined=pred_refined1, 
							target_clusters=target_clusters1, 
							c_m=c_m1)

cv2.imwrite('./results/'+sample+'/'+cell_type+'/annotation/'+sample+'_'+cell_type+'_res='+str(res)+'_num='+str(num)+'_umi='+str(min_UMI)+'.jpg', ret_img1)

#cell type ratio
cell_type_ind=pd.Series(pred_refined1).astype(int).isin(target_clusters1).values
ratio=np.sum(cell_type_ind)/np.sum(binary_resized_scanx)
cell_type_ratio[cell_type]=np.round(ratio*100,1)
print(cell_type_ratio)

#save the results
d_r["pred"]=pred_refined1
d_r["target_clusters"]=target_clusters1
d_r["c_m"]=c_m1
res_path=data_dir+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_pred_results.pkl"
with open(res_path,'wb') as f: pickle.dump(d_r,f)


#---------- subtype ----------
subtype_ratio={}
genes=d_s_g[cell_type]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
print(genes)
min_UMI=0

for g in genes:
	pred_refined2, target_clusters2, c_m2=tesla.annotation(img=img,
														#binary=binary_scanx,
														binary=binary_detect,
														sudo_adata=enhanced_exp_adata,
														genes=[g],
														resize_factor=resize_factor,
														num_required=1,
														target_size="small",
														min_UMI=min_UMI)
	pred_refined2[cell_type_ind==False]=-1
	ratio=np.sum(pd.Series(pred_refined2).astype(int).isin(target_clusters2))/np.sum(cell_type_ind)
	subtype_ratio[g]=np.round(ratio*100,1)
	ret_img2=tesla.visualize_annotation(img=ret_img1,
									binary=binary_resized_cv2,
									resize_factor=1,
									pred_refined=pred_refined2,
									target_clusters=target_clusters2,
									c_m=c_m2,
									cnt_color=clr.LinearSegmentedColormap.from_list('blue',["#EAE7CC", '#0000BA'], N=256))
	text=sample+": "+cell_type+" "+g+" ("+str(subtype_ratio[g])+"%)"
	ret_img2_upd=cv2.putText(ret_img2,text,(50,50),cv2.FONT_ITALIC,1.2,(0,0,0),2,cv2.LINE_AA)
	cv2.imwrite("./results/"+sample+"/"+cell_type+"/annotation/"+sample+"_"+cell_type+"_subtype_"+g+".jpg", ret_img2_upd)

print(subtype_ratio)
overall_ratio[cell_type]=subtype_ratio

#========== Summarize all the ratios ==========
#overall_ratio[sample]=cell_type_ratio
#ratio_path="./results/"+sample+"/"+sample+"_overall_ratios.txt"
#with open(ratio_path,'w') as file:
#	for key,value in overall_ratio.items():
#		file.write(f'{key}: {value}\n')


#========================================================================================
#---------------------------------------- Tumor ----------------------------------------
#========================================================================================
cell_type="Tumor"
genes=d_g[cell_type]
genes=[i for i in genes if i in enhanced_exp_adata.var.index.tolist()]
num=d_num[cell_type]
print(genes)
print(num)
min_UMI=0

#========================== Step 1. Detect tumor regions using EPCAM ==========================
d_r={}
#target_size can be set to "small" or "large"
binary_detect=binary_scanx
binary_visual=binary_cv2
pred_refined1, target_clusters1, c_m1=tesla.annotation(img=img, 
													#binary=binary_scanx,
													binary=binary_detect,
													sudo_adata=enhanced_exp_adata, 
													genes=genes, 
													resize_factor=resize_factor,
													num_required=num, 
													target_size="small",
													min_UMI=min_UMI)

#check target clusters
for cluster in target_clusters1:
	ret_img_cluster=tesla.visualize_annotation(img=img, 
										binary=binary_visual, #correct the boundary issue
										resize_factor=resize_factor,
										pred_refined=pred_refined1, 
										target_clusters=[cluster], 
										c_m=c_m1)
	cv2.imwrite("./results/"+sample+"/"+cell_type+"/annotation/"+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_cluster="+str(cluster)+".jpg", ret_img_cluster)

#update the detected clusters
target_clusters1=[]

#plot clusters
ret_img1=tesla.visualize_annotation(img=img, 
							#binary=binary_cv2, #correct the boundary issue
							binary=binary_visual,
							resize_factor=resize_factor,
							pred_refined=pred_refined1, 
							target_clusters=target_clusters1, 
							c_m=c_m1)

cv2.imwrite('./results/'+sample+'/'+cell_type+'/annotation/'+sample+'_'+cell_type+'_res='+str(res)+'_num='+str(num)+'_umi='+str(min_UMI)+'.jpg', ret_img1)

#cell type ratio
cell_type_ind=pd.Series(pred_refined1).astype(int).isin(target_clusters1).values
ratio=np.sum(cell_type_ind)/np.sum(binary_resized_scanx)
cell_type_ratio[cell_type]=np.round(ratio*100,1)
print(cell_type_ratio)

#save the results
d_r["pred"]=pred_refined1
d_r["target_clusters"]=target_clusters1
d_r["c_m"]=c_m1
res_path=data_dir+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_pred_results.pkl"
with open(res_path,'wb') as f: pickle.dump(d_r,f)

#========== Summarize all the ratios ==========
overall_ratio[sample]=cell_type_ratio
ratio_path="./results/"+sample+"/"+sample+"_overall_ratios.txt"
with open(ratio_path,'w') as file:
	for key,value in overall_ratio.items():
		file.write(f'{key}: {value}\n')


#========================== Step 2. Overlay tumor regions with T_CD4, T_CD8, and Macrophages ==========================
#read in data
sample_list=["P12","PF12","P18","PF18","XP_upd","XPF"]
sample="P12"
cnt_name="cnt_scanx" #or "cnt_cv2" for XP_upd samples ("cnt_scanx" for most cases)
res=50
min_UMI=0
data_dir="./data/"+sample+"/"
#enhanced_exp_adata=sc.read(data_dir+sample+"_sudo_log_s=1_res="+str(res)+"_nbr=10_k=2_"+cnt_name+".h5ad")

img=cv2.imread(data_dir+"he.tif")
resize_factor=1000/np.min(img.shape[0:2])
resize_width=int(img.shape[1]*resize_factor)
resize_height=int(img.shape[0]*resize_factor)
img_resized=cv2.resize(img, ((resize_width, resize_height)))
binary_cv2=np.load(data_dir+sample+"_"+cnt_name+"_binary.npy")
binary_resized_cv2=cv2.resize(binary_cv2, ((resize_width, resize_height)))

#overlay
#==================== overlay ====================
overlay_cell_types=["Tumor","T_CD4","T_CD8","Macrophages","Monocytes"]
#overlay_cell_types=["Tumor","Monocytes"]

target_clusters_dic={}
pred_refined_dic={}
#read in saved results
for cell_type in overlay_cell_types:
	num=d_num[cell_type]
	res_path=data_dir+sample+"_"+cell_type+"_res="+str(res)+"_num="+str(num)+"_umi="+str(min_UMI)+"_pred_results.pkl"
	with open(res_path, 'rb') as f: d_r=pickle.load(f)
	target_clusters=d_r["target_clusters"]
	pred_refined=d_r["pred"]
	target_clusters_dic[cell_type]=target_clusters
	pred_refined_dic[cell_type]=pred_refined

cnt_color_tumor=clr.LinearSegmentedColormap.from_list('red', ["#EAE7CC", '#BA0000'], N=256)
cnt_color_immune=clr.LinearSegmentedColormap.from_list('blue', ["#EAE7CC", '#08038B'], N=256)
cnt_color_macro=clr.LinearSegmentedColormap.from_list('green', ["#EAE7CC", '#006600'], N=256)
cnt_color_mono=clr.LinearSegmentedColormap.from_list('green', ["#EAE7CC", '#006600'], N=256)

#-------------------- A. Tumor + T_CD4 --------------------
overlay_pairs=["Tumor","T_CD4"]
target_clusters0=target_clusters_dic[overlay_pairs[0]]
target_clusters1=target_clusters_dic[overlay_pairs[1]]
pred_refined0=pred_refined_dic[overlay_pairs[0]]
pred_refined1=pred_refined_dic[overlay_pairs[1]]

ret_img=img_resized.copy()

#whiten
white_ratio=0.5
ret_img=ret_img*(1-white_ratio)+np.array([255, 255, 255])*(white_ratio)
#ret_img[binary_resized_cv2!=0]=ret_img[binary_resized_cv2!=0]*(1-white_ratio)+np.array([255, 255, 255])*(white_ratio)
#cv2.imwrite("./results/"+sample+"/Overlay/"+sample+"_"+overlay_pairs[0]+"+"+overlay_pairs[1]+"_whiteratio="+str(white_ratio)+".jpg",ret_img)

alpha=0.8

#Tumor regions
for i in range(len(target_clusters0)):
	color=((np.array(cnt_color_tumor(int((len(target_clusters0)-i)/len(target_clusters0)*255)))[0:3])*255).astype(int)[::-1]
	target_img=(1*(pred_refined0==target_clusters0[i])).reshape(resize_height, resize_width)
	target_img[binary_resized_cv2==0]=0
	ret_img[target_img!=0]=ret_img[target_img!=0]*(1-alpha)+np.array(color)*(alpha)

#T_CD4 regions
for i in range(len(target_clusters1)):
	color=((np.array(cnt_color_immune(int((len(target_clusters1)-i)/len(target_clusters1)*255)))[0:3])*255).astype(int)[::-1]
	target_img=(1*(pred_refined1==target_clusters1[i])).reshape(resize_height, resize_width)
	target_img[binary_resized_cv2==0]=0
	ret_img[target_img!=0]=ret_img[target_img!=0]*(1-alpha)+np.array(color)*(alpha)

cv2.imwrite("./results/"+sample+"/Overlay/"+sample+"_"+overlay_pairs[0]+"+"+overlay_pairs[1]+"_white="+str(white_ratio)+"_alpha="+str(alpha)+".jpg",ret_img)

#-------------------- B. Tumor + T_CD8 --------------------
overlay_pairs=["Tumor","T_CD8"]
target_clusters0=target_clusters_dic[overlay_pairs[0]]
target_clusters1=target_clusters_dic[overlay_pairs[1]]
pred_refined0=pred_refined_dic[overlay_pairs[0]]
pred_refined1=pred_refined_dic[overlay_pairs[1]]

ret_img=img_resized.copy()

#whiten
white_ratio=0.5
ret_img=ret_img*(1-white_ratio)+np.array([255, 255, 255])*(white_ratio)
#ret_img[binary_resized_cv2!=0]=ret_img[binary_resized_cv2!=0]*(1-white_ratio)+np.array([255, 255, 255])*(white_ratio)
#cv2.imwrite("./results/"+sample+"/Overlay/"+sample+"_"+overlay_pairs[0]+"+"+overlay_pairs[1]+"_whiteratio="+str(white_ratio)+".jpg",ret_img)

alpha=0.8

#Tumor regions
for i in range(len(target_clusters0)):
	color=((np.array(cnt_color_tumor(int((len(target_clusters0)-i)/len(target_clusters0)*255)))[0:3])*255).astype(int)[::-1]
	target_img=(1*(pred_refined0==target_clusters0[i])).reshape(resize_height, resize_width)
	target_img[binary_resized_cv2==0]=0
	ret_img[target_img!=0]=ret_img[target_img!=0]*(1-alpha)+np.array(color)*(alpha)

#T_CD4 regions
for i in range(len(target_clusters1)):
	color=((np.array(cnt_color_immune(int((len(target_clusters1)-i)/len(target_clusters1)*255)))[0:3])*255).astype(int)[::-1]
	target_img=(1*(pred_refined1==target_clusters1[i])).reshape(resize_height, resize_width)
	target_img[binary_resized_cv2==0]=0
	ret_img[target_img!=0]=ret_img[target_img!=0]*(1-alpha)+np.array(color)*(alpha)

cv2.imwrite("./results/"+sample+"/Overlay/"+sample+"_"+overlay_pairs[0]+"+"+overlay_pairs[1]+"_white="+str(white_ratio)+"_alpha="+str(alpha)+".jpg",ret_img)

#-------------------- C. Tumor + Macrophages --------------------
overlay_pairs=["Tumor","Macrophages"]
target_clusters0=target_clusters_dic[overlay_pairs[0]]
target_clusters1=target_clusters_dic[overlay_pairs[1]]
pred_refined0=pred_refined_dic[overlay_pairs[0]]
pred_refined1=pred_refined_dic[overlay_pairs[1]]

ret_img=img_resized.copy()

#whiten
white_ratio=0.5
ret_img=ret_img*(1-white_ratio)+np.array([255, 255, 255])*(white_ratio)
#ret_img[binary_resized_cv2!=0]=ret_img[binary_resized_cv2!=0]*(1-white_ratio)+np.array([255, 255, 255])*(white_ratio)
#cv2.imwrite("./results/"+sample+"/Overlay/"+sample+"_"+overlay_pairs[0]+"+"+overlay_pairs[1]+"_whiteratio="+str(white_ratio)+".jpg",ret_img)

alpha=0.8

#Tumor regions
for i in range(len(target_clusters0)):
	color=((np.array(cnt_color_tumor(int((len(target_clusters0)-i)/len(target_clusters0)*255)))[0:3])*255).astype(int)[::-1]
	target_img=(1*(pred_refined0==target_clusters0[i])).reshape(resize_height, resize_width)
	target_img[binary_resized_cv2==0]=0
	ret_img[target_img!=0]=ret_img[target_img!=0]*(1-alpha)+np.array(color)*(alpha)

#T_CD4 regions
for i in range(len(target_clusters1)):
	color=((np.array(cnt_color_macro(int((len(target_clusters1)-i)/len(target_clusters1)*255)))[0:3])*255).astype(int)[::-1]
	target_img=(1*(pred_refined1==target_clusters1[i])).reshape(resize_height, resize_width)
	target_img[binary_resized_cv2==0]=0
	ret_img[target_img!=0]=ret_img[target_img!=0]*(1-alpha)+np.array(color)*(alpha)

cv2.imwrite("./results/"+sample+"/Overlay/"+sample+"_"+overlay_pairs[0]+"+"+overlay_pairs[1]+"_white="+str(white_ratio)+"_alpha="+str(alpha)+".jpg",ret_img)

#-------------------- D. Tumor + Monocytes --------------------
overlay_pairs=["Tumor","Monocytes"]
target_clusters0=target_clusters_dic[overlay_pairs[0]]
target_clusters1=target_clusters_dic[overlay_pairs[1]]
pred_refined0=pred_refined_dic[overlay_pairs[0]]
pred_refined1=pred_refined_dic[overlay_pairs[1]]

ret_img=img_resized.copy()

#whiten
white_ratio=0.5
ret_img=ret_img*(1-white_ratio)+np.array([255, 255, 255])*(white_ratio)
#ret_img[binary_resized_cv2!=0]=ret_img[binary_resized_cv2!=0]*(1-white_ratio)+np.array([255, 255, 255])*(white_ratio)
#cv2.imwrite("./results/"+sample+"/Overlay/"+sample+"_"+overlay_pairs[0]+"+"+overlay_pairs[1]+"_whiteratio="+str(white_ratio)+".jpg",ret_img)

alpha=0.8

#Tumor regions
for i in range(len(target_clusters0)):
	color=((np.array(cnt_color_tumor(int((len(target_clusters0)-i)/len(target_clusters0)*255)))[0:3])*255).astype(int)[::-1]
	target_img=(1*(pred_refined0==target_clusters0[i])).reshape(resize_height, resize_width)
	target_img[binary_resized_cv2==0]=0
	ret_img[target_img!=0]=ret_img[target_img!=0]*(1-alpha)+np.array(color)*(alpha)

#T_CD4 regions
for i in range(len(target_clusters1)):
	color=((np.array(cnt_color_mono(int((len(target_clusters1)-i)/len(target_clusters1)*255)))[0:3])*255).astype(int)[::-1]
	target_img=(1*(pred_refined1==target_clusters1[i])).reshape(resize_height, resize_width)
	target_img[binary_resized_cv2==0]=0
	ret_img[target_img!=0]=ret_img[target_img!=0]*(1-alpha)+np.array(color)*(alpha)

cv2.imwrite("./results/"+sample+"/Overlay/"+sample+"_"+overlay_pairs[0]+"+"+overlay_pairs[1]+"_white="+str(white_ratio)+"_alpha="+str(alpha)+".jpg",ret_img)

