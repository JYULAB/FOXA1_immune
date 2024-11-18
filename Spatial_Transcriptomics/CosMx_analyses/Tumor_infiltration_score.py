# Tumor Immune Infiltration Score -- Split each fov into patches
import os,csv
import warnings
warnings.filterwarnings('ignore')
import slideio
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import cv2
from scipy.spatial.distance import cdist

cat_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
cat_color1=["#F56867","#FEB915","#C798EE","#59BE86"] #['Epithelial', 'Club Cells', 'NEPC', 'Tumor']
cat_color2=["#6D1A9C","#15821E","#3A84E6"] #['T-cell', 'B-cell', 'macrophage']
cat_color3=["#1167A9", "#FF7F0E", "#2C7D2C"] #['immune_cell', 'infiltrating_immune_cell', 'tumor_cell']
#['#1F77B4', '#FF7F0E', '#2CA02C'] (default cat color)
color_map1={
	'Epithelial': "#F56867",
	'Club Cells': "#FEB915",
	'NEPC': "#C798EE",
	'Tumor': "#59BE86"
}
color_map2={
	'T-cell': "#6D1A9C",
	'B-cell': "#15821E",
	'macrophage': "#3A84E6"
}
color_map3={
	'immune_cell': '#1167A9',
	'infiltrating_immune_cell': '#FF7F0E',
	'tumor_cell': '#2C7D2C'
}
cnt_color=clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)

os.chdir("/Users/jinghuang/Library/CloudStorage/Dropbox/Jian_Jing/Jindan_collaboration/data/CosMx")

#------------------ Step 1. Read in preprocessed adata ------------------
adata=sc.read("./adata_cell_type_upd.h5ad")
cell_type_upd_key="cell_types_upd"

#-------------------------------------------------------------------------
# annotate immune subtypes (no need to run again)
# CD4: Treg, T CD4 naive, T CD4 memory
# CD8: T CD8 naive, T CD8 memory
# macrophage

adata.obs["immune_subtype"]=adata.obs["cell_types_upd"].copy()
adata.obs["immune_subtype"]=adata.obs["immune_subtype"].astype(str)

# CD4:
cd4_indices=adata[adata.obs["cell_types"].isin(["Treg", "T CD4 memory", "T CD4 naive"])].obs.index.tolist()
adata.obs.loc[cd4_indices,"immune_subtype"]="T_CD4"
# CD8:
cd8_indices=adata[adata.obs["cell_types"].isin(["T CD8 memory", "T CD8 naive"])].obs.index.tolist()
adata.obs.loc[cd8_indices, "immune_subtype"]="T_CD8"
# macrophage:
macro_indices=adata[adata.obs["cell_types"]=="macrophage"].obs.index.tolist()
adata.obs.loc[macro_indices, "immune_subtype"]="macrophage"

adata.obs["immune_subtype"]=adata.obs["immune_subtype"].astype("category")
adata.obs["immune_subtype"].value_counts()

# save the updated adata
adata.write_h5ad("./adata_cell_type_upd.h5ad")

#------------------ Step 2. Determine the radius threshold (include all fovs) ------------------
# fov_list=cell_exp_dat["fov"].value_counts().index.tolist()
# fov_list=sorted(fov_list)
fov_list=[1, 2, 3, 4, 6, 7, 9, 10, 13, 14, 15, 16, 17, 18, 25] # 15 fovs

#Tumor cells: ["Tumor", "PSA High", "FASN High", "NEPC"]
#Immune cells: ["T-cell", "B-cell", "macrophage"]
target_tumor_cells=['Tumor', 'Club Cells', 'PSA High', 'FASN High', 'NEPC']
target_immune_cells=['T-cell', 'B-cell', 'macrophage']
total_cells=adata.obs[cell_type_upd_key].value_counts().index.tolist()
non_tumor_cells=[i for i in total_cells if i not in target_tumor_cells]
nn_small=5
nn_median=10
nn_large=20
dists_thres_small=[]
dists_thres_median=[]
dists_thres_large=[]

for fov_num in fov_list:
	fov_df=adata.obs.loc[adata.obs["fov"]==fov_num,["global_x","global_y",cell_type_upd_key]].copy()
	tumor_cells_num=fov_df.loc[fov_df[cell_type_upd_key].isin(target_tumor_cells),:].shape[0]
	non_tumor_cells_num=fov_df.loc[fov_df[cell_type_upd_key].isin(non_tumor_cells),:].shape[0]
	print("==================== fov"+str(fov_num)+" ====================")
	print("The number of tumor cells: "+str(tumor_cells_num))
	print("The number of non tumor cells: "+str(non_tumor_cells_num))
	if ((tumor_cells_num)>0) & ((non_tumor_cells_num)>0):
		tumor_df=fov_df.loc[fov_df[cell_type_upd_key].isin(target_tumor_cells),["global_x","global_y"]].copy()
		non_tumor_df=fov_df.loc[fov_df[cell_type_upd_key].isin(non_tumor_cells),["global_x","global_y"]].copy()
		dists=cdist(tumor_df.values,non_tumor_df.values)
		dists_df=pd.DataFrame(dists, index=tumor_df.index, columns=non_tumor_df.index)
		for _, row in dists_df.iterrows():
			sorted_row=sorted(row)
			if len(sorted_row)>=nn_large:
				dists_thres_small.append(sorted_row[nn_small-1])
				dists_thres_median.append(sorted_row[nn_median-1])
				dists_thres_large.append(sorted_row[nn_large-1])
			elif len(sorted_row)>=nn_median:
				dists_thres_small.append(sorted_row[nn_small-1])
				dists_thres_median.append(sorted_row[nn_median-1])
			elif len(sorted_row)>=nn_small:
				dists_thres_small.append(sorted_row[nn_small-1])


#determine the threshold
thres_small=np.mean(dists_thres_small)
thres_median=np.mean(dists_thres_median)
thres_large=np.mean(dists_thres_large)
print(thres_small)
print(thres_median)
print(thres_large)

# >>> print(thres_small)
# 295.52383062118884
# >>> print(thres_median)
# 385.3947562548894
# >>> print(thres_large)
# 517.2949241711834

dists_small_quantiles=np.quantile(dists_thres_small, np.arange(0,1.1,0.1))
dists_median_quantiles=np.quantile(dists_thres_median, np.arange(0,1.1,0.1))
dists_large_quantiles=np.quantile(dists_thres_large, np.arange(0,1.1,0.1))
print(dists_small_quantiles)
print(dists_median_quantiles)
print(dists_large_quantiles)

# >>> print(dists_small_quantiles)
# [  47.29693436  128.75752392  154.9225613   180.22763384  207.06762181
#   239.70919882  282.70833026  341.83768073  429.51600669  552.05117512
#  1244.50994371]
# >>> print(dists_median_quantiles)
# [  73.          186.13167382  219.08217636  249.01305182  281.51731741
#   321.27869522  369.74856322  435.33090862  538.86176335  690.5964813
#  1469.99217685]
# >>> print(dists_large_quantiles)
# [ 113.50770899  267.29571633  310.38685539  348.76066273  391.76013069
#   439.18731764  496.00403224  577.63352971  696.21907472  883.34081746
#  2068.29809264]

# final radius
nn=nn_small # nn = 5 nearest neighbors
dists_thres=dists_small_quantiles[5] # 50% percentile (239.70919881997082) | 50% percentile for fov 6, 7, 9, 10, 18, 25 is around 186
dists_thres=np.round(dists_thres) # 240
print(dists_thres)

#================================== check radius patterns ==================================
label_key="label"
d_props={}

for fov in fov_list:
	adata_sub=adata[(adata.obs["fov"]==fov) & (adata.obs[cell_type_upd_key].isin(target_tumor_cells+target_immune_cells)),:].copy()
	print("===== "+str(fov)+" =====")
	adata_sub.obs[cell_type_upd_key].value_counts()
	tumor_indices=adata_sub[adata_sub.obs[cell_type_upd_key].isin(target_tumor_cells)].obs.index.tolist()
	immune_indices=adata_sub[adata_sub.obs[cell_type_upd_key].isin(target_immune_cells)].obs.index.tolist()
	adata_sub.obs[label_key]="nan"
	adata_sub.obs[label_key]=adata_sub.obs[label_key].astype(str)
	adata_sub.obs.loc[tumor_indices,label_key]="tumor_cell"
	adata_sub.obs.loc[immune_indices,label_key]="immune_cell"
	if (len(tumor_indices)>0) & (len(immune_indices)>0):
		tumor_df=adata_sub.obs.loc[(adata_sub.obs[cell_type_upd_key].isin(target_tumor_cells)),["global_x","global_y"]].copy()
		immune_df=adata_sub.obs.loc[(adata_sub.obs[cell_type_upd_key].isin(target_immune_cells)),["global_x","global_y"]].copy()
		#identify immune infiltrating cells
		dists=cdist(tumor_df.values,immune_df.values)
		dists_df=pd.DataFrame(dists, index=tumor_df.index, columns=immune_df.index)
		index=(dists_df<=dists_thres)*1
		index_counts=np.sum(index, axis=0)
		infiltrating_indices=[index_counts.index[i] for i in range(len(index_counts)) if index_counts[i]>0]
		#calculate the proportion of infiltrating immune cells / all cells
		infiltrating_prop=len(infiltrating_indices)/adata[adata.obs["fov"]==fov].shape[0]
		d_props[fov]=infiltrating_prop
		#update the columns for immune infiltrating
		adata_sub.obs.loc[infiltrating_indices,label_key]="infiltrating_immune_cell"
		adata_sub.obs[label_key]=adata_sub.obs[label_key].astype("category")
		final_cell_types=adata_sub.obs[label_key].value_counts().index.tolist()
		cat_color_upd=[color_map3[i] for i in list(color_map3.keys()) if i in final_cell_types]
		#generate the plots to testify the radius threshold
		fig_title="fov "+str(fov)+": immune infiltrating cells (radius = "+str(int(dists_thres))+" | nn = "+str(nn)+")"
		fig_path="./immune_percentages/fov"+str(fov)+"_immune_infiltrating_cells_radius="+str(int(dists_thres))+"_nn="+str(nn)+".png"
		fig=sc.pl.scatter(adata_sub,alpha=1,x="global_y",y="global_x",color=label_key, palette=cat_color_upd, show=False,size=40)
		fig.set_aspect("equal","box")
		fig.invert_xaxis()
		fig.set_title(fig_title)
		fig.figure.savefig(fig_path,dpi=100)
		del adata_sub.uns[label_key+"_colors"]
		plt.clf()
		plt.close()


#==============================================================================================================
# results
# immune infiltrating proportion
d_props={1: 0.0012033694344163659, 
		 2: 0.0029875067897881585, 
		 3: 0.0056772100567721, 
		 4: 0.037900874635568516, 
		 6: 0.07955088389870998, 
		 7: 0.1320921985815603, 
		 9: 0.02301736765013601, 
		 10: 0.05704563650920737, 
		 13: 0.008609444299504304, 
		 14: 0.007062146892655367, 
		 15: 0.02115965924704589, 
		 16: 0.020276100086281276, 
		 17: 0.017672557519173057, 
		 18: 0.11805026656511805, 
		 25: 0.1311098332374928}


