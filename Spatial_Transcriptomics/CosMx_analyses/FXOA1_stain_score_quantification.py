# FOXA1 stain intensity -- Split each fov into patches
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
from scipy.stats import pearsonr, spearmanr
from sklearn.cluster import KMeans

os.chdir("/Users/jinghuang/Library/CloudStorage/Dropbox/Jian_Jing/Jindan_collaboration/data/CosMx")

#========================================= Step 1. Read in datasets =========================================
# read in the whole whole image data
FOXA1_slide = slideio.open_slide("FOXA1_IHC_update.tif","GDAL")
scene=FOXA1_slide.get_scene(0)
scene.size #(26471, 18328)
img=scene.read_block()
img=cv2.cvtColor(img,cv2.COLOR_RGB2BGR) #RGB -> BGR
img.shape #(18328, 26471, 3)

# read in cell expresion data
cell_exp_dat=pd.read_csv("cell_exp_dat.csv", header=0, index_col=0, na_filter=False)
cell_exp_dat["global_x"]=round(cell_exp_dat["global_x"]).astype("int").tolist()
cell_exp_dat["global_y"]=round(cell_exp_dat["global_y"]).astype("int").tolist()
# fov_list=cell_exp_dat["fov"].value_counts().index.tolist()
# fov_list=sorted(fov_list)

# read in IF spatial coordinates
IF_coords=pd.read_csv("FOXA1_IHCstain_FOV_coordinates.csv", header=0, index_col=None, na_filter=False)

# read in gene expression adata
adata=sc.read("./adata_cell_type_upd.h5ad")

#========================================= Step 2. Read in masks and gray images =========================================
fov_list=[6, 7, 9, 10, 18] # exclude 25

# read in tumor masks that Luly extract
tumor_mask_dic={}
for fov_num in fov_list:
	# read in tumor mask
	tumor_mask=np.load("./FOXA1_quantification/FOV"+str(fov_num)+"_tumor_mask.npy")
	tumor_mask_dic[fov_num]=tumor_mask
	# print(tumor_mask.shape)
	# print(np.unique(tumor_mask))


# read in FOXA1 positive regions
stain_positive_dic={}
for fov_num in fov_list:
	stain_positive_mask=np.load("./FOXA1_quantification/FOV"+str(fov_num)+"_FOXA1_positive_mask.npy")
	stain_positive_dic[fov_num]=stain_positive_mask
	# print(stain_positive_mask.shape)
	# print(np.unique(stain_positive_mask))


# read in FOXA1 gray images
stain_gray_dic={}
for fov_num in fov_list:
	gray_img=cv2.imread("./FOXA1_quantification/FOV"+str(fov_num)+"_FOXA1_positive_gray.png")
	stain_gray_dic[fov_num]=gray_img
	# print(gray_img.shape)
	# print(np.unique(gray_img))


#========================================= Step 3. Split fovs into patches =========================================
#----------------------------------------- Patches: 2x2 -----------------------------------------
patch_num=2

# patch names
patch_names_list=[]
for fov_num in fov_list:
	for row_idx in range(patch_num):
		for col_idx in range(patch_num):
			patch_name="fov"+str(fov_num)+"_"+str(row_idx)+"_"+str(col_idx)
			print(patch_name)
			patch_names_list.append(patch_name)


patch_names_list=['fov6_0_0', 'fov6_0_1', 'fov6_1_0', 'fov6_1_1', 'fov7_0_0', 'fov7_0_1', 'fov7_1_0', 'fov7_1_1', 'fov9_0_0', 'fov9_0_1', 'fov9_1_0', 'fov9_1_1', 'fov10_0_0', 'fov10_0_1', 'fov10_1_0', 'fov10_1_1', 'fov18_0_0', 'fov18_0_1', 'fov18_1_0', 'fov18_1_1']
metrics=["tumor_mask_prop","FOXA1_positive_prop","FOXA1_over_tumor_prop","FOXA1_intensity","stain_score","immune_prop"]
quan_df=pd.DataFrame(np.zeros((len(fov_list*patch_num*patch_num), len(metrics))), index=patch_names_list, columns=metrics)

#----------------------------------------- Immune percentages -----------------------------------------
# first add a patch label column
adata.obs["patch_name"]="nan"

for fov_num in np.unique(cell_exp_dat["fov"]).tolist():
	# CosMX coords
	fov_gene_coords=adata.obs.loc[adata.obs["fov"]==fov_num, ["global_x", "global_y"]].copy()
	# rotate 180 degrees to map gene exp with stain image coords
	fov_gene_coords["x_tmp"]=np.max(fov_gene_coords["global_x"])-fov_gene_coords["global_x"]
	fov_gene_coords["y_tmp"]=np.max(fov_gene_coords["global_y"])-fov_gene_coords["global_y"]
	# split each fov into patches
	fov_gene_width=fov_gene_coords["global_x"].max()-fov_gene_coords["global_x"].min()
	fov_gene_height=fov_gene_coords["global_y"].max()-fov_gene_coords["global_y"].min()
	print("fov"+str(fov_num)+" gene exp height = "+str(fov_gene_height)+" | width = "+str(fov_gene_width))
	patch_width=int(fov_gene_width/patch_num)
	patch_height=int(fov_gene_height/patch_num)
	print("patch height = "+str(patch_height)+" | patch width = "+str(patch_width))
	for row_idx in range(patch_num):
		for col_idx in range(patch_num):
			patch_name="fov"+str(fov_num)+"_"+str(row_idx)+"_"+str(col_idx)
			patch_gene_coords=fov_gene_coords[(fov_gene_coords["x_tmp"]>=(row_idx*patch_height)) &
											  (fov_gene_coords["x_tmp"]<=((row_idx+1)*patch_height)) &
											  (fov_gene_coords["y_tmp"]>=(col_idx*patch_width)) &
											  (fov_gene_coords["y_tmp"]<=((col_idx+1)*patch_width))]
			patch_indices=patch_gene_coords.index.tolist()
			adata.obs.loc[patch_indices, "patch_name"]=patch_name
			# test patch coords (match -- fov 1)
			# adata_sub=adata[adata.obs["patch_name"]==patch_name].copy()
			# fig_title=patch_name+" cell type annotations"
			# fig_path="./immune_percentages/"+patch_name+"_cell_type_annotations.png"
			# fig=sc.pl.scatter(adata_sub,alpha=1,x="global_y",y="global_x",color=cell_type_upd_key, palette=cat_color, show=False,size=40)
			# fig.set_aspect("equal","box")
			# fig.invert_xaxis()
			# fig.set_title(fig_title)
			# fig.figure.savefig(fig_path,dpi=100)
			# del adata_sub.uns[cell_type_upd_key+"_colors"]
			# plt.clf()
			# plt.close()


# calculate immune percentages
cell_type_upd_key="cell_types_upd"
label_key="label"
target_tumor_cells=['Tumor', 'Club Cells', 'PSA High', 'FASN High', 'NEPC']
target_immune_cells=['T-cell', 'B-cell', 'macrophage']
color_map3={
	'immune_cell': '#1167A9', # blue
	'infiltrating_immune_cell': '#FF7F0E', # orange
	'tumor_cell': '#2C7D2C' # green
}

# prespecified radius
nn=5
dists_thres=240

d_props={}
for fov_num in fov_list:
	for row_idx in range(patch_num):
		for col_idx in range(patch_num):
			patch_name="fov"+str(fov_num)+"_"+str(row_idx)+"_"+str(col_idx)
			adata_sub=adata[(adata.obs["patch_name"]==patch_name) & (adata.obs[cell_type_upd_key].isin(target_tumor_cells+target_immune_cells))].copy()
			print("===== "+str(patch_name)+" =====")
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
				infiltrating_prop=len(infiltrating_indices)/adata[adata.obs["patch_name"]==patch_name].shape[0]
				d_props[patch_name]=infiltrating_prop
				quan_df.loc[patch_name, "immune_prop"]=infiltrating_prop
				#update the columns for immune infiltrating
				adata_sub.obs.loc[infiltrating_indices,label_key]="infiltrating_immune_cell"
				adata_sub.obs[label_key]=adata_sub.obs[label_key].astype("category")
				final_cell_types=adata_sub.obs[label_key].value_counts().index.tolist()
				cat_color_upd=[color_map3[i] for i in list(color_map3.keys()) if i in final_cell_types]
				#generate the plots to visualize immune infiltration
				fig_title=patch_name+": immune infiltrating cells (radius = "+str(int(dists_thres))+" | nn = "+str(nn)+")"
				fig_path="./immune_percentages/patch_num=2/"+patch_name+"_immune_infiltrating_cells_radius="+str(int(dists_thres))+"_nn="+str(nn)+".png"
				fig=sc.pl.scatter(adata_sub,alpha=1,x="global_y",y="global_x",color=label_key, palette=cat_color_upd, show=False,size=65)
				fig.set_aspect("equal","box")
				fig.invert_xaxis()
				fig.set_title(fig_title)
				fig.figure.savefig(fig_path,dpi=100)
				del adata_sub.uns[label_key+"_colors"]
				plt.clf()
				plt.close()


print(d_props)
d_props={'fov6_0_0': 0.10540788267644363, 
		 'fov6_0_1': 0.0769980506822612, 
		 'fov6_1_0': 0.059443911792905084, 
		 'fov6_1_1': 0.07472959685349066, 
		 'fov7_0_0': 0.19983753046303818, 
		 'fov7_0_1': 0.127826941986234, 
		 'fov7_1_0': 0.1054421768707483, 
		 'fov7_1_1': 0.08886810102899906, 
		 'fov9_0_0': 0.02358887952822241, 
		 'fov9_0_1': 0.01972062448644207, 
		 'fov9_1_0': 0.02545743834526651, 
		 'fov9_1_1': 0.022401433691756272, 
		 'fov10_0_0': 0.04114189756507137, 
		 'fov10_0_1': 0.028148148148148148, 
		 'fov10_1_0': 0.10336341263330599, 
		 'fov10_1_1': 0.05867970660146699, 
		 'fov18_0_0': 0.19808612440191387, 
		 'fov18_0_1': 0.052534562211981564, 
		 'fov18_1_0': 0.14594594594594595, 
		 'fov18_1_1': 0.07394766780432309}

immune_prop=[]
for patch_name, prop in d_props.items():
	immune_prop.append(prop)


immune_prop=[0.10540788267644363, 0.0769980506822612, 0.059443911792905084, 0.07472959685349066, 0.19983753046303818, 0.127826941986234, 0.1054421768707483, 0.08886810102899906, 0.02358887952822241, 0.01972062448644207, 0.02545743834526651, 0.022401433691756272, 0.04114189756507137, 0.028148148148148148, 0.10336341263330599, 0.05867970660146699, 0.19808612440191387, 0.052534562211981564, 0.14594594594594595, 0.07394766780432309]

#----------------------------------------- FOXA1 intensity score -----------------------------------------
# split tumor mask
split2_tumor_mask_dic={}
for fov_num in fov_list:
	tumor_mask=tumor_mask_dic[fov_num]
	print("fov"+str(fov_num)+" tumor mask height = "+str(tumor_mask.shape[0])+" | width = "+str(tumor_mask.shape[1]))
	patch_height=int(tumor_mask.shape[0]/patch_num)
	patch_width=int(tumor_mask.shape[1]/patch_num)
	print("patch height = "+str(patch_height)+" | patch width = "+str(patch_width))
	for row_idx in range(patch_num):
		for col_idx in range(patch_num):
			patch_name="fov"+str(fov_num)+"_"+str(row_idx)+"_"+str(col_idx)
			patch_fov=tumor_mask[(row_idx*patch_height):((row_idx+1)*patch_height),(col_idx*patch_width):((col_idx+1)*patch_width)].copy() # 2 channels
			split2_tumor_mask_dic[patch_name]=patch_fov
			# print(patch_fov.shape)
			cv2.imwrite("./FOXA1_quantification/patches/patch_num=2/"+patch_name+"_tumor_mask.png", patch_fov*255)


# split FOXA1 positive mask
split2_stain_positive_dic={}
for fov_num in fov_list:
	stain_positive_mask=stain_positive_dic[fov_num]
	print("fov"+str(fov_num)+" stain positive mask height = "+str(stain_positive_mask.shape[0])+" | width = "+str(stain_positive_mask.shape[1]))
	patch_height=int(stain_positive_mask.shape[0]/patch_num)
	patch_width=int(stain_positive_mask.shape[1]/patch_num)
	print("patch height = "+str(patch_height)+" | patch width = "+str(patch_width))
	for row_idx in range(patch_num):
		for col_idx in range(patch_num):
			patch_name="fov"+str(fov_num)+"_"+str(row_idx)+"_"+str(col_idx)
			patch_fov=stain_positive_mask[(row_idx*patch_height):((row_idx+1)*patch_height),(col_idx*patch_width):((col_idx+1)*patch_width)].copy() # 2 channels
			split2_stain_positive_dic[patch_name]=patch_fov
			# print(patch_fov.shape)
			cv2.imwrite("./FOXA1_quantification/patches/patch_num=2/"+patch_name+"_FOXA1_positive_mask.png", patch_fov*255)


# split FOXA1 positive gray image -> exp intensity
split2_stain_gray_dic={}
for fov_num in fov_list:
	gray_img=stain_gray_dic[fov_num]
	print("fov"+str(fov_num)+" stain gray image height = "+str(gray_img.shape[0])+" | width = "+str(gray_img.shape[1]))
	patch_height=int(gray_img.shape[0]/patch_num)
	patch_width=int(gray_img.shape[1]/patch_num)
	print("patch height = "+str(patch_height)+" | patch width = "+str(patch_width))
	for row_idx in range(patch_num):
		for col_idx in range(patch_num):
			patch_name="fov"+str(fov_num)+"_"+str(row_idx)+"_"+str(col_idx)
			patch_fov=gray_img[(row_idx*patch_height):((row_idx+1)*patch_height),(col_idx*patch_width):((col_idx+1)*patch_width),:].copy() # 3 channels
			split2_stain_gray_dic[patch_name]=patch_fov
			# print(patch_fov.shape)
			cv2.imwrite("./FOXA1_quantification/patches/patch_num=2/"+patch_name+"_FOXA1_positive_gray.png", patch_fov)


# evaluate intensity score
# quantification 
for fov_num in fov_list:
	for row_idx in range(patch_num):
		for col_idx in range(patch_num):
			# read in tumor mask, stain positive mask, and stain gray image
			patch_name="fov"+str(fov_num)+"_"+str(row_idx)+"_"+str(col_idx)
			tumor_mask=split2_tumor_mask_dic[patch_name]
			stain_positive_mask=split2_stain_positive_dic[patch_name]
			gray_img=split2_stain_gray_dic[patch_name]
			# calculate FOXA1 percentages
			size_ratio=(tumor_mask.shape[0]*tumor_mask.shape[1])/(stain_positive_mask.shape[0]*stain_positive_mask.shape[1])
			tumor_prop=np.sum(tumor_mask)/(tumor_mask.shape[0]*tumor_mask.shape[1])
			FOXA1_positive_prop=np.sum(stain_positive_mask)/(stain_positive_mask.shape[0]*stain_positive_mask.shape[1])
			FOXA1_over_tumor_prop=(size_ratio*np.sum(stain_positive_mask))/np.sum(tumor_mask)
			quan_df.loc[patch_name, "tumor_mask_prop"]=tumor_prop
			quan_df.loc[patch_name, "FOXA1_positive_prop"]=FOXA1_positive_prop
			quan_df.loc[patch_name, "FOXA1_over_tumor_prop"]=FOXA1_over_tumor_prop
			# FOXA1 exp intensity
			quan_df.loc[patch_name, "FOXA1_intensity"]=np.mean(gray_img[gray_img!=255])/255 # scale to 0-1


quan_df["stain_score"]=quan_df["FOXA1_over_tumor_prop"]*quan_df["FOXA1_intensity"]


#----------------------------------------- Correlation analyses -----------------------------------------
y=quan_df["stain_score"].copy()
x=quan_df["immune_prop"].copy()

# correlation tests
pearsonr(x, y) # pearson: normal
spearmanr(x, y) # spearman: rank-based

PearsonRResult(statistic=-0.7370596287487247, pvalue=0.00020942848292534446)
SignificanceResult(statistic=-0.5203007518796992, pvalue=0.018683221143120664)

# generate a scatter plot
plt.scatter(x, y, color='blue')
plt.title("FOXA1 stain score = FOXA1 proportion * FOXA1 intensity")
plt.xlabel("immune infiltrating proportion")
plt.ylabel("FOXA1 stain score")
plt.show()
plt.clf()
plt.close()

# save results
quan_df.to_csv("./immune_percentages/patch_num=2/quan_res_patch=2by2.csv", index=True)

