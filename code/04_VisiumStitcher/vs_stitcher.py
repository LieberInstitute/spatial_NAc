import scanpy as sc
import visium_stitcher as vs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import scipy.sparse as sparse
import scipy.io as sio
import scipy.stats as stats

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# Read in the .csv with sample key
homeDir = "/data/abattle4/prashanthi/spatial_NAc/"
sample_key = pd.read_csv(homeDir + "raw-data/sample_key_spatial_NAc.csv")
transforms_all = pd.read_csv(homeDir + "processed-data/02_image_stitching/adjusted_transforms_highres.csv")
sample_key = sample_key[sample_key['In analysis'] == "Yes"]
sample_key = sample_key[sample_key['Refined transforms'] == "Yes"]

brain_list = []
for brain in sample_key.Brain.unique():
    print(brain)
    sample_key_brain = sample_key.loc[sample_key.Brain == brain].copy()
    sample_key_brain = sample_key_brain.reset_index(drop = True)
    adatas = []
    # Read in the data from all capture areas correspnding to the brain and specify the corrected transforms
    for i in range(sample_key_brain.shape[0]):
        adata = sc.read_visium(homeDir + "processed-data/01_spaceranger_reorg/" +
                       sample_key_brain.Slide[i] + "/outs", count_file="filtered_feature_bc_matrix.h5")
        adata.var_names_make_unique()
        adata.obs_names = [sample_key_brain.Slide[i]+j for j in adata.obs_names]
        adata.obs['sample'] = sample_key_brain.Slide[i]
        if(sample_key_brain.loc[i, 'Refined transforms'] == 'No'):
            transform = vs.transform_finder(homeDir + sample_key_brain['XML file name'][i])
            adata.uns['transform'] = transform[sample_key_brain['Transform_index'][i]]
        else:
            M = np.array([[1, 0, 0], [0, 1, 0]], np.float64)
            M[0,2] = transforms_all.x[transforms_all.sample_id == sample_key_brain.Slide[i]]
            M[1,2] = transforms_all.y[transforms_all.sample_id == sample_key_brain.Slide[i]]
            M[0,0] = np.cos(transforms_all.theta[transforms_all.sample_id == sample_key_brain.Slide[i]])
            M[0,1] = -np.sin(transforms_all.theta[transforms_all.sample_id == sample_key_brain.Slide[i]])
            M[1,0] = np.sin(transforms_all.theta[transforms_all.sample_id == sample_key_brain.Slide[i]])
            M[1,1] = np.cos(transforms_all.theta[transforms_all.sample_id == sample_key_brain.Slide[i]])
            adata.uns['transform'] = M
        adatas.append(adata)
    print("Completed reading the annData")
    # Get pairwise overlap for all possible regions
    for i in range(sample_key_brain.shape[0]):
        for j in range(sample_key_brain.shape[0]):
            if(j != i):
                adata_stitched = vs.stitch([adatas[j], adatas[i]])
                df = adata_stitched.obs
                df = df.loc[df['sample'] == sample_key_brain.Slide[i]]
                df = df.reindex(index = adatas[i].obs.index)
                adatas[i].obs['Overlap_' + sample_key_brain.Slide[j]] = df.overlap
    # Summarize the overlap in each capture area
    for i in range(sample_key_brain.shape[0]):
        final_overlap = adatas[i].obs.filter(regex='Overlap').astype(int).replace(0, np.nan).idxmax(axis=1).replace(np.nan,'None',regex=True).to_list()
        final_overlap = [sub.replace('Overlap_', '') for sub in final_overlap]
        adatas[i].obs['Overlap'] = final_overlap
        #sc.pl.spatial(adatas[i], color = 'Overlap')
    # Get one merged object
    if len(adatas) == 2:
        adata = adatas[0].concatenate(adatas[1])
    elif len(adatas) == 3:
        adata = adatas[0].concatenate(adatas[1], adatas[2])
    elif len(adatas) == 4:
        adata = adatas[0].concatenate(adatas[1], adatas[2], adatas[3])
    elif len(adatas) == 5: 
        adata = adatas[0].concatenate(adatas[1], adatas[2], adatas[3], adatas[4])
    else:
        print("We need to specify an option to merge these capture areas")

    # mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-') 
    # ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
    adata.obs['sample_overlap'] = adata.obs['sample'].astype(str) + "_Overlap_" + adata.obs['Overlap'].astype(str)
    #sc.pl.violin(adata, 'n_genes_by_counts',jitter=0.4, groupby = 'sample_overlap', rotation= 90)
    #sc.pl.violin(adata, 'total_counts',jitter=0.4, groupby = 'sample_overlap', rotation= 90)
    #sc.pl.violin(adata, 'pct_counts_mt',jitter=0.4, groupby = 'sample_overlap', rotation= 90)
    #sc.pl.violin(adata, 'pct_counts_ribo',jitter=0.4, groupby = 'sample_overlap', rotation= 90)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)

    sample_order = sample_key_brain.Slide.values.copy()
    for i in sample_key_brain.Slide.values:
        for j in sample_key_brain.Slide.values:
            if(i != j):
                N_ibj = sum(adata.obs.sample_overlap == (i + "_Overlap_" + j))
                N_jbi = sum(adata.obs.sample_overlap == (j + "_Overlap_" + i))
                if(N_ibj >= N_jbi):
                    ind_i = np.where(sample_order == i)[0][0]
                    ind_j = np.where(sample_order == j)[0][0]
                    if(ind_j < ind_i):
                        tmp = sample_order[ind_i]
                        sample_order[ind_i] = sample_order[ind_j]
                        sample_order[ind_j] = tmp
                else:
                    ind_i = np.where(sample_order == i)[0][0]
                    ind_j = np.where(sample_order == j)[0][0]
                    if(ind_j > ind_i):
                        tmp = sample_order[ind_i]
                        sample_order[ind_i] = sample_order[ind_j]
                        sample_order[ind_j] = tmp  
    # Get one stitched object
    if len(adatas) == 2:
        brain_ad = vs.stitch([adatas[np.where(sample_key_brain.Slide.values == sample_order[0])[0][0]],
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[1])[0][0]]],
                          image = homeDir + sample_key_brain['Stitched image file'][0])
    elif len(adatas) == 3:
        brain_ad = vs.stitch([adatas[np.where(sample_key_brain.Slide.values == sample_order[0])[0][0]], 
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[1])[0][0]], 
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[2])[0][0]]], 
                           image = homeDir + sample_key_brain['Stitched image file'][0])
    elif len(adatas) == 4:
        brain_ad = vs.stitch([adatas[np.where(sample_key_brain.Slide.values == sample_order[0])[0][0]], 
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[1])[0][0]], 
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[2])[0][0]], 
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[3])[0][0]]], 
                           image = homeDir + sample_key_brain['Stitched image file'][0])
    elif len(adatas) == 5: 
        brain_ad = vs.stitch([adatas[np.where(sample_key_brain.Slide.values == sample_order[0])[0][0]], 
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[1])[0][0]], 
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[2])[0][0]], 
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[3])[0][0]], 
                           adatas[np.where(sample_key_brain.Slide.values == sample_order[4])[0][0]]], 
                           image = homeDir + sample_key_brain['Stitched image file'][0])
    else:
        print("We need to specify an option to merge these capture areas")       
    brain_ad.var = adatas[0].var
    sc.pl.spatial(brain_ad, color="sample", show = False)
    plt.savefig("/data/abattle4/prashanthi/spatial_NAc/plots/04_VisiumStitcher/"+ brain + "/all_spots.pdf")
    sc.pl.spatial(brain_ad[~brain_ad.obs["overlap"]], color="sample", show = False)
    plt.savefig("/data/abattle4/prashanthi/spatial_NAc/plots/04_VisiumStitcher/"+ brain + "/non-overlapping_spots.pdf")
    # Write the data
    saveDir = "/data/abattle4/prashanthi/spatial_NAc/processed-data/04_VisiumStitcher/" + brain + "/"
    brain_ad.obs.rename(columns={'Overlap': 'overlap_slide', 'overlap': 'exclude_overlapping'}, inplace=True)
    brain_ad.obs.drop(['in_tissue','array_row', 'array_col'], axis=1).to_csv(saveDir + "vs_stitcher.csv")
    # and delete individual datasets to save space
    brain_list.append(brain_ad)
    del(adatas, adata, brain_ad)
