#conda activate mouse_retro
# wget -O  CEMBA.epiretro.mcds.tar.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230782/suppl/GSE230782%5FCEMBA%2Eepiretro%2Emcds%2Etar%2Egz"
# gunzip CEMBA.epiretro.mcds.tar.gz
# tar -xvf CEMBA.epiretro.mcds.tar

import pandas as pd

import anndata
import scanpy as sc
from ALLCools.mcds import MCDS
from ALLCools.clustering import log_scale
var_dim_er ='geneslop2k'
chrom_to_remove = ['chrX', 'chrY', 'chrM', 'chrL']
excluded_L1_annot = ['FALSE']
mc_type='CHN'
region_to_subregion = {'CTX': ['MOp','SSp','ACA','AI','RSP','AUD','PTLp','VIS'],
                       'HIP': ['CAa','CAp','DGa','DGp'],
                       'RHP': ['ENT'],
                       'OLF': ['MOB'],
                       'PIR': ['PIRa','PIRp'],
                       'STR': ['STR'],
                       'PAL': ['PAL'],
                       'AMY': ['AMY'],
                       'TH': ['THm','THl','THp'],
                       'HY': ['HY'],
                       'MB': ['SC','MRN','VTA','PAG','IC'],
                       'HB': ['P','MY']
                      }
rs2_subregion_dict = {
    'MOp': ['MOp'],
    'SSp': ['SSp'],
    'ACA': ['ACA'],
    'AI': ['AI'],
    'RSP': ['RSP'],
    'AUD': ['AUD'],
    # Entries from here on follow the same pattern as the first six
    'PTLp': ['PTLp'],
    'VIS': ['VIS'],
    'ENT': ['ENT'],
    'CAa': ['CAa'],
    'CAp': ['CAp'],
    'DGa': ['DGa'],
    'DGp': ['DGp'],
    'PIRa': ['PIRa'],
    'PIRp': ['PIRp'],
    'MOB': ['MOB'],
    'PAL': ['PAL'],
    'STR': ['STR'],
    'AMY': ['AMY'],
    'THl': ['THl'], 
    'THm': ['THm'], 
    'THp': ['THp'],
    'HY': ['HY'],
    'SC': ['SC'],
    'MRN': ['MRN'],
    'VTA': ['VTA'],
    'PAG': ['PAG'],
    'IC': ['IC'],
    'P': ['P'],
    'MY': ['MY']
}
selected_ER_slice = ['STR']


meta_path = '/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_revisions/03-projections_retroSeq/cell_48032_RS2_meta_nooutlier.csv'
meta = pd.read_csv(meta_path, index_col=0, header=0)
meta

selc = meta.index[meta['Source'].isin(selected_ER_slice) & ~meta['PassTargetFilter'].isin(excluded_L1_annot)]
print(len(selc))

mcds = MCDS.open('/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_revisions/03-projections_retroSeq/CEMBA.epiretro.mcds', var_dim=var_dim_er, use_obs=selc)
mcds

er = mcds.get_adata(mc_type='CHN')
er

# Subset the meta DataFrame using the indices in selc
# This ensures that only the metadata for the selected cells is included
filtered_meta = meta.loc[selc].copy()

# Align the indices of 'er' with 'filtered_meta' and add the metadata to 'er.obs'
# Assuming 'er' is already subsetted to match the cells in 'selc'
er.obs = filtered_meta

# Check the first few rows to ensure that the metadata has been added correctly
print(er.obs.head())
# Drop Cocluster and Subclass because they are empty
er.obs = er.obs.drop(columns=["Cocluster", "Subclass"])

log_scale(er, with_mean=True)
er.X = -er.X

path_2 = '/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/25_revisions/03-projections_retroSeq/rs2_mch_matrix.h5ad'
er.write_h5ad(path_2)