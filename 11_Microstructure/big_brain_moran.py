import pandas as pd
import numpy as np
import nibabel as nib
import sklearn.metrics
from neuromaps import stats
from neuromaps.datasets import fetch_atlas
from brainspace.null_models.moran import MoranRandomization
import pickle

datapath = 'insuladev/HCPDbold_1/BigBrain/'

n_perm = 10000
seed = 6789

outpath = 'insuladev/dropbox/'
import os
os.makedirs(outpath, exist_ok=True)

# Load fsLR atlas
fslr = fetch_atlas('fsLR', '32k')
vertices = {
    'lh': nib.load(fslr['inflated'].L).agg_data('NIFTI_INTENT_POINTSET'),
    'rh': nib.load(fslr['inflated'].R).agg_data('NIFTI_INTENT_POINTSET')
}

# Define file paths for different histology types
hist_types = {
    'Micro G1': {
        'lh': datapath + 
        'tpl-fs_LR_hemi-L_den-32k_desc-Micro_G1.shape.gii',
        'rh': datapath + 
        'tpl-fs_LR_hemi-R_den-32k_desc-Micro_G1.shape.gii'
    },
    'Micro G2': {
        'lh': datapath + 
        'tpl-fs_LR_hemi-L_den-32k_desc-Micro_G2.shape.gii',
        'rh': datapath + 
        'tpl-fs_LR_hemi-R_den-32k_desc-Micro_G2.shape.gii'
    }
}

# PC1 and insula indices
pc1_files = {
    'lh': datapath + 'vertex_PC1_scores_L.gii',
    'rh': datapath + 'vertex_PC1_scores_R.gii'
}

label_files = {
    'lh': datapath + 'left_hemisphere_insula_labels_mw.csv',
    'rh': datapath + 'right_hemisphere_insula_labels_mw.csv'
}

# Load PC1 and label indices
pc1_data = {
    hemi: nib.load(pc1_files[hemi]).darrays[0].data for hemi in ['lh', 'rh']
}
insula_indices = {
    hemi: pd.read_csv(label_files[hemi])['Vertex'].values - 1
    for hemi in ['lh', 'rh']
}

# Compute distance matrices
distance_matrices = {}
for hemi in ['lh', 'rh']:
    verts = vertices[hemi][insula_indices[hemi], :]
    dist = sklearn.metrics.pairwise_distances(verts)
    np.fill_diagonal(dist, 1)
    dist **= -1
    distance_matrices[hemi] = dist.astype('float64')

# Run spin tests for all hist types
results = {}
mrs = MoranRandomization(joint=True, tol=1e-6, n_rep=n_perm, random_state=seed)

for hist_label, paths in hist_types.items():
    results[hist_label] = {}
    for hemi in ['lh', 'rh']:
        data = nib.load(paths[hemi]).darrays[0].data
        data_insula = data[insula_indices[hemi]]
        pc1_insula = pc1_data[hemi][insula_indices[hemi]]
        dist = distance_matrices[hemi]

        hsurr = mrs.fit(dist).randomize(data_insula).T
        
        corr, pval, null_dist = stats.compare_images(
            data_insula, pc1_insula,
            metric='spearmanr',
            nulls=hsurr,
            ignore_zero=False,
            return_nulls=True
        )
        null_vec = np.asarray(null_dist).squeeze()
        results[hist_label][hemi] = {
            'observed': float(corr),
            'p_value': float(pval),
            'null': null_vec
        }
        
        print(f'{hist_label} ({hemi}): r = {corr:.3f}, p = {pval:.3f}')

        safe_label = hist_label.replace(' ', '_')  # e.g. "Micro G1" -> "Micro_G1"
        out_file = os.path.join(outpath, f'{safe_label}_{hemi}_moran.pkl')
        with open(out_file, 'wb') as f:
            pickle.dump(results[hist_label][hemi], f)
