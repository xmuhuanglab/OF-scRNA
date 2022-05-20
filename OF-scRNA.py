import click
import configparser

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

@click.version_option(version='1.0.0')
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
# Extract cell
def extract_cells(parameters_path):
    
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg

    cp = configparser.ConfigParser()
    cp.read(parameters_path)

    # workdir
    identification_dir = cp.get('workdir', 'identification_dir')
    
    anno_col = cp.get('Extract', 'anno_col')
    celltype = cp.get('Extract', 'celltype_interest')
    
    batch_col = cp.get('Extract', 'batch_col')
    batch = eval(cp.get('Extract', 'batch'))
    
    
    total_data_path = cp.get('workdir', 'data_path') + '/' + cp.get('Extract', 'total_data_path')
    
    
    if not os.path.exists(identification_dir):
        os.makedirs(identification_dir)
    
    os.chdir(identification_dir)

    data = sc.read(total_data_path)

    if batch is None and celltype is None:
        logg.error("Please pass the batch or celltype parameter")
    elif batch is not None and celltype is None:
        print("Extract " + batch + ' cells')
        adata_sub = data[data.obs[batch_col].isin([batch])]
        filename = batch + '.h5ad'
        adata_sub.write_h5ad(filename, compression='gzip')
    elif batch is None and celltype is not None:
        print("Extract " + celltype)
        adata_sub = data[data.obs[anno_col].isin([celltype])]
        filename = celltype + '.h5ad'
        adata_sub.write_h5ad(filename, compression='gzip')
    elif batch is not None and celltype is not None:
        print("Extract " + batch + ' ' + celltype)
        adata_sub = data[(data.obs[batch_col].isin([batch])) & (data.obs[anno_col].isin([celltype])) ]
        filename = batch + '_' + celltype + '.h5ad'
        adata_sub.write_h5ad(filename, compression='gzip')


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def scanpy_workflow(parameters_path):
    '''
    min_mean
    max_mean
    min_disp
    n_comps
    show_var
    show_pc
    save_dir
    '''
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg
    
    
    
    cp = configparser.ConfigParser()
    cp.read(parameters_path)
    
    
    
    data_path = cp.get('scanpy', 'data_path')
    
    min_mean = eval(cp.get('scanpy', 'min_mean'))
    max_mean = eval(cp.get('scanpy', 'max_mean'))
    min_disp = eval(cp.get('scanpy', 'min_disp'))
    n_comps = eval(cp.get('scanpy', 'n_comps'))
    show_var = eval(cp.get('scanpy', 'show_var'))
    show_pc = eval(cp.get('scanpy', 'show_pc'))
    
    
    identification_dir = cp.get('workdir', 'identification_dir')
    
    if not os.path.exists(identification_dir):
        os.makedirs(identification_dir)
    
    os.chdir(identification_dir)
    
    
    data = sc.read(data_path)
    
    # scanpy basic workflow
    sc.pp.highly_variable_genes(data, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
    
    data.raw = data
    data = data[:, data.var.highly_variable]
    sc.pp.scale(data, max_value=10)
    sc.tl.pca(data, svd_solver='arpack', n_comps = n_comps)
    sc.pp.neighbors(data, n_neighbors=20, n_pcs=25)
    if show_var:
        sc.pl.pca_variance_ratio(data, log=True, n_pcs = show_pc, save = 'pca_var.pdf', show = False)
    
    celltype = cp.get('Extract', 'celltype_interest')
    csv_path = cp.get('scanpy', 'csv_path')
    
    data.write_h5ad(celltype+'_nor'+'.h5ad', compression='gzip')
    data.raw.to_adata().write_csvs('. ' + csv_path, skip_data=False)

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def CSS(parameters_path):
    
    """
    input_dir = os.path.abspath('./scanpy2seurat/')
    output_filename = os.path.abspath('tnf_css.h5ad')
    work_dir = os.path.abspath('./')
    css_script = os.path.abspath('../../code/1_1.CSS.R')
    batch_tag = 'batch'
    css_cluster_resolution = 0.5
    seurat_resolution = 0.55
    k_neighbor = 15
    """
    
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg
    
    
    
    cp = configparser.ConfigParser()
    cp.read(parameters_path)
    
    identification_dir = cp.get('workdir', 'identification_dir')
    
    if not os.path.exists(identification_dir):
        os.makedirs(identification_dir)
    
    os.chdir(identification_dir)
    os.getcwd()
    
    work_dir = cp.get('workdir', 'work_dir')
    
    css_script = work_dir + cp.get('CSS', 'css_script')
    input_dir = '.' + cp.get('scanpy', 'csv_path')
    batch_tag = cp.get('CSS', 'batch_tag')
    
    celltype = cp.get('Extract', 'celltype_interest')
    
    output_filename = celltype + '_css' + '.h5ad'
    
    
    css_cluster_resolution = eval(cp.get('CSS', 'css_cluster_resolution'))
    seurat_resolution = eval(cp.get('CSS', 'seurat_resolution'))
    k_neighbor = eval(cp.get('CSS', 'k_neighbor'))
    
    
    os.system("Rscript {} -i {} -b {} -o {} -w {} -c {} -s {} -k {}".format(css_script, input_dir, batch_tag, 
                                                                            output_filename, identification_dir, css_cluster_resolution, seurat_resolution, k_neighbor))
    
    tmp = sc.read(output_filename)
    
    csv_path = cp.get('scanpy', 'csv_path')
    
    var_path = '.' + csv_path + '/' + 'var.csv'
    var = pd.read_csv(var_path, index_col=0)
    
    if tmp.raw is None:
        tmp.raw = tmp
        tmp = tmp[:, var.highly_variable]
    
    
    tmp.write_h5ad(output_filename, compression='gzip')

    

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def findOncoFetal(parameters_path):
    """
    batch = 'batch', cell_anno = 'subtype', var_path = './scanpy2seurat/var.csv', width = 8, height = 4, path2save = './results/1.Identification/'
    """
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg
    
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    
    cp = configparser.ConfigParser()
    cp.read(parameters_path)
    
    
    identification_dir = cp.get('workdir', 'identification_dir')
    
    if not os.path.exists(identification_dir):
        os.makedirs(identification_dir)
    
    os.chdir(identification_dir)
    
    
    
    batch = cp.get('findOncofetal', 'batch_col')
    anno_col = cp.get('findOncofetal', 'anno_col')
    
    celltype = cp.get('Extract', 'celltype_interest')

    width = eval(cp.get('findOncofetal', 'width'))
    height = eval(cp.get('findOncofetal', 'height') )
    
    save = '.' + cp.get('findOncofetal', 'save')
    
    data_path = cp.get('findOncofetal', 'data_path')
    tmp = sc.read(data_path)

    
    # Ro/e
    a = stats.chi2_contingency(pd.crosstab(tmp.obs[anno_col], tmp.obs[batch]))

    ct = pd.crosstab(tmp.obs[anno_col], tmp.obs[batch])

    ct = ct.divide(a[3])
    
    ct.to_csv('relative_enrichment.csv')


    sc.pp.neighbors(tmp, n_neighbors=20, n_pcs=tmp.obsm['X_css'].shape[1], use_rep='X_css')    

    sc.tl.paga(tmp, groups=anno_col)
    
    df = scv.get_df(tmp, 'paga/connectivities', precision=3, columns=tmp.obs[anno_col].values.categories.to_list(), index = tmp.obs[anno_col].values.categories.to_list()).T
    
    # df.style.background_gradient(cmap='Blues').format('{:.2g}').to_excel('./paga_total.xlsx')


    df_plot = df[ct[(ct['Fetal'] > 1) & (ct['Tumor'] < 1) & (ct['Adj Normal'] < 1)].index.to_list()].loc[ct[ct['Fetal'] <= 1].index.to_list()]
    df_plot.to_csv('paga_subset.csv')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (width, height))

    g1 = sns.heatmap(ax=ax1, data = ct, annot=True, cbar_kws={'label': 'Ro/e'}, cmap = 'Blues_r')
    g1.set_xlabel('')
    g1.set_ylabel('')

    g2 = sns.heatmap(ax=ax2, data = df_plot, annot=True, cbar_kws={'label': 'Connectivity'}, fmt='.2f')
    fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    
    plt.savefig('FindOncoFetal.pdf', bbox_inches='tight')
    tmp.write_h5ad(celltype + '_final.h5ad', compression='gzip')

    tmp.write_csvs(save, skip_data =False)
    
    celltype = cp.get('Extract', 'celltype_interest')
    
    convert_script = cp.get('workdir', 'script_dir') + '/' + cp.get('findOncofetal', 'convert_script')
    
    os.system("Rscript {} -f {} -c {} -w {}".format(convert_script, save, celltype, identification_dir) )


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def findgene(parameters_path):
    """
    batch = 'batch', cell_anno = 'subtype', case_gene = [''], top_n = 50, tumor_query = None, fetal_ref = None, work_dir = '2.Gene'
    """
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg
    
    
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    
    cp = configparser.ConfigParser()
    cp.read(parameters_path)
    
    
    
    work_dir = cp.get('workdir', 'work_dir')
    tn_markers = cp.get('genesimilarity', 'tn_markers')
    fetal_markers = cp.get('genesimilarity', 'fetal_markers')
    
    full_tn_markers = cp.get('genesimilarity', 'full_tn_markers')
    full_fetal_markers = cp.get('genesimilarity', 'full_fetal_markers')
    
    tumor_query = cp.get('genesimilarity', 'tumor_query')
    fetal_ref = cp.get('genesimilarity', 'fetal_ref')
    
    batch = cp.get('findOncofetal', 'batch_col')
    anno_col = cp.get('findOncofetal', 'anno_col')
    print(fetal_ref)
    
    top_n = eval(cp.get('genesimilarity', 'top_n'))
    
    genesimi_script = work_dir + cp.get('genesimilarity', 'genesimi_script')
    
    genesimi_dir = cp.get('workdir', 'genesimi_dir')
    
    if not os.path.exists(genesimi_dir):
        os.makedirs(genesimi_dir)
    
    os.chdir(genesimi_dir)
   
    identification_dir = cp.get('workdir', 'identification_dir')
    
    celltype = cp.get('Extract', 'celltype_interest')
    
    tmp = sc.read(identification_dir + '/' + celltype + '_final.h5ad')
    
    ct = pd.read_csv(identification_dir + '/relative_enrichment.csv', index_col=0)

    t_type = ct[(ct['Tumor']>1) & (ct['Adj Normal']<=1) & (ct['Fetal']<=1)].index.values.tolist()

    print(t_type)
    n_type = ct[(ct['Tumor']<=1) & (ct['Adj Normal']>1) & (ct['Fetal']<=1)].index.values.tolist()

    print(n_type)
    tn_type = ct[(ct['Tumor']>1) & (ct['Adj Normal']>1) & (ct['Fetal']<=1)].index.values.tolist()

    print(tn_type)
    fetal_type = ct[(ct['Tumor']<=1) & (ct['Adj Normal']<=1) & (ct['Fetal']>1)].index.values.tolist()

    print(fetal_type)
    if len(t_type) != 0 & len(n_type) != 0 & len(tn_type) != 0 & len(fetal_type) != 0:
        

        
        
        tn = tmp[tmp.obs[batch].isin(['Tumor']) & \
             tmp.obs[anno_col].isin(t_type) | tmp.obs.batch.isin(['Adj Normal']) & tmp.obs[anno_col].isin(n_type) | tmp.obs[batch].isin(['Tumor', 'Adj Normal'])\
                 & tmp.obs[anno_col].isin(tn_type) ]

        sc.tl.rank_genes_groups(tn, anno_col, method='wilcoxon')
        pd.DataFrame(tn.uns['rank_genes_groups']['names']).head(top_n).to_csv(tn_markers)

        tn_df = pd.DataFrame(tn.uns['rank_genes_groups']['names']).head(top_n)
        
        fetal = tmp[tmp.obs[batch].isin(['Fetal']) & tmp.obs[anno_col].isin(fetal_type) ]

        sc.tl.rank_genes_groups(fetal, anno_col, method='wilcoxon')
        pd.DataFrame(fetal.uns['rank_genes_groups']['names']).head(top_n).to_csv(fetal_markers)

        fetal_df = pd.DataFrame(fetal.uns['rank_genes_groups']['names']).head(top_n)
        

    elif (len(t_type) != 0) & (len(n_type) != 0) & (len(tn_type) == 0) & (len(fetal_type) != 0):
        
        
        print('do this')
        tn = tmp[tmp.obs[batch].isin(['Tumor']) & \
             tmp.obs[anno_col].isin(t_type) | tmp.obs.batch.isin(['Adj Normal']) & tmp.obs[anno_col].isin(n_type) ]

        sc.tl.rank_genes_groups(tn, anno_col, method='wilcoxon')
        tn_df = pd.DataFrame(tn.uns['rank_genes_groups']['names']).head(top_n)
        pd.DataFrame(tn.uns['rank_genes_groups']['names']).head(top_n).to_csv(tn_markers)

        
        
        
        result = tn.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges', 'scores']}).head(top_n).to_csv(full_tn_markers)
        
        

        fetal = tmp[tmp.obs[batch].isin(['Fetal']) & tmp.obs[anno_col].isin(fetal_type) ]

        sc.tl.rank_genes_groups(fetal, anno_col, method='wilcoxon')
        pd.DataFrame(fetal.uns['rank_genes_groups']['names']).head(top_n).to_csv(fetal_markers)

        fetal_df = pd.DataFrame(fetal.uns['rank_genes_groups']['names']).head(top_n)


        fetal_result = fetal.uns['rank_genes_groups']
        groups = fetal_result['names'].dtype.names
        pd.DataFrame(
            {group + '_' + key[:1]: fetal_result[key][group]
            for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges', 'scores']}).head(top_n).to_csv(full_fetal_markers)

        
        
    
    elif not all( [len(t_type), len(n_type), len(fetal_type)] ):
        logg.error("Can't compare gene similarity")

    tn_df = pd.read_csv(tn_markers, index_col = 0)
    fetal_df = pd.read_csv(fetal_markers, index_col = 0)
    
    Overlap_genes = list(set(tn_df[tumor_query]).intersection(set(fetal_df[fetal_ref])))

    
    pd.DataFrame(Overlap_genes).to_csv('overlap_genes.csv')
    
    
    
    os.system("Rscript {} -w {} -t {} -f {} -m {} -n {} -q {} -r {}".format(genesimi_script, os.path.abspath(genesimi_dir), tn_markers, fetal_markers, full_tn_markers, full_fetal_markers, tumor_query, fetal_ref) )

    
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def findgrn(parameters_path):
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg
    
    
    
    cp = configparser.ConfigParser()
    cp.read(parameters_path)
    
    grnsimi_dir = cp.get('workdir', 'grnsimi_dir')
    if not os.path.exists(grnsimi_dir):
        os.makedirs(grnsimi_dir)
    
    os.chdir(grnsimi_dir)
    
    script_dir = cp.get('workdir', 'script_dir')
    grnsimi_script = script_dir + '/' + cp.get('grnsimilarity', 'grnsimi_script')
    
    tumor_query = cp.get('genesimilarity', 'tumor_query')
    fetal_ref = cp.get('genesimilarity', 'fetal_ref')

    run_SCENIC = cp.get('grnsimilarity', 'run_SCENIC')
    
    n_cores = eval(cp.get('grnsimilarity', 'n_cores'))
    
    minCountsPerGene_rate = eval(cp.get('grnsimilarity', 'minCountsPerGene_rate'))
    
    minSamples_rate = eval(cp.get('grnsimilarity', 'minSamples_rate'))
    
    input_file = cp.get('workdir', 'identification_dir') + '/' + cp.get('grnsimilarity', 'input')
    
    width = cp.get('grnsimilarity', 'width')
    height = cp.get('grnsimilarity', 'height')
    
    
    os.system("Rscript {} -w {} -i {} -l {} -n {} -c {} -s {} -q {} -r {} -y {} -h {}"\
              .format(grnsimi_script, grnsimi_dir, input_file, run_SCENIC, n_cores, minCountsPerGene_rate, minSamples_rate, tumor_query, fetal_ref, width, height) )
    
    
    
    
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def findCCC(parameters_path):
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg
    
    
    
    cp = configparser.ConfigParser()
    cp.read(parameters_path)
    
    ccc_dir = cp.get('workdir', 'ccc_dir')
    if not os.path.exists(ccc_dir):
        os.makedirs(ccc_dir)
    
    os.chdir(ccc_dir)
    
    
    tumor_data = cp.get('CCC', 'tumor_data')
    normal_data = cp.get('CCC', 'normal_data')
    fetal_data = cp.get('CCC', 'fetal_data')
    
    tumor_type = cp.get('CCC', 'tumor_type')
    normal_type = cp.get('CCC', 'normal_type')
    fetal_type = cp.get('CCC', 'fetal_type')
    
    CCCscript = cp.get('workdir', 'script_dir') + '/' + cp.get('CCC', 'CCCscript')
    
    group_by = cp.get('CCC', 'group_by')
    
    cellchat_type = cp.get('CCC', 'cellchat_type')
    
    trim = eval(cp.get('CCC', 'trim'))
    run_Cellchat = cp.get('CCC', 'run_Cellchat')

    
    os.system("Rscript {} -w {} -t {} -n {} -f {} -g {} -m {} -r {} -q {} -p {} -e {} -l {}"\
              .format(CCCscript, ccc_dir, tumor_data, normal_data, fetal_data, group_by, cellchat_type, trim, tumor_type, normal_type, fetal_type, run_Cellchat) )

    
    
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def drive(parameters_path):
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg
    
    
    
    cp = configparser.ConfigParser()
    cp.read(parameters_path)
    
    ccc_dir = cp.get('workdir', 'ccc_dir')
    if not os.path.exists(ccc_dir):
        os.makedirs(ccc_dir)
    
    os.chdir(ccc_dir)
    
    data = cp.get('workdir', 'data_path') + '/' + cp.get('drive', 'data')
    
    sender = cp.get('drive', 'sender')
    receiver = cp.get('drive', 'receiver')
    
    ligand_target_matrix = cp.get('workdir', 'data_path') + '/' + cp.get('drive', 'ligand_target_matrix')
    lr_network = cp.get('workdir', 'data_path') + '/' + cp.get('drive', 'lr_network')
    weighted_networks = cp.get('workdir', 'data_path') + '/' + cp.get('drive', 'weighted_networks')
    
    label = cp.get('drive', 'label')
    
    lfc_cutoff = eval(cp.get('drive', 'lfc_cutoff'))
    condition_colname = cp.get('drive', 'condition_colname')
    
    condition_oi = cp.get('drive', 'condition_oi')
    condition_reference = cp.get('drive', 'condition_reference')
    
    organism = cp.get('drive', 'organism')
    
    width = cp.get('drive', 'width')
    height = eval(cp.get('drive', 'height'))
    drive_script = cp.get('workdir', 'script_dir') + '/' + cp.get('drive', 'drive_script')
    
    os.system("Rscript {} -w {} -d {} -s {} -r {} -l {} -n {} -c {} -e {} -o {} -i {} -f {} -m {} -u {} -h {} -b {}"\
              .format(drive_script, ccc_dir, data, sender, receiver, ligand_target_matrix, lr_network, lfc_cutoff, weighted_networks, condition_colname \
                      , condition_oi, condition_reference, organism, width, height, label) )


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def trajectory(parameters_path):
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg
    
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    
    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.set_figure_params('scvelo',dpi=150,dpi_save=300, format='pdf')
    
    cp = configparser.ConfigParser()
    cp.read(parameters_path)
    
    trajectory_dir = cp.get('workdir', 'trajectory_dir')
    if not os.path.exists(trajectory_dir):
        os.makedirs(trajectory_dir)
    
    os.chdir(trajectory_dir)
    
    
    s = scv.read(cp.get('trajectory', 'spliced_path'))    
    u = scv.read(cp.get('trajectory', 'unspliced_path'))    
    
    adata = s
    adata.layers['spliced'] = s.X
    adata.layers['unspliced'] = u.X
    
    m = pd.read_csv(cp.get('trajectory', 'meta_path'), index_col=0)
    adata.obs = m
    
    basis = cp.get('trajectory', 'basis')
    
    umap = pd.read_csv(cp.get('trajectory', 'embedding_path'), index_col=0)
    adata.obsm['X_' + basis] = umap.values
    
    adata.write(cp.get('Extract', 'celltype_interest') + '_orig.h5ad', compression='gzip')    
    
    min_shared_counts = eval(cp.get('trajectory', 'min_shared_counts'))
    n_top_genes = eval(cp.get('trajectory', 'n_top_genes'))
    n_neighbors = eval(cp.get('trajectory', 'n_neighbors'))
    n_pcs = eval(cp.get('trajectory', 'n_pcs'))
    mode = cp.get('trajectory', 'mode')
    color = cp.get('trajectory', 'color')
    alpha = eval(cp.get('trajectory', 'alpha'))
    
    title = cp.get('trajectory', 'title')
    
    min_corr = eval(cp.get('trajectory', 'min_corr'))
    
    save = cp.get('trajectory', 'save')
    
    scv.pp.filter_and_normalize(adata, min_shared_counts=min_shared_counts, n_top_genes=n_top_genes)
    scv.pp.moments(adata, n_neighbors=n_neighbors, n_pcs = n_pcs)
    scv.tl.velocity(adata,mode=mode)
    print(adata.var.velocity_genes.sum())
    scv.tl.velocity_graph(adata)
    
    scv.pl.velocity_embedding_stream(adata, basis=basis,color=color, legend_loc='right margin', alpha=alpha
                                 , title=title, save=save, show=False)
    
    scv.tl.rank_velocity_genes(adata, groupby=color, min_corr=min_corr)

    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    df.to_csv('scvelo_diff_genes.csv')
    
    adata.write(cp.get('Extract', 'celltype_interest') + '_final.h5ad', compression='gzip') 
    
    


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--parameters_path', '-p', required=True, type=click.Path(), help='parameters file')
def clinical(parameters_path):
    
    import os
    import numpy as np
    import pandas as pd

    import scanpy as sc
    import scvelo as scv
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import stats

    from matplotlib.backends.backend_pdf import PdfPages
    import logging as logg
    
    
    
    cp = configparser.ConfigParser()
    cp.read(parameters_path)
    
    clinical_dir = cp.get('workdir', 'clinical_dir')
    if not os.path.exists(clinical_dir):
        os.makedirs(clinical_dir)
    
    os.chdir(clinical_dir)
    
    clinical_script = cp.get('workdir', 'script_dir') + '/' + cp.get('Clinical', 'clinical_script')
    
    oncofetal = cp.get('workdir', 'genesimi_dir') + '/' + cp.get('Clinical', 'oncofetal')
    
    data_rna = cp.get('Clinical', 'data_rna')
    data_surv = cp.get('Clinical', 'data_surv')
    
    data_tag = cp.get('Clinical', 'data_tag')
    
    os.system("Rscript {} -w {} -g {} -r {} -c {} -t {}"\
              .format(clinical_script, clinical_dir, oncofetal, data_rna, data_surv, data_tag) )
    
    
    
    
    
    
    
    

cli.add_command(extract_cells)
cli.add_command(scanpy_workflow)
cli.add_command(CSS)
cli.add_command(findOncoFetal)
cli.add_command(findgene)
cli.add_command(findgrn)
cli.add_command(findCCC)
cli.add_command(trajectory)
cli.add_command(clinical)
cli.add_command(drive)


if __name__ == '__main__':

    cli()


