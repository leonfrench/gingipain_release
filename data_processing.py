import glob
from pathlib import Path
import os
import pandas as pd
from pandas.api.types import is_numeric_dtype
from scipy.stats import zscore


adult_data_path = Path('./data/raw/allen_HBA')
gene_universe_path = Path('./data/gene_lists/gene_universe.txt')


def read_expression_file(file_name):
    expression_df = pd.read_csv(file_name, index_col=0, header=None)
    expression_df.index.rename('probe_id', inplace=True)
    return expression_df


def read_samples_file(samples_file):
    sample_df = pd.read_csv(samples_file)
    sample_df.set_index(sample_df.index + 1, inplace=True)
    sample_df.index.rename('sample_id', inplace=True)
    return sample_df


def filter_gene_universe(probes_df): #, gene_list):
    # drop probes with no annotated gene symbol
    probes_to_drop = probes_df[(probes_df.gene_symbol.str.startswith('A_')) | (probes_df.gene_symbol.str.startswith('CUST_'))]
    probes_df = probes_df[~probes_df.index.isin(probes_to_drop.index)]
    
    updated_annotations = pd.read_csv('./data/probe_annotations/updated_gene_symbol_annotations.csv')
    probes_df = probes_df.drop(['gene_symbol'], axis=1).merge(updated_annotations, on='probe_name')
    
    return probes_df


def get_probes_data(probes_file, probes_strategy='default'):
    strats = ['default', 'reannotator', 'qc_filter', 'qc_scale']
    assert probes_strategy in strats

    # depending on strategy, may merge in diff tables to update probes info
    probes_df = pd.read_csv(probes_file)

    # rename columns for consistency between adult and fetal brain datasets
    if 'probeset_name' in probes_df.columns:
        probes_df.rename(columns={'probeset_name': 'probe_name',
                                  'probeset_id': 'probe_id'}, inplace=True)
    cols = ['probe_id', 'probe_name', 'gene_symbol']
    print(f'probes_file: {probes_file}')
    probes_df = probes_df.loc[:, cols]
    print(f'cols: {probes_df.columns}')

    if probes_strategy == 'reannotator':
        reannotations = get_probe_reannotations(
            './data/raw/gene_symbol_annotations/AllenInstitute_custom_Agilent_Array.txt')
        # drop the original gene_symbol column
        probes_df.drop('gene_symbol', axis=1, inplace=True)
        # merge in the reannotated gene_symbols
        probes_df = probes_df.merge(reannotations, on='probe_name')

    elif probes_strategy in ['qc_filter', 'qc_scale']:
        # update gene_symbols using same method as was used in Rscript spatial transcriptomic analysis
        qc_filt = get_probe_qc_filter(
            './data/probe_annotations/Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx')
        probes_df = probes_df.merge(
            qc_filt, left_on='probe_name', right_on='probe')
        probes_df = probes_df[probes_df.qc_filter == True]
        assert is_numeric_dtype(probes_df.m)
        assert is_numeric_dtype(probes_df.b)

        print('After getting probes_df which merged qc data, shape is {}'.format(
            probes_df.shape))

    
    #gene_universe = pd.read_table(gene_universe_path, header=None)
    #gene_universe = gene_universe.iloc[:, 0]
    print(f'probes shape: {probes_df.shape}')
    probes_df = filter_gene_universe(probes_df)
    probes_df.set_index('probe_id', inplace=True)
    print(f'probes shape after filter: {probes_df.shape}')
    print('-----')
    #updated_annotations = pd.read_csv('./data/probe_annotations/updated_gene_symbol_annotations.csv')
    #probes_df.drop(['gene_symbol'], axis=1).merge(updated_annotations, on='probe_name')
    
    return probes_df


def get_probe_reannotations(re_annotations_file):
    # pre-processing function to prepare the probe reannotations file to be merge with probes df
    re_annotations = pd.read_table(re_annotations_file,
                                   usecols=['#PROBE_ID', 'Gene_symbol'])
    re_annotations.rename(columns={'#PROBE_ID': 'probe_name'}, inplace=True)
    re_annotations.dropna(inplace=True)
    re_annotations.set_index('probe_name', inplace=True)
    # split gene_symbols which have multiple genes associated with a single
    # probe_name creates a new row for each of the gene_symbols
    re_annotations = (re_annotations.Gene_symbol.str.split(';', expand=True)
                                    .stack()
                                    .reset_index()
                                    .drop('level_1', axis=1)
                                    .rename(columns={0: 'gene_symbol'}))

    return re_annotations


def get_probe_qc_filter(qc_file):
    # pre-processes the Allen qc filter file
    df = pd.read_excel(qc_file)
    df = df.iloc[:, [0, 1, 2, 3, 7]].drop(0)  # .rename(columns=[col_names])
    col_names = ['probe', 'gene_symbol', 'm', 'b', 'qc_filter']
    df.columns = col_names
    df[['m', 'b']] = df[['m', 'b']].apply(pd.to_numeric)
    return df.drop('gene_symbol', axis=1)


def scale_expression_vals(exp_df, probes_qc):
    print('Expression df shape before scaling: {}'.format(exp_df.shape))
    globScale = exp_df.values.flatten().mean()
    exp_df = exp_df - globScale

    df = exp_df.merge(probes_qc.loc[:, ['m', 'b']],
                      left_index=True, right_index=True)
    df.m = df.m.astype(float)
    df.b = df.b.astype(float)
    try:
        df = df.apply(lambda x: df.m * x + df.b)
    except Exception as e:
        print('Exception occurred on row when applying scale fxn, skip row')
        print(e)
    df = df.drop(['m', 'b'], axis=1)

    df.dropna(inplace=True)
    print('Expression df shape after scaling: {}'.format(df.shape))
    return df + globScale


def get_donor_data(donor_file_list, probes_strategy):
    # to work with both fetal and adult metadata
    probe_file_strings = ['Probes', 'rows_meta']
    samples_file_strings = ['Sample', 'columns_meta']
    expression_file_strings = ['Expression', 'expression']

    for file in donor_file_list:
        if any(string in file.stem for string in probe_file_strings):
            probes_df = get_probes_data(file, probes_strategy=probes_strategy)
        elif any(string in file.stem for string in samples_file_strings):
            samples_df = read_samples_file(file)
        elif any(string in file.stem for string in expression_file_strings):
            exp_df = read_expression_file(file)
        else:
            continue

    return exp_df, samples_df, probes_df


def get_mean_expression_by_brain_area(exp_df, samples_df):
    assert(exp_df.T.shape[0] == samples_df.shape[0])

    # merge in metadata (brain area of sample)
    annotated_df = exp_df.T.merge(samples_df[['structure_name']],
                                  left_index=True, right_index=True)

    # get mean expression level for samples within a brain area
    expression_by_structure = annotated_df.groupby('structure_name').mean()
    expression_by_structure.T.index.rename('gene_symbol', inplace=True)

    return expression_by_structure.T


def get_exp_by_genes(exp_df, probes_df):
    """
    input is exp_df and probes metadata
    output is exp_df grouped by gene_symbols and averaged
    """
    annotated_df = exp_df.merge(probes_df[['gene_symbol']],
                                left_index=True, right_index=True)

    exp_by_genes = (annotated_df.groupby('gene_symbol')
                                .mean())

    return exp_by_genes


def strip_left_right(structure_name):
    brain_area_fragments = structure_name.split(',')
    clean_fragments = []

    for frag in brain_area_fragments:
        if frag.strip() not in ['left', 'right', 'Left', 'Right']:
            clean_fragments.append(frag)

    clean_structure_name = ','.join(clean_fragments)
    return clean_structure_name


def get_single_donor_tidy_df(exp_df, samples_df, probes_df, donor_id, probes_strategy):
    # remove left/right from brain structure_names
    samples_df.structure_name = samples_df.structure_name.apply(
        strip_left_right)

    if probes_strategy == 'qc_scale':
        exp_df = scale_expression_vals(exp_df, probes_df)
        print('size of df after scaling: {}'.format(exp_df.shape))

    # merge in probes metadata and get expression by gene_symbol
    expression_by_genes = get_exp_by_genes(exp_df, probes_df)
    print('size of df after merging probe info, grouping by gene: {}'.format(
        expression_by_genes.shape))
    
    # merge in sample metadata and aggregate expression by brain area
    exp_brain_area_by_genes = get_mean_expression_by_brain_area(expression_by_genes, samples_df)
    

    ranked_exp_by_area = exp_brain_area_by_genes.rank(ascending=True)
    zscored_exp_by_area = pd.DataFrame(zscore(ranked_exp_by_area, axis=1), 
                                       index=ranked_exp_by_area.index,
                                       columns= ranked_exp_by_area.columns)
    #zscored_exp_by_area = ranked_exp_by_area.apply(zscore, axis=1, result_type='broadcast')

    melted = pd.melt(zscored_exp_by_area.reset_index(),
                     id_vars='gene_symbol',
                     var_name='brain_area')
    melted['donor_id'] = donor_id

    return melted


def generate_HBA_dataset(dataset, probes_strategy, selected_donor_ids=['10021', '9861', '14380', '15697', '15496', '12876']):
    assert dataset in ['adult', 'fetal']
    if dataset == 'adult':
        print('--- Generating adult human HBA dataset ---')
        #data_path = Path('./data/raw/allen_HBA')
        data_path = adult_data_path
    elif dataset == 'fetal':
        print('--- Generating fetal human HBA dataset ---')
        #data_path = Path('./data/raw/allen_human_fetal_brain')
        data_path = fetal_data_path


    selected_donor_paths = []
    for donor_id in selected_donor_ids:
        selected_donor_paths.append(list(data_path.glob(f'normalized_microarray*{donor_id}*')))

    all_donors = []
    for i, donor_folder in enumerate(selected_donor_paths):
        print(f'Processing donor #{i+1}')
        #donor_id = donor_folder.split('/')[-1].split('_')[-1]
        #donor_folder = list(donor_folder)
        donor_id = donor_folder[0].stem.split('_')[-1]
        print(f'Donor ID: {donor_id}')
        #donor_files = glob.glob(os.path.join(donor_folder, '*'))
        donor_files = list(donor_folder[0].glob('*.csv'))
        expression, samples, probes = get_donor_data(donor_files,
                                                     probes_strategy)

        tidy_donor = get_single_donor_tidy_df(expression,
                                              samples,
                                              probes,
                                              donor_id=donor_id,
                                              probes_strategy=probes_strategy)
        all_donors.append(tidy_donor)

    all_donors_long = pd.concat(all_donors)

    structure_genes_exp_matrix = pd.pivot_table(all_donors_long,
                                                values='value',
                                                index='gene_symbol',
                                                columns='brain_area')
    return structure_genes_exp_matrix


def get_dataset(dataset, probes_strategy='default', selected_donor_ids=['10021', '9861', '14380', '15697', '15496', '12876']):
    # to reduce code duplication between get_HBA_dataset and
    #  get_fetal_HBA_dataset
    strats = ['reannotator', 'qc_filter', 'qc_scale', 'default']
    datasets = ['adult', 'fetal']
    assert probes_strategy in strats
    assert dataset in datasets

    donor_str = '-'.join(selected_donor_ids)
    filename = '{0}_brainarea_vs_genes_exp_{1}_donors_{2}.tsv'.format(dataset,
                                                                      probes_strategy, donor_str)
    HBA_data_out_path = os.path.join('data', 'processed_HBA', filename)

    if os.path.exists(HBA_data_out_path):
        print(('Processed HBA brain dataset found locally. ' +
               'Loading from {}'.format(HBA_data_out_path)))
        structure_genes_exp_matrix = pd.read_csv(HBA_data_out_path,
                                                 index_col='gene_symbol',
                                                 sep='\t')

    else:
        os.makedirs('./data/processed_HBA', exist_ok=True)
        structure_genes_exp_matrix = generate_HBA_dataset(dataset,
                                                          probes_strategy,
                                                          selected_donor_ids=selected_donor_ids)
        print('-- Writing data to ' + HBA_data_out_path + ' -- ')
        structure_genes_exp_matrix.to_csv(HBA_data_out_path, sep='\t')
    return structure_genes_exp_matrix


def generate_aggregate_data(probes_strategy='default', 
                             donor_ids=['9861','10021', '12876', '14380', '15496', '15697'], 
                             min_donors=1):
    """
    """
    data_path = adult_data_path
    hba_donor_folders = [donor_dir for donor_dir in data_path.iterdir() if (donor_dir.is_dir() & ('microarray' in donor_dir.stem))]
    
    
    selected_donor_paths = []
    for path in hba_donor_folders:
        if any(donor_id in path.stem for donor_id in donor_ids):
            selected_donor_paths.append(path)
    
    all_donors = []
    for i, donor_folder in enumerate(selected_donor_paths):
        print(f'Processing donor #{i+1}')
        print(f'Donor directory: {donor_folder}')
        donor_id = donor_folder.stem
        print(f'Donor ID: {donor_id}')
        donor_files = list(donor_folder.glob('*.csv')) 
        
        expression, samples, probes = get_donor_data(donor_files, probes_strategy=probes_strategy)
        donor_gene_exp = get_exp_by_genes(expression, probes)
        samples.structure_name = samples.structure_name.apply(strip_left_right)
        # this takes mean of gene expression from different brain structures
        # an alternative way to try doing it is by concatenating all samples (from all donors) 
        # and average the final values
        #exp_by_brain_area = merge_samples_w_exp(donor_gene_exp, samples)
        exp_by_brain_area = get_mean_expression_by_brain_area(donor_gene_exp, samples)
        #exp_by_brain_area['donor_id'] = donor_id
        donor_expression  = exp_by_brain_area.T
        donor_expression['donor_id'] = donor_id
        
        all_donors.append(donor_expression)
        
    concat_data = pd.concat(all_donors)
    n_donors_per_structure = concat_data.groupby('structure_name').donor_id.nunique().reset_index()
    selected_structures = n_donors_per_structure[n_donors_per_structure.donor_id >= min_donors].structure_name
    # the index of concatenated data is the sampled structure_name
    selected_samples = concat_data[concat_data.index.isin(selected_structures)]
    final_matrix = selected_samples.groupby(['donor_id', 'structure_name']).mean().groupby('structure_name').mean()
    
    return final_matrix.T
