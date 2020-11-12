import math
import os
import json
import xml.etree.ElementTree as ET
import requests
from bs4 import BeautifulSoup
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import rgb2hex
import numpy as np
import data_processing as data
import analysis as hba


def get_ontology(json_file='./data/ontology.json', atlas_id=265297125, graph_id=10):
    if os.path.exists(json_file):
        with open(json_file) as data_file:
            ontology = json.load(data_file)
    else:
        query_url = "http://api.brain-map.org/api/v2/structure_graph_download/{}.json".format(graph_id)
        r = requests.get(query_url)
        response = r.json()
        ontology = response['msg'][0]
        # dump out the ontology to a file
        with open(json_file, 'w') as outfile:
            json.dump(ontology, outfile)

    return ontology


def find_structure(self, attr, value):
    # function to find structure in ontology by attribute
    if self[attr] == value:
        return self
    else:
        for child in self['children']:
            match = find_structure(child, attr, value)
            if match:
                return match


def get_sIDs_in_SVG(svg_filename, ontology):
    # function to parse SVG and return a dict with sID as key and value is structure_name
    ns = {'svg':'http://www.w3.org/2000/svg'}
    structures_names = {}

    with open(svg_filename, 'r') as svg:
            doc = ET.parse(svg)
            root = doc.getroot()

            for elem in root.iter():
                # if this is path element and it has a structure_id
                # add structure_id to structures dictionary if it doesn't already exist
                if elem.tag == '{%s}%s' % (ns['svg'], 'path'):
                    if 'structure_id' in elem.attrib:
                        s_id = int(elem.attrib['structure_id'])
                        # is there a need for structures_names dict? would a list of sIDs suffice?
                        # could remove the need for ontology in the function if simply append s_id to a list
                        if s_id not in structures_names:
                            srec = find_structure(ontology, 'id', int(s_id))
                            structures_names[s_id] = srec['name']
    return structures_names


def find_children_sID(self, attr, value):
    # function to find structure IDs of the children of a given structure
    children_sID = []
    reduced_onto = find_structure(self, attr, value)
    for child in reduced_onto['children']:
        children_sID.append(child['id'])
    return children_sID


def roll_up(sID_list, ontology, HBA_lookup, AUC_vals=None):
    if AUC_vals is None:
        AUC_vals = []
    for sID in sID_list:
        # if the structure with sID is found in HBA_lookup
        if HBA_lookup[HBA_lookup.id == float(sID)].shape[0] > 0:
            value = float(HBA_lookup[HBA_lookup.id == float(sID)].loc[:,'AUROC'])
            AUC_vals.append(value)
        elif len(find_children_sID(ontology, 'id', sID)) > 0:
            # the sID was not found in the results table
            # get a list of sIDs of the children of the structure of interest
            children_sIDs = find_children_sID(ontology, 'id', sID)
            roll_up(children_sIDs, ontology, HBA_lookup, AUC_vals)
            #print(f'sID ({sID}) not found in table, rolling up values from children structures {children_sIDs}')
        else:
            # the structure was not found in HBA_lookup, nor were any of its children
            # get the values of parent structure

            # get parent sID
            structure_info = find_structure(ontology, 'id', sID)
            parent_sID = structure_info['parent_structure_id']
            #print(f'sID: {sID} not found in lookup, getting values from parent: {parent_sID}')
            # get AUC val for that parent sID
            try:
                value = float(HBA_lookup[HBA_lookup.id == float(parent_sID)].loc[:,'AUROC'])
                AUC_vals.append(value)
            except TypeError:
                pass
                # print('No children to roll-up AUC values or parent to inherit AUC val from for sID: {}'.format(sID))
    return AUC_vals


def get_AUC_vals_for_sIDs(sIDs, ontology, HBA_lookup):
    AUC_vals = {}

    for sID in sIDs:
        children_auc_vals = roll_up([sID], ontology, HBA_lookup)
        AUC_vals[sID] = np.mean(children_auc_vals)

    return AUC_vals


def convert_AUCvals_to_hex(sID_AUC_map, cmap='RdBu_r', generate_cbar=False): #'bwr'
    centered_vals = np.abs(np.asarray(list(sID_AUC_map.values())) - 0.5)
    #drop nans
    centered_vals = centered_vals[~np.isnan(centered_vals)]
    extreme_val = max(centered_vals)
    #extreme_val = max(np.abs(np.asarray(list(sID_AUC_map.values())) - 0.5))

    norm = mpl.colors.Normalize(vmin=0.5 - extreme_val , vmax=0.5 + extreme_val)
    color_map = plt.cm.get_cmap('{}'.format(cmap))
    sID_hex_map = {}
    for key, AUC_value in sID_AUC_map.items():
        #print('AUC: {} and normed: {}'.format(AUC_value, norm(AUC_value)))
        #normalized_value = norm(AUC_value)
        if math.isnan(AUC_value):
            colour = '#FFFFFF'
        else:
            colour = color_map(norm(AUC_value))

        sID_hex_map[key] = rgb2hex(colour)

    if generate_cbar:
        fig, ax = plt.subplots(figsize=(1, 8))
        cbar = mpl.colorbar.ColorbarBase(ax, cmap=color_map,norm=norm)
                                 #orientation='horizontal')
        cbar.set_label('AUC', size=12)

        return sID_hex_map, fig

    return sID_hex_map


def modify_structure_color(input_svg_file, sID_AUC_map, output_file, cbar_file=None):
    # function to modify the colour of a structure in the SVG
    with open(output_file, 'w') as outfile, open(input_svg_file, 'r') as infile:
        # parse the svg file
        soup = BeautifulSoup(infile, 'xml')

        if cbar_file:
            sID_hex_map, cbar = convert_AUCvals_to_hex(sID_AUC_map, generate_cbar=True)
            #fig.savefig(fname='./figures/colorbar.svg', format="svg")
            print(f'saving colourbar to: {cbar_file}')
            cbar.savefig(fname=cbar_file)

        else:
            sID_hex_map = convert_AUCvals_to_hex(sID_AUC_map)

        # searches svg file for descendants with the structure_id of interest
        for sID in sID_hex_map:

            for desc in soup.descendants:
                try:
                    attributes = desc.attrs
                    try:
                        # need a special rule for white matter of forebrain
                        if int(attributes['structure_id']) == 9219: # 9219 = sID of telencephalic whitematter
                            attributes['style'] = 'stroke:black;fill:{}'.format('#FFFFFF') #FFC0CB
                        if int(attributes['structure_id']) == sID:
                            #print("found structure of interest: ", sID)
                            #print('structure_id: ', desc['structure_id'])
                            #print(attributes['style'])
                            attributes['style'] = 'stroke:black;fill:{}'.format(sID_hex_map[sID])
                            #print('colouring in with: ', attributes['style'])
                            #print()
                        else:
                            continue
                    except (KeyError, AttributeError):
                        #print('no attributes for this descendant:', desc)
                        continue
                except AttributeError:
                    #print('no stucture attributes for this descendant:', attributes)
                    continue
        print(f'Writing modified svg to {output_file}')
        outfile.write(str(soup))
        #return outfile


def create_auc_lookup(exp_df, gene_list, ontology):
    if ontology == 'adult':
        # modify the default ontology to remove left/right from structure names
        ontology_df = pd.read_csv('./data/raw/allen_HBA/normalized_microarray_donor10021/Ontology.csv')
        # name column to be the same as fetal ontology
        ontology_df['structure_name'] = ontology_df.name.apply(data.strip_left_right)

    elif ontology == 'fetal':
        # not all structures are contained in any single ontology.csv for fetal brains
        # concatenate all and drop duplicated rows
        fetal_ontology1 = pd.read_csv('./data/raw/allen_human_fetal_brain/lmd_matrix_12566/columns_metadata.csv')
        fetal_ontology2 = pd.read_csv('./data/raw/allen_human_fetal_brain/lmd_matrix_12690/columns_metadata.csv')
        fetal_ontology3 = pd.read_csv('./data/raw/allen_human_fetal_brain/lmd_matrix_12840/columns_metadata.csv')
        fetal_ontology4 = pd.read_csv('./data/raw/allen_human_fetal_brain/lmd_matrix_14751/columns_metadata.csv')
        ontology_df = pd.concat([fetal_ontology1, fetal_ontology2, fetal_ontology3, fetal_ontology4])
        ontology_df.drop_duplicates(subset=['structure_id'], inplace=True)
        # rename id column to be same as with adult ontology
        ontology_df.rename(columns={'structure_id': 'id'}, inplace=True)

    # create dicts that map from sID to structure names and inverse
    # this isn't used?
    sID_to_brainstructure = ontology_df.set_index('id').loc[:, 'structure_name'].to_dict()
    brainstructure_to_sID = {}
    for k, v in sID_to_brainstructure.items():
        brainstructure_to_sID[v] = brainstructure_to_sID.get(v, [])
        brainstructure_to_sID[v].append(k)

    # map sIDs to results table
    results = hba.generate_stats_table(exp_df, gene_list)
    results = results.reset_index()
    results = results.rename(columns={'index': 'structure'})
    results['id'] = results.structure.map(brainstructure_to_sID)

    # expand rows where there are multiple sIDs for each structure, then merge results
    HBA_lookup = results.set_index('structure').id.apply(pd.Series).stack().reset_index(level=1, drop=True)
    #HBA_lookup = HBA_lookup.reset_index(name='id').merge(ontology_df[['id', 'parent_structure_id']], on='id')
    HBA_lookup = HBA_lookup.reset_index(name='id').merge(ontology_df[['id']], on='id')
    HBA_lookup = HBA_lookup.merge(results.drop('id', axis=1), on='structure')

    return HBA_lookup


def clean_SVG(svg_filename, output_filename):
    with open(svg_filename, 'r') as svg, open(output_filename, 'w') as outfile:
        soup = BeautifulSoup(svg, 'xml')
        to_remove = soup.find_all(attrs={"graphic_group_label": ["Atlas - Human Sulci", "Atlas - Human Hotspots"]})
        for r in to_remove:
            r.decompose()
        
        print(f'Writing cleaned svg to {output_filename}')
        outfile.write(str(soup))


def modify_svg(svg_file, svg_output, graph_id, lookup_table, cbar_file=None):
    if graph_id == 'adult':
        ontology = get_ontology(json_file='./data/ontology.json', atlas_id=265297125, graph_id=10)
    elif graph_id == 'fetal':
        ontology = get_ontology(json_file='./data/fetal21ontology.json', atlas_id=3, graph_id=16)
    elif graph_id == 'fetal_brainstem':
        ontology = get_ontology(json_file='./data/fetal_brainstem_ontology.json', atlas_id=287730656, graph_id=16)

    svg_sIDs = get_sIDs_in_SVG(svg_file, ontology)

    sID_auc_map = get_AUC_vals_for_sIDs(svg_sIDs.keys(), ontology, lookup_table)

    svg_output = modify_structure_color(svg_file, sID_auc_map, svg_output, cbar_file=cbar_file)

    #return svg_output
    

def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result
