import copy
import os
import logging
import shutil
import string
from Bio.PDB import PDBParser, Selection

from libs import arcimboldo_air, bioutils, features, hhsearch, utils
from typing import Dict, List
import pandas as pd


class Template:

    def __init__ (self, parameters_list: List, output_dir: str, input_dir: str, num_of_copies: int):

        self.pdb_path: str
        self.pdb_id: str
        self.template_path: str
        self.chain: str = ''
        self.polyala_res_dict: Dict = None
        self.add_to_msa: bool = False
        self.add_to_templates: bool = False
        self.sum_prob: bool = False
        self.aligned: bool = False
        self.template_features: Dict = None
        self.merged_templates_list: List = None
        
        self.pdb_path = self.check_pdb(parameters_list['pdb'], input_dir)
        self.pdb_id = utils.get_file_name(self.pdb_path)
        self.add_to_msa = parameters_list.get('add_to_msa', self.add_to_msa)
        self.add_to_templates = parameters_list.get('add_to_templates', self.add_to_templates)
        self.sum_prob = parameters_list.get('sum_prob', self.sum_prob)
        self.aligned = parameters_list.get('aligned', self.aligned)
        self.template_path = f'{output_dir}/{self.pdb_id}_template.pdb'
        self.chain = parameters_list.get('chain', self.chain)

        if num_of_copies == 1:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(self.pdb_id, self.pdb_path)
            chains = [chain.get_id() for chain in structure.get_chains()]
            if self.chain == '':
                if len(chains) == 1:
                    self.chain = chains.pop()
                else:
                    raise Exception('There is more than one chain in the structure. Select one in the configuration file.')
            else:
                if not self.chain in chains:
                    raise Exception('The choosen chain does not exist in the structure.')
    
        polyala_res = parameters_list.get('polyala_res', self.polyala_res_dict)
        if polyala_res is not None:
            pdb_out_path = os.path.join(output_dir, self.pdb_id +"_polyala.pdb")
            self.polyala_res_dict = bioutils.convert_template_to_polyala(pdb_in_path=self.pdb_path, pdb_out_path=pdb_out_path, polyala_res=polyala_res)
            self.pdb_path = pdb_out_path

    def check_pdb(self, pdb: str, output_dir: str) -> str:

        if not os.path.exists(pdb):
            bioutils.download_pdb(pdb_id=pdb, output_dir=output_dir)
            pdb = f'{output_dir}/{pdb}.pdb'
        else:
            pdb_aux = f'{output_dir}/{os.path.basename(pdb)}'
            if pdb != pdb_aux:
                shutil.copy2(pdb, pdb_aux)
                pdb = pdb_aux                

        return pdb

    def choose_best_offset(self, target, merged_list: List, output_dir: str) -> Dict:
        
        self.merged_list = merged_list
        if self is target:
            return

        results_pdist = []
        results_pdist.append(['reference'] + [utils.get_file_name(file) for file in target.merged_list])

        results_superpose = []
        results_superpose.append(['reference'] + [utils.get_file_name(file) for file in target.merged_list])

        results_algorithm = []

        for x, query_pdb in enumerate(merged_list):
            reference_pdist_list = []
            reference_superpose_list = []
            reference_algorithm = []
            res_list_length = len([res for res in Selection.unfold_entities(PDBParser().get_structure('test', query_pdb), 'R')])

            for y, target_pdb in enumerate(target.merged_list):
                rmsd, nalign, quality_q, aligned_res_list = bioutils.superpose_pdbs(query_pdb=query_pdb,
                                                                       target_pdb=target_pdb, output_superposition = False)
                reference_superpose_list.append(f'{rmsd}, {nalign} ({res_list_length})')
                reference_pdist_list.append(bioutils.pdist(query_pdb=query_pdb, target_pdb=target_pdb))
                reference_algorithm.append((x, y, bioutils.pdist(query_pdb=query_pdb, target_pdb=target_pdb)))
            
            results_pdist.append([utils.get_file_name(query_pdb)] + reference_pdist_list)
            results_superpose.append([utils.get_file_name(query_pdb)] + reference_superpose_list)
            results_algorithm.append(reference_algorithm)

        best_offset_list = bioutils.calculate_best_offset(results_algorithm)
        return_offset_dict = {}
        for x,y,_ in best_offset_list:
            return_offset_dict[merged_list[x]] = y

        output_txt = f'{output_dir}/best_offset.txt'

        with open(output_txt, 'w') as f_in:

            f_in.write('Calculated with distances:\n')

            rows = []
            for x in results_pdist[1:]:
                rows.append(x)
            df = pd.DataFrame(rows, columns=results_pdist[0])
            f_in.write(df.to_markdown())        
        
            f_in.write('\n\n')

            f_in.write('Calculated with superpose:\n')

            rows = []
            for x in results_superpose[1:]:
                rows.append(x)
            df = pd.DataFrame(rows, columns=results_superpose[0])
            f_in.write(df.to_markdown())

            f_in.write('\n\n')
            f_in.write('Best combination is (using pdist):')
            f_in.write(f'{str(return_offset_dict)}')
        
        return return_offset_dict
    
    def generate_features(self, a_air):

        output_hhr = f'{a_air.run_dir}/output.hhr'
        pdb70_path = f'{a_air.run_dir}/pdb70'

        if not self.aligned:
            bioutils.pdb2mmcif(output_dir=a_air.run_dir, pdb_in_path=self.pdb_path, cif_out_path=f'{a_air.run_dir}/{self.pdb_id}.cif')
            hhsearch.generate_hhsearch_db(template_cif_path=f'{a_air.run_dir}/{self.pdb_id}.cif', output_dir=a_air.run_dir)
            hhsearch.run_hhsearch(fasta_path=a_air.fasta_path, pdb70_db=pdb70_path,output_path=output_hhr)

            if a_air.num_of_copies > 1:
                chain_list, transformations_list = bioutils.read_remark_350(pdb_path=f'{a_air.run_dir}/{self.pdb_id}.cif')
            else:
                chain_list = [self.chain]
                transformations_list = [[
                    ['1.000000', '0.000000', '0.000000'], 
                    ['0.000000', '1.000000', '0.000000'], 
                    ['0.000000', '0.000000', '1.000000'], 
                    ['0.00000', '0.00000', '0.00000']
                ]]

            new_chain_list = list(string.ascii_uppercase)[:len(transformations_list) * len(chain_list)]

            logging.info('Assembly can be build using chain(s) '+ str(chain_list) + ' by applying the following transformations:')
            for matrix in transformations_list:
                logging.info(str(matrix))

            if len(new_chain_list) != a_air.num_of_copies:
                raise Exception(f'Assembly description from REMARK 350 contains {len(new_chain_list)} subunits. Please, try to '
                    f'generate a new REMARK 350 (manually or using e.g. PISA) for considering a new assembly with '
                    f'{a_air.num_of_copies} subunits')

            query_seq_length = len(a_air.query_sequence)
            list_of_paths_of_pdbs_trans = []
            list_of_paths_of_pdbs_to_merge = []

            for chain in chain_list:
                template_features = features.extract_template_features_from_pdb(
                    query_sequence=a_air.query_sequence,
                    hhr_path=output_hhr,
                    pdb_id=self.pdb_id,
                    chain_id=chain,
                    mmcif_db=a_air.run_dir)

                g = features.Features(query_sequence=a_air.query_sequence)
                g.append_new_template_features(new_template_features=template_features, custom_sum_prob=self.sum_prob)
                templates_dict = g.write_all_templates_in_features(output_dir=a_air.run_dir)
                template_path = list(templates_dict.values())[0]
                pdb_path = f'{a_air.run_dir}/{self.pdb_id}_{chain}'

                for i, transformation in enumerate(transformations_list):
                    bioutils.change_chain(pdb_in_path=template_path, 
                                pdb_out_path=f'{pdb_path}_{i+1}.pdb',
                                rot_tra_matrix=transformation)
                    list_of_paths_of_pdbs_trans.append(f'{pdb_path}_{i+1}.pdb')

            offset_dict = self.choose_best_offset(a_air.templates[0], list_of_paths_of_pdbs_trans, a_air.output_dir)
            for i, pdb_path in enumerate(list_of_paths_of_pdbs_trans):
                if bool(offset_dict):
                    offset_value = offset_dict[pdb_path]
                else:
                    offset_value = i
                offset = query_seq_length * (offset_value) + a_air.glycines * (offset_value)
                new_pdb_path = f'{a_air.run_dir}/{self.pdb_id}_{new_chain_list[offset_value]}.pdb'
                bioutils.change_chain(pdb_in_path=pdb_path,
                                pdb_out_path=new_pdb_path,
                                offset=offset, chain='A')
                list_of_paths_of_pdbs_to_merge.append(new_pdb_path)
            bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=sorted(list_of_paths_of_pdbs_to_merge),
                    merged_pdb_path=self.template_path)

        else:
            shutil.copy2(self.pdb_path, self.template_path)

        template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
            query_sequence=a_air.query_sequence_assembled,
            pdb_path=self.template_path,
            pdb_id=self.pdb_id,
            chain_id='A')
            
        self.template_features = copy.deepcopy(template_features)
        
    def __repr__(self):
        
        return f' \
        pdb_path: {self.pdb_path} \n \
        pdb_id: {self.pdb_id} \n \
        template_path: {self.template_path} \n \
        chain: {self.chain} \n \
        polyala_res_dict: {self.polyala_res_dict} \n \
        add_to_msa: {self.add_to_msa} \n \
        add_to_templates: {self.add_to_templates} \n \
        sum_prob: {self.sum_prob} \n \
        aligned: {self.aligned} \n \
        template_features: {self.template_features}'