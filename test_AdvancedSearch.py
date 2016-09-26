from unittest import TestCase
from PDBsearch.AdvancedSearch import *


class TestAdvancedSearch(TestCase):

    def test_all_experimental_type_molecule_type(self):
        a = AdvancedSearch()
        a.all_experimental_type_molecule_type(ExperimentalMethod.Electron_Microscopy, MoleculeType.Nucleic_Acid)
        a.send()
        self.assertIn('1C2W', a.result)

    def test_pdb_id_s(self):
        self.fail()

    def test_entity_id_s(self):
        self.fail()

    def test_chain_id_s(self):
        self.fail()

    def test_pubmed_id_s(self):
        self.fail()

    def test_uniprotkb_accession_number_s(self):
        self.fail()

    def test_text_search(self):
        self.fail()

    def test_mmcif_keyword_search_classification(self):
        self.fail()

    def test_pfam_accession_number_s(self):
        self.fail()

    def test_uniprot_gene_name_s(self):
        self.fail()

    def test_sequence_cluster_name_s(self):
        self.fail()

    def test_structure_title(self):
        self.fail()

    def test_structure_description(self):
        self.fail()

    def test_macromolecule_name_s(self):
        self.fail()

    def test_large_structures(self):
        self.fail()

    def test_author_name(self):
        self.fail()

    def test_deposit_date(self):
        self.fail()

    def test_release_date(self):
        self.fail()

    def test_revise_date(self):
        self.fail()

    def test_latest_released_date(self):
        self.fail()

    def test_latest_modifies_structures(self):
        self.fail()

    def test_structural_genomics_project(self):
        self.fail()

    def test_macromolecule_type(self):
        self.fail()

    def test_number_of_chains_au(self):
        self.fail()

    def test_number_of_chains_ba(self):
        self.fail()

    def test_number_of_entities(self):
        self.fail()

    def test_protein_stoichiometry(self):
        self.fail()

    def test_protein_symmetry(self):
        self.fail()

    def test_number_of_models(self):
        self.fail()

    def test_number_of_disulfide_bonds(self):
        self.fail()

    def test_link_records(self):
        self.fail()

    def test_molecular_weight_struct(self):
        self.fail()

    def test_secondary_structure_content(self):
        self.fail()

    def test_secondary_structure_length(self):
        self.fail()

    def test_sequence_blast_fasta_psi_blast(self):
        self.fail()

    def test_wild_type_protein(self):
        self.fail()

    def test_translated_nucleotide_sequence(self):
        self.fail()

    def test_sequence_motif(self):
        self.fail()

    def test_chain_length(self):
        self.fail()

    def test_chemical_name(self):
        self.fail()

    def test_chemical_id_s(self):
        self.fail()

    def test_inchi_descriptor(self):
        self.fail()

    def test_chemical_structure_smiles(self):
        self.fail()

    def test_molecular_weight_chem_comp(self):
        self.fail()

    def test_chemical_formula(self):
        self.fail()

    def test_chemical_component_type(self):
        self.fail()

    def test_binding_affinity(self):
        self.fail()

    def test_has_ligand_s(self):
        self.fail()

    def test_has_modified_residue_s(self):
        self.fail()

    def test_sub_component_s(self):
        self.fail()

    def test_biologically_interesting_molecules(self):
        self.fail()

    def test_expression_organism(self):
        self.fail()

    def test_enzyme_classification(self):
        self.fail()

    def test_experimental_method(self):
        self.fail()

    def test_x_ray_resolution(self):
        self.fail()

    def test_x_ray_average_b_factor(self):
        self.fail()

    def test_refinement_r_factors(self):
        self.fail()

    def test_diffraction_source(self):
        self.fail()

    def test_structure_determination_method(self):
        self.fail()

    def test_reflections(self):
        self.fail()

    def test_cell_dimensions(self):
        self.fail()

    def test_software(self):
        self.fail()

    def test_space_group(self):
        self.fail()

    def test_crystal_properties(self):
        self.fail()

    def test_em_assembly(self):
        self.fail()

    def test_detector(self):
        self.fail()

    def test_citation(self):
        self.fail()

    def test_pubmed_abstract(self):
        self.fail()

    def test_has_external_links(self):
        self.fail()

    def test_search_by_mixed_criteria(self):
        """
        On 11.06.2016 given criteria gave result of 4 structures ['1JBN', '1PAR', '2DNJ', '2K2N'] so test expects them
        as result.
        """
        as_prop = Bio.PDB.AdvancedSearch.AdvancedSearch()
        as_prop.structure_title(Bio.PDB.AdvancedSearch.StructComparator.Starts_with, 'DNA')
        as_prop.number_of_chains_ba('5', '6')
        as_prop.send()

        print(as_prop.search)
        self.assertEqual(4, len(as_prop.search))