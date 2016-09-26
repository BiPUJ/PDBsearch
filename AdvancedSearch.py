# PDBsearch, version 0.1
#
# Copyright 2016 by Wiktoria Karwicka & Jacek Smietanski
# (http://jaceksmietanski.net)
#
# Freely available at https://github.com/BiPUJ/PDBsearch.git
#
# This code is governed by biopython license. Please see the LICENSE file
# that should have been included as part of this package


""" Advanced search on RCSB Protein Data Bank server. """

from xml.etree.ElementTree import Element, SubElement, tostring
from enum import Enum
import requests

# Enum classes for search purposes
# universal - used in more then one function
class YesOrNo(Enum):
    """ To build query parameters in function wild_type_protein() """
    Yes = 'Y'
    No = 'N'


class YesNo(Enum):
    """ To build query parameters in functions: sequence_blast_fasta_psi_blast(), has_ligand_s(),
    has_modified_residue_s()
    """
    Yes = 'yes'
    No = 'no'


class CountComparator(Enum):
    """ To build query parameters in functions: number_of_models(), citation() """
    Equals = '='
    More = '>'
    Less = '<'
    Different = '!='


class ComparatorChoice(Enum):
    """ To build query parameters in functions: diffraction_source(), citation() """
    Contains = 'contains'
    Starts_with = 'startswith'
    Ends_with = 'endswith'


class StructComparator(Enum):
    """ To build query parameters in functions: mmcif_keyword_search_classification(), structure_title(),
    structure_description(), expression_organism(), structure_determination_method(), reflections(), software()
    """
    Contains = 'contains'
    Equals = 'equals'
    Starts_with = 'startswith'
    Ends_with = 'endswith'
    Does_not_contain = '!contains'
    Does_not_start_with = '!startswith'
    Does_not_end_with = '!endswith'


class PolymericType(Enum):
    """ To build query parameters in functions: chemical_name(), chemical_id_s(), chemical_structure_smiles()
    """
    Any = 'Any'
    Free = 'Free'
    Polymeric = 'Polymeric'


# Quick Search
class ExperimentalMethod(Enum):
    """ To build query parameters in function all_experimental_type_molecule_type() """
    X_ray = 'X-RAY'
    NMR = 'NMR'
    Electron_Microscopy = 'ELECTRON MICROSCOPY'
    Hybrid = 'HYBRID'
    Other = 'other'
    Ignore = 'ignore'


class MoleculeType(Enum):
    """ To build query parameters in function all_experimental_type_molecule_type() """
    Protein = 'protein'
    Nucleic_Acid = 'nucleic'
    Protein_NA_Complex = 'complex'
    Contains_Nucleic_Acid = 'nucleicAcidContaining'
    Other = 'other'
    Ignore = 'ignore'


# ID(s) and Keywords
class SequenceStructureName(Enum):
    """ To build query parameters in function sequence_cluster_name_s() """
    Contains = 'contains'
    Equals = 'equals'


# Structure Annotation
class LargeStructures(Enum):
    """ To build query parameters in function large_structures() """
    Yes = 'yes'
    Omit_large_structures = 'Omit Large Structures'


# Deposition
class AuthorName(Enum):
    """ To build query parameters in function author_name() """
    All_authors = 'All Authors'
    Structure_author = 'Structure Author'
    Citation_author = 'Citation Author'


class TrueOrFalse(Enum):
    """ To build query parameters in function author_name() """
    true = 'true'
    false = 'false'


class ProjectName(Enum):
    """ To build query parameters in function structural_genomics_project() """
    All = 'AllSGProjectName'
    Not_assigned = 'NaSGProjectName'
    Enzyme_Function_Initiative = 'Enzyme Function Initiative'
    NIAID = 'NIAID, National Institute of Allergy and Infectious Diseases'
    NPPSFA = 'NPPSFA, National Project on Protein Structural and Functional Analyses'
    PSI = 'PSI, Protein Structure Initiative'
    PSI_biology = 'PSI:Biology'


class CenterInitials(Enum):
    """ To build query parameters in function structural_genomics_project() """
    AllSGProjectInitials = 'AllSGProjectInitials'
    NaSGProjectInitials = 'NaSGProjectInitials'
    ATCG3D = 'ATCG3D'
    BIGS = 'BIGS'
    BSGC = 'BSGC'
    BSGI = 'BSGI'
    CEBS = 'CEBS'
    CELLMAT = 'CELLMAT'
    CESG = 'CESG'
    CHTSB = 'CHTSB'
    CSGID = 'CSGID'
    CSMP = 'CSMP'
    GPCR = 'GPCR'
    IFN = 'IFN'
    ISFI = 'ISFI'
    ISPC = 'ISPC'
    JCSG = 'JCSG'
    MCSG = 'MCSG'
    MPP = 'MPP'
    MPSBC = 'MPSBC'
    MPSbyNMR = 'MPSbyNMR'
    MSGP = 'MSGP'
    MSGPP = 'MSGPP'
    MTBI = 'MTBI'
    NESG = 'NESG'
    NHRS = 'NHRS'
    NHRs = 'NHRs'
    NPCXstals = 'NPCXstals'
    NYCOMPS = 'NYCOMPS'
    NYSGRC = 'NYSGRC'
    NYSGXRC = 'NYSGXRC'
    NatPro = 'NatPro'
    OCSP = 'OCSP'
    OPPF = 'OPPF'
    PCSEP = 'PCSEP'
    PSF = 'PSF'
    RSGI = 'RSGI'
    S2F = 'S2F'
    SECSG = 'SECSG'
    SGC = 'SGC'
    SGCGES = 'SGCGES'
    SGPP = 'SGPP'
    SPINE = 'SPINE'
    SPINE_2 = 'SPINE-2'
    SSGCID = 'SSGCID'
    SSPF = 'SSPF'
    STEMCELL = 'STEMCELL'
    TBSGC = 'TBSGC'
    TCELL = 'TCELL'
    TEMIMPS = 'TEMIMPS'
    UC4CDI = 'UC4CDI'
    XMTB = 'XMTB'
    YSG = 'YSG'


class CenterName(Enum):
    """ To build query parameters in function structural_genomics_project() """
    All = 'AllSGProjects'
    Not_assigned = 'NaSGProjects'
    Accelerated_technologies_center = 'Accelerated Technologies Center for Gene to 3D Structure'
    Assembly_dynamics_and_evolution = 'Assembly, Dynamics and Evolution of Cell-Cell and Cell-Matrix Adhesions'
    Atoms_to_animals = 'Atoms-to-Animals: The Immune Function Network'
    Bacterial_targets_at_IGS_CNRS = 'Bacterial targets at IGS-CNRS, France'
    Berkeley_structural_genomics_center = 'Berkeley Structural Genomics Center'
    Center_for_eukaryotic_structural_genomics = 'Center for Eukaryotic Structural Genomics'
    Center_for_high_throughput_structural_biology = 'Center for High-Throughput Structural Biology'
    Center_for_structural_genomics_of_ID = 'Center for Structural Genomics of Infectious Diseases'
    Center_for_structural_of_membrane_proteins = 'Center for Structures of Membrane Proteins'
    Chaperone_enabled_studies = 'Chaperone-Enabled Studies of Epigenetic Regulation Enzymes'
    Enzyme_discovery = 'Enzyme Discovery for Natural Product Biosynthesis'
    GPCR_network = 'GPCR Network'
    Integrated_center = 'Integrated Center for Structure and Function Innovation'
    Israel_structural_proteomics_center = 'Israel Structural Proteomics Center'
    Joint_center = 'Joint Center for Structural Genomics'
    Marseilles_structural_genomics_program = 'Marseilles Structural Genomics Program @ AFMB'
    Medical_structural_genomisc = 'Medical Structural Genomics of Pathogenic Protozoa'
    Membrane_protein_structural_biology_consortium = 'Membrane Protein Structural Biology Consortium'
    Membrane_protein_structural_by_solution = 'Membrane Protein Structures by Solution NMR'
    Midwest_center_for_structural_genomics = 'Midwest Center for Structural Genomics'
    Mitochondrial_protein_partnership = 'Mitochondrial Protein Partnership'
    Montreal_kingston_bacterial_structural_genomics_initiative = \
        'Montreal-Kingston Bacterial Structural Genomics Initiative'
    Mycobacterium_tuberculosis_structural_proteomics_project = \
        'Mycobacterium Tuberculosis Structural Proteomics Project'
    New_york_consortium_on_membrane_protein_structure = 'New York Consortium on Membrane Protein Structure'
    New_york_SGX_research_center_for_structural_genomics = 'New York SGX Research Center for Structural Genomics'
    New_york_structural_genomics_research_consortium = 'New York Structural Genomics Research Consortium'
    Northeast_structural_genomics_consortium = 'Northeast Structural Genomics Consortium'
    Nucleocytoplasmic_transport_a_target_for_cellular_control = \
        'Nucleocytoplasmic Transport: a Target for Cellular Control'
    Ontario_centre_for_structural_proteomics = 'Ontario Centre for Structural Proteomics'
    Oxford_protein_production_facility = 'Oxford Protein Production Facility'
    Paris_sud_ueast_structural_genomics = 'Paris-Sud Yeast Structural Genomics'
    Partnership_for_nuclear_receptor_signaling_code_biology = \
        'Partnership for Nuclear Receptor Signaling Code Biology'
    Partnership_for_stem_cell_biology = 'Partnership for Stem Cell Biology'
    Partnership_for_T_cell_biology = 'Partnership for T-Cell Biology'
    Program_for_the_characterization_of_secreted_effector_proteins = \
        'Program for the Characterization of Secreted Effector Proteins'
    Protein_structure_factory = 'Protein Structure Factory'
    RIKEN_structural_genomics_proteomics_initiative = 'RIKEN Structural Genomics/Proteomics Initiative'
    Scottish_structural_proteomics_facility = 'Scottish Structural Proteomics Facility'
    Seattle_structural_genomics_center_for_infectious_disease = \
        'Seattle Structural Genomics Center for Infectious Disease'
    Southeast_collaboratory_for_structural_genomics = 'Southeast Collaboratory for Structural Genomics'
    Structural_genomics_consortium = 'Structural Genomics Consortium'
    Structural_genomics_consortium_for_research_on_gene_expression = \
        'Structural Genomics Consortium for Research on Gene Expression'
    Structural_genomics_of_pathogenic_protozoa_consortium = 'Structural Genomics of Pathogenic Protozoa Consortium'
    Structural_proteomics_in_Europe = 'Structural Proteomics in Europe'
    Structural_proteomics_in_Europe_2 = 'Structural Proteomics in Europe 2'
    Structure_2_function_project = 'Structure 2 Function Project'
    Structure_function_analysis = 'Structure-Function Analysis of Polymorphic CDI Toxin-Immunity Protein Complexes'
    Structures_of_Mtb_proteins = 'Structures of Mtb Proteins Conferring Susceptibility to Known Mtb Inhibitors'
    TB_structural_genomics_consortium = 'TB Structural Genomics Consortium'
    Transcontinental_EM_initiative_for_membrane_protein_structure = \
        'Transcontinental EM Initiative for Membrane Protein Structure'


# Structure Features
class YesNoIgnore(Enum):
    """ To build query parameters in function macromolecule_type() """
    Yes = 'Y'
    No = 'N'
    Ignore = '?'


class SearchTool(Enum):
    """ To build query parameters in function sequence_blast_fasta_psi_blast() """
    BLAST = 'blast'
    FASTA = 'fasta'
    PSI_BLAST = 'psiblast'


class EntityType(Enum):
    """ To build query parameters in function number_of_entities() """
    Any = ' '
    Protein = 'p'
    RNA = 'r'
    DNA = 'd'
    Ligand = 'n'
    Other = '?'


class ProteinSymmetry(Enum):
    """ To build query parameters in function protein_symmetry() """
    Asymmetric = 'A1'
    Cyclic_C2 = 'C2'
    Cyclic_C3 = 'C3'
    Cyclic_C4 = 'C4'
    Cyclic_C5 = 'C5'
    Cyclic_C6 = 'C6'
    Cyclic_C7 = 'C7'
    Cyclic_C8 = 'C8'
    Cyclic_C9 = 'C9'
    Cyclic_C10 = 'C10'
    Cyclic_C11 = 'C11'
    Cyclic_C12 = 'C12'
    Cyclic_C13 = 'C13'
    Cyclic_C14 = 'C14'
    Cyclic_C17 = 'C17'
    Cyclic_C22 = 'C22'
    Cyclic_C24 = 'C24'
    Cyclic_C31 = 'C31'
    Cyclic_C38 = 'C38'
    Cyclic_C39 = 'C39'
    Dihedral_D2 = 'D2'
    Dihedral_D3 = 'D3'
    Dihedral_D4 = 'D4'
    Dihedral_D5 = 'D5'
    Dihedral_D6 = 'D6'
    Dihedral_D7 = 'D7'
    Dihedral_D8 = 'D8'
    Dihedral_D9 = 'D9'
    Dihedral_D11 = 'D11'
    Dihedral_D12 = 'D12'
    Dihedral_D17 = 'D17'
    Dihedral_D48 = 'D48'
    Tetrahedral = 'T'
    Octahedral = 'O'
    Icosahedral = 'I'
    Helical = 'H'


class LinkConnectionType(Enum):
    """ To build query parameters in function link_records() """
    Any = 'any'
    Covalent_bond = 'covalent_bond'
    Disulfide_bridge = 'disulfide_bridge'
    Hydrogen_bond = 'hydrogen_bond'
    Metal_coordination = 'metal_coordination'
    Mismatches_base_pairs = 'mismatched_base_pairs'
    Ionic_interaction_saltbridge = 'ionic_interaction_saltbridge'


class SecondaryStructureMethod(Enum):
    """ To build query parameters in function secondary_structure_length() """
    DSSP = 'DSSP'
    Stride = 'Stride'


class HelixSubtype(Enum):
    """ To build query parameters in function secondary_structure_length() """
    Three_to_ten_helix = '3-10 helix'
    Alpha_helix = 'alpha helix'
    Pi_helix = 'pi helix'


class BetaStrandSubtype(Enum):
    """ To build query parameters in function secondary_structure_length() """
    Beta_bridge = 'beta bridge'
    Beta_strand = 'beta_strand'


# Sequence Features
class PercentCoverageOfWTPSequence(Enum):
    """ To build query parameters in function wild_type_protein() """
    Hundred = '100'
    Ninety_five = '95'
    Ninety = '90'
    Eighty_five = '85'
    Eighty = '80'
    Seventy_five = '75'
    Seventy = '70'
    Sixty_five = '65'
    Sixty = '60'


# Chemical Components
class ChemicalNameComparator(Enum):
    """ To build query parameters in function chemical_name() """
    Contains = 'contains'
    Equals = 'equals'
    Sounds_like = 'sounds_like'


class InChIDescriptorType(Enum):
    """ To build query parameters in inchi_descriptor() """
    InChI = 'InChI'
    InChI_key = 'InChI_key'


class StructureSearchType(Enum):
    """ To build query parameters in chemical_structure_smiles() """
    Substructure = 'Substructure'
    Exact = 'Exact'
    Similar = 'Similar'
    Superstructure = 'Superstructure'


class ChemicalComponentType(Enum):
    """ To build query parameters in chemical_component_type() """
    non_polymer = 'non-polymer'
    peptide_linking = 'peptide linking'
    L_peptide_linking = 'L-peptide linking'
    D_peptide_linking = 'D-peptide linking'
    L_peptide_NH3_amino_terminus = 'L-peptide NH3 amino terminus'
    D_peptide_NH3_amino_terminus = 'D-peptide NH3 amino terminus'
    L_peptide_NH3_COOH_carboxy_terminus = 'L-peptide COOH carboxy terminus'
    D_peptide_NH3_COOH_carboxy_terminus = 'D-peptide COOH carboxy terminus'
    L_gamma_peptide_C_delta_linking = 'L-gamma-peptide, C-delta linking'
    D_gamma_peptide_C_delta_linking = 'D-gamma-peptide, C-delta linking'
    L_beta_peptide_C_gamma_linking = 'L-beta-peptide, C-gamma linking'
    D_beta_peptide_C_gamma_linking = 'D-beta-peptide, C-gamma linking'
    peptide_like = 'peptide-like'
    DNA_linking = 'DNA linking'
    L_DNA_linking = 'L-DNA linking'
    DNA_OH_5_prime_terminus = 'DNA OH 5 prime terminus'
    DNA_OH_3_prime_terminus = 'DNA OH 3 prime terminus'
    RNA_linking = 'RNA linking'
    L_RNA_linking = 'L-RNA linking'
    RNA_OH_5_prime_terminus = 'RNA OH 5 prime terminus'
    RNA_OH_3_prime_terminus = 'RNA OH 5 prime terminus'
    Saccharide = 'saccharide'
    L_saccharide = 'L-saccharide'
    D_saccharide = 'D-saccharide'
    L_saccharide_1_4_and_1_4_linking = 'L-saccharide 1,4 and 1,4 linking'
    D_saccharide_1_4_and_1_4_linking = 'D-saccharide 1,4 and 1,4 linking'
    L_saccharide_1_4_and_1_6_linking = 'L-saccharide 1,4 and 1,6 linking'
    D_saccharide_1_4_and_1_6_linking = 'D-saccharide 1,4 and 1,6 linking'
    Other = 'other'


class AffinityType(Enum):
    """ To build query parameters in binding_affinity() """
    Ki = 'Ki'
    Kd = 'Kd'
    EC50 = 'EC50'
    IC50 = 'IC50'
    deltaG = '&Delta;G'
    deltaH = '&Delta;H'
    minusTdeltaS = '-T&Delta;S'
    Ka = 'Ka'


# Biologically Interesting Molecules (from BIRD)
class BIRDNameComparator(Enum):
    """ To build query parameters in biologically_interesting_molecules() """
    Contains = 'contains'
    Equals = 'equals'
    Sounds_like = 'sounds_like'
    Starts_with = 'starts_with'


class BIRDType(Enum):
    """ To build query parameters in biologically_interesting_molecules() """
    Any = 'Any'
    Chalkophore = 'Chalkophore'
    Cyclic_depsipeptide = 'Cyclic Depsipeptide'
    Cyclic_lipopeptide = 'Cyclic Lipopeptide'
    Cyclic_peptide = 'Cyclic Peptide'
    Glycopeptide = 'Glycopeptide'
    Lipoglycopeptide = 'Lipoglycopeptide'
    Lipopeptide = 'Lipopeptide'
    Macrolide = 'Macrolide'
    Non_polymer = 'Non-polymer'
    Oligopeptide = 'Oligopeptide'
    Oligosaccharide = 'Oligosaccharide'
    Peptide_like = 'Peptide-like'
    Polypeptide = 'Polypeptide'
    Thiopeptide = 'Thiopeptide'


class BIRDClass(Enum):
    """ To build query parameters in biologically_interesting_molecules() """
    Any = 'Any'
    Antagonist = 'Antagonist'
    Antibiotic = 'Antibiotic'
    Anticancer = 'Anticancer'
    Anticoagulant = 'Anticoagulant'
    Antiinflammatory = 'Antiinflammatory'
    Antimicrobial = 'Antimicrobial'
    Antiretroviral = 'Antiretroviral'
    Antithombotic = 'Antithombotic'
    Antitumor = 'Antitumor'
    Antiviral = 'Antiviral'
    Capase_inhibitor = 'Capase Inhibitor'
    Enzyme_inhibitor = 'Enzyme Inhibitor'
    Immunosuppressant = 'Immunosuppressant'
    Inhibitor = 'Inhibitor'
    Lantibiotic = 'Lantibiotic'
    Metabolism = 'Metabolism'
    Metal_transport = 'Metal Transport'
    Thrombin_inhibitor = 'Thrombin Inhibitor'
    Toxin = 'Toxin'
    Transport_activator = 'Transport Activator'
    Trypsin_inhibitor = 'Trypsin Inhibitor'
    Unknown = 'Unknown'


# Methods
class ExperimentalMethods(Enum):
    """ To build query parameters in experimental_method() """
    X_RAY = 'X-RAY'
    Solution_NMR = 'SOLUTION NMR'
    Solid_state_NMR = 'SOLID-STATE NMR'
    Electron_microscopy = 'ELECTRON MICROSCOPY'
    Electron_crystallography = 'ELECTRON CRYSTALLOGRAPHY'
    Fiber_diffraction = 'FIBER DIFFRACTION'
    Neutron_diffraction = 'NEUTRON DIFFRACTION'
    Solution_scattering = 'SOLUTION SCATTERING'
    Other = 'OTHER'
    Hybrid = 'HYBRID'


class YesNoOrIgnore(Enum):
    """ To build query parameters in function experimental_method() """
    Yes = 'Y'
    No = 'N'
    Ignore = 'Ignore'


class SpaceGroup(Enum):
    """ To build query parameters in function space_group() """
    A_1 = 'A 1'
    A_2 = 'A 2'
    B_1_1_2 = 'B 1 1 2'
    B_2 = 'B 2'
    B_2_21_2 = 'B 2 21 2'
    C_1_2_1 = 'C 1 2 1'
    C_1_2c_1 = 'C 1 2/c 1'
    C_1_21_1 = 'C 1 21 1'
    C_2_2_2 = 'C 2 2 2'
    C_2_2_21 = 'C 2 2 21'
    C_21 = 'C 21'
    C_4_21_2 = 'C 4 21 2'
    F_2_2_2 = 'F 2 2 2'
    F_2_3 = 'F 2 3'
    F_4_2_2 = 'F 4 2 2'
    F_4_3_2 = 'F 4 3 2'
    F_41_3_2 = 'F 41 3 2'
    H_minus3 = 'H -3'
    H_3 = 'H 3'
    H_3_2 = 'H 3 2'
    I_m4_2_d = 'I -4 2 d'
    I_m4_c_2 = 'I -4 c 2'
    I_1_2_1 = 'I 1 2 1'
    I_1_21_1 = 'I 1 21 1'
    I_2_2_2 = 'I 2 2 2'
    I_2_3 = 'I 2 3'
    I_21 = 'I 21'
    I_21_21_21 = 'I 21 21 21'
    I_21_3 = 'I 21 3'
    I_4 = 'I 4'
    I_4_2_2 = 'I 4 2 2'
    I_4_3_2 = 'I 4 3 2'
    I_41 = 'I 41'
    I_41_2_2 = 'I 41 2 2'
    I_41_3_2 = 'I 41 3 2'
    I_41a = 'I 41/a'
    P_m1 = 'P -1'
    P_m3 = 'P -3'
    P_1 = 'P 1'
    P_1_1_2 = 'P 1 1 2'
    P_1_1_21 = 'P 1 1 21'
    P_1_2_1 = 'P 1 2 1'
    P_1_21_1 = 'P 1 21 1'
    P_1_21c_1 = 'P 1 21/c 1'
    P_1_21n_1 = 'P 1 21/n 1'
    P_2_2_2 = 'P 2 2 2'
    P_2_2_21 = 'P 2 2 21'
    P_2_21_21 = 'P 2 21 21'
    P_2_3 = 'P 2 3'
    P_21_2_2 = 'P 21 2 2'
    P_21_2_21 = 'P 21 2 21'
    P_21_21_2 = 'P 21 21 2'
    P_21_21_2_A = 'P 21 21 2 A'
    P_21_3 = 'P 21 3'
    P_3 = 'P 3'
    P_3_1_2 = 'P 3 1 2'
    P_3_2_1 = 'P 3 2 1'
    P_31 = 'P 31'
    P_31_1_2 = 'P 31 1 2'
    P_31_2_1 = 'P 31 2 1'
    P_32 = 'P 32'
    P_32_1_2 = 'P 32 1 2'
    P_32_2_1 = 'P 32 2 1'
    P_4 = 'P 4'
    P_4_2_2 = 'P 4 2 2'
    P_4_21_2 = 'P 4 21 2'
    P_41 = 'P 41'
    P_41_2_2 = 'P 41 2 2'
    P_41_21_2 = 'P 41 21 2'
    P_41_3_2 = 'P 41 3 2'
    P_42 = 'P 42'
    P_42_2_2 = 'P 42 2 2'
    P_42_21_2 = 'P 42 21 2'
    P_42_3_2 = 'P 42 3 2'
    P_43 = 'P 43'
    P_43_2_2 = 'P 43 2 2'
    P_43_21_2 = 'P 43 21 2'
    P_43_3_2 = 'P 43 3 2'
    P_6 = 'P 6'
    P_6_2_2 = 'P 6 2 2'
    P_61 = 'P 61'
    P_61_2_2 = 'P 61 2 2'
    P_62 = 'P 62'
    P_62_2_2 = 'P 62 2 2'
    P_63 = 'P 63'
    P_63_2_2 = 'P 63 2 2'
    P_64 = 'P 64'
    P_64_2_2 = 'P 64 2 2'
    P_65 = 'P 65'
    P_65_2_2 = 'P 65 2 2'
    P_b_c_a = 'P b c a'
    R_3 = 'R 3'
    R_3_2 = 'R 3 2'


class EMAssembly(Enum):
    """ To build query parameters in function em_assembly() """
    TwoD_crystal = ''
    ThreeD_crystal = ''
    Cell = 'CELL'
    Filament = 'FILAMENT'
    Scosahedral = 'ICOSAHEDRAL'
    Particle = 'PARTICLE'
    Single_particle = 'SINGLE PARTICLE'
    Tissue = 'TISSUE'


class Detector(Enum):
    """ To build query parameters in function detector() """
    ADSC = 'ADSC'
    ADSC_QUANTUM_1 = 'ADSC QUANTUM 1'
    ADSC_QUANTUM_4 = 'ADSC QUANTUM 4'
    ADSC_QUANTUM_4r = 'ADSC QUANTUM 4r'
    ADSC_QUANTUM_210 = 'ADSC QUANTUM 210'
    ADSC_QUANTUM_270 = 'ADSC QUANTUM 270'
    ADSC_QUANTUM_315 = 'ADSC QUANTUM 315'
    ADSC_QUANTUM_210r = 'ADSC QUANTUM 210r'
    ADSC_QUANTUM_315r = 'ADSC QUANTUM 315r'
    APEX_II_CCD = 'APEX II CCD'
    APS = 'APS'
    BRANDEIS = 'BRANDEIS'
    BRUKER = 'BRUKER'
    Bruker_AXIOM_200 = 'Bruker AXIOM 200'
    Bruker_DIP_6040 = 'Bruker DIP-6040'
    Bruker_Platinum_135 = 'Bruker Platinum 135'
    BRUKER_SMART_2000 = 'BRUKER SMART 2000'
    BRUKER_SMART_6000 = 'BRUKER SMART 6000'
    BRUKER_SMART_6500 = 'BRUKER SMART 6500'
    CUSTOM_MADE = 'CUSTOM-MADE'
    DECTRIS_EIGER_R_1M = 'DECTRIS EIGER R 1M'
    DECTRIS_EIGER_R_4M = 'DECTRIS EIGER R 4M'
    DECTRIS_EIGER_X_1M = 'DECTRIS EIGER X 1M'
    DECTRIS_EIGER_X_4M = 'DECTRIS EIGER X 4M'
    DECTRIS_EIGER_X_9M = 'DECTRIS EIGER X 9M'
    DECTRIS_EIGER_X_16M = 'DECTRIS EIGER X 16M'
    DECTRIS_MYTHEN2_R_1K = 'DECTRIS MYTHEN2 R 1K'
    DECTRIS_MYTHEN2_R_1D = 'DECTRIS MYTHEN2 R 1D'
    DECTRIS_PILATUS_200K = 'DECTRIS PILATUS 200K'
    DECTRIS_PILATUS_300K = 'DECTRIS PILATUS 300K'
    DECTRIS_PILATUS_2M = 'DECTRIS PILATUS 2M'
    DECTRIS_PILATUS_2M_F = 'DECTRIS PILATUS 2M-F'
    DECTRIS_PILATUS_6M = 'DECTRIS PILATUS 6M'
    DECTRIS_PILATUS_6M_F = 'DECTRIS PILATUS 6M-F'
    DECTRIS_PILATUS_12M = 'DECTRIS PILATUS 12M'
    DECTRIS_PILATUS3_1M = 'DECTRIS PILATUS 1M'
    DECTRIS_PILATUS3_2M = 'DECTRIS PILATUS 2M'
    DECTRIS_PILATUS3_6M = 'DECTRIS PILATUS 6M'
    DECTRIS_PILATUS3_R_100K_A = 'DECTRIS PILATUS3 R 100K-A'
    DECTRIS_PILATUS3_X_100K_A = 'DECTRIS PILATUS3 X 100K-A'
    ENRAF_NONIUS = 'ENRAF-NONIUS'
    ENRAF_NONIUS_CAD4 = 'ENRAF-NONIUS CAD4'
    ENRAF_NONIUS_FAST = 'ENRAF-NONIUS FAST'
    FUJI = 'FUJI'
    HENDRIX_LENTFER = 'HENDRIX-LENTFER'
    KODAK = 'KODAK'
    MACSCIENCE = 'MACSCIENCE'
    MACSCIENCE_DIP100 = 'MACSCIENCE DIP100'
    MACSCIENCE_DIP100S = 'MACSCIENCE DIP100S'
    MAC_Science = 'MAC Science'
    MAC_Science_DIP_2000 = 'MAC Science DIP-2000'
    MAC_Science_DIP_2030 = 'MAC Science DIP-2030'
    MAC_Science_DIP_3000 = 'MAC Science DIP-3000'
    MAC_Science_DIP_320 = 'MAC Science DIP-320'
    MAR555_FLAT_PANEL = 'MAR555 FLAT PANEL'
    MAR_CCD_165_mm = 'MAR CCD 165 mm'
    MAR_CCD_130_mm = 'MAR CCD 130 mm'
    MARMOSAIC_225_mm_CCD = 'MARMOSAIC 225 mm CCD'
    MARMOSAIC_300_mm_CCD = 'MARMOSAIC 300 mm CCD'
    MARMOSAIC_325_mm_CCD = 'MARMOSAIC 325 mm CCD'
    MARRESEARCH = 'MARRESEARCH'
    MAR_scanner_180_mm_plate = 'MAR scanner 180 mm plate'
    MAR_scanner_300_mm_plate = 'MAR scanner 300 mm plate'
    MAR_scanner_345_mm_plate = 'MAR scanner 345 mm plate'
    NICOLET = 'NICOLET'
    NICOLET_P3 = 'NICOLET P3'
    NICOLET_P3X = 'NICOLET P3X'
    NOIR_1 = 'NOIR-1'
    NONIUS_CAD4 = 'NONIUS CAD4'
    Nonius_Kappa_CCD = 'Nonius Kappa CCD'
    Other = 'OTHER'
    OXFORD = 'OXFORD'
    OXFORD_ONYX_CCD = 'OXFORD ONYX CCD'
    OXFORD_RUBY_CCD = 'OXFORD RUBY CCD'
    OXFORD_SAPPHIRE_CCD = 'OXFORD SAPPHIRE CCD'
    OXFORD_TITAN_CCD = 'OXFORD TITAN CCD'
    PHILLIPS = 'PHILLIPS'
    PILATUS = 'PILATUS'
    PRINCETON_2K = 'PRINCETON 2K'
    PSI_PILATUS_6M = 'PSI PILATUS 6M'
    RAYONIX_MX_225 = 'RAYONIX MX-225'
    RAYONIX_MX_300 = 'RAYONIX MX-300'
    RAYONIX_MX_325 = 'RAYONIX MX-325'
    RAYONIX_MX225HE = 'RAYONIX MX225HE'
    RAYONIX_MX300HE = 'RAYONIX MX300HE'
    RAYONIX_MX325HE = 'RAYONIX MX325HE'
    RAYONIX_SX_165mm = 'RAYONIX SX-165mm'
    RIGAKU = 'RIGAKU'
    RIGAKU_AFC_5R = 'RIGAKU AFC-5R'
    RIGAKU_AFC_6R = 'RIGAKU AFC-6R'
    RIGAKU_AFC_6S = 'RIGAKU AFC-6S'
    RIGAKU_AFC11 = 'RIGAKU AFC11'
    RIGAKU_AFC11_KAPPA = 'RIGAKU AFC11-KAPPA'
    RIGAKU_AFC9 = 'RIGAKU AFC9'
    RIGAKU_JUPITER_140 = 'RIGAKU JUPITER 140'
    RIGAKU_JUPITER_210 = 'RIGAKU JUPITER 210'
    RIGAKU_MERCURY = 'RIGAKU MERCURY'
    RIGAKU_RAXIS = 'RIGAKU RAXIS'
    RIGAKU_RAXIS_RUH2R = 'RIGAKU RAXIS RUH2R'
    RIGAKU_RAXIS_HTC = 'RIGAKU RAXIS HTC'
    RIGAKU_RAXIS_II = 'RIGAKU RAXIS II'
    RIGAKU_RAXIS_IIC = 'RIGAKU RAXIS IIC'
    RIGAKU_RAXIS_IV = 'RIGAKU RAXIS IV'
    RIGAKU_RAXIS_IVpp = 'RIGAKU RAXIS IV++'
    RIGAKU_RAXIS_V = 'RIGAKU RAXIS V'
    RIGAKU_RAXIS_VII = 'RIGAKU RAXIS VII'
    RIGAKU_SATURN_70 = 'RIGAKU SATURN 70'
    RIGAKU_SATURN_724 = 'RIGAKU SATURN 724'
    RIGAKU_SATURN_92 = 'RIGAKU SATURN 92'
    RIGAKU_SATURN_944 = 'RIGAKU SATURN 944'
    RIGAKU_SATURN_944p = 'RIGAKU SATURN 944+'
    RIGAKU_SATURN_A200 = 'RIGAKU SATURN A200'
    SBC = 'SBC'
    SBC_2 = 'SBC-2'
    SBC_3 = 'SBC-3'
    SDMS = 'SDMS'
    SIEMENS = 'SIEMENS'
    SIEMENS_AED2 = 'SIEMENS AED2'
    SIEMENS_FOUR_CIRCLE = 'SIEMENS FOUR-CIRCLE'
    SIEMENS_HI_STAR = 'SIEMENS HI-STAR'
    SIEMENS_P4 = 'SIEMENS P4'
    SIEMENS_NICOLET = 'SIEMENS-NICOLET'
    SIEMENS_NICOLET_X100 = 'SIEMENS-NICOLET X100'
    STOE = 'STOE'
    STOE_SIEMENS_AED2 = 'STOE-SIEMENS AED2'
    SYNTEX = 'SYNTEX'
    UCSD_MARK_II = 'UCSD MARK II'
    UCSD_MARK_III = 'UCSD MARK III'
    WEISSENBERG = 'WEISSENBERG'
    XENTRONICS = 'XENTRONICS'
    XUONG_HAMLIN = 'XUONG-HAMLIN'


# Publication
class Journals(Enum):
    """ To build query parameters in function citation() """
    All = 'all'
    Acs_Chem_Biol = 'Acs Chem.Biol.'
    Acta_Crystallogr_Sect_D = 'Acta Crystallogr.,Sect.D'
    Acta_Crystallogr_Sect_F = 'Acta Crystallogr.,Sect.F'
    Angew_Chem_Int_Ed_Engl = 'Angew.Chem.Int.Ed.Engl.'
    Biochem_Biophys_Res_Commun = 'Biochem.Biophys.Res.Commun.'
    Biochem_J = 'Biochem.J.'
    Biochemistry = 'Biochemistry'
    Bioorg_Med_Chem_Lett = 'Bioorg.Med.Chem.Lett.'
    Cell = 'Cell'
    Chem_Biol = 'Chem.Biol.'
    EMBO_J = 'EMBO J.'
    FEBS_J = 'FEBS J.'
    FEBS_Lett = 'FEBS Lett.'
    J_Am_Chem_Soc = 'J.Am.Chem.Soc.'
    J_Biol_Chem = 'J.Biol.Chem.'
    J_Med_Chem = 'J.Med.Chem.'
    J_Mol_Biol = 'J.Mol.Biol.'
    J_Struct_Biol = 'J.Struct.Biol.'
    J_Virol = 'J.Virol.'
    Mol_Cell = 'Mol.Cell'
    Nat_Commun = 'Nat Commun'
    Nat_Struct_Mol_Biol = 'Nat.Struct.Mol.Biol.'
    Nature = 'Nature'
    Nucleic_Acids_Res = 'Nucleic Acids Res.'
    Plos_One = 'Plos One'
    Proc_Natl_Acad_Sci_USA = 'Proc.Natl.Acad.Sci.USA'
    Protein_Sci = 'Protein Sci.'
    Proteins = 'Proteins'
    Science = 'Science'
    Structure = 'Structure'
    Other = 'other'


# Misc
class ExternalLink(Enum):
    """ To build query parameters in function has_external_links() """
    EDS_IDS = 'EDS_IDS'
    RECOORD_PDB_ID = 'RECOORD_PDB_ID'
    OLDERADO_PDB_ID = 'OLDERADO_PDB_ID'
    HIVTOOLBOX_PDB_ID = 'HIVTOOLBOX_PDB_ID'
    SBGRID_DOI_ID = 'SBGRID_DOI_ID'


# A class with advanced search functions
class AdvancedSearch():
    # Variables for main query and refinement counter
    def __init__(self, url='http://www.rcsb.org/pdb/rest/search'):
        self.query = Element('orgPdbCompositeQuery', version="1.0")
        self.counter = 0  # counter for refinement level
        self.pdb_url = url  # remote pdb server
        self.result = ""  # search results

    # Quick Search
    def all_experimental_type_molecule_type(self, experimental_method, molecule_type):
        """Builds a query for Quick search purposes with specification All/Experimental Type/Molecule Type.

        @param experimental_method: experimental method
        @type experimental_method: ExperimentalMethod

        @param molecule_type: molecule type
        @type molecule_type: MoleculeType
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.HoldingsQuery'

        description = SubElement(title, 'description')
        description.text = 'Holdings Search : Molecule Type=' + molecule_type.value + '  Experimental Method=' \
                           + experimental_method.value + '  '

        experimental_method_query = SubElement(title, 'experimentalMethod')
        experimental_method_query.text = experimental_method.value

        molecule_type_query = SubElement(title, 'moleculeType')
        molecule_type_query.text = molecule_type.value

        self.query.append(query_refinement)
        print(self.query)
        self.counter += 1


    # ID(s) and Keywords
    def pdb_id_s(self, pdb_ids):
        """Builds a query for ID(s) and Keywords search purposes with specification PDB ID(s).

        @param pdb_ids: structure ID(s) from PDB (e.g. 4HHB, 2CPK)
        @type pdb_ids: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.StructureIdQuery'

        description = SubElement(title, 'description')
        description.text = 'Simple query for a list of PDB IDs: ' + str(pdb_ids)

        structure_id_list_query = SubElement(title, 'structureIdList')
        structure_id_list_query.text = str(pdb_ids)

        self.query.append(query_refinement)
        self.counter += 1

    def entity_id_s(self, entity_ids):
        """Builds a query for ID(s) and Keywords search purposes with specification Entity ID(s).

        @param entity_ids: entity ID(s) from PDB (e.g. 4HHB:1, 2CPK:1)
        @type entity_ids: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.EntityIdQuery'

        description = SubElement(title, 'description')
        description.text = 'Simple query for a list of PDB Entity IDs: ' + str(entity_ids)

        entity_id_list_query = SubElement(title, 'entityIdList')
        entity_id_list_query.text = str(entity_ids)

        self.query.append(query_refinement)
        self.counter += 1

    def chain_id_s(self, chain_ids):
        """Builds a query for ID(s) and Keywords search purposes with specification Chain ID(s).

        @param chain_ids: chain ID(s) from PDB (e.g. 4HHB.A, 2CPK.E); the chain IDs are case sensitive
        @type chain_ids: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.ChainIdQuery'

        description = SubElement(title, 'description')
        description.text = 'Simple query for a list of PDB Chain IDs: ' + str(chain_ids)

        entity_id_list_query = SubElement(title, 'entityIdList')
        entity_id_list_query.text = str(chain_ids)

        self.query.append(query_refinement)
        self.counter += 1

    def pubmed_id_s(self, pubmed_ids):
        """Builds a query for ID(s) and Keywords search purposes with specification PubMed ID(s).

        @param pubmed_ids: PubMed ID(s) from PDB (e.g. 6726807, 10490104)
        @type pubmed_ids: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.PubmedIdQuery'

        description = SubElement(title, 'description')
        description.text = 'Simple query for a list of PubMed IDs: ' + str(pubmed_ids)

        pubmed_id_list_query = SubElement(title, 'pubmedIdList')
        pubmed_id_list_query.text = str(pubmed_ids)

        self.query.append(query_refinement)
        self.counter += 1

    def uniprotkb_accession_number_s(self, up_accession_ids):
        """Builds a query for ID(s) and Keywords search purposes with specification UniProtKB Accession Number(s).

        @param up_accession_ids: UniProtKB Accession Number(s) (e.g. P69905)
        @type up_accession_ids: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.UpAccessionIdQuery'

        description = SubElement(title, 'description')
        description.text = 'Simple query for a list of Uniprot Accession IDs: ' + str(up_accession_ids)

        accession_id_list_query = SubElement(title, 'accessionIdList')
        accession_id_list_query.text = str(up_accession_ids)

        self.query.append(query_refinement)
        self.counter += 1

    def text_search(self, text_search):
        """Builds a query for ID(s) and Keywords search purposes with specification Text Search.

        @param advanced_keywords: the full text of the mmCIF coordinate file (e.g. actin).
        @type advanced_keywords: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.AdvancedKeywordQuery'

        description = SubElement(title, 'description')
        description.text = 'Text Search for: ' + str(text_search)

        keywords_query = SubElement(title, 'keywords')
        keywords_query.text = str(text_search)

        self.query.append(query_refinement)
        self.counter += 1

    def mmcif_keyword_search_classification(self, struct_keywords_comparator, struct_keywords):
        """Builds a query for ID(s) and Keywords search purposes with specification mmCIF Keyword Search
        (Classificaton).

        @param struct_keywords_comparator: Comparator
        @type struct_keywords_comparator: StructComparator

        @param struct_keywords: Search mmCIF item (e.g. Unknown Function).
        @type struct_keywords: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.TokenKeywordQuery'

        description = SubElement(title, 'description')
        description.text = 'TokenKeywordQuery: struct_keywords.pdbx_keywords.comparator=' \
                           + struct_keywords_comparator.value + ' struct_keywords.pdbx_keywords.value=' \
                           + str(struct_keywords)

        struct_keywords_comparator_query = SubElement(title, 'struct_keywords.pdbx_keywords.comparator')
        struct_keywords_comparator_query.text = struct_keywords_comparator.value

        struct_keywords_query = SubElement(title, 'struct_keywords')
        struct_keywords_query.text = str(struct_keywords)

        self.query.append(query_refinement)
        self.counter += 1

    def pfam_accession_number_s(self, pfam_accession_ids):
        """Builds a query for ID(s) and Keywords search purposes with specification Pfam Accession Number(s).

        @param pfam_accession_ids: Pfam ID(s) (e.g. PF00104).
        @type struct_keywords_comparator: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.PfamIdQuery'

        description = SubElement(title, 'description')
        description.text = 'Pfam Accession Numbers: ' + str(pfam_accession_ids)

        pfam_id_query = SubElement(title, 'pfamID')
        pfam_id_query.text = str(pfam_accession_ids)

        self.query.append(query_refinement)
        self.counter += 1

    def uniprot_gene_name_s(self, uniprot_gene_name):
        """Builds a query for ID(s) and Keywords search purposes with specification UniProt Gene Name(s).

        @param uniprot_gene_name: UniProt Gene Name(s) (e.g. HBA1).
        @type uniprot_gene_name: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.UniprotGeneNameQuery'

        description = SubElement(title, 'description')
        description.text = 'UniProt Gene Name: ' + str(uniprot_gene_name)

        uniprot_gene_name_query = SubElement(title, 'pfamID')
        uniprot_gene_name_query.text = str(uniprot_gene_name)

        self.query.append(query_refinement)
        self.counter += 1

    def sequence_cluster_name_s(self, sequence_cluster_name, sequence_cluster_name_comparator):
        """Builds a query for ID(s) and Keywords search purposes with specification Sequence Cluster Name(s).

        @param sequence_cluster_name_comparator: Comparator.
        @type uniprot_gene_name: SequenceStructureName

        @param sequence_cluster_name: Sequence Cluster Name(s) (e.g. Arginase).
        @type uniprot_gene_name: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.SequenceClusterQuery'

        description = SubElement(title, 'description')
        description.text = 'Sequence Cluster Name Search : Name=' + str(sequence_cluster_name)

        sequence_cluster_name_query = SubElement(title, 'sequenceClusterName')
        sequence_cluster_name_query.text = str(sequence_cluster_name)

        sequence_cluster_name_comparator_query = SubElement(title, 'comparator')
        sequence_cluster_name_comparator_query.text = sequence_cluster_name_comparator.value

        self.query.append(query_refinement)
        self.counter += 1

    # Structure Annotation
    def structure_title(self, struct_title_comparator, struct_title):
        """Builds a query for Structure Annotation search purposes with specification Structure Title.

        @param struct_title_comparator: Comparator.
        @type struct_title_comparator: StructComparator

        @param struct_title: PDB 'TITLE' record or mmCIF title (e.g. Solution NMR structure of protein).
        @type struct_title: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.StructTitleQuery'

        description = SubElement(title, 'description')
        description.text = 'StructTitleQuery: struct.title.comparator=' + struct_title_comparator.value \
                           + ' struct.title.value=' + str(struct_title)

        struct_title_comparator_query = SubElement(title, 'struct.title.comparator')
        struct_title_comparator_query.text = struct_title_comparator.value

        struct_title_query = SubElement(title, 'struct.title.value')
        struct_title_query.text = str(struct_title)

        self.query.append(query_refinement)
        self.counter += 1

    def structure_description(self, struct_description_comparator, struct_description):
        """Builds a query for Structure Annotation search purposes with specification Structure Description.

        @param struct_description_comparator: Comparator.
        @type struct_description_comparator: StructComparator

        @param struct_description: PDB 'COMPND' record or mmCIF description (e.g. Tyrosine-protein).
        @type struct_description: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.StructDescQuery'

        description = SubElement(title, 'description')
        description.text = 'StructDescQuery: entity.pdbx_description.comparator=' \
                           + struct_description_comparator.value + ' entity.pdbx_description.value=' \
                           + str(struct_description)

        struct_description_comparator_query = SubElement(title, 'entity.pdbx_description.comparator')
        struct_description_comparator_query.text = struct_description_comparator.value

        struct_description_query = SubElement(title, 'entity.pdbx_description.value')
        struct_description_query.text = str(struct_description)

        self.query.append(query_refinement)
        self.counter += 1

    def macromolecule_name_s(self, macromolecule_name):
        """Builds a query for Structure Annotation search purposes with specification Macromolecule Name.

        @param macromolecule_name:  macromolecule name (e.g. mCherry).
        @type macromolecule_name: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.MoleculeNameQuery'

        description = SubElement(title, 'description')
        description.text = 'Molecule Name Search : Molecule Name=' + str(macromolecule_name)

        macromolecule_name_query = SubElement(title, 'macromoleculeName')
        macromolecule_name_query.text = str(macromolecule_name)

        self.query.append(query_refinement)
        self.counter += 1

    def large_structures(self, large_structures):
        """Builds a query for Structure Annotation search purposes with specification Large Structures.

        @param large_structures: _pdbx_database_status.pdb_format_compatible.
        @type large_structures: LargeStructures
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.LargeStructureQuery'

        description = SubElement(title, 'description')
        description.text = 'Search for Large Structures: search option=' + large_structures.value

        large_structures_query = SubElement(title, 'macromoleculeName')
        large_structures_query.text = large_structures.value

        self.query.append(query_refinement)
        self.counter += 1

    # Deposition
    def author_name(self, search_type, author_name, exact_match):
        """Builds a query for Deposition search purposes with specification Author Name.

        @param search_type: type of authors.
        @type search_type: AuthorName

        @param author_name: structure author name (audit author) and/or primary citation author name
        (e.g. Perutz, M.F.).
        @type author_name: str

        @param exact_match: exact match or not.
        @type exact_match: TrueOrFalse
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.AdvancedAuthorQuery'

        description = SubElement(title, 'description')
        description.text = 'Author Name: Search type is ' + search_type.value + ' and Author is ' + str(author_name) \
                           + ' and Exact match is ' + exact_match.value

        search_type_query = SubElement(title, 'searchType')
        search_type_query.text = search_type.value

        author_name_query = SubElement(title, 'audit_author.name')
        author_name_query.text = str(author_name)

        exact_match_query = SubElement(title, 'exactMatch')
        exact_match_query.text = exact_match.value

        self.query.append(query_refinement)
        self.counter += 1

    def deposit_date(self, date_min, date_max):
        """Builds a query for Deposition search purposes with specification Deposit Date.

        @param date_min: start date (e.g. 2008-01-01).
        @type date_min: str

        @param date_max: stop date (e.g. 2008-12-31).
        @type date_max: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.DepositDateQuery'

        original_comparator_query = SubElement(title, 'database_PDB_rev.date_original.comparator')
        original_comparator_query.text = 'between'

        date_min_query = SubElement(title, 'database_PDB_rev.date_original.min')
        date_min_query.text = str(date_min)

        date_max_query = SubElement(title, 'database_PDB_rev.date_original.max')
        date_max_query.text = str(date_max)

        mode_type_query = SubElement(title, 'database_PDB_rev.mod_type.value')
        mode_type_query.text = '1'

        self.query.append(query_refinement)
        self.counter += 1

    def release_date(self, date_min, date_max):
        """Builds a query for Deposition search purposes with specification Release Date.

        @param date_min: start date (e.g. 2008-01-01).
        @type date_min: str

        @param date_max: stop date (e.g. 2008-12-31).
        @type date_max: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.ReleaseDateQuery'

        comparator_query = SubElement(title, 'database_PDB_rev.date.comparator')
        comparator_query.text = 'between'

        date_min_query = SubElement(title, 'database_PDB_rev.date.min')
        date_min_query.text = str(date_min)

        date_max_query = SubElement(title, 'database_PDB_rev.date.max')
        date_max_query.text = str(date_max)

        mode_type_query = SubElement(title, 'database_PDB_rev.mod_type.value')
        mode_type_query.text = '1'

        self.query.append(query_refinement)
        self.counter += 1

    def revise_date(self, date_min, date_max):
        """Builds a query for Deposition search purposes with specification Revise Date.

        @param date_min: start date (e.g. 2008-01-01).
        @type date_min: str

        @param date_max: stop date (e.g. 2008-12-31).
        @type date_max: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.ReviseDateQuery'

        comparator_query = SubElement(title, 'database_PDB_rev.date.comparator')
        comparator_query.text = 'between'

        date_min_query = SubElement(title, 'database_PDB_rev.date.min')
        date_min_query.text = str(date_min)

        date_max_query = SubElement(title, 'database_PDB_rev.date.max')
        date_max_query.text = str(date_max)

        mode_type_query = SubElement(title, 'database_PDB_rev.mod_type.value')
        mode_type_query.text = '1'

        self.query.append(query_refinement)
        self.counter += 1

    def latest_released_date(self):
        """Builds a query for Deposition search purposes with specification Latest Released Date.
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.LastLoadQuery'

        self.query.append(query_refinement)
        self.counter += 1

    def latest_modifies_structures(self):
        """Builds a query for Deposition search purposes with specification Latest Released Date.
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.ModifiedStructuresQuery'

        self.query.append(query_refinement)
        self.counter += 1

    def structural_genomics_project(self, project_name, center_initials, center_name):
        """Builds a query for Deposition search purposes with specification Latest Released Date.

        @param project_name: project name.
        @type project_name: ProjectName

        @param center_initials: center initials.
        @type center_initials: CenterInitials

        @param center_name: center name.
        @type center_name: CenterName
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.SGProjectQuery'

        description = SubElement(title, 'description')
        description.text = 'Structural Genomics Search : Center Initials=' + center_initials.value + '  Center Name=' \
                           + center_name.value

        project_name_query = SubElement(title, 'pdbx_SG_project.project_name.value')
        project_name_query.text = project_name.value

        center_initials_query = SubElement(title, 'pdbx_SG_project.initial_of_center.value')
        center_initials_query.text = center_initials.value

        center_name_query = SubElement(title, 'pdbx_SG_project.full_name_of_center.value')
        center_name_query.text = center_name.value

        project_id_query = SubElement(title, 'pdbx_SG_project.id.value')
        project_id_query.text = '-1'

        project_id_comparator_query = SubElement(title, 'pdbx_SG_project.id.comparator')
        project_id_comparator_query.text = '<![CDATA[>]]>'

        self.query.append(query_refinement)
        self.counter += 1

    # Structure Features
    def macromolecule_type(self, contains_protein, contains_dna, contains_rna, contains_dna_rna_hybrid):
        """Builds a query for Structure Features search purposes with specification Macromolecule Type.

        @param contains_protein: contains protein.
        @type contains_protein: YesNoIgnore

        @param contains_dna: contains DNA.
        @type contains_dna: YesNoIgnore

        @param contains_rna: contains RNA.
        @type contains_rna: YesNoIgnore

        @param contains_dna_rna_hybrid: contains RNA.
        @type contains_dna_rna_hybrid: YesNoIgnore
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.ChainTypeQuery'

        description = SubElement(title, 'description')
        description.text = 'Chain Type Search : Contains Protein=' + contains_protein.value + ', Contains DNA=' \
                           + contains_dna.value + ', Contains RNA=' + contains_rna.value + ', Contains DNA/RNA Hybrid=' \
                           + contains_dna_rna_hybrid.value

        contains_protein_query = SubElement(title, 'containsProtein')
        contains_protein_query.text = contains_protein.value

        contains_dna_query = SubElement(title, 'containsDna')
        contains_dna_query.text = contains_dna.value

        contains_rna_query = SubElement(title, 'containsRna')
        contains_rna_query.text = contains_rna.value

        contains_dna_rna_hybrid_query = SubElement(title, 'containsHybrid')
        contains_dna_rna_hybrid_query.text = contains_dna_rna_hybrid.value

        self.query.append(query_refinement)
        self.counter += 1

    def number_of_chains_au(self, min_number, max_number):
        """Builds a query for Structure Features search purposes with specification Number of Chains (Asymmetric Unit).

        @param min_number: minimal number of chains (e.g. 5).
        @type min_number: str

        @param max_number: maximal number of chains (e.g. 6).
        @type max_number: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.NumberOfChainsQuery'

        description = SubElement(title, 'description')
        description.text = 'Number of Chains Search : Min Number of Chains=' + str(min_number) \
                           + ' Max Number of Chains=' + str(max_number)

        min_number_query = SubElement(title, 'struct_asym.numChains.min')
        min_number_query.text = str(min_number)

        max_number_query = SubElement(title, 'struct_asym.numChains.max')
        max_number_query.text = str(max_number)

        self.query.append(query_refinement)
        self.counter += 1

    def number_of_chains_ba(self, min_number, max_number):
        """Builds a query for Structure Features search purposes with specification Number of Chains
        (Biological Assembly).

        @param min_number: minimal number of chains (e.g. 5).
        @type min_number: str

        @param max_number: maximal number of chains (e.g. 6).
        @type max_number: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.simple.BiolUnitQuery'

        description = SubElement(title, 'description')
        description.text = 'Oligomeric state Search : Min Number of oligomeric state=' + str(min_number) \
                           + ' Max Number of oligomeric state=' + str(max_number)

        min_number_query = SubElement(title, 'oligomeric_statemin')
        min_number_query.text = str(min_number)

        max_number_query = SubElement(title, 'oligomeric_statemax')
        max_number_query.text = str(max_number)

        self.query.append(query_refinement)
        self.counter += 1

    def number_of_entities(self, entity_type, min_number, max_number):
        """Builds a query for Structure Features search purposes with specification Number of Entities.

        @param entity_type: entity type.
        @type entity_type: EntityType

        @param min_number: minimal number of chains (e.g. 5).
        @type min_number: str

        @param max_number: maximal number of chains (e.g. 6).
        @type max_number: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.NumberOfEntitiesQuery'

        description = SubElement(title, 'description')
        description.text = 'Number of Entities Search : Entity Type=' + entity_type.value + ' Min Number of Entities=' \
                           + str(min_number) + ' Max Number of Entities=' + str(max_number)

        entity_type_query = SubElement(title, 'struct_asym.numEntities.min')
        entity_type_query.text = entity_type.value

        min_number_query = SubElement(title, 'struct_asym.numEntities.min')
        min_number_query.text = str(min_number)

        max_number_query = SubElement(title, 'struct_asym.numEntities.max')
        max_number_query.text = str(max_number)

        self.query.append(query_refinement)
        self.counter += 1

    def protein_stoichiometry(self, stoichiometry):
        """Builds a query for Structure Features search purposes with specification Protein Stoichiometry.

        @param stoichiometry: composition of biological assembly (e.g. A2B2).
        @type stoichiometry: str

        @param min_number: minimal number of chains (e.g. 4).
        @type min_number: str

        @param max_number: maximal number of chains (e.g. 6).
        @type max_number: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.StoichiometryQuery'

        description = SubElement(title, 'description')
        description.text = 'Stoichiometry in biological assembly: Stoichiometry is ' + str(stoichiometry)

        stoichiometry_query = SubElement(title, 'stoichiometry')
        stoichiometry_query.text = str(stoichiometry)

        self.query.append(query_refinement)
        self.counter += 1

    def protein_symmetry(self, symmetry, rmsd_min, rmsd_max):
        """Builds a query for Structure Features search purposes with specification Protein Symmetry.

        @param symmetry: Comparator.
        @type symmetry: ProteinSymmetry

        @param rmsd_min: minimal RMSD (e.g. 0.5).
        @type rmsd_min: str

        @param rmsd_max: maximal RMSD (e.g. 1.0).
        @type rmsd_max: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.PointGroupQuery'

        description = SubElement(title, 'description')
        description.text = 'Finds PDB entries based on symmetry: Protein Symmetry is ' + symmetry.value \
                           + ' and R m s d min is ' + str(rmsd_min) + ' and R m s d max is ' + str(rmsd_max)

        symmetry_query = SubElement(title, 'pointGroup')
        symmetry_query.text = symmetry.values

        rmsd_comparator_query = SubElement(title, 'rMSDComparator')
        rmsd_comparator_query.text = 'between'

        rmsd_min_query = SubElement(title, 'rMSDMin')
        rmsd_min_query.text = str(rmsd_min)

        rmsd_max_query = SubElement(title, 'rMSDMax')
        rmsd_max_query.text = str(rmsd_max)

        self.query.append(query_refinement)
        self.counter += 1

    def protein_symmetry_browser(self):   # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    def number_of_models(self, model_count_comparator, model_count):
        """Builds a query for Structure Features search purposes with specification Number of Model(s).

        @param model_count_comparator: Comparator.
        @type model_count_comparator: CountComparator

        @param model_count: number of MODEL records in the PDB file (e.g. 2).
        @type model_count: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ModelCountQuery'

        description = SubElement(title, 'description')
        description.text = 'ModelCountQuery: mvStructure.modelCount.comparator=' + model_count_comparator.value \
                           + ' mvStructure.modelCount.value=' + str(model_count)

        model_count_comparator_query = SubElement(title, 'mvStructure.modelCount.comparator')
        model_count_comparator_query.text = model_count_comparator.value

        model_count_query = SubElement(title, 'mvStructure.modelCount.value')
        model_count_query.text = str(model_count)

        self.query.append(query_refinement)
        self.counter += 1

    def number_of_disulfide_bonds(self, min_number, max_number):
        """Builds a query for Structure Features search purposes with specification Number of Disulfide Bond(s).

        @param min_number: minimal number od disulfide bonds (e.g. 5).
        @type min_number: str

        @param max_number: maximal number od disulfide bonds (e.g. 6).
        @type max_number: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.CloseContactsQuery'

        description = SubElement(title, 'description')
        description.text = '# of Disulfide link records: Min is ' + str(min_number) + ' and Max is ' + str(max_number)

        min_number_query = SubElement(title, 'min')
        min_number_query.text = str(min_number)

        max_number_query = SubElement(title, 'max')
        max_number_query.text = str(max_number)

        self.query.append(query_refinement)
        self.counter += 1

    def link_records(self, link_connection_type, component, atom_label, connected_component, connected_atom_label):
        """Builds a query for Structure Features search purposes with specification Link Records(s).

        @param link_connection_type: Camparator.
        @type link_connection_type: LinkConnectionType

        @param component: component (e.g. HEM).
        @type component: str

        @param atom_label: atom label (e.g. FE).
        @type atom_label: str

        @param connected_component: connected component (e.g. NBN).
        @type connected_component: str

        @param connected_atom_label: connected atom label (e.g. C).
        @type connected_atom_label: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.LinkConnectionQuery'

        description = SubElement(title, 'description')
        description.text = 'Link records: Link connection type is ' + link_connection_type.value \
                           + ' and Component is ' + str(component) + ' and Atom label is ' + str(atom_label) \
                           + ' and Connected component is ' + str(connected_component) \
                           + ' and Connected atom label is ' + str(connected_atom_label)

        link_connection_type_query = SubElement(title, 'linkConnectionType')
        link_connection_type_query.text = link_connection_type.value

        component_query = SubElement(title, 'component')
        component_query.text = str(component)

        atom_label_query = SubElement(title, 'atomLabel')
        atom_label_query.text = str(atom_label)

        connected_component_query = SubElement(title, 'connectedComponent')
        connected_component_query.text = str(connected_component)

        connected_atom_label_query = SubElement(title, 'connectedAtomLabel')
        connected_atom_label_query.text = str(connected_atom_label)

        self.query.append(query_refinement)
        self.counter += 1

    def molecular_weight_struct(self, min_weight, max_weight):
        """Builds a query for Structure Features search purposes with specification Number of Molecular Weight
        (Structure).

        @param min_number: minimal weight (e.g. 26700.0).
        @type min_number: str

        @param max_number: maximal number weight (e.g. 26800.0).
        @type max_number: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.MolecularWeightQuery'

        description = SubElement(title, 'description')
        description.text = 'Molecular Weight Search : Min Molecular Weight=' + str(
            min_weight) + ' Max Molecular Weight=' \
                           + str(max_weight)

        min_weigh_query = SubElement(title, 'entity.formula_weight.min')
        min_weigh_query.text = str(min_weight)

        max_weigh_query = SubElement(title, 'entity.formula_weight.max')
        max_weigh_query.text = str(max_weight)

        self.query.append(query_refinement)
        self.counter += 1

    def secondary_structure_content(self, helix_min_per, helix_max_per, helix_min_num, helix_max_num,
                                    sheets_min_per, sheets_max_per, sheets_min_num, sheets_max_num):
        """Builds a query for Structure Features search purposes with specification Secondary Structure Content.

        @param helix_min_per: minimal percent content of alpha helices (e.g. 12).
        @type helix_min_per: str

        @param helix_max_per: maximal percent content of alpha helices (e.g. 100).
        @type helix_max_per: str

        @param helix_min_num: minimal number of alpha helices (e.g. 44).
        @type helix_min_num: str

        @param helix_min_num: maximal number of alpha helices (e.g. 100).
        @type helix_min_num: str

        @param sheets_min_per: minimal percent content of beta sheets (e.g. 12).
        @type sheets_min_per: str

        @param sheets_max_per: maximal percent content of beta sheets (e.g. 100).
        @type sheets_max_per: str

        @param sheets_min_num: minimal number of beta sheets (e.g. 22).
        @type sheets_min_num: str

        @param sheets_min_num: maximal number of beta sheets (e.g. 100).
        @type sheets_min_num: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.SecondaryStructureQuery'

        description = SubElement(title, 'description')
        description.text = 'SecondaryStructureQuery: polyStats.helixPercent.comparator=between' \
                           + ' polyStats.helixCount.comparator=between' \
                           + ' polyStats.sheetPercent.comparator=between polyStats.sheetCount.comparator=between' \
                           + ' polyStats.helixPercent.min=' + str(helix_min_per) + ' polyStats.helixPercent.max=' \
                           + str(helix_max_per) + ' polyStats.helixCount.min=' + str(helix_min_num) \
                           + ' polyStats.helixCount.max=' + str(helix_max_num) + ' polyStats.sheetPercent.min=' \
                           + str(sheets_min_per) + ' polyStats.sheetPercent.max=' + sheets_max_per \
                           + ' polyStats.sheetCount.min=' + str(sheets_min_num) + 'polyStats.sheetCount.max=' \
                           + str(sheets_max_num)

        helix_per_comparator_query = SubElement(title, 'polyStats.helixPercent.comparator')
        helix_per_comparator_query.text = 'between'

        helix_num_comparator_query = SubElement(title, 'polyStats.helixCount.comparator')
        helix_num_comparator_query.text = 'between'

        sheets_per_comparator_query = SubElement(title, 'polyStats.sheetPercent.comparator')
        sheets_per_comparator_query.text = 'between'

        sheets_num_comparator_query = SubElement(title, 'polyStats.sheetCount.comparator')
        sheets_num_comparator_query.text = 'between'

        helix_min_per_query = SubElement(title, 'polyStats.helixPercent.min')
        helix_min_per_query.text = str(helix_min_per)

        helix_max_per_query = SubElement(title, 'polyStats.helixPercent.max')
        helix_max_per_query.text = str(helix_max_per)

        helix_min_num_query = SubElement(title, 'polyStats.helixCount.min')
        helix_min_num_query.text = str(helix_min_num)

        helix_max_per_query = SubElement(title, 'polyStats.helixCount.max')
        helix_max_per_query.text = str(helix_max_per)

        sheets_min_per_query = SubElement(title, 'polyStats.sheetPercent.min')
        sheets_min_per_query.text = str(sheets_min_num)

        sheets_max_per_query = SubElement(title, 'polyStats.sheetPercent.max')
        sheets_max_per_query.text = str(sheets_max_per)

        sheets_min_num_query = SubElement(title, 'polyStats.sheetCount.min')
        sheets_min_num_query.text = str(sheets_min_num)

        sheets_max_num_query = SubElement(title, 'polyStats.sheetCount.max')
        sheets_max_num_query.text = str(sheets_max_num)

        self.query.append(query_refinement)
        self.counter += 1

    def secondary_structure_length(self, secondary_structure_method, helix_subtype, helix_length_min, helix_length_max,
                                   beta_subtype, beta_strand_length_min, beta_strand_length_max):
        """Builds a query for Structure Features search purposes with specification Secondary Structure Length.

        @param secondary_structure_method: secondary structure method.
        @type secondary_structure_method: SecondaryStructureMethod

        @param helix_subtype: subtype of helix
        @type helix_subtype: HelixSubtype

        @param helix_length_min: minimal helix length (e.g. 22).
        @type helix_length_min: str

        @param helix_length_max: maximal helix length (e.g. 50).
        @type helix_length_max: str

        @param beta_subtype: subtype of beta strand
        @type beta_subtype: BetaStrandSubtype

        @param beta_strand_length_min: minimal beta strand length (e.g. 12).
        @type beta_strand_length_min: str

        @param beta_strand_length_max: maximal beta strand length (e.g. 100).
        @type beta_strand_length_max:
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.SecondaryStructureLengthQuery'

        description = SubElement(title, 'description')
        description.text = 'SecondaryStructureLengthQuery: ss.method.value=' + secondary_structure_method.value \
                           + ' helix.subtype.value=' + helix_subtype.value \
                           + ' helix.length.comparator=between beta_strand.subtype.value=' + beta_subtype.value \
                           + ' helix.length.min=' + str(helix_length_min) + ' helix.length.max=' \
                           + str(helix_length_max) + ' beta_strand.length.comparator=between beta_strand.length.min=' \
                           + str(beta_strand_length_min) + 'beta_strand.length.max=' + str(beta_strand_length_max)

        secondary_structure_method_query = SubElement(title, 'ss.method.value')
        secondary_structure_method_query.text = secondary_structure_method.value

        helix_subtype_query = SubElement(title, 'helix.subtype.value')
        helix_subtype_query.text = helix_subtype.value

        helix_length_comparator_query = SubElement(title, 'helix.length.comparator')
        helix_length_comparator_query.text = 'between'

        helix_length_min_query = SubElement(title, 'helix.length.min')
        helix_length_min_query.text = str(helix_length_min)

        helix_length_max_query = SubElement(title, 'helix.length.max')
        helix_length_max_query.text = str(helix_length_max)

        beta_subtype_query = SubElement(title, 'beta_strand.subtype.value')
        beta_subtype_query.text = beta_subtype.value

        beta_strand_comparator_query = SubElement(title, 'beta_strand.length.comparator')
        beta_strand_comparator_query.text = 'between'

        beta_strand_length_min_query = SubElement(title, 'beta_strand.length.min')
        beta_strand_length_min_query.text = str(beta_strand_length_min)

        beta_strand_length_max_query = SubElement(title, 'beta_strand.length.max')
        beta_strand_length_max_query.text = str(beta_strand_length_max)

        self.query.append(query_refinement)
        self.counter += 1

    def scop_classification_browser(self):   # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    def cath_classification_browser(self):    # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    def taxon_browser(self):   # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    # Sequence features
    def sequence_blast_fasta_psi_blast(self, structure_id, chain_id, sequence, search_tool, mask_low_complexity,
                                       e_value_cutoff, seq_identity_cutoff):
        """Builds a query for Sequence Features search purposes with specification Sequence (BLAST/FASTA/PSI_BLAST).

        @param structure_id: structure ID from PDB (e.g. 4HHB, 2CPK).
        @type structure_id: str

        @param chain_id: chain ID(s) from PDB (e.g. 4HHB.A, 2CPK.E; the chain IDs are case sensitive)
        @type chain_id: str

        @param sequence: sequence (e.g. VLSPADKTNVKAAWGKVGAHAGEYGAEAL).
        @type sequence: str

        @param search_tool: search tool.
        @type search_tool: SearchTool

        @param mask_low_complexity: mask low complexity regions or not
        @type mask_low_complexity: YesNo

        @param e_value_cutoff: expectation value cutoff (e.g. 10.0).
        @type e_value_cutoff: str

        @param seq_identity_cutoff: sequence identity cutoff (e.g. 30).
        @type seq_identity_cutoff: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.SequenceQuery'

        description = SubElement(title, 'description')
        description.text = 'Sequence Search (Structure:Chain = ' + str(structure_id) + ':' + str(chain_id) \
                           + ', Expectation Value = ' + str(e_value_cutoff) + ', Sequence Identity = ' \
                           + str(seq_identity_cutoff) + '%, Search Tool = ' + search_tool.value \
                           + ', Mask Low Complexity=' + mask_low_complexity.value + ')'

        structure_id_query = SubElement(title, 'structureId')
        structure_id_query.text = str(structure_id)

        chain_id_query = SubElement(title, 'chainId')
        chain_id_query.text = str(chain_id)

        sequence_query = SubElement(title, 'sequence')
        sequence_query.text = str(sequence)

        search_tool_query = SubElement(title, 'searchTool')
        search_tool_query.text = search_tool.value

        mask_low_complexity_query = SubElement(title, 'searchTool')
        mask_low_complexity_query.text = mask_low_complexity.value

        e_value_cutoff_query = SubElement(title, 'eCutOff')
        e_value_cutoff_query.text = str(e_value_cutoff)

        seq_identity_cutoff_query = SubElement(title, 'sequenceIdentityCutoff')
        seq_identity_cutoff_query.text = seq_identity_cutoff

        self.query.append(query_refinement)
        self.counter += 1

    def wild_type_protein(self, include_expression_tags, per_coverage_of_wtp_seq):
        """Builds a query for Sequence Features search purposes with specification Wild Type Protein.

        @param include_expression_tags: include expression tags or no.
        @type include_expression_tags: YesOrNo

        @param per_coverage_of_wtp_seq: percent of coverage of wild type protein sequence
        @type per_coverage_of_wtp_seq: PercentCoverageOfWTPSequence
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.WildTypeProteinQuery'

        description = SubElement(title, 'description')
        description.text = 'Simple query for a list of PDB Entity IDs'

        include_expression_tags_query = SubElement(title, 'includeExprTag')
        include_expression_tags_query.text = include_expression_tags.value

        per_coverage_of_wtp_seq_query = SubElement(title, 'percentSeqAlignment')
        per_coverage_of_wtp_seq_query.text = per_coverage_of_wtp_seq.value

        self.query.append(query_refinement)
        self.counter += 1

    def translated_nucleotide_sequence(self, sequence, e_value_cutoff):
        """Builds a query for Sequence Features search purposes with specification Translated Nucleotide Sequence
        (BLASTX).

        @param sequence: sequence (e.g. ATGCCAGGGGCAGT).
        @type sequence: str

        @param e_value_cutoff: expectation value cutoff (e.g. 10.0).
        @type e_value_cutoff: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.BlastXQuery'

        description = SubElement(title, 'description')
        description.text = 'BLASTX Search (Sequence = ' + sequence + '), Expectation Value = ' + e_value_cutoff

        sequence_query = SubElement(title, 'includeExprTag')
        sequence_query.text = sequence

        e_value_cutoff_query = SubElement(title, 'eCutOff')
        e_value_cutoff_query.text = e_value_cutoff

        self.query.append(query_refinement)
        self.counter += 1

    def sequence_motif(self, motif):
        """Builds a query for Sequence Features search purposes with specification Translated Nucleotide Sequence
        (BLASTX).

        @param motif: sequence (e.g. T[AG]Y).
        @type motif: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.MotifQuer'

        description = SubElement(title, 'description')
        description.text = 'Motif Query For: ' + motif

        motif_query = SubElement(title, 'motif')
        motif_query.text = motif

        self.query.append(query_refinement)
        self.counter += 1

    def chain_length(self, min_length, max_length):
        """Builds a query for Sequence Features search purposes with specification Chain Length.

        @param min_length: minimal chain length (e.g. 4).
        @type min_length: str

        @param max_length: maximal chain length (e.g. 5).
        @type max_length: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.SequenceLengthQuery'

        description = SubElement(title, 'description')
        description.text = 'Sequence Length Search : Min Sequence Length=' + str(min_length) + ' Max Sequence Length=' \
                           + (max_length)

        min_length_query = SubElement(title, 'v_sequence.chainLength.min')
        min_length_query.text = (min_length)

        max_length_query = SubElement(title, 'v_sequence.chainLength.max')
        max_length_query.text = (max_length)

        self.query.append(query_refinement)
        self.counter += 1

    def genome_location_browser(self):  # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    # Chemical Components
    def chemical_name(self, chemical_name_comparator, name, polymeric_type):
        """Builds a query for Chemical Components search purposes with specification Chemical Name.

        @param chemical_name_comparator: chemical name comparator.
        @type chemical_name_comparator: ChemicalNameComparator

        @param name: name of the chemical (e.g. biotin).
        @type name: str

        @param polymeric_type: polymeric type of the chemical.
        @type polymeric_type: PolymericType
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ChemCompNameQuery'

        description = SubElement(title, 'description')
        description.text = 'Chemical Name: Name ' + chemical_name_comparator.value + ' ' + str(name) \
                           + ' and Polymeric type is ' + polymeric_type.value

        chemical_name_comparator_query = SubElement(title, 'comparator')
        chemical_name_comparator_query.text = chemical_name_comparator.value

        name_query = SubElement(title, 'name')
        name_query.text = str(name)

        polymeric_type_query = SubElement(title, 'polymericType')
        polymeric_type_query.text = polymeric_type.value

        self.query.append(query_refinement)
        self.counter += 1

    def chemical_id_s(self, chemical_compound_id, polymeric_type):
        """Builds a query for Chemical Components search purposes with specification Chemical ID(s).

        @param chemical_compound_id: chemical compound ID (e.g CL).
        @type chemical_compound_id: ChemicalNameComparator

        @param polymeric_type: polymeric type of the chemical.
        @type polymeric_type: PolymericType
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ChemCompIdQuery'

        description = SubElement(title, 'description')
        description.text = 'Chemical ID(s):  ' + str(chemical_compound_id) + ' and Polymeric type is ' \
                           + polymeric_type.value

        chemical_compound_id_query = SubElement(title, 'chemCompId')
        chemical_compound_id_query.text = str(chemical_compound_id)

        polymeric_type_query = SubElement(title, 'polymericType')
        polymeric_type_query.text = polymeric_type.value

        self.query.append(query_refinement)
        self.counter += 1

    def inchi_descriptor(self, inchi, descriptor_type):
        """Builds a query for Chemical Components search purposes with specification InChI Descriptor.

        @param inchi: InChI string or its InChI key (e.g InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H).
        @type inchi: ChemicalNameComparator

        @param descriptor_type: descriptor type of InChI.
        @type descriptor_type: InChIDescriptorType
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ChemCompDescriptorQuery'

        description = SubElement(title, 'description')
        description.text = 'InChI Descriptor: Descriptor is InChI=' + str(inchi) + ' and Descriptor type is ' \
                           + descriptor_type.value

        inchi_query = SubElement(title, 'descriptor')
        inchi_query.text = str(inchi)

        descriptor_type_query = SubElement(title, 'descriptorType')
        descriptor_type_query.text = descriptor_type.value

        self.query.append(query_refinement)
        self.counter += 1

    def chemical_structure_smiles(self, smiles, search_type, similarity, polymeric_type):
        """Builds a query for Chemical Components search purposes with specification Chemical Structure.

        @param smiles: chemical structure SMILES (e.g CC(=O)NC1C(C(C(OC1O)CO)O)O).
        @type smiles: str

        @param search_type: search type
        @type search_type: StructureSearchType

        @param similarity: similarity (e.g. 0.7)
        @type similarity: str

        @param polymeric_type: polymeric type
        @type polymeric_type: PolymericType
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ChemSmilesQuery'

        description = SubElement(title, 'description')
        description.text = 'Chemical structure (SMILES): Structure is  and Smiles is ' + str(smiles) \
                           + ' and Search type is ' + search_type.value + ' and Similarity is ' + str(similarity) \
                           + ' and Polymeric type is ' + polymeric_type.value

        structure_query = SubElement(title, 'structure')
        structure_query.text = ''

        smiles_query = SubElement(title, 'smiles')
        smiles_query.text = str(smiles)

        search_type_query = SubElement(title, 'searchType')
        search_type_query.text = search_type.value

        similarity_query = SubElement(title, 'similarity')
        similarity_query.text = str(similarity)

        polymeric_type_query = SubElement(title, 'polymericType')
        polymeric_type_query.text = polymeric_type.value

        self.query.append(query_refinement)
        self.counter += 1

    def molecular_weight_chem_comp(self, min_weight, max_weight):
        """Builds a query for Chemical Components search purposes with specification Molecular Weight
        (Chemical component).

        @param min_weight: minimal molecular weight (e.g 0).
        @type min_weight: str

        @param max_weight: maximal molecular weight (e.g. 1000)
        @type max_weight: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ChemFormulaWeightQuery'

        description = SubElement(title, 'description')
        description.text = 'Molecular Weight (Chemical component): Molecular weight min is ' + str(min_weight) \
                           + ' and Molecular weight max is ' + str(max_weight)

        weight_comparator_query = SubElement(title, 'molecularWeightComparator')
        weight_comparator_query.text = 'between'

        min_weight_query = SubElement(title, 'molecularWeightMin')
        min_weight_query.text = str(min_weight)

        max_weight_query = SubElement(title, 'molecularWeightMax')
        max_weight_query.text = str(max_weight)

        self.query.append(query_refinement)
        self.counter += 1

    def chemical_formula(self, formula):
        """Builds a query for Chemical Components search purposes with specification Chemical Formula.

        @param formula: chemical formula (C10-15 N4 O* P0).
        @type min_weight: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ChemCompFormulaQuery'

        description = SubElement(title, 'description')
        description.text = 'Chemical Formula: Formula is ' + str(formula)

        formula_query = SubElement(title, 'formula')
        formula_query.text = str(formula)

        self.query.append(query_refinement)
        self.counter += 1

    def chemical_component_type(self, chem_component_type, polymeric_type):
        """Builds a query for Chemical Components search purposes with specification Chemical Component Type.

        @param chem_component_type: type of chemical component.
        @type type: ChemicalComponentType

        @param polymeric_type: polymeric type (on PDB site is always = Any).
        @type polymeric_type: PolymericType
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ChemCompTypeQuery'

        description = SubElement(title, 'description')
        description.text = 'Chemical Component Type: Type is ' + chem_component_type.value + ' and Polymeric type is' \
                           + polymeric_type.value

        chem_component_type_query = SubElement(title, 'type')
        chem_component_type_query.text = chem_component_type.value

        polymeric_type_query = SubElement(title, 'polymericType')
        polymeric_type_query.text = polymeric_type.value

        self.query.append(query_refinement)
        self.counter += 1

    def binding_affinity(self, binding_affinity_min, binding_affinity_max, affinity_type):
        """Builds a query for Chemical Components search purposes with specification Binding Afinity.

        @param binding_affinity_min: minimal binding affinity (e.g. -100).
        @type binding_affinity_min: str

        @param binding_affinity_max: maximal binding affinity (e.g. 100).
        @type binding_affinity_max: str

        @param affinity_type: afinity type.
        @type affinity_type: AffinityType
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.BindingAffinityQuery'

        description = SubElement(title, 'description')
        description.text = 'Binding Affinity: Binging affinity min is ' + str(binding_affinity_min) \
                           + ' and Binding affinity max is ' + str(binding_affinity_max) + ' and Affinity Type is ' \
                           + affinity_type.value

        binding_affinity_comparator_query = SubElement(title, 'bindingAffinityComparator')
        binding_affinity_comparator_query.text = 'between'

        binding_affinity_min_query = SubElement(title, 'bingingAffinityMin')
        binding_affinity_min_query.text = str(binding_affinity_min)

        binding_affinity_max_query = SubElement(title, 'bingingAffinityMax')
        binding_affinity_max_query.text = str(binding_affinity_max)

        affinity_type_query = SubElement(title, 'affinityType')
        affinity_type_query.text = affinity_type.value

        self.query.append(query_refinement)
        self.counter += 1

    def has_ligand_s(self, has_ligands):
        """Builds a query for Chemical Components search purposes with specification Has Ligand(s).

        @param has_ligands: structure contains any free ligands or not.
        @type has_ligands: YesNo
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.NoLigandQuery'

        description = SubElement(title, 'description')
        description.text = 'Ligand Search : Has free Ligands=' + has_ligands.value

        have_ligands_query = SubElement(title, 'haveLigands')
        have_ligands_query.text = has_ligands.value

        self.query.append(query_refinement)
        self.counter += 1

    def has_modified_residue_s(self, has_modified_residues):
        """Builds a query for Chemical Components search purposes with specification Has Modified Residue(s).

        @param has_modified_residues: structure contains any modifies residues or not.
        @type has_modified_residues: YesNo
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.NoModResQuery'

        description = SubElement(title, 'description')
        description.text = 'Ligand Search : Has modified polymeric residues=' + has_modified_residues.value

        have_modified_residues_query = SubElement(title, 'haveModifiedResidues')
        have_modified_residues_query.text = has_modified_residues.value

        self.query.append(query_refinement)
        self.counter += 1

    def sub_component_s(self, sub_components):
        """Builds a query for Chemical Components search purposes with specification Has Modified Residue(s).

        @param sub_components: sub-components of chemical component (e.g. 0QE, DPN PHE).
        @type sub_components: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ChemCompSubCompQuery'

        description = SubElement(title, 'description')
        description.text = 'Sub-components:  ' + str(sub_components)

        sub_components_query = SubElement(title, 'sub-components')
        sub_components_query.text = str(sub_components)

        self.query.append(query_refinement)
        self.counter += 1

    # Biologically Interesting Molecules (from BIRD)
    def biologically_interesting_molecules(self, bird_comparator, name_or_id, bird_type, bird_class):
        """Builds a query for Biologically Interesting Molecules (from BIRD) search purposes with specification
        Biologically Interesting Molecules (from BIRD).

        @param bird_comparator: Comparator.
        @type bird_comparator: BIRDNameComparator

        @param name_or_id: name or ID of BIRD molecule (e.g. Antibiotic).
        @type name_or_id: str

        @param bird_type: BIRD molecule type.
        @type bird_type: BIRDType

        @param bird_class: BIRD molecule class.
        @type bird_class: BIRDClass
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.BirdQuery'

        description = SubElement(title, 'description')
        description.text = 'Biologically Interesting Molecules (from BIRD): Name or ID ' + bird_comparator.value + ' ' \
                           + str(name_or_id) + ' and Type is ' + bird_type.value + ' and Class is ' + bird_class.value

        bird_comparator_query = SubElement(title, 'comparator')
        bird_comparator_query.text = bird_comparator.value

        name_or_id_query = SubElement(title, 'name')
        name_or_id_query.text = str(name_or_id)

        type_query = SubElement(title, 'birdType')
        type_query.text = bird_type.value

        bird_class_query = SubElement(title, 'birdClass')
        bird_class_query.text = bird_class.value

        self.query.append(query_refinement)
        self.counter += 1

    # Biology
    def source_organism_browser(self):  # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    def expression_organism(self, expression_organism_comparator, organism_name):
        """Builds a query for Biology search purposes with specification Expression Organism.

        @param expression_organism_comparator: Comparator.
        @type expression_organism_comparator: StructComparator

        @param organism_name: scientific name of the expression host organism (e.g. Saccharomyces).
        @type organism_name: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ExpressionOrganismQuery'

        description = SubElement(title, 'description')
        description.text = 'ExpressionOrganismQuery: entity_src_gen.pdbx_host_org_scientific_name.comparator=' \
                           + expression_organism_comparator.value \
                           + ' entity_src_gen.pdbx_host_org_scientific_name.value=' + str(organism_name)

        expression_organism_comparator_query = SubElement(title,
                                                          'entity_src_gen.pdbx_host_org_scientific_name.comparator')
        expression_organism_comparator_query.text = expression_organism_comparator.value

        organism_name_query = SubElement(title, 'entity_src_gen.pdbx_host_org_scientific_name.value')
        organism_name_query.text = str(organism_name)

        self.query.append(query_refinement)
        self.counter += 1

    def enzyme_classification_browser(self):  # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    def enzyme_classification(self, enzyme_classification):
        """Builds a query for Biology search purposes with specification Enzyme Classification.

        @param enzyme_classification: enzyme classification number (e.g. 2.7.11); wildcards are supported
        *for all protein-serine/threonine kinase.
        @type enzyme_classification: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.EnzymeClassificationQuery'

        description = SubElement(title, 'description')
        description.text = 'Enzyme Classification Search : EC=' + str(enzyme_classification)

        enzyme_classification_query = SubElement(title, 'Enzyme_Classification')
        enzyme_classification_query.text = str(enzyme_classification)

        self.query.append(query_refinement)
        self.counter += 1

    def biological_process_browser(self):  # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    def cell_component_browser(self):  # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    def molecular_function_browser(self):  # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    def transporter_classifiction_browser(self):  # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    # Methods
    def experimental_method(self, experimental_method, has_experimental_data):
        """Builds a query for Methods search purposes with specification Experimental Method.

        @param experimental_method: Comparator
        @type experimental_method: ExperimentalMethods

        @param has_experimental_data: Comparator
        @type has_experimental_data: YesNoOrIgnore
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ExpTypeQuery'

        description = SubElement(title, 'description')
        description.text = 'Experimental Method Search : Experimental Method=' + experimental_method.value \
                           + ', Has Experimental Data=' + has_experimental_data.value

        experimental_method_query = SubElement(title, 'mvStructure.expMethod.value')
        experimental_method_query.text = experimental_method.value

        has_experimental_data_query = SubElement(title, 'mvStructure.hasExperimentalData.value')
        has_experimental_data_query.text = has_experimental_data.value

        self.query.append(query_refinement)
        self.counter += 1

    def x_ray_resolution(self, min_resolution, max_resolution):
        """Builds a query for Methods search purposes with specification X-ray Resolution.

        @param min_resolution: minimal X-RAY resolution (e.g. 2.3).
        @type min_resolution: str

        @param max_resolution: maximal X-RAY resolution (e.g. 2.5)
        @type max_resolution: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ResolutionQuery'

        description = SubElement(title, 'description')
        description.text = 'ResolutionQuery: refine.ls_d_res_high.comparator=between refine.ls_d_res_high.min=' \
                           + str(min_resolution) + ' refine.ls_d_res_high.max=' + str(max_resolution)

        resolution_comparator_query = SubElement(title, 'refine.ls_d_res_high.comparator')
        resolution_comparator_query.text = 'between'

        min_resolution_query = SubElement(title, 'refine.ls_d_res_high.min')
        min_resolution_query.text = str(min_resolution)

        max_resolution_query = SubElement(title, 'refine.ls_d_res_high.max')
        max_resolution_query.text = str(max_resolution)

        self.query.append(query_refinement)
        self.counter += 1

    def x_ray_average_b_factor(self, min_mean, max_mean):
        """Builds a query for Methods search purposes with specification X-ray Average B Factor.

        @param min_mean: minimal X-RAY average B factor (e.g. 2.3).
        @type min_mean: str

        @param max_mean: maximal X-RAY average B factor (e.g. 2.5)
        @type max_mean: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.AverageBFactorQuery'

        description = SubElement(title, 'description')
        description.text = 'AverageBFactorQuery: refine.B_iso_mean.comparator=between refine.B_iso_mean.min=' \
                           + str(min_mean) + ' refine.B_iso_mean.max=' + str(max_mean)

        mean_comparator_query = SubElement(title, 'refine.ls_d_res_high.comparator')
        mean_comparator_query.text = 'between'

        min_mean_query = SubElement(title, 'refine.ls_d_res_high.min')
        min_mean_query.text = str(min_mean)

        max_mean_query = SubElement(title, 'refine.ls_d_res_high.max')
        max_mean_query.text = str(max_mean)

        self.query.append(query_refinement)
        self.counter += 1

    def refinement_r_factors(self, observed_min, observed_max, all_min, all_max, r_work_min, r_work_max, r_free_min,
                             r_free_max):
        """Builds a query for Methods search purposes with specification X-ray Average B Factor.

        @param observed_min: minimal observed R factor (e.g. 11).
        @type observed_min: str

        @param observed_max: maximal observed R factor (e.g. 77).
        @type observed_max: str

        @param all_min: minimal all R factor (e.g. 11).
        @type all_min: str

        @param all_max: maximal all R factor (e.g. 77).
        @type all_max: str

        @param r_work_min: minimal R Work factor (e.g. 12).
        @type r_work_min: str

        @param r_work_max: maximal R Work factor (e.g. 44).
        @type r_work_max: str

        @param r_free_min: minimal R Free factor (e.g. 14).
        @type r_free_min: str

        @param r_free_max: maximal R Free factor (e.g. 55).
        @type r_free_max: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.XrayRefinementQuer'

        description = SubElement(title, 'description')
        description.text = 'XrayRefinementQuery: refine.ls_R_factor_obs.comparator=between refine.ls_R_factor_obs.min='\
                           + str(observed_min) + ' refine.ls_R_factor_obs.max=' + str(observed_max) \
                           + ' refine.ls_R_factor_all.comparator=between refine.ls_R_factor_all.min=' + str(all_min) \
                           + ' refine.ls_R_factor_all.max=' + str(all_max) \
                           + ' refine.ls_R_factor_R_work.comparator=between refine.ls_R_factor_R_work.min=' \
                           + str(r_work_min) + ' refine.ls_R_factor_R_work.max=' + str(r_work_max) \
                           + ' refine.ls_R_factor_R_free.comparator=between refine.ls_R_factor_R_free.min=' \
                           + str(r_free_min) + 'refine.ls_R_factor_R_free.max=' + str(r_free_max)

        observed_comparator_query = SubElement(title, 'refine.ls_R_factor_obs.comparator')
        observed_comparator_query.text = 'between'

        observed_min_query = SubElement(title, 'refine.ls_R_factor_obs.min')
        observed_min_query.text = str(observed_min)

        observed_max_query = SubElement(title, 'refine.ls_R_factor_obs.max')
        observed_max_query.text = str(observed_max)

        all_comparator_query = SubElement(title, 'refine.ls_R_factor_all.comparator')
        all_comparator_query.text = 'between'

        all_min_query = SubElement(title, 'refine.ls_R_factor_all.min')
        all_min_query.text = str(all_min)

        all_max_query = SubElement(title, 'refine.ls_R_factor_all.max')
        all_max_query.text = str(all_max)

        r_work_comparator_query = SubElement(title, 'refine.ls_R_factor_R_work.comparator')
        r_work_comparator_query.text = 'between'

        r_work_min_query = SubElement(title, 'refine.ls_R_factor_R_work.min')
        r_work_min_query.text = str(r_work_min)

        r_work_max_query = SubElement(title, 'refine.ls_R_factor_R_work.max')
        r_work_max_query.text = str(r_work_max)

        r_free_comparator_query = SubElement(title, 'refine.ls_R_factor_R_free.comparator')
        r_free_comparator_query.text = 'between'

        r_free_min_query = SubElement(title, 'refine.ls_R_factor_R_free.min')
        r_free_min_query.text = str(r_free_min)

        r_free_max_query = SubElement(title, 'refine.ls_R_factor_R_free.max')
        r_free_max_query.text = str(r_free_max)

        self.query.append(query_refinement)
        self.counter += 1

    def diffraction_source(self, site_comparator, synchrotron_site, beamline_comparator, synchrotron_beamline,
                           min_temperature, max_temperature):
        """Builds a query for Methods search purposes with specification X-ray Average B Factor.

        @param site_comparator: Comparator.
        @type site_comparator: ComparatorChoice

        @param synchrotron_site: synchroton site (e.g. PHOTON FACTORY).
        @type synchrotron_site: str

        @param beamline_comparator: Comparator.
        @type beamline_comparator: ComparatorChoice

        @param synchrotron_beamline: synchroton beamline (e.g. AR-NW14A).
        @type synchrotron_beamline: str

        @param min_temperature: minimal temperature (e.g. 22).
        @type min_temperature: str

        @param max_temperature: maximal temperature (e.g. 77).
        @type max_temperature: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.XrayDiffrnSourceQuery'

        description = SubElement(title, 'description')
        description.text = 'XrayDiffrnSourceQuery: diffrn_source.source.comparator=like' \
                           + ' diffrn_source.pdbx_synchrotron_site.comparator=' + site_comparator.value \
                           + ' diffrn_source.pdbx_synchrotron_site.value=' + str(synchrotron_site) \
                           + ' diffrn_source.pdbx_synchrotron_beamline.comparator=' + beamline_comparator.value \
                           + ' diffrn_source.pdbx_synchrotron_beamline.value=' + str(synchrotron_beamline) \
                           + ' diffrn.ambient_temp.comparator=between diffrn.ambient_temp.min=' + str(min_temperature) \
                           + ' diffrn.ambient_temp.max= ' + str(max_temperature)

        diffraction_source_comparator_query = SubElement(title, 'diffrn_source.source.comparator')
        diffraction_source_comparator_query.text = 'like'

        site_comparator_query = SubElement(title, 'diffrn_source.pdbx_synchrotron_site.comparator')
        site_comparator_query.text = site_comparator.value

        synchrotron_site_query = SubElement(title, 'diffrn_source.pdbx_synchrotron_site.value')
        synchrotron_site_query.text = str(synchrotron_site)

        beamline_comparator_query = SubElement(title, 'diffrn_source.pdbx_synchrotron_beamline.comparator')
        beamline_comparator_query.text = beamline_comparator.value

        synchrotron_beamline_query = SubElement(title, 'diffrn_source.pdbx_synchrotron_beamline.value')
        synchrotron_beamline_query.text = str(synchrotron_beamline)

        temperature_comparator_query = SubElement(title, 'diffrn.ambient_temp.comparator')
        temperature_comparator_query.text = 'between'

        min_temperature_query = SubElement(title, 'diffrn.ambient_temp.min')
        min_temperature_query.text = str(min_temperature)

        max_temperature_query = SubElement(title, 'diffrn.ambient_temp.max')
        max_temperature_query.text = str(max_temperature)

        self.query.append(query_refinement)
        self.counter += 1

    def structure_determination_method(self, method_comparator, method_to_determine):
        """Builds a query for Methods search purposes with specification Structure Determination Method.

        @param method_comparator: Comparator.
        @type method_comparator: StructComparator

        @param method_to_determine: method to determine structure (e.g. molecular replacement).
        @type method_to_determine: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.MethodToDetermineStructQuery'

        description = SubElement(title, 'description')
        description.text = 'MethodToDetermineStructQuery: '

        method_comparator_query = SubElement(title, 'refine.pdbx_method_to_determine_struct.comparator')
        method_comparator_query.text = method_comparator.value

        method_to_determine_query = SubElement(title, 'refine.pdbx_method_to_determine_struct.value')
        method_to_determine_query.text = str(method_to_determine)

        self.query.append(query_refinement)
        self.counter += 1

    def reflections(self, d_resolution_high_min, d_resolution_high_max, b_iso_wilson_estimate_min,
                    b_iso_wilson_estimate_max, overall_redundancy_min, overall_redundancy_max,
                    per_of_possible_reflections_min, per_of_possible_reflections_max, r_value_min, r_value_max):
        """Builds a query for Methods search purposes with specification Reflections.

        @param d_resolution_high_min: minimal D resolution high (angstroms) (e.g. 1.5).
        @type d_resolution_high_min: str

        @param d_resolution_high_max: maximal D resolution high (angstroms) (e.g. 2.8).
        @type d_resolution_high_max: str

        @param b_iso_wilson_estimate_min: minimal B iso Wilson estimate (angstroms^2) (e.g. 1).
        @type b_iso_wilson_estimate_min: str

        @param b_iso_wilson_estimate_max: maximal B iso Wilson estimate (angstroms^2) (e.g. 9).
        @type b_iso_wilson_estimate_max: str

        @param overall_redundancy_min: minimal overall redundancy (e.g. 2).
        @type overall_redundancy_min: str

        @param overall_redundancy_max: maximal overall redundancy (e.g. 50).
        @type overall_redundancy_max: str

        @param per_of_possible_reflections_min: minimal percentage of possible reflections (%) (e.g. 2).
        @type per_of_possible_reflections_min: str

        @param per_of_possible_reflections_max: maximal percentage of possible reflections (%) (e.g. 8).
        @type per_of_possible_reflections_max: str

        @param r_value_min: R Value for Merging Intensities (Observed) (e.g. 1).
        @type r_value_min: str

        @param r_value_max: R Value for Merging Intensities (Observed) (e.g. 2).
        @type r_value_max: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.XrayReflnsQuery'

        description = SubElement(title, 'description')
        description.text = 'XrayReflnsQuery: reflns.d_resolution_low.comparator=between' \
                           + ' reflns.d_resolution_high.comparator=between reflns.d_resolution_high.min=' \
                           + str(d_resolution_high_min) + ' reflns.d_resolution_high.max=' + str(d_resolution_high_max)\
                           + ' reflns.B_iso_Wilson_estimate.comparator=between reflns.B_iso_Wilson_estimate.min=' \
                           + str(b_iso_wilson_estimate_min) + ' reflns.B_iso_Wilson_estimate.max=' \
                           + str(b_iso_wilson_estimate_max) + ' reflns.pdbx_netI_over_av_sigmaI.comparator=between' \
                           + ' reflns.pdbx_redundancy.comparator=between reflns.pdbx_redundancy.min=' \
                           + str(overall_redundancy_min) + ' reflns.pdbx_redundancy.max=' + str(overall_redundancy_max)\
                           + ' reflns.percent_possible_obs.comparator=between reflns.percent_possible_obs.min=' \
                           + str(per_of_possible_reflections_min) + ' reflns.percent_possible_obs.max=' \
                           + str(per_of_possible_reflections_max) + ' reflns.pdbx_Rmerge_I_obs.comparator=between ' \
                           + 'reflns.pdbx_Rmerge_I_obs.min=' + r_value_min + ' reflns.pdbx_Rmerge_I_obs.max=' \
                           + r_value_max + ' '

        resolution_low_comparator_query = SubElement(title, 'reflns.d_resolution_low.comparator')
        resolution_low_comparator_query.text = 'between'

        resolution_high_comparator_query = SubElement(title, 'reflns.d_resolution_high.comparator')
        resolution_high_comparator_query.text = 'between'

        d_resolution_high_min_query = SubElement(title, 'reflns.d_resolution_high.min')
        d_resolution_high_min_query.text = str(d_resolution_high_min)

        d_resolution_high_max_query = SubElement(title, 'reflns.d_resolution_high.max')
        d_resolution_high_max_query.text = str(d_resolution_high_max)

        b_iso_wilson_estimate_comparator_query = SubElement(title, 'reflns.B_iso_Wilson_estimate.comparator')
        b_iso_wilson_estimate_comparator_query.text = 'between'

        b_iso_wilson_estimate_min_query = SubElement(title, 'reflns.B_iso_Wilson_estimate.min')
        b_iso_wilson_estimate_min_query.text = str(b_iso_wilson_estimate_min)

        b_iso_wilson_estimate_max_query = SubElement(title, 'reflns.B_iso_Wilson_estimate.max')
        b_iso_wilson_estimate_max_query.text = str(b_iso_wilson_estimate_max)

        overall_redundancy_comparator_query = SubElement(title, 'reflns.pdbx_redundancy.comparator')
        overall_redundancy_comparator_query.text = 'between'

        overall_redundancy_min_query = SubElement(title, 'reflns.pdbx_redundancy.min')
        overall_redundancy_min_query.text = str(overall_redundancy_min)

        overall_redundancy_max_query = SubElement(title, 'reflns.pdbx_redundancy.max')
        overall_redundancy_max_query.text = str(overall_redundancy_max)

        per_of_possible_reflections_query = SubElement(title, 'reflns.percent_possible_obs.comparator')
        per_of_possible_reflections_query.text = 'between'

        per_of_possible_reflections_min_query = SubElement(title, 'reflns.percent_possible_obs.min')
        per_of_possible_reflections_min_query.text = str(per_of_possible_reflections_min)

        per_of_possible_reflections_max_query = SubElement(title, 'reflns.percent_possible_obs.min')
        per_of_possible_reflections_max_query.text = str(per_of_possible_reflections_max)

        r_value_query = SubElement(title, 'reflns.pdbx_Rmerge_I_obs.comparator')
        r_value_query.text = 'between'

        r_value_min_query = SubElement(title, 'reflns.pdbx_Rmerge_I_obs.min')
        r_value_min_query.text = str(r_value_min)

        r_value_max_query = SubElement(title, 'reflns.pdbx_Rmerge_I_obs.max')
        r_value_max_query.text = str(r_value_max)

        self.query.append(query_refinement)
        self.counter += 1

    def cell_dimensions(self, a_length_min, a_length_max, b_length_min, b_length_max, c_length_min, c_length_max,
                        angle_alpha_min, angle_alpha_max, angle_beta_min, angle_beta_max, angle_gamma_min,
                        angle_gamma_max):
        """Builds a query for Methods search purposes with specification Cell Dimensions.

        @param a_length_min: minimal a length (angstroms) (e.g. 1).
        @type a_length_min: str

        @param a_length_max: maximal a length (angstroms) (e.g. 3).
        @type a_length_max: str

        @param b_length_min: minimal b length (angstroms) (e.g. 1).
        @type b_length_min: str

        @param b_length_max: maximal b length (angstroms) (e.g. 4).
        @type b_length_max: str

        @param c_length_min: minimal c length (angstroms) (e.g. 1).
        @type c_length_min: str

        @param c_length_max: maximal c length (angstroms) (e.g. 5).
        @type c_length_max: str

        @param angle_alpha_min: minimal angle alpha (degrees) (e.g. 90).
        @type angle_alpha_min: str

        @param angle_alpha_max: maximal angle alpha (degrees) (e.g. 140).
        @type angle_alpha_max: str

        @param angle_beta_min: minimal angle beta (degrees) (e.g. 90).
        @type angle_beta_min: str

        @param angle_beta_max: maximal angle beta (degrees) (e.g. 130).
        @type angle_beta_max: str

        @param angle_gamma_min: minimal angle gamma (degrees) (e.g. 90).
        @type angle_gamma_min: str

        @param angle_gamma_max: maximal angle gamma (degrees) (e.g. 120).
        @type angle_gamma_max: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.EnzymeClassificationQuery'

        description = SubElement(title, 'description')
        description.text = 'XrayCellQuery: cell.length_a.comparator=between cell.length_a.min=' + str(a_length_min) \
                           + ' cell.length_a.max=' + str(a_length_max) \
                           + ' cell.length_b.comparator=between cell.length_b.min=' + str(b_length_min) \
                           + ' cell.length_b.max=' + str(b_length_max) \
                           + ' cell.length_c.comparator=between cell.length_c.min=' + str(c_length_min) \
                           + ' cell.length_c.max=' + str(c_length_max) + ' cell.angle_alpha.comparator=between' \
                           + ' cell.angle_alpha.min=' + str(angle_alpha_min) + ' cell.angle_alpha.max=' \
                           + str(angle_alpha_max) + ' cell.angle_beta.comparator=between cell.angle_beta.min=' \
                           + str(angle_beta_min) + ' cell.angle_beta.max=' + str(angle_beta_max) \
                           + ' cell.angle_gamma.comparator=between cell.angle_gamma.min=' + str(angle_gamma_min) \
                           + ' cell.angle_gamma.max=' + str(angle_gamma_max) + ' '

        a_length_comparator_query = SubElement(title, 'cell.length_a.comparator')
        a_length_comparator_query.text = 'between'

        a_length_min_query = SubElement(title, 'cell.length_a.min')
        a_length_min_query.text = str(a_length_min)

        a_length_max_query = SubElement(title, 'cell.length_a.max')
        a_length_max_query.text = str(a_length_max)

        b_length_comparator_query = SubElement(title, 'cell.length_b.comparator')
        b_length_comparator_query.text = 'between'

        b_length_min_query = SubElement(title, 'cell.length_b.min')
        b_length_min_query.text = str(b_length_min)

        b_length_max_query = SubElement(title, 'cell.length_b.max')
        b_length_max_query.text = str(b_length_max)

        c_length_comparator_query = SubElement(title, 'cell.length_c.comparator')
        c_length_comparator_query.text = 'between'

        c_length_min_query = SubElement(title, 'cell.length_c.min')
        c_length_min_query.text = str(c_length_min)

        c_length_max_query = SubElement(title, 'cell.length_c.max')
        c_length_max_query.text = str(c_length_max)

        angle_alpha_comparator_query = SubElement(title, 'cell.angle_alpha.comparator')
        angle_alpha_comparator_query.text = 'between'

        angle_alpha_min_query = SubElement(title, 'cell.angle_alpha.min')
        angle_alpha_min_query.text = str(angle_alpha_min)

        angle_alpha_max_query = SubElement(title, 'cell.angle_alpha.max')
        angle_alpha_max_query.text = str(angle_alpha_max)

        angle_beta_comparator_query = SubElement(title, 'cell.angle_beta.comparator')
        angle_beta_comparator_query.text = 'between'

        angle_beta_min_query = SubElement(title, 'cell.angle_beta.min')
        angle_beta_min_query.text = str(angle_beta_min)

        angle_beta_max_query = SubElement(title, 'cell.angle_beta.max')
        angle_beta_max_query.text = str(angle_beta_max)

        angle_gamma_comparator_query = SubElement(title, 'cell.angle_gamma.comparator')
        angle_gamma_comparator_query.text = 'between'

        angle_gamma_min_query = SubElement(title, 'cell.angle_gamma.min')
        angle_gamma_min_query.text = str(angle_gamma_min)

        angle_gamma_max_query = SubElement(title, 'cell.angle_gamma.max')
        angle_gamma_max_query.text = str(angle_gamma_max)

        self.query.append(query_refinement)
        self.counter += 1

    def software(self, name_comparator, name):
        """Builds a query for Methods search purposes with specification Software.

        @param name_comparator: Comparator.
        @type name_comparator: StructComparator

        @param name_to_determine: name of software programs used in various structures (e.g. Xengen).
        @type name_to_determine: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.SoftwareQuery'

        description = SubElement(title, 'description')
        description.text = 'SoftwareQuery: software.name.comparator=' + name_comparator.value + ' software.name.value=' \
                           + str(name)

        name_comparator_query = SubElement(title, 'vsoftware.name.comparator')
        name_comparator_query.text = name_comparator.value

        name_query = SubElement(title, 'vsoftware.name.value')
        name_query.text = str(name)

        self.query.append(query_refinement)
        self.counter += 1

    def space_group(self, space_group):
        """Builds a query for Methods search purposes with specification Space Group.

        @param space_group: Comparator.
        @type space_group: SpaceGroup
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.SpaceGroupQuery'

        description = SubElement(title, 'description')
        description.text = 'SpaceGroupQuery: symmetry.space_group_name_H-M.value=' + space_group.value

        space_group_query = SubElement(title, 'symmetry.space_group_name_H-M.value')
        space_group_query.text = space_group.value

        self.query.append(query_refinement)
        self.counter += 1

    def crystal_properties(self, density_metthews_min, density_matthews_max, density_solvent_min, density_solvent_max,
                           temperature_min, temperature_max, ph_min, ph_max):
        """Builds a query for Methods search purposes with specification Crystal Properties.

        @param density_metthews_min: minimal density (Matthews) (e.g. 1).
        @type density_metthews_min: str

        @param density_matthews_max: maximal density (Matthews) (e.g. 2).
        @type density_matthews_max: str

        @param density_solvent_min: minimal density (% solvent) (e.g. 33).
        @type density_solvent_min: str

        @param density_solvent_max: maximal density (% solvent) (e.g. 88).
        @type density_solvent_max: str

        @param temperature_min: minimal temperature (K) (e.g. 12).
        @type temperature_min: str

        @param temperature_max: maximal temperature (K) (e.g. 55).
        @type temperature_max: str

        @param ph_min: minimal pH (e.g. 4.3).
        @type ph_min: str

        @param ph_max: maximal pH (e.g. 7.4).
        @type ph_max: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.CrystalQuery'

        description = SubElement(title, 'description')
        description.text = 'CrystalQuery: exptl_crystal.density_Matthews.comparator=between' \
                           + ' exptl_crystal.density_Matthews.min=' + str(density_metthews_min) \
                           + ' exptl_crystal.density_Matthews.max=' + str(density_matthews_max) \
                           + ' exptl_crystal.density_percent_sol.comparator=between' \
                           + ' exptl_crystal.density_percent_sol.min=' + str(density_solvent_min) \
                           + ' exptl_crystal.density_percent_sol.max=' + str(density_solvent_max) \
                           + ' exptl_crystal_grow.temp.comparator=between exptl_crystal_grow.temp.min=' \
                           + str(temperature_min) + ' exptl_crystal_grow.temp.max=' + str(temperature_max) \
                           + ' exptl_crystal_grow.pH.comparator=between exptl_crystal_grow.pH.min=' + str(ph_min) \
                           + ' exptl_crystal_grow.pH.max=' + str(ph_max) + ' '

        density_matthews_query = SubElement(title, 'exptl_crystal.density_Matthews.comparator')
        density_matthews_query.text = 'between'

        density_metthews_min_query = SubElement(title, 'exptl_crystal.density_Matthews.min')
        density_metthews_min_query.text = str(density_metthews_min)

        density_matthews_max_query = SubElement(title, 'exptl_crystal.density_Matthews.max')
        density_matthews_max_query.text = str(density_matthews_max)

        density_solvent_query = SubElement(title, 'exptl_crystal.density_percent_sol.comparator')
        density_solvent_query.text = 'between'

        density_solvent_min_query = SubElement(title, 'exptl_crystal.density_percent_sol.min')
        density_solvent_min_query.text = str(density_solvent_min)

        density_solvent_max_query = SubElement(title, 'exptl_crystal.density_percent_sol.max')
        density_solvent_max_query.text = str(density_solvent_max)

        temperature_query = SubElement(title, 'exptl_crystal_grow.temp.comparator')
        temperature_query.text = 'between'

        temperature_min_query = SubElement(title, 'exptl_crystal_grow.temp.min')
        temperature_min_query.text = str(temperature_min)

        temperature_max_query = SubElement(title, 'exptl_crystal_grow.temp.max')
        temperature_max_query.text = str(temperature_max)

        ph_query = SubElement(title, 'exptl_crystal_grow.pH.comparator')
        ph_query.text = 'between'

        ph_min_query = SubElement(title, 'exptl_crystal_grow.pH.min')
        ph_min_query.text = str(ph_min)

        ph_max_query = SubElement(title, 'exptl_crystal_grow.pH.max')
        ph_max_query.text = str(ph_max)

        self.query.append(query_refinement)
        self.counter += 1

    def em_assembly(self, em_assembly_comparator):
        """Builds a query for Methods search purposes with specification EM Assembly.

        @param em_assembly_comparator: Comparator.
        @type em_assembly_comparator: EMAssembly
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.EmAssemblyQuery'

        description = SubElement(title, 'description')
        description.text = 'EmAssemblyQuery: em_assembly.aggregation_state.value=' + em_assembly_comparator.value

        em_assembly_comparator_query = SubElement(title, 'em_assembly.aggregation_state.value')
        em_assembly_comparator_query.text = em_assembly_comparator.value

        self.query.append(query_refinement)
        self.counter += 1

    def detector(self, detector_type):
        """Builds a query for Methods search purposes with specification EM Assembly.

        @param detector_type: Comparator.
        @type detector_type: Detector
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.DetectorTypeQuery'

        description = SubElement(title, 'description')
        description.text = 'DetectorTypeQuery: diffrn_detector.type.value=' + detector_type.value + ' '

        detector_type_query = SubElement(title, 'diffrn_detector.type.value')
        detector_type_query.text = detector_type.value

        self.query.append(query_refinement)
        self.counter += 1

    # Publication
    def citation(self, author_name, title_comparator, citation_title, primary_citation_only, journal,
                 year_comparator, year):
        """Builds a query for Methods search purposes with specification EM Assembly.

        @param author_name: author name (e.g. Macek, P.).
        @type author_name: str

        @param title_comparator: Comparator.
        @type title_comparator: StructComparator

        @param citation_title: citation title (e.g. Structure).
        @type citation_title: str

        @param primary_citation_only: primary citation only (must equal 'primary' or '').
        @type primary_citation_only: str

        @param journal: Comparator.
        @type journal:

        @param year_comparator: Comparator.
        @type year_comparator:

        @param year: year - yyyy (e.g. 2009).
        @type year:
        """

        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.CitationTwoQuery'

        description = SubElement(title, 'description')
        description.text = 'CitationTwoQuery: ' + 'citation_author.citation_id.value=' + str(primary_citation_only) \
                           + 'citation.title.comparator=' + title_comparator.value + ' citation.id.value=primary' \
                           + ' citation.journal_abbrev.comparator=like citation.journal_abbrev.value=' \
                           + journal.value + ' citation.year.comparator=' + year_comparator.value \
                           + ' citation.year.value=' + str(year) + ' '

        author_name_query = SubElement(title, 'citation_author.name.value')
        author_name_query.text = str(author_name)

        title_comparator_query = SubElement(title, 'citation.title.comparator')
        title_comparator_query.text = title_comparator.value

        citation_title_query = SubElement(title, 'citation.title.value')
        citation_title_query.text = str(citation_title)

        primary_citation_only_query = SubElement(title, 'citation.id.value')
        primary_citation_only_query.text = primary_citation_only

        journal_comparator_query = SubElement(title, 'citation.journal_abbrev.comparator')
        journal_comparator_query.text = 'like'

        journal_query = SubElement(title, 'citation.journal_abbrev.value')
        journal_query.text = journal.value

        year_comparator_query = SubElement(title, 'citation.year.comparator')
        year_comparator_query.text = year_comparator.value

        year_query = SubElement(title, 'citation.year.value')
        year_query.text = str(year)

        self.query.append(query_refinement)
        self.counter += 1

    def medical_subject_headings_browser(self):  # TODO: not implemented yet. Will need communication to the external browser.
        print("This function is not supported yet")
        return None

    def pubmed_abstract(self, input_pubmed):
        """Builds a query for Methods search purposes with specification PubMed Abstract.

        @param input_pubmed: PubMed Abstract (e.g. ACTIN).
        @type input_pubmed: str
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.PubMedQuery'

        description = SubElement(title, 'description')
        description.text = 'PubMed Abstract Text Search for: ' + str(input_pubmed)

        input_pubmed_query = SubElement(title, 'inputPubMed_Value')
        input_pubmed_query.text = str(input_pubmed)

        self.query.append(query_refinement)
        self.counter += 1

    # Misc
    def has_external_links(self, has_external_link):
        """Builds a query for Methods search purposes with specification PubMed Abstract.

        @param has_external_link: Comparator..
        @type has_external_link: ExternalLink
        """
        query_refinement = Element('queryRefinement')
        query_refinement_level = SubElement(query_refinement, 'queryRefinementLevel')
        query_refinement_level.text = str(self.counter)

        title = SubElement(query_refinement, 'orgPdbQuery')

        query_type = SubElement(title, 'queryType')
        query_type.text = 'org.pdb.query.ExternalLinkQuery'

        description = SubElement(title, 'description')
        description.text = 'External Link Search : External Link=' + has_external_link.value

        external_link_query = SubElement(title, 'externalLink')
        external_link_query.text = has_external_link.value

        self.query.append(query_refinement)
        self.counter += 1

    # Perform Advanced Search - send to server function
    def send(self):
        """Sends request to a sever.
        """
        headers = {'Content-Type':'application/x-www-form-urlencoded'}
        url = 'http://www.rcsb.org/pdb/rest/search'
        req = requests.post(url, headers=headers, data=tostring(self.query).decode())
        self.result = req.text.split('\n')
        self.result = self.result[:-1]

    # Clear all query and results
    def clear(self):
        """
        Clear query and results.
        """
        self.query.clear()
        self.result = {}
        self.counter = 0

if __name__ == '__main__':
    # TODO call search from the command line
    pass
    '''
    import sys
    doc = """AdvancedSearch.py
    Usage:
    AdvancedSearch.py send() - send your query to search server.
    """
    print(doc)

    a = AdvancedSearch()
    a.structure_title(StructComparator.Starts_with, 'DNA')
    a.number_of_chains_ba('5', '6')
    a.send()
    result = a.result
    '''

    # TODO search history management
    # TODO unit tests