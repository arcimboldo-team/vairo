ID_TO_HHBLITS_AA_3LETTER_CODE = {0: 'ALA', 1: 'CYS', 2: 'ASP', 3: 'GLU', 4: 'PHE', 5: 'GLY', 6: 'HIS',
                                 7: 'ILE', 8: 'LYS', 9: 'LEU', 10: 'MET', 11: 'ASN', 12: 'PRO', 13: 'GLN',
                                 14: 'ARG', 15: 'SER', 16: 'THR', 17: 'VAL', 18: 'TRP', 19: 'TYR', 20: 'X',
                                 21: '-'}



ATOM_TYPES = ['N', 'CA', 'C', 'CB', 'O', 'CG', 'CG1', 'CG2', 'OG', 'OG1', 'SG', 'CD', 'CD1', 'CD2',
              'ND1', 'ND2', 'OD1', 'OD2', 'SD', 'CE', 'CE1', 'CE2', 'CE3', 'NE', 'NE1', 'NE2', 'OE1',
              'OE2', 'CH2', 'NH1', 'NH2', 'OH', 'CZ', 'CZ2', 'CZ3', 'NZ', 'OXT']

ATOM_ORDER = {atom_type: i for i, atom_type in enumerate(ATOM_TYPES)}
ORDER_ATOM = {v: k for k, v in ATOM_ORDER.items()}

RAMACHANDRAN_TABLE = \
    [
        [16, 45, 64, 51, 40, 39, 32, 38, 43, 63, 35, 4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 16],
        [1, 10, 23, 25, 18, 16, 12, 15, 19, 21, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 1],
        [0, 1, 7, 5, 2, 7, 3, 5, 3, 3, 4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 3, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [1, 5, 2, 1, 3, 3, 4, 4, 3, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 11, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1],
        [1, 1, 2, 0, 1, 2, 8, 2, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 2, 0, 3, 2, 5, 8, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 6, 11, 6, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0],
        [0, 0, 0, 0, 1, 1, 3, 4, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 4, 7, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 1, 3, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 2, 0, 2, 4, 2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 2, 0, 2, 4, 2, 6, 5, 2, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [1, 0, 0, 0, 1, 7, 5, 2, 2, 4, 1, 2, 0, 3, 2, 0, 2, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1],
        [0, 0, 4, 2, 5, 5, 7, 6, 7, 9, 6, 4, 4, 5, 7, 1, 2, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 5, 4, 5, 13, 22, 12, 16, 22, 30, 49, 52, 31, 18, 7, 2, 0, 0, 0, 0, 0, 0, 0, 0, 6, 1, 0, 0, 0, 1, 1, 0, 1,
         0, 0, 0],
        [1, 1, 1, 4, 10, 13, 38, 34, 49, 72, 160, 709, 1013, 205, 45, 1, 2, 0, 1, 0, 1, 0, 0, 0, 4, 7, 1, 1, 1, 0, 0, 0,
         0, 0, 1, 0, 1],
        [0, 0, 0, 1, 10, 22, 42, 67, 101, 229, 1035, 6328, 3530, 357, 18, 6, 2, 0, 0, 1, 0, 1, 0, 0, 6, 2, 2, 1, 0, 2,
         0, 1, 2, 0, 0, 1, 0],
        [0, 0, 1, 3, 8, 25, 51, 91, 172, 309, 1589, 4895, 1978, 146, 10, 1, 0, 0, 0, 0, 0, 0, 0, 1, 7, 5, 1, 0, 1, 0, 0,
         1, 1, 1, 0, 1, 0],
        [1, 0, 2, 4, 19, 46, 97, 125, 199, 387, 966, 1888, 714, 29, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 1, 0, 0, 1, 0,
         0, 1, 0, 0, 0, 1],
        [0, 0, 2, 7, 23, 78, 145, 180, 298, 521, 872, 873, 140, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 3, 1, 1, 0, 0, 1,
         1, 0, 0, 0, 0],
        [1, 0, 1, 9, 25, 93, 127, 251, 438, 625, 517, 172, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 11, 6, 3, 0, 0, 0, 1,
         0, 0, 0, 0, 1],
        [0, 0, 3, 17, 38, 117, 167, 379, 456, 374, 124, 12, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 24, 31, 6, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0],
        [0, 2, 6, 14, 42, 105, 205, 282, 247, 96, 20, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 76, 40, 4, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0],
        [0, 0, 4, 12, 32, 78, 143, 133, 62, 18, 11, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 89, 109, 19, 3, 1, 0, 0, 0, 0,
         0, 0, 0, 1, 0],
        [0, 1, 4, 14, 30, 56, 59, 30, 15, 13, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 31, 206, 129, 10, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0],
        [0, 1, 1, 13, 18, 30, 27, 14, 13, 13, 7, 3, 1, 0, 0, 1, 0, 0, 0, 1, 0, 7, 111, 263, 74, 2, 3, 1, 0, 0, 0, 0, 0,
         0, 0, 0, 0],
        [0, 0, 5, 12, 28, 25, 7, 9, 13, 40, 21, 3, 2, 0, 0, 0, 0, 0, 0, 0, 1, 10, 84, 103, 17, 0, 3, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0],
        [0, 4, 12, 27, 41, 47, 15, 5, 39, 65, 32, 1, 1, 0, 0, 0, 0, 0, 0, 0, 3, 10, 33, 29, 4, 1, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 1, 0],
        [2, 1, 14, 30, 82, 59, 16, 15, 53, 106, 43, 1, 1, 0, 2, 0, 0, 0, 0, 0, 1, 5, 4, 9, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0,
         1, 0, 2],
        [1, 9, 15, 28, 59, 69, 44, 28, 57, 104, 46, 4, 0, 0, 0, 2, 0, 0, 0, 0, 1, 2, 3, 1, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 1],
        [1, 7, 19, 35, 46, 82, 85, 73, 100, 115, 36, 10, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 2,
         0, 0, 1, 1],
        [3, 13, 27, 58, 114, 179, 237, 264, 239, 180, 85, 26, 5, 4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 1, 1, 0, 1, 0,
         0, 0, 1, 0, 0, 3],
        [3, 15, 27, 101, 219, 360, 458, 528, 443, 331, 206, 120, 58, 22, 7, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 2, 1, 1, 0,
         0, 0, 1, 0, 0, 0, 3],
        [4, 27, 69, 202, 397, 798, 940, 762, 632, 499, 436, 442, 336, 110, 15, 4, 0, 1, 0, 0, 0, 0, 0, 1, 2, 0, 1, 2, 0,
         0, 0, 0, 0, 0, 0, 0, 4],
        [12, 56, 155, 293, 617, 883, 852, 717, 565, 559, 668, 801, 587, 115, 2, 2, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0,
         0, 0, 0, 0, 1, 0, 0, 1, 12],
        [13, 94, 231, 462, 690, 734, 567, 480, 384, 457, 685, 856, 434, 36, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0,
         0, 0, 0, 0, 0, 0, 4, 13],
        [20, 188, 423, 486, 737, 612, 440, 300, 303, 401, 594, 549, 150, 14, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0,
         0, 1, 0, 0, 0, 0, 1, 2, 20],
        [57, 189, 362, 392, 468, 342, 239, 201, 199, 313, 394, 209, 21, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 5, 57],
        [32, 100, 162, 169, 160, 118, 99, 94, 102, 117, 137, 43, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 2, 1, 0, 0, 0,
         0, 0, 1, 0, 0, 2, 32],
        [16, 45, 64, 51, 40, 39, 32, 38, 43, 63, 35, 4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 16]
    ]

INPUT_PARAMETERS = {
    'global_input': {
        'mode': {
            'required': True,
            'mode': ['naive', 'guided']
        },
        'output_dir': {
            'required': True,
            'mode': ['naive', 'guided']
        },
        'af2_dbs_path': {
            'required': True,
            'mode': ['naive', 'guided']
        },
        'run_dir': {
            'required': False,
            'mode': ['naive', 'guided'],
            'default': None
        },
        'small_bfd': {
            'required': False,
            'mode': ['naive'],
            'default': False
        },
        'glycines': {
            'required': False,
            'mode': ['naive', 'guided'],
            'default': 50
        },
        'run_af2': {
            'required': False,
            'mode': ['naive', 'guided'],
            'default': True
        },
        'stop_after_msa': {
            'required': False,
            'mode': ['naive'],
            'default': False
        },
        'reference': {
            'required': False,
            'mode': ['guided'],
            'default': None
        },
        'experimental_pdbs': {
            'required': False,
            'mode': ['naive', 'guided'],
            'default': ''
        },
        'mosaic': {
            'required': False,
            'mode': ['naive', 'guided'],
            'default': 1
        },
        'mosaic_partition': {
            'required': False,
            'mode': ['naive', 'guided'],
            'default': []
        },
        'mosaic_seq_partition': {
            'required': False,
            'mode': ['naive', 'guided'],
            'default': []
        },
        'cluster_templates': {
            'required': False,
            'mode': ['naive', 'guided'],
            'default': False
        },
        'cluster_templates_msa': {
            'required': False,
            'mode': ['naive'],
            'default': -1
        },
        'cluster_templates_msa_mask': {
            'required': False,
            'mode': ['naive'],
            'default': ''
        },
        'cluster_templates_sequence': {
            'required': False,
            'mode': ['naive'],
            'default': None
        },
        'sequences': {
            'required': True,
            'mode': ['naive', 'guided'],
        },
        'features': {
            'required': False,
            'mode': ['guided'],
            'default': []
        },
        'templates': {
            'required': False,
            'mode': ['guided'],
            'default': []
        },
        'append_library': {
            'required': False,
            'mode': ['guided'],
            'default': []
        },
        'show_pymol': {
            'required': False,
            'mode': ['naive', 'guided'],
            'default': ''
        }
    },
    'sequence_input': {
        'fasta_path': {
            'required': True
        },
        'num_of_copies': {
            'required': False,
            'default': 1
        },
        'positions': {
            'required': False,
            'default': None
        },
        'name': {
            'required': False,
            'default': None
        },
        'mutations': {
            'required': False,
            'default': None
        },
        'predict_region': {
            'required': False,
            'default': []
        }
    },
    'template_input': {
        'pdb': {
            'required': True
        },
        'add_to_msa': {
            'required': False,
            'default': False
        },
        'add_to_templates': {
            'required': False,
            'default': True
        },
        'generate_multimer': {
            'required': False,
            'default': True
        },
        'sum_prob': {
            'required': False,
            'default': False
        },
        'strict': {
            'required': False,
            'default': True
        },
        'aligned': {
            'required': False,
            'default': False
        },
        'legacy': {
            'required': False,
            'default': False
        },
        'modifications': {
            'required': False,
            'default': []
        },
    },
    'features_input': {
        'path': {
            'required': True
        },
        'keep_msa': {
            'required': False,
            'default': -1
        },
        'keep_templates': {
            'required': False,
            'default': -1
        },
        'msa_mask': {
            'required': False,
            'default': ''
        },
        'sequence': {
            'required': False,
            'default': None
        },
        'numbering_features': {
            'required': False,
            'default': None
        },
        'numbering_query': {
            'required': False,
            'default': None
        },
        'positions': {
            'required': False,
            'default': None
        },
        'mutations': {
            'required': False,
            'default': None
        }
    },
    'append_library_input': {
      'path': {
          'required': True
      },
      'aligned': {
          'required': False,
          'default': True
      },
      'add_to_msa': {
          'required': False,
          'default': True
      },
      'add_to_templates': {
            'required': False,
            'default': False
      },
      'numbering_library': {
            'required': False,
            'default': None
      },
      'numbering_query': {
            'required': False,
            'default': None
      },
    },
    'modifications_input': {
        'chain': {
            'required': True,
        },
        'position': {
            'required': False,
            'default': -1,
        },
        'maintain_residues': {
            'required': False,
            'default': [],
        },
        'delete_residues': {
            'required': False,
            'default': [],
        },
        'when': {
            'required': False,
            'default': 'after_alignment',
        },
        'mutations': {
            'required': False,
            'default': None,
        },
    },
    'mutations_input': {
        'numbering_residues': {
            'required': True,
        },
        'mutate_with': {
            'required': True,
        },
    },
}
