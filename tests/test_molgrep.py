import molgrep

from rdkit import Chem


def test_get_mol_records_from_smiles():
    recs = list(molgrep.get_mol_records(['C', 'CC']))

    assert len(recs) == 2

    assert recs[0].index == 1
    assert recs[0].file == 'arg1'
    assert Chem.MolToSmiles(recs[0].mol) == 'C'

    # This is the "per-file" index, which is always 1 for SMILES
    assert recs[1].index == 1

    assert recs[1].file == 'arg2'
    assert Chem.MolToSmiles(recs[1].mol) == 'CC'
