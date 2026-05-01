import molgrep

import pytest
from rdkit import Chem


@pytest.fixture()
def record():
    mol = Chem.MolFromSmiles('O')
    mol.SetProp('_Name', 'water')
    mol.SetProp('x', 'wet')
    return molgrep.MolRecord(mol, 'test.smi', 42)


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


@pytest.mark.parametrize('fmt, match, count, expected', [
    ('', [], 0, ''),
    ('{file} {index} {name} {count}', [], 2, 'test.smi 42 water 2'),
    ('{smiles} {match}', [1, 2], 2, 'O 1,2'),
    ('{p_x} {does_not_exist}', [], 0, 'wet '),
])
def test_format_result(record, fmt, match, count, expected):
    got = molgrep.format_result(fmt, record, match=match, count=count)
    assert got == expected


def test_format_result_no_name(record):
    record.mol.ClearProp('_Name')
    got = molgrep.format_result('{name}', record)
    assert got == ''
