"""
Find molecules matching a SMARTS. Results may be shown as a table, as images, or
written to structure files.
"""

import argparse
import itertools
import os
from dataclasses import dataclass
from typing import Any

import molcat
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

__version__ = '0.1.0'

DEFAULT_FMT = '{file} {index} {match}'
DEFAULT_COUNT_FMT = '{file} {index} {count}'


@dataclass
class MolRecord:
    mol: Chem.Mol
    file: str
    index: int


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('smarts', help='substructure to search')
    parser.add_argument('files_or_smiles', nargs='+')
    parser.add_argument('--max-per-mol',
                        '-m',
                        type=int,
                        help='maximum matches to report per structure')
    # parser.add_argument('--max-per-file', '-M', type=int,
    #     help='maximum matches to report per file')  # TODO
    parser.add_argument('--matched',
                        metavar='<filename>',
                        help='file to write structures matching the SMARTS')
    parser.add_argument('--not_matched',
                        metavar='<filename>',
                        help='file to write structure that did not match')
    parser.add_argument('--count-per-mol',
                        '-c',
                        action='store_true',
                        help='only count the number of matches per molecule')
    # parser.add_argument('--count-per-file', '-C', action='store_true',
    #     help='only count the number of matches per file')  # TODO
    parser.add_argument('--format',
                        '-f',
                        help='formatting string for text output')
    parser.add_argument(
        '--union',
        '-u',
        action='store_true',
        help='for each matching molecule, show the union of all the matches')
    parser.add_argument(
        '--image',
        '-i',
        action='store_true',
        help='show results as images (requires graphics-supporting terminal)')
    parser.add_argument('--keep-h',
                        '-H',
                        action='store_true',
                        help='keep hydrogens')
    parser.add_argument('--size-x',
                        '-x',
                        type=int,
                        help='for image output, X dimension in pixels; '
                        'default: determine automatically')
    return parser.parse_args(argv)


def get_mol_records(files_or_smiles, removeHs=True):
    for arg_idx, file_or_smiles in enumerate(files_or_smiles, 1):
        if os.path.isfile(file_or_smiles):
            reader = molcat.get_reader(file_or_smiles, removeHs=removeHs)
            for mol_idx, mol in enumerate(reader, 1):
                yield MolRecord(mol, file_or_smiles, mol_idx)
        else:
            if mol := Chem.MolFromSmiles(file_or_smiles):
                yield MolRecord(mol, f'arg{arg_idx}', index=1)


def get_png_with_match(mol: Chem.Mol, match,
                       size: tuple[int, int] = (500, 300)) -> bytes:
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    opts: Any = d.drawOptions()
    opts.addStereoAnnotation = True
    d.DrawMolecule(mol, highlightAtoms=match)
    d.FinishDrawing()
    return d.GetDrawingText()


def main():
    args = parse_args()
    size = molcat.determine_size(args.size_x)
    query = Chem.MolFromSmarts(args.smarts)
    for rec in get_mol_records(args.files_or_smiles):
        matches = rec.mol.GetSubstructMatches(query)
        try:
            name = rec.mol.GetProp('_Name')
        except KeyError:
            name = ''

        if args.union and matches:
            matches = [sorted(set(itertools.chain.from_iterable(matches)))]

        if matches and args.count_per_mol:
            fmt = args.format or DEFAULT_COUNT_FMT
            print(
                fmt.format(file=rec.file,
                           index=rec.index,
                           count=len(matches),
                           name=name))
            continue

        for match_idx, match in enumerate(matches, 1):
            fmt = args.format or DEFAULT_FMT
            print(
                fmt.format(file=rec.file,
                           index=rec.index,
                           match=match,
                           name=name))
            if args.image:
                mol2d = molcat.to_2d(rec.mol)
                png_data = get_png_with_match(mol2d, match, size)
                molcat.show_image(png_data)
            if args.max_per_mol and match_idx == args.max_per_mol:
                break
