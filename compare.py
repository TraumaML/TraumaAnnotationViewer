"""Comparing annotation directories.

This is only for the directory in

    NLP-Trauma Project (R21)/brat_backup/8_2_22

Results are written to the "out" directory.

"""

import os
import typing
import brat

# data directories
DATA = '/Users/marc/Dropbox/Shared/NLP-Trauma Project (R21)/brat_backup/8_2_22'
COMBINED_DATA = os.path.join(DATA, 'combined')
ADJUDICATED_DATA = os.path.join(DATA, 'adjudicated')
BRAT_ADJUDICATED_DATA = os.path.join(DATA, 'brat_adjudicated')


def compare(dir1: str, dir2: str):
    files1 = [f for f in os.listdir(dir1) if f.endswith('.ann')]
    files2 = [f for f in os.listdir(dir2) if f.endswith('.ann')]
    files = {}
    for f in files1:
        files.setdefault(f, ['', ''])[0] = os.path.join(dir1, f)
    for f in files2:
        files.setdefault(f, ['', ''])[1] = os.path.join(dir2, f)
    write_file_list(dir1, dir2, files)
    compare_files(dir1, dir2, files)


def write_file_list(dir1: str, dir2: str, files: dict):
    dir1 = os.path.basename(dir1)
    dir2 = os.path.basename(dir2)
    with open(f'out/files-{dir1}-{dir2}.tab', 'w') as fh:
        fh.write(f'column1: {dir1}\n')
        fh.write(f'column2: {dir2}\n\n')
        files_in_dir1 = 0
        files_in_dir2 = 0
        for f in sorted(files):
            dir1_has_file = 1 if files[f][0] else 0
            dir2_has_file = 1 if files[f][1] else 0
            if dir1_has_file:
                files_in_dir1 += 1
            if dir2_has_file:
                files_in_dir2 += 1
            fh.write(f'{os.path.basename(f)}  {dir1_has_file}  {dir2_has_file}\n')
        fh.write(f'\nfiles in {dir1} --> {files_in_dir1}')
        fh.write(f'\nfiles in {dir2} --> {files_in_dir2}\n')


def compare_files(dir1: str, dir2: str, files: dict):
    d1 = os.path.basename(dir1)
    d2 = os.path.basename(dir2)
    with open(f'out/extents-{d1}-{d2}.txt', 'w') as fh:
        print(type(fh), fh)
        for fname in sorted(files):
            f1, f2 = files[fname]
            if f1 and f2:
                compare_annotations(fh, files[fname])


def compare_annotations(fh: typing.TextIO, files: dict):
    file1, file2 = files
    annotations1 = brat.BratAnnotations(file1)
    annotations2 = brat.BratAnnotations(file2)
    print(annotations1)
    print(annotations2)
    annotations1.compare(fh, annotations2)


if __name__ == '__main__':

    compare(ADJUDICATED_DATA, BRAT_ADJUDICATED_DATA)
    compare(COMBINED_DATA, ADJUDICATED_DATA)
    compare(COMBINED_DATA, BRAT_ADJUDICATED_DATA)
