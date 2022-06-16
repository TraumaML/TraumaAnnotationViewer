import os
import sys
import collections


BRAT_BACKUP = '/Users/marc/Dropbox/Shared/NLP-Trauma Project (R21)/brat_backup'

ANNOTATORS = set(['Ann', 'Mei', 'Phil'])


def check_files(directory):
    print(directory)
    for fname in os.listdir(directory):
        if not fname.endswith('.ann'):
            continue
        check_file(directory, fname)


def check_file(directory, fname):
    fpath = os.path.join(directory, fname)
    lines = open(fpath).readlines()
    identifiers = []
    for line in lines:
        fields = line.split('\t')
        identifiers.append(fields[0])
    #print('   ',  len(identifiers), len(set(identifiers)))
    if len(identifiers) != len(set(identifiers)):
        print('    WARNING: duplicate identifiers in', fname, len(identifiers), len(set(identifiers)))


if __name__ == '__main__':
    
    backup_date = sys.argv[1]
    for annotator in ANNOTATORS:
        check_files(os.path.join(BRAT_BACKUP, backup_date, 'EHR', annotator))
