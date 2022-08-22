"""Simple annotation viewer.

Show all the annotations in a series of concordances for each file and show some
global analyses on the left context.

$ python show.py <BACKUP>

This runs on directory BACKUP inside BRAT_BACKUP, which has all Brat backups.

This works with backups where the directory has a subdirectory named EHR, which
has subdirectories named for the annotators (Ann, Mei and Phil). The annotator
directories are flat and have files like

    102607056_0.ann
    102607056_0.txt
    103902665_0.ann
    103902665_0.txt

This basically includes all backups starting from 3_28_22, although the earlier
ones (before 5_12_22) often have one text file missing. The code is not guaranteed
to work on all backup directories, even if they have the EHR subdirectory, but it
should work on the most recent ones.

"""

# TODO: figure out the problem with unit 105663850_0

import os, sys, glob
from math import log
from typing import Union
from collections import Counter

from utils import split, index_characters
from html import Tag, Text, Span, P, TR, TD, H1, H3
from html import Href, Anchor, NL, SPACE

# location of the annotation files
BRAT_BACKUP = '/Users/marc/Dropbox/Shared/NLP-Trauma Project (R21)/brat_backup'

# tags that we are interested in
EXTENT_TAGS = {'Event', 'Symptom', 'Substance', 'Temporal_Frame', 'Perpetrator'}
RELATION_TAGS = {'Grounded_To'}

# number of characters to print to the left and right
CONTEXT = 50

# settings to determine what left contexts to print in the analysis
MINIMUM_FREQUENCY = 2
MINIMUM_PMI = 5

# Set of annotators that we are interested in
ANNOTATORS = {'Ann', 'Mei', 'Phil'}


class Corpus(object):

    """Based on a directory with files from Brat annotations and relates sets of
    annotations and text files with annotation units.

    directory    -  path to the corpus directory
    name         -  basename of the directory
    fnames       -  list of files in the corpus directory
    units        -  { unit-name => AnnotationUnit }
    vocabulary   -  Vocabulary
    annotations  -  Annotations
    contexts     -  Contexts

    """

    def __init__(self, directory: str):
        self.directory = directory
        self.name = os.path.basename(self.directory)
        self.fnames = self._get_filenames()
        self.units = self._get_units()
        self.vocabulary = Vocabulary(self)
        self.annotations = Annotations()
        self.annotations.add_corpus_annotations(self)
        self.contexts = Contexts(self)

    def __str__(self):
        return ("<Corpus '%s' with %d units and %d files>"
                % (self.name, len(self.units), len(self.fnames)))

    def _get_filenames(self):
        """Return a list of filenames for the corpus, but only include those
        where the annotator in the ANNOTATORS set."""
        fnames = glob.glob("%s/EHR/*/[0-9]*_[0-9]*" % self.directory)
        filtered_fnames = []
        for f in fnames:
            annotator = os.path.normpath(f).split(os.sep)[-2]
            if annotator in ANNOTATORS:
                filtered_fnames.append(f)
        return filtered_fnames

    def _get_units(self):
        units = {}
        for unit_name, file_names in self.split_files().items():
            try:
                units[unit_name] = AnnotationUnit(self, unit_name, file_names)
            except KeyError:
                print("Warning: error loading AnnotationUnit %s" % unit_name)
        return units

    def split_files(self):
        """Split the files into groups with the same file identifier. Returns
        a dictionary indexed on unit names where the values are lists of file
        names."""
        files = {}
        for fname in self.fnames:
            base = os.path.basename(fname)[:11]
            files.setdefault(base, []).append(fname)
        return files

    def get_unit_names(self):
        """Return a sorted list of unit names in the corpus."""
        return sorted(k for k in list(self.units.keys()))

    def write_reports(self):
        """Write all annotations from all units to files as well as all context
        observations. Also write an index that ties all files together."""
        html_dir = os.path.join('html', self.name)
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)
        units_written = []
        for unit_name in self.get_unit_names():
            unit = self.units[unit_name]
            ### THIS NEEDS TO BE CHANGED AFTER ALL DEBUGGING
            # if unit_name != '100370074_1':
            #     continue
            UnitFile(self, unit).write()
            units_written.append(unit)
        for tag in self.contexts.data:
            ExtentsFile(self, tag).write()
            ContextFile(self, tag).write()
        # breakpoint()
        ExtentsFile(self, 'Event', event_type='other').write()
        IndexFile(self, units_written).write()

    def print_units(self):
        """Print all units in the corpus."""
        for unit_name in sorted(self.units):
            unit = self.units[unit_name]
            print(unit)


class Vocabulary(object):

    """Keeps track of counts of tokens and bigrams as well as the PMI between all
    token pairs.

    corpus             -  the Corpus that the vocabulary is for
    tokens             -  list of all tokens from the vocabulary
    bigrams            -  list of all bigrams in the corpus
    token_counter      -  Counter of all tokens
    bigram_counter     -  Counter of all bigrams
    number_of_tokens   -  number of tokens
    number_of_types    -  number of types (size of the vocabulary)
    number_of_bigrams  -  number of bigrams
    pmis               -  { (x,y) => (P(x,y), P(x), P(y), PMI(x,y)) }

    The word pairs used to calculate the pointwise mutual information score are
    taken from the bigram counter.

    """

    def __init__(self, corpus: Corpus):
        """Initialize from a Corpus. Create the tokens and bigrams lists and
        Counters and calculate all PMIs. This requires the Corpus to have an
        instantiated dictionary of AnnotationUnits."""
        self.corpus = corpus
        self.tokens = []
        for unit in corpus.units.values():
            if unit.text is not None:
                tokens = [t.lower() for t in unit.text.split()]
                self.tokens.extend(tokens)
        self.token_counter = Counter(self.tokens)
        self.number_of_tokens = len(self.tokens)
        self.number_of_types = len(self.token_counter)
        self.bigrams = self._get_bigrams()
        self.number_of_bigrams = len(self.bigrams)
        self.bigram_counter = Counter(self.bigrams)
        self.pmis = self._calculate_pmi_scores()

    def _calculate_pmi_scores(self) -> dict:
        """Calculate the pointwise mutual information for all word pairs. This
        is restricted to word that are adjacent. For each pair <x,y> we store
        P(x,y), P(x), P(y) and PMI(x,y)."""
        pmi_scores = {}
        for pair, count in self.bigram_counter.items():
            p_x_y = count / self.number_of_bigrams
            p_x = self.token_counter[pair[0]] / self.number_of_tokens
            p_y = self.token_counter[pair[1]] / self.number_of_tokens
            pmi = log(p_x_y / (p_x * p_y))
            pmi_scores[pair] = (p_x_y, p_x, p_y, pmi)
        return pmi_scores

    def __len__(self):
        return len(self.token_counter)

    def __str__(self):
        return ("<Vocabulary on %s tokens with %d elements>"
                % (self.number_of_tokens, len(self)))
    
    def _get_bigrams(self):
        bigrams = []
        for i in range(self.number_of_tokens - 1):
            bigrams.append((self.tokens[i], self.tokens[i+1]))
        return bigrams

    def get_pmi(self, x: str, y: str) -> Union[float, None]:
        """Return PMI(x,y) for a pair of words. If the PMI of a pair is not
        defined return None."""
        _P_x_y, _p_x, _p_y, pmi = self.pmis.get((x, y), (None, None, None, None))
        return pmi

    def print_pmi(self, threshold: float, limit=None):
        """Print the PMI for all pairs as long as the PMI > threshold. Stop
        printing when the number of scores printed exceeds the limit."""
        if limit is None:
            limit = sys.maxsize
        printed = 0
        for pair, (p_x_y, p_x, p_y, pmi) in self.pmis.items():
            if printed <= limit:
                c_x_y = self.bigram_counter[pair]
                c_x = self.token_counter[pair[0]]
                c_y = self.token_counter[pair[1]]
                if pmi > threshold and c_x_y > 1:
                    print("%2d %f  %3d %f  %3d %f  %7.4f  %s"
                          % (c_x_y, p_x_y, c_x, p_x, c_y, p_y, pmi, pair))
                    printed += 1


class AnnotationUnit(object):

    """All information associated with an annotated file, including all the
    annotations and the text of the file.

    corpus       -  Corpus
    name         -  string
    fnames       -  list of files in the unit
    text         -  text of the unit (the primary source)
    char2sent    -  { character-offset => sentence-number }
    annotators   -  set of annotator names
    annotations  -  Annotations

    The name of the unit is the part of the filename that all files in the unit
    share. For example, if the files in the unit are

      .../brat_backup/5_18_22/EHR/Ann/104054131_0.ann
      .../brat_backup/5_18_22/EHR/Ann/104054131_0.txt
      .../brat_backup/5_18_22/EHR/Mei/104054131_0.ann
      .../brat_backup/5_18_22/EHR/Mei/104054131_0.txt

    Then the name is 104054131_0.

    """

    def __init__(self, corpus: Corpus, unit_name: str, file_names: list):
        """Create an AnnotationUnit from a list of files in a Corpus. The name of
        the unit was derived from the list of files."""
        print("Loading %s" % unit_name)
        self.corpus = corpus
        self.name = unit_name
        self.fnames = file_names
        self.text = self._get_text()
        self.char2sent = index_characters(split(self.text))
        self.annotators = set(self.get_annotator(fname) for fname in self.fnames)
        # Initialize an Annotations instance and then add the annotations from
        # the unit into that instance
        self.annotations = Annotations()
        self.annotations.add_unit_annotations(self)

    def __str__(self):
        return ("<AnnotationUnit %s files=%s extents=%d relations=%d>"
                % (self.name, len(self.fnames),
                   len(self.annotations.extents), len(self.annotations.relations)))

    def __getitem__(self, i: int):
        return self.fnames[i]

    def _get_text(self):
        """Pull the text out of the text file, if there are more text files pick the
        first one after sorting them on length of the file name."""
        text_files = [f for f in self.fnames if f.endswith('.txt')]
        text_files = sorted(text_files, key=lambda fname: len(fname))
        if text_files:
            return open(text_files[0]).read()
        else:
            print("Warning: no text file for %s" % self)
            return None

    def text_size(self):
        """Return the length of the text. Return 0 if the text is None."""
        return 0 if self.text is None else len(self.text)

    def get_extent(self, annotator: str, identifier: str):
        return self.annotations.get_extent(annotator, identifier)

    def get_extents(self):
        return self.annotations.get_extents()

    def get_relations(self):
        return self.annotations.get_relations()

    def get_number_of_extents(self):
        return self.annotations.number_of_extents
    
    def get_number_of_relations(self):
        return self.annotations.number_of_relations

    @staticmethod
    def get_annotator(file_name: str):
        """Get the annotator for the file. We expect file paths like
        brat_backup/5_12_22/EHR/Phil/105541561_2.ann."""
        return os.path.split(os.path.split(file_name)[0])[1]

    def get_annotation_files(self):
        """Get all files with annotations."""
        return [f for f in sorted(self.fnames) if f.endswith('ann')]

    def print_files(self):
        for fname in self.fnames:
            print(os.path.basename(fname))

    def print_annotations(self):
        for annotator, tag, offsets, text in self.annotations.extents:
            print("%s\t%s\t%s\t%s" % (annotator, tag, offsets, text))

    def print_annotations_index(self):
        for annotator in self.annotations.extents_idx.keys():
            print(annotator, end=': ')
            for identifier in self.annotations.extents_idx[annotator].keys():
                print(identifier, end=' ')
            print()


class Annotation(object):

    """All information relevant to an annotation. Abstract class that is only
    used when initializing an ExtentAnnotation or a RelationAnnotation.

    unit        -  AnnotationUnit
    text        -  unit.text
    char2sent   -  unit.char2sent
    annotator   -  annotator (string)
    identifier  -  annotation identifier (for example T154)
    tag         -  annotation tag ("Event"|"Symptom"|...)

    """
    
    def __init__(self, unit: AnnotationUnit, annotator: str, fields: list):
        """Initialization takes the unit an annotation is in, the annotator name
        and the fields from the annotation file. For extent annotations there are
        three fields and for relation annotations there are two:

            T21     Temporal_Frame   3363 3369            Recent
            R5      Grounded_To      Arg1:T19 Arg2:T21
        """
        self.unit = unit
        self.text = unit.text
        self.char2sent = unit.char2sent
        self.annotator = annotator
        self.identifier = fields[0]
        self.tag = fields[1].split()[0]


class ExtentAnnotation(Annotation):

    """Instance variables (not including the ones inherited from Annotation):

    extent         -  annotation extent (string, for example 'made threats')
    offsets        -  start and end offsets ("8284:8295", "8284:8295;8304:8309") 
    positions      -  offsets as a list ([[828,8295],[8304,8309]])
    start          -  start of first fragment (integer)
    end            -  end of last fragment (integer)
    sentence_no    -  sentence number of the first fragment (integer)
    left_context   -  N=CONTEXT characters to the left of start position
    right_context  -  N=CONTEXT characters to the right of end position
    attributes     -  dictionary of attributes: {'attr' -> (attr_id, attr_value)}

    The more complicated case for offsets and positions is for discontinuous
    annotations (that is, annotations with more than one fragment). There is
    currently no middle context, but there may be a future need for it.

    """

    def __init__(self, unit: AnnotationUnit, annotator: str, fields: list):
        super().__init__(unit, annotator, fields)
        self.extent = fields[2]
        self.offsets = ':'.join(fields[1].split()[1:])
        pairs = [pair for pair in self.offsets.split(';')]
        self.positions = [[int(pos) for pos in pair.split(':')] for pair in pairs]
        self.start = self.positions[0][0]
        self.end = self.positions[-1][-1]
        self.sentence_no = self.char2sent[self.start] + 1
        self.left_context = self.get_left_context()
        self.right_context = self.get_right_context()
        self.attributes = {}

    def __str__(self):
        return ("<%s %s %s %s '%s' %s>"
                % (self.tag, self.unit.name, self.identifier, self.offsets,
                   self.extent, self.annotator))

    def __eq__(self, other):
        # for sorting purposes, we do not care about the identifier or annotator
        return self.positions == other.positions

    def __lt__(self, other):
        return self.positions[0][0] < other.positions[0][0]

    def add_attribute(self, attribute_id, attribute, value):
        self.attributes[attribute] = (attribute_id, value)

    def get_attribute(self, name):
        val = self.attributes.get(name)
        return None if val is None else val[1]

    def get_key_string(self):
        # This just returns the string from the text
        return self.text[self.start:self.end]

    def first_token(self):
        """Return the first token of the key phrase."""
        try:
            return self.get_key_string().split()[0]
        except IndexError:
            # for those cases where we have truncated text files where the
            # annotation offsets are not in the text file
            return None

    def previous_token(self):
        """Return the token just before the key phrase."""
        try:
            return self.left_context.split()[-1]
        except IndexError:
            # for those cases where we have truncated text files where the
            # annotation offsets are not in the text file
            return None

    def get_left_context(self):
        return protect(self.text[self.start-CONTEXT:self.start])

    def get_right_context(self):
        return protect(self.text[self.end:self.end+CONTEXT])


class RelationAnnotation(Annotation):

    def __init__(self, unit: AnnotationUnit, annotator: str, fields: list):
        super().__init__(unit, annotator, fields)
        self.args = {}
        for arg in fields[1].split()[1:]:
            if arg.startswith('Arg') and ':' in arg:
                argname, value = arg.split(':')
                self.args[argname] = value
        self.arg1 = unit.get_extent(annotator, self.args.get('Arg1'))
        self.arg2 = unit.get_extent(annotator, self.args.get('Arg2'))

    def __str__(self):
        return ("<%s %s %s %s %s>"
                % (self.tag, self.identifier, self.annotator,
                   self.arg1.identifier, self.arg2.identifier))

    def locations(self):
        """Return the starting positions of the two arguments with the arguments
        ordered on text position."""
        return (self.sorted_args()[0].positions[0][0],
                self.sorted_args()[1].positions[0][0])

    def sorted_args(self):
        """Returns the arguments sorted on text position."""
        return sorted([self.arg1, self.arg2], key=lambda x: x.positions[0][0])

    def text_between_args(self):
        args = self.sorted_args()
        p1 = args[0].positions[-1][-1]
        p2 = args[1].positions[0][0]
        return self.text[p1:p2]

    def pp(self):
        print(self)
        print('   ', self.arg1)
        print('   ', self.arg2)


class Annotations(object):

    """Maintains the extent and relation annotations for an AnnotationUnit or a Corpus.

    source               -  the source of the annotations: Corpus or AnnotationUnit
    extents              -  ExtentAnnotation list
    extents_idx          -  { annotator => { identifier => ExtentAnnotation }}
    grouped_extents      -  { tag => ExtentAnnotation list }
    grouped_extents_idx  -  { tag => { extent => ExtentAnnotation list }}
    relations            -  RelationAnnotation list
    number_of_extents    -  integer, length of self.extents
    number_of_relations  -  integer, length of self.relations

    """

    def __init__(self):
        self.source = None
        self.extents = []
        self.extents_idx = {}
        self.grouped_extents = None
        self.grouped_extents_idx = None
        self.relations = []
        self.number_of_extents = 0
        self.number_of_relations = 0

    def __str__(self):
        return "<Annotations %d %d>" % (self.number_of_extents,
                                        self.number_of_relations)

    def add_unit_annotations(self, unit: AnnotationUnit):
        """Add extent and relation annotations from the annotation unit."""
        if unit.text is None:
            print('Warning: not adding annotations for', self)
        else:
            for file_name in unit.get_annotation_files():
                short_file_name = os.path.basename(file_name)
                annotator = unit.get_annotator(file_name)
                with open(file_name) as fh:
                    for line in fh:
                        fields = line.strip().split('\t')
                        if is_brat_extent(fields):
                            adjust_tag(fields)
                            a = ExtentAnnotation(unit, annotator, fields)
                            self.add_extent(a)
                        elif is_brat_attribute(fields):
                            attribute_id = fields[0]
                            rest = fields[1].split()
                            attribute = rest[0]
                            extent_id = rest[1]
                            value = True if len(rest) == 2 else rest[2]
                            extent = self.get_extent(annotator, extent_id)
                            extent.add_attribute(attribute_id, attribute, value)
                        elif is_brat_relation(fields):
                            adjust_tag(fields)
                            a = RelationAnnotation(unit, annotator, fields)
                            self.add_relation(a)
                        else:
                            print('Warning, unexpected line in %s:' % short_file_name)
                            print('---', '\t'.join(fields), '\n---', file_name)

    def add_corpus_annotations(self, corpus: Corpus):
        """Add extent and relation annotations from the corpus."""
        # populate the list of extents and create the sorted extents dictionary
        self.grouped_extents = {tag: [] for tag in EXTENT_TAGS}
        for unit in corpus.units.values():
            for extent in unit.get_extents():
                self.add_extent(extent)
                self.grouped_extents[extent.tag].append(extent)
            for relation in unit.get_relations():
                self.add_relation(relation)
        # group the extents on tag and extent
        self.grouped_extents_idx = {tag: {} for tag in EXTENT_TAGS}
        for extent in self.extents:
            # ignoring the difficult cases with discontinuous annotations
            if len(extent.positions) > 1:
                continue
            text = extent.get_key_string()
            self.grouped_extents_idx[extent.tag].setdefault(text, []).append(extent)

    def add_extent(self, a: ExtentAnnotation):
        self.extents.append(a)
        self.extents_idx.setdefault(a.annotator, {})
        self.extents_idx[a.annotator][a.identifier] = a
        self.number_of_extents += 1

    def add_relation(self, a: RelationAnnotation):
        self.relations.append(a)
        self.number_of_relations += 1

    def get_extent(self, annotator: str, identifier: str):
        return self.extents_idx[annotator][identifier]

    def get_extents(self):
        return self.extents

    def get_relations(self):
        return self.relations

    def get_grouped_extents(self, tag: str, extent=None):
        """Get all the extents of a particular tag with a particular extent. If
        no extent is given this returns a dictionary indexed on extents, if an
        extent is given this returns a list of instances of ExtentAnnotation."""
        if extent is None:
            return self.grouped_extents_idx[tag]
        else:
            return self.grouped_extents_idx[tag].get(extent, [])

    def print_grouped_extents_summary(self):
        """Print count for each tag."""
        for tag in self.grouped_extents:
            print("%5d  %s" % (len(self.grouped_extents[tag]), tag))

    def print_grouped_extents(self):
        for tag in self.grouped_extents_idx:
            for key in self.grouped_extents_idx[tag]:
                annotations = self.grouped_extents_idx[tag][key]
                print(tag, len(annotations), key)


class Contexts(object):

    """Keeps track of the contexts of annotations.

    corpus       -  Corpus
    vocabulary   -  Vocabulary (same as corpus.vocabulary)
    annotations  -  Annotations (same as corpus.annotations)
    data         -  { tag => { extent => list of (term1, term2, pmi) }}

    """

    def __init__(self, corpus: Corpus):
        """Collect data on left contexts of extents. In particular, get the
        last word of the left context and the first word of the key phrase and
        the pmi of those pairs."""
        self.corpus = corpus
        self.vocabulary = corpus.vocabulary
        self.annotations = corpus.annotations
        self.data = {tag: {} for tag in EXTENT_TAGS}
        self.pmis_found = {}
        self.pmis_missed = {}
        for tag in self.annotations.grouped_extents_idx:
            self.pmis_found[tag] = 0
            self.pmis_missed[tag] = 0
            for keyphrase in self.annotations.grouped_extents_idx[tag]:
                self.data[tag][keyphrase] = {'count': 0, 'data': {}}
                for extent in self.annotations.grouped_extents_idx[tag][keyphrase]:
                    self._add_extent(tag, keyphrase, extent)
        print('PMIS FOUND: ', self.pmis_found)
        print('PMIS MISSED:', self.pmis_missed)
        # breakpoint()
        # self.print_data()
        # self.print_significant_data()

    def _add_extent(self, tag: str, keyphrase: str, extent: ExtentAnnotation):
        t1, t2 = extent.previous_token(), extent.first_token()
        pmi = self.vocabulary.get_pmi(t1, t2)
        if None not in (t1, t2, pmi):
            d = self.data[tag][keyphrase]
            d['count'] += 1
            d['data'].setdefault((t1, t2), {'count': 0, 'pmi': pmi, 'data': []})
            d['data'][(t1, t2)]['count'] += 1
            d['data'][(t1, t2)]['data'].append(extent)
            self.pmis_found[tag] += 1
        else:
            # TODO: this happens way too often, find out why
            # print('--- MISSED PMI', extent, t1, t2, pmi)
            self.pmis_missed[tag] += 1

    def get_significant_contexts(self, tag: str):
        """Returns a subset of the contexts for the tag, only including those
        where the frequency and PMI are above a certain threshold."""
        contexts = {}
        for keyphrase in list(self.data[tag].keys()):
            count = self.data[tag][keyphrase]['count']
            contexts[keyphrase] = { 'count': count, 'data': {} }
            phrase_data = self.data[tag][keyphrase]['data']
            for pair in phrase_data:
                pair_data = phrase_data[pair]
                pair_count = pair_data['count']
                pair_pmi = pair_data['pmi']
                if pair_pmi >= MINIMUM_PMI and pair_count >= MINIMUM_FREQUENCY:
                    contexts[keyphrase]['data'][pair] = pair_data
        return contexts

    def print_data(self):
        for tag in list(self.data.keys())[:10]:
            print(tag)
            for keyphrase in list(self.data[tag].keys())[:10]:
                count = self.data[tag][keyphrase]['count']
                phrase_data = self.data[tag][keyphrase]['data']
                print('    %3d  %s' % (count, protect(keyphrase)))
                for pair in phrase_data:
                    pair_data = phrase_data[pair]
                    pair_count = pair_data['count']
                    pair_pmi = pair_data['pmi']
                    print('       %3d  %.2f  %s ' % (pair_count, pair_pmi, pair[0]))

    def print_significant_data(self):
        for tag in list(self.data.keys()):
            print(tag)
            for keyphrase in list(self.data[tag].keys()):
                count = self.data[tag][keyphrase]['count']
                phrase_data = self.data[tag][keyphrase]['data']
                collected_data = []
                for pair in phrase_data:
                    pair_data = phrase_data[pair]
                    pair_count = pair_data['count']
                    pair_pmi = pair_data['pmi']
                    if pair_pmi >= MINIMUM_PMI and pair_count >= MINIMUM_FREQUENCY:
                        collected_data.append('       %3d  %.2f  %s ' % (pair_count, pair_pmi, pair[0]))
                if collected_data:
                    print('    %3d  %s' % (count, protect(keyphrase)))
                    for c in collected_data:
                        print(c)


class ExtentKwic(object):

    """Implements a context for an ExtentAnnotation. The main task of this class
    is to generate a line in the KWIC for the extent annotation.

    annotation     -  ExtentAnnotation that the KWIC is being created for
    annotator      -  the annotator (str)
    identifier     -  identifier copied from the annotation 
    left_context   -  left context string copied from the annotation
    right_context  -  right context string copied from the annotation
    keyphrase      -  list of pairs for the keyphrase

    The keyphrase usually is something like

        [('KEY', 'emotional abuse')]

    But for discontinuous tags it is more complex since it includes a key for
    each fragment as well as the internal contexts spans between the fragments:
    
        [('KEY', 'car'), ('CONTEXT', 'behind'), ('KEY', 'hit her')]

    """
    
    def __init__(self, annotation: ExtentAnnotation):
        self.annotation = annotation
        self.annotator = annotation.annotator
        self.identifier = annotation.identifier
        self.left_context = annotation.left_context
        self.right_context = annotation.right_context
        self.keyphrase = []
        for i, (p1, p2) in enumerate(annotation.positions):
            self.keyphrase.append(('KEY', annotation.text[p1:p2]))
            # insert the text between fragments (the key-internal context)
            if i+1 < len(annotation.positions):
                text = annotation.text[p2:annotation.positions[i+1][0]]
                self.keyphrase.append(('CONTEXT', text))

    def as_html(self, add_location=False):
        """Returns self as an instance of html.TR."""
        td_source = TD(attrs={'class': 'source'}, dtrs=[Text(self.annotator)])
        span_id = Span(attrs={'class': 'tagid'}, dtrs=[Text('(%s)' % self.identifier)])
        td_left = TD(attrs={'align': 'right'}, dtrs=[Text(self.left_context)])
        td_key = self.keyphrase_as_html()
        td_key.add(Text('&nbsp;&nbsp;' + self.right_context))
        td_key.add(SPACE, span_id)
        dtrs = [td_source, td_left, td_key]
        if add_location:
            name = self.annotation.unit.name
            fname = Href("kwic-%s.html" % name, name)
            s = self.annotation.sentence_no
            p = self.annotation.start
            dtrs.insert(0, TD(dtrs=[Text('%s : %s : %s' % (fname, s, p))]))
        return TR(dtrs=dtrs)

    def keyphrase_as_html(self):
        """Return the keyphrase as an HTML <Span> tag."""
        result = TD()
        for _type, span in self.keyphrase:
            if _type == 'KEY':
                result.add(Span(attrs={'class': 'keyword'},
                                dtrs=[Text('[%s]' % protect(span))]))
            elif _type == 'CONTEXT':
                result.add(Span(dtrs=[Text(protect(span))]))
        return result


class HtmlFile(object):

    """Some shared functionality for all classes that create HTML files.

    html_obj  -  an instance of html.Tag, the top-level tag
    content   -  another instance of html.Tag, used for the main content
    fname     -  the file that the HTML instance will be written to

    """

    def __init__(self, title: str):
        self.fname = None          # will be set in subclasses
        self.content = Tag('div')
        dtrs = [Tag('head', dtrs=[self.style_tag(), NL]), NL,
                Tag('body', dtrs=[H1(title), self.content, NL]), NL]
        self.html_obj = Tag('html', dtrs=dtrs)

    def style_tag(self):
        """Return a <styles> tag for the head portion of the HTML."""
        return Tag('style', dtrs=[Text(self.styles())])

    @staticmethod
    def styles():
        """Returns all style definitions to be wrapped in a <styles> tag."""
        return (
            ".keyword { color: blue; }\n"
            ".signature { color: darkred; }\n"
            ".source { color: darkgreen; }\n"
            ".tagid { font-size: 12px; color: darkgrey; }\n"
            ".green { color: green; }\n"
            ".darkred { color: darkred; }\n"
            ".extent { font-size: 24px; }\n"
            ".count { text-align: right; width: 25px; color: darkgreen; }\n"
            ".pmi { text-align: right; width: 45px; color: darkgreen; }\n"
            "table.fancy { margin-top: 10px; padding: 5px; border: thin solid lightgray; }\n")

    def write(self):
        """Write the HTML file, using the fname variable as the destination."""
        with open(self.fname, 'w') as fh:
            fh.write(str(self.html_obj))


class IndexFile(HtmlFile):

    """Creates an HTML file for the index."""
    
    def __init__(self, corpus, units):
        title = "Analysis of annotations in %s" % corpus.name
        super().__init__(title)
        self.units = units
        self.fname = 'html/%s/index.html' % corpus.name
        self._add_extents_links()
        self._add_left_context_links()
        self._add_annotation_units_table()

    def _add_extents_links(self):
        header = "All extents by tag:\n"
        hrefs = Tag('ul')
        self.content.add(P(header), hrefs, NL)
        # the extents for each tag
        for tag in sorted(EXTENT_TAGS):
            hrefs.add(Tag('li', dtrs=[Href('extents-%s.html' % tag, tag)]))
            if tag == 'Event':
                hrefs.add(Tag('li', dtrs=[Href('extents-%s-other.html' % tag, tag + '-other')]))

    def _add_left_context_links(self):
        header = "Global observations on left context:\n"
        hrefs = Tag('ul')
        self.content.add(P(header), hrefs, NL)
        # the contexts for each tag
        for tag in sorted(EXTENT_TAGS):
            hrefs.add(Tag('li', dtrs=[Href('context-%s.html' % tag, tag)]))

    def _add_annotation_units_table(self):
        header = "Annotations per file:\n"
        table = Tag('table', attrs={'cellspacing': 0, 'cellpadding': 5, 'border': 1})
        self.content.add(P(header), NL, Tag('blockquote', dtrs=[table, NL]), NL)
        tds = [SPACE, Tag('td', dtrs=[Text('')]), NL,
               SPACE, Tag('td', dtrs=[Text('file')]), NL,
               SPACE, Tag('td', dtrs=[Text('text size')]), NL,
               SPACE, Tag('td', dtrs=[Text('extents')]), NL,
               SPACE, Tag('td', dtrs=[Text('relations')]), NL,
               SPACE, Tag('td', dtrs=[Text('annotators')]), NL]
        tr = Tag('tr', dtrs=tds)
        table.add(tr, NL)
        for i, unit in enumerate(self.units):
            name = unit.name
            annotators = ', '.join(unit.annotators)
            count = Text(str(i))
            size = Text(unit.text_size())
            extents = Text(unit.get_number_of_extents())
            relations = Text(unit.get_number_of_relations())
            tds = [SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[count]), NL,
                   SPACE, Tag('td', dtrs=[Href('kwic-%s.html' % name, name)]), NL,
                   SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[size]), NL,
                   SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[extents]), NL,
                   SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[relations]), NL,
                   SPACE, Tag('td', dtrs=[Text(annotators)]), NL]
            table.add(TR(dtrs=tds), NL)


class UnitFile(HtmlFile):

    """Creates an HTML file for an AnnotationUnit, which is basically a
   concordance of all annotations."""
    
    def __init__(self, corpus: Corpus, unit: AnnotationUnit):
        super().__init__(unit.name)
        self.unit = unit
        self.fname = 'html/%s/kwic-%s.html' % (corpus.name, self.unit.name)
        # TODO: should really get these from the Annotations on the unit
        extents = { tag: [] for tag in EXTENT_TAGS }
        for extent in sorted(self.unit.get_extents()):
            extents[extent.tag].append(extent)
        self._add_navigation(extents)
        self._add_extents(extents)
        self._add_relations()

    def _add_navigation(self, extents: dict):
        # only add those tags for which there are examples
        annotations_filtered = {k: v for k, v in extents.items() if v}
        relations = self.unit.get_relations()
        div = Tag('div', attrs={'class': 'navigation'})
        for i, tag in enumerate(sorted(annotations_filtered)):
            prefix = '\n| ' if i > 0 else '[ '
            div.add(Text(prefix), Href('#' + tag, tag + 's'))
        if relations:
            div.add(Text('\n| '), Href('#Relation', 'Relations'), NL)
        div.add(Text(']\n'))
        self.content.add(Text('\n\n'), div, NL)

    def _add_extents(self, extents: dict):
        for tag in sorted(extents):
            if extents[tag]:
                table = Tag('table', nl=True, attrs={'class': 'fancy', 'cellpadding': 5})
                self.content.add(NL, Anchor(tag), H3('%ss' % tag), NL, table)
                for a in remove_duplicates(extents[tag]):
                    table.add(ExtentKwic(a).as_html())

    def _add_relations(self):
        relations = self.unit.get_relations()
        if relations:
            table = Tag('table', nl=True, attrs={'class': 'fancy', 'cellpadding': 5})
            for rel in sorted(relations, key=lambda x: x.locations()):
                arg1_id = rel.arg1.identifier
                arg2_id = rel.arg2.identifier
                args = rel.sorted_args()
                # signature = "%s(%s, %s)" % (rel.identifier, arg1_id, arg2_id)
                signature = "%s&nbsp;&rightarrow;&nbsp;%s" % (arg1_id, arg2_id)
                td_ann = TD(attrs={'class': 'source', 'valign': 'top'},
                            dtrs=[Text(rel.annotator)])
                td_signature = TD(attrs={'class': 'signature', 'valign': 'top'},
                                  dtrs=[Text(signature)])
                td_context = self._context_as_html(rel, args)
                table.add(TR(dtrs=[td_ann, td_signature, td_context]))
            self.content.add(NL, Anchor('Relation'), H3('%ss' % 'Relation'), NL, table)

    @staticmethod
    def _context_as_html(rel, args):
        return TD(dtrs=([Text(args[0].left_context)]
                        + [green('{')]
                        + ExtentKwic(args[0]).keyphrase_as_html().dtrs
                        + [green('}'), red_subscript(args[0].identifier)]
                        + [Text(rel.text_between_args())]
                        + [green('{')]
                        + ExtentKwic(args[1]).keyphrase_as_html().dtrs
                        + [green('}'), red_subscript(args[1].identifier)]
                        + [Text(args[1].right_context)]))


class ExtentsFile(HtmlFile):

    """Used to print all extents of a particular tag with all its occurrences."""

    def __init__(self, corpus: Corpus, tag: str, event_type=None):
        super().__init__("%s instances" % tag)
        self.tag = tag
        self.corpus = corpus
        if event_type is None:
            self.fname = 'html/%s/extents-%s.html' % (corpus.name, tag)
        else:
            self.fname = 'html/%s/extents-%s-%s.html' % (corpus.name, tag, event_type)
        annotations_dict = self.corpus.annotations.get_grouped_extents(tag)
        for extent in sorted(annotations_dict):
            annotations = annotations_dict[extent]
            if tag == 'Event' and event_type == 'other':
                event_types = [a.get_attribute('Event_Type') for a in annotations]
                if 'Other' not in event_types:
                    continue
            s = "&lt;%s>" % extent
            self.content.add(P(s, nl=True, attrs={'class': 'extent'}))
            self._add_extents(annotations)

    def _add_extents(self, annotations: list):
        table = Tag('table', nl=True, attrs={'cellpadding': 8})
        self.content.add(Tag('blockquote', dtrs=[table]))
        for extent_annotation in annotations:
            table.add(ExtentKwic(extent_annotation).as_html(add_location=True))


class ContextFile(HtmlFile):

    """Contains all that is needed to print the context observations for a
    particular tag.

    tag          -  "Event" | "Perpetrator" | "Substance" | "Symptom" | "Temporal_Frame"
    fname        -  html/<CORPUS_NAME>/context-<tag>.html
    corpus       -  Corpus
    contexts     -  Contexts (same as corpus.contexts)
    annotations  -  Annotations (same as corpus.annotations)
    extents      -  { extent: str => [ ExtentAnnotation, ... ] }

    """

    def __init__(self, corpus: Corpus, tag: str):
        super().__init__("%s left context observations" % tag)
        self.tag = tag
        self.fname = 'html/%s/context-%s.html' % (corpus.name, tag)
        self.corpus = corpus
        self.contexts = corpus.contexts
        self.annotations = corpus.annotations
        self.extents = self.annotations.get_grouped_extents(tag)
        contexts_for_tag = self.contexts.get_significant_contexts(tag)
        # print('>>>', contexts_for_tag.keys())
        # print('|||', contexts_for_tag.get('bullied'))
        for extent in sorted(contexts_for_tag):
            self.contexts_for_extent = contexts_for_tag[extent]
            if self.contexts_for_extent['data']:
                total_count = str(self.contexts_for_extent['count'])
                for pair, pair_data in self.contexts_for_extent['data'].items():
                    self._add_summary(total_count, pair, pair_data, extent)
                    self._add_contexts(pair_data)

    def _add_summary(self, total_count, pair, pair_data, extent):
        table = Tag('table', nl=True, attrs={'class': 'fancy', 'cellpadding': 8})
        self.content.add(table)
        pmi = "%.2f" % pair_data['pmi']
        count = str(pair_data['count'])
        text_left = Span(dtrs=[Text(protect(pair[0]))])
        text_key = Span(class_='keyword', dtrs=[Text(protect('[%s]' % extent))])
        wider_extent = '%s %s' % (pair[0], extent)
        wider_extent_annos = self.extents.get(wider_extent, [])
        wider_extent_count = len(wider_extent_annos)
        table.add(TR(dtrs=[TD(attrs=self.attrs(25), dtrs=[Text(total_count)]),
                           TD(attrs=self.attrs(20), dtrs=[Text(count)]),
                           TD(attrs=self.attrs(20), dtrs=[Text(wider_extent_count)]),
                           TD(attrs=self.attrs(40), dtrs=[Text(pmi)]),
                           TD(dtrs=[text_left, SPACE, text_key])]))

    def _add_contexts(self, pair_data: dict):
        table = Tag('table', nl=True, attrs={'cellpadding': 8})
        self.content.add(Tag('blockquote', dtrs=[table]))
        lextents = pair_data['data']
        for lextent in lextents:
            table.add(ExtentKwic(lextent).as_html(add_location=True))

    @staticmethod
    def attrs(width):
        return {'width': width, 'align': 'right'}


def adjust_tag(fields: list) -> bool:
    """Adjust the tag in the fields by stripping a trailing _red or _yellow
    suffix. This is done because for adjudication we have introduced variants
    of tags so we can distinguish between annotations of two annotators. This
    is only relevant for extents and relations, not for attributes. Returns
    True or False depending on whether the fields argument was adjusted."""
    (tag, *rest) = fields[1].split(' ')
    adjusted = False
    if '_red' in tag:
        adjusted = True
        tag = tag[:-4]
    elif '_yellow' in tag:
        adjusted = True
        tag = tag[:-7]
    if adjusted:
        fields[1] = "%s %s" % (tag, ' '.join(rest))
    return adjusted


def is_brat_extent(fields: list) -> bool:
    return fields[0].startswith('T') and len(fields) >= 3


def is_brat_attribute(fields: list) -> bool:
    return fields[0].startswith('A') and len(fields) >= 2


def is_brat_relation(fields: list) -> bool:
    return fields[0].startswith('R') and len(fields) >= 2


def remove_duplicates(annotations: list) -> list:
    """Remove the duplicates from a list of annotations. This is done for the
    pre-annotations which appear to all have multiple instances with the same
    annotator, tag and offset, but different identifiers."""
    seen = set()
    filtered_annotations = []
    for anno in annotations:
        tup = (anno.annotator, anno.tag, anno.offsets)
        if tup not in seen:
            seen.add(tup)
            filtered_annotations.append(anno)
    return filtered_annotations


def protect(text: str) -> str:
    """Replace new line with a space and replace & and < with html entities."""
    return text.replace('\n', ' ').replace('&', '&amp;').replace('<', '&lt;')


def green(text: str) -> Tag:
    """Return an XML <span> tag with class=green."""
    return Tag('span', attrs={'class': 'green'}, dtrs=[Text(text)])


def green_subscript(text: str) -> Tag:
    """Return an XML <sub> tag with class=green."""
    return Tag('sub', attrs={'class': 'green'}, dtrs=[Text(text)])


def red_subscript(text: str) -> Tag:
    """Return an XML <sub> tag with class=darkred."""
    return Tag('sub', attrs={'class': 'darkred'}, dtrs=[Text(text)])


def run():
    """The main function to create all HTML files."""
    corpus = Corpus(os.path.join(BRAT_BACKUP, sys.argv[1]))
    corpus.write_reports()


if __name__ == '__main__':
    
    run()
