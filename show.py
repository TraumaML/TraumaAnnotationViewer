"""show.py

Show all the annotations in a series of concordances for each file and show some
global analyses on the left context.

$ python run.py <BACKUP>

This runs on a directory with all Brat backups as defined in the BRAT_BACKUP
global variable, the argument identifies the backup we want.

This works with backups where the directory has a subdirectory named EHR, which
has subdirectories named for the annotatotors (Ann, Mei and Phil). THe annotator
directories are flat and have files like

    102607056_0.ann
    102607056_0.txt
    103902665_0.ann
    103902665_0.txt

This basically includes all backups starting from 3_28_22, although the earlier
ones (before 5_12_22) often have one text file missing.

"""


import os, sys, glob
from math import log
from collections import Counter
from html import Tag, Text, Span, P, TR, TD, H1, H3, H4
from html import Href, Anchor, NL, BR, SPACE

# location of the annotation files
BRAT_BACKUP = '/Users/marc/Dropbox/Shared/NLP-Trauma Project (R21)/brat_backup'

# tags that we are interested in
EXTENT_TAGS = {'Event', 'Symptom', 'Substance', 'Temporal_Frame'}
RELATION_TAGS = {'Grounded_To'}

# number of characters to print to the left and right
CONTEXT = 50

# settings to determine what left contexts to print in the global analysis
MINIMUM_FREQUENCY = 2
MINIMUM_PMI = 5

# Set of annotators that we are interested in
ANNOTATORS = set(['Ann', 'Mei', 'Phil'])


class Corpus(object):

    """Based on a directory with files from Brat annotations and relates sets of
    annotation and text files with annotation unit. All annotation units and the
    overal corpus are asscoiated with annotations.

    directory    -  path to the corpus directory
    name         -  basename of the directory
    fnames       -  list of files in the corpus directory
    units        -  { unit-name => AnnotationUnit }
    vocabulary   -  Vocabulary
    annotations  -  Annotations
    contexts     -  Contexts

    """

    def __init__(self, directory):
        self.directory = directory
        self.name = os.path.basename(self.directory)
        self.fnames = self._get_filenames()
        self.units = self._get_units()
        self.vocabulary = Vocabulary(self)
        self.annotations = Annotations(corpus=self)
        self.annotations.add_corpus_annotations()
        self.contexts = Contexts(self)
        print('>>>', self)

    def __str__(self):
        return ("<Corpus '%s' with %d units and %d files>"
                % (self.name, len(self.units), len(self.fnames)))

    def _get_filenames(self):
        """Return a list of filenames for the corpus, but only include those
        where the annotator in in the ANNOTATORS set."""
        fnames = glob.glob("%s/EHR/*/[0-9]*_[0-9]*" % self.directory)
        filtered_fnames = []
        for f in fnames:
            annotator = os.path.normpath(f).split(os.sep)[-2]
            if annotator in ANNOTATORS:
                filtered_fnames.append(f)
        return filtered_fnames

    def _get_units(self):
        units = {}
        for unit_name, fnames in self.split_files().items():
            if False:
                print('>>>', unit_name)
                for fname in fnames:
                    print('        ', os.sep.join(fname.split(os.sep)[-2:]))
            units[unit_name] = AnnotationUnit(self, unit_name, fnames)
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
            UnitFile(self, unit).write()
            units_written.append(unit)
            break # >>>
        contexts = self.contexts.get_significant_contexts()
        #print(contexts['Event']['bullied'])
        for tag in self.contexts.data:
            extents =  self.annotations.get_grouped_extents(tag)
            significant_contexts = self.contexts.get_significant_contexts2(tag)
            ContextFile(self, tag, significant_contexts, extents).write()
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

    The word pairs used to calculate the pointwise mutual are taken from the
    bigram counter.

    """

    def __init__(self, corpus):
        """Initialize from a Corpus. Create all token/bigram lists and Counters
        and calculate all PMIs."""
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

    def _calculate_pmi_scores(self):
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

    def get_pmi(self, pair):
        """Return <P(x,y), P(x), P(y), PMI(x,y)> for a pair. If the PMI of a
        pair is not defined return <None, None, None, None>."""
        return self.pmis.get(pair, (None, None, None, None))[3]

    def print_pmi(self, threshold, limit=None):
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
    annotators   -  set of annotator names
    annotations  -  Annotations

    The name of the unit is the part of the filename that all files in the unit
    share, for example, if the two files in the unit are 107176272_0.ann and
    107176272_0.txt then the name is 107176272_0.

    """

    def __init__(self, corpus: Corpus, name: str, fnames: list):
        """Create an AnnotationUnit from a list of files in a Corpus. The name of
        the unit was derived from the list of files."""
        self.corpus = corpus
        self.name = name
        self.fnames = fnames
        self.text = None
        self._initialize_text()
        self.annotators = set(self.get_annotator(fname) for fname in self.fnames)
        self.annotations = Annotations(unit=self)
        self.annotations.add_unit_annotations()

    def __str__(self):
        return ("<AnnotationUnit %s files=%s extents=%d relations=%d>"
                % (self.name, len(self.fnames),
                   self.number_of_annotations, self.number_of_relations))

    def __getitem__(self, i):
        return self.fnames[i]

    def _initialize_text(self):
        """Pull the text out of the text file, if there are more text files pick the
        first after sorting them on length of the file name."""
        text_files = [f for f in self.fnames if f.endswith('.txt')]
        text_files = sorted(text_files, key=lambda fname: len(fname))
        if text_files:
            self.text = open(text_files[0]).read()
        else:
            print("Warning: no text file for %s" % self)

    def text_size(self):
        return 0 if self.text is None else len(self.text)

    def get_extent(self, annotator, identifier):
        return self.annotations.get_extent(annotator, identifier)

    def get_extents(self):
        return self.annotations.get_extents()

    def get_relations(self):
        return self.annotations.get_relations()

    def get_number_of_extents(self):
        return self.annotations.number_of_extents
    
    def get_number_of_relations(self):
        return self.annotations.number_of_relations
    
    def get_annotator(self, fname):
        """Get the annotator for the file. We expect paths like
        brat_backup/5_12_22/EHR/Phil/105541561_2.ann."""
        return os.path.split(os.path.split(fname)[0])[1]

    def get_annotation_files(self):
        """Get all files with annotations."""
        return [f for f in sorted(self.fnames) if f.endswith('ann')]

    def get_grouped_extents(self):
        # TODO: shouldn't this be in the Annotations object?
        extents_sorted = { tag: [] for tag in EXTENT_TAGS }
        for extent in sorted(self.get_extents()):
            extents_sorted[extent.tag].append(extent)
        return extents_sorted

    def print_files(self):
        for fname in self.fnames:
            print(os.path.basename(fname))

    def print_annotations(annotations):
        for annotator, tag, offsets, text in self.annotations.extents:
            print("%s\t%s\t%s\t%s" % (annotator, tag, offsets, text))

    def print_annotations_index(self):
        for annotator in self.annotations.extents_idx.keys():
            print(annotator, end=': ')
            for identifier in self.annotations.extents_idx[annotator].keys():
                print(identifier, end=' ')
            print()


class Contexts(object):

    """Keeps track of the contexts of annotations.

    corpus       -  Corpus
    vocabulary   -  Vocabulary (shared with the corpus)
    annotations  -  Annotations (shared with the corpus)
    data         -  { tag => { extent => list of (term1, term2, pmi) }}

    """

    def __init__(self, corpus):
        """Collect triples of the last word of the left context and the first
        word of the key phrase and the pmi of those pairs."""
        self.corpus = corpus
        self.vocabulary = corpus.vocabulary
        self.annotations = corpus.annotations
        self.data = {tag: {} for tag in EXTENT_TAGS}
        self.data2 = {tag: {} for tag in EXTENT_TAGS}
        for tag in self.annotations.grouped_extents_idx:
            for keyphrase in self.annotations.grouped_extents_idx[tag]:
                self.data[tag][keyphrase] = []
                self.data2[tag][keyphrase] = { 'count': 0, 'data': {} }
                for extent in self.annotations.grouped_extents_idx[tag][keyphrase]:
                    t1, t2 = extent.previous_token(), extent.first_token()
                    pmi = self.vocabulary.get_pmi((t1, t2))
                    if None not in (t1, t2, pmi):
                        self.data[tag][keyphrase].append((t1, t2, pmi))
                        d = self.data2[tag][keyphrase]
                        d['count'] += 1
                        d['data'].setdefault((t1, t2), { 'count': 0, 'pmi': pmi, 'data': [] })
                        d['data'][(t1, t2)]['count'] += 1
                        d['data'][(t1, t2)]['data'].append((extent, t1, t2, pmi))
                    else:
                        # TODO: this happens way too often, find out why
                        pass
        #self.print_data()
        #self.print_significant_data()

    def get_significant_contexts(self):
        """Returns a subset of the contexts, only including those where the
        frequency and PMI are above a certain threshold."""
        significant_contexts =  {tag: {} for tag in EXTENT_TAGS}
        for tag in self.data:
            for keyphrase in sorted(self.data[tag]):
                triples = sorted(self.data[tag][keyphrase])
                triples_bag = Counter(triples)
                significant = [(t, c) for (t, c) in triples_bag.most_common()
                               if c >= MINIMUM_FREQUENCY and t[2] >= MINIMUM_PMI]
                if significant:
                    significant_contexts[tag][keyphrase] = significant
        return significant_contexts

    def print_data(self):
        for tag in list(self.data2.keys())[:10]:
            print(tag)
            for keyphrase in list(self.data2[tag].keys())[:10]:
                count = self.data2[tag][keyphrase]['count']
                phrase_data = self.data2[tag][keyphrase]['data']
                print('    %3d  %s' % (count, protect(keyphrase)))
                for pair in phrase_data:
                    pair_data = phrase_data[pair]
                    pair_count = pair_data['count']
                    pair_pmi = pair_data['pmi']
                    print('       %3d  %.2f  %s ' % (pair_count, pair_pmi, pair[0]))

    def print_significant_data(self):
        for tag in list(self.data2.keys()):
            print(tag)
            for keyphrase in list(self.data2[tag].keys()):
                count = self.data2[tag][keyphrase]['count']
                phrase_data = self.data2[tag][keyphrase]['data']
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

    def get_significant_contexts2(self, tag):
        contexts = {}
        for keyphrase in list(self.data2[tag].keys()):
            count = self.data2[tag][keyphrase]['count']
            contexts[keyphrase] = { 'count': count, 'data': {} }
            phrase_data = self.data2[tag][keyphrase]['data']
            for pair in phrase_data:
                pair_data = phrase_data[pair]
                pair_count = pair_data['count']
                pair_pmi = pair_data['pmi']
                if pair_pmi >= MINIMUM_PMI and pair_count >= MINIMUM_FREQUENCY:
                    contexts[keyphrase]['data'][pair] = pair_data
        return contexts

    
class Annotations(object):

    """Maintains the Annotations for an AnnotationUnit or a Corpus.

    extents              -  ExtentAnnotations
    extents_idx          -  { annotator => { identifier => ExtentAnnotation }}
    grouped_extents      -  { tag => ExtentAnnotations }
    grouped_extents_idx  -  { tag => { extent => ExtentAnnotations }}
    relations            -  RelationAnnotations

    """

    def __init__(self, unit=None, corpus=None):
        self.unit = unit
        self.corpus = corpus
        self.extents = []
        self.extents_idx = {}
        self.grouped_extents = None
        self.grouped_extents_idx = None
        self.relations = []
        self.number_of_extents = 0
        self.number_of_relations = 0

    def add_unit_annotations(self):
        if self.unit.text is None:
            print('Warning: not adding annotations for', self)
        else:
            for fname in self.unit.get_annotation_files():
                annotator = self.unit.get_annotator(fname)
                with open(fname) as fh:
                    for line in fh:
                        fields = line.strip().split('\t')
                        if len(fields) == 1:
                            print('Warning, unexpected line:')
                            print('---', '\t'.join(fields), '\n---', fname)
                            continue
                        tag = fields[1].split()[0]
                        if tag in EXTENT_TAGS:
                            a = ExtentAnnotation(self.unit, annotator, fields)
                            self.add_extent(a)
                        elif tag in RELATION_TAGS:
                            a = RelationAnnotation(self.unit, annotator, fields)
                            self.add_relation(a)

    def add_corpus_annotations(self):
        # populate the list of extents and create the sorted extents dictionary
        self.grouped_extents = {tag: [] for tag in EXTENT_TAGS}
        for unit in self.corpus.units.values():
            for extent in unit.get_extents():
                self.extents.append(extent)
                self.grouped_extents[extent.tag].append(extent)
        # group the extents on tag and extent
        self.grouped_extents_idx = {tag: {} for tag in EXTENT_TAGS}
        for extent in self.extents:
            # ignoring the difficult cases with discontinuous annotations
            if len(extent.positions) > 1:
                continue
            text = extent.get_key_string()
            self.grouped_extents_idx[extent.tag].setdefault(text, []).append(extent)

    def add_extent(self, a):
        self.extents.append(a)
        self.extents_idx.setdefault(a.annotator, {})
        self.extents_idx[a.annotator][a.identifier] = a
        self.number_of_extents += 1

    def add_relation(self, a):
        self.relations.append(a)
        self.number_of_relations += 1

    def get_extent(self, annotator, identifier):
        return self.extents_idx[annotator][identifier]

    def get_extents(self):
        return self.extents

    def get_relations(self):
        return self.relations

    def get_grouped_extents(self, tag):
        return self.grouped_extents_idx[tag]

    def print_grouped_extents_summary(self):
        """Print count for each tag."""
        for tag in self.grouped_extents:
            print("%5d  %s" % (len(self.grouped_extents[tag]), tag))

    def print_grouped_extents(self):
        for tag in self.grouped_extents_idx:
            for key in self.grouped_extents_idx[tag]:
                annotations = self.grouped_extents_idx[tag][key]
                print(tag, len(annotations), key)


class Annotation(object):

    """All information relevant to an annotation. Abstract class that is only
    used when initializing an ExtentAnnotation or a RelationAnnotation.

    unit        -  AnnotationUnit
    text        -  unit.text
    annotator   -  annotator (string)
    identifier  -  annotation identifier (for example T154)
    tag         -  annotation tag (Event|Symptom|...)

    Initialization takes the unit an annotation is in, the annotator name and
    the fields from the annotation file. For extent annotations there are three
    fields and for relation annotations there are two:

    T21     Temporal_Frame 3363 3369        Recent
    R5      Grounded_To Arg1:T19 Arg2:T21

    """
    
    def __init__(self, unit, annotator, fields):
        self.unit = unit
        self.text = unit.text
        self.annotator = annotator
        self.identifier = fields[0]
        self.tag = fields[1].split()[0]


class ExtentAnnotation(Annotation):

    """Instance variables (not including the ones ingerited from Annotation):

    extent         -  annotation extent (string, for example 'made threats')
    offsets        -  start and end offsets ("8284:8295", "8284:8295;8304:8309") 
    positions      -  offsets as a list ([[828,8295],[8304,8309]])
    start          -  start of first fragment (integer)
    end            -  end of last fragment (integer)
    left_context   -  N=CONTEXT characters to the left of start position
    right_context  -  N=CONTEXT characters to the right of end position

    The more complicated case for offsets and positions is for discontinuous
    annotations (that is, annotaitons with more than one fragment). There is
    currently no middle context, but there may be a future need for it.

    The keyphrase is a bit of an odd duck here, want to get rid of it.
    
    """
    
    def __init__(self, unit, annotator, fields):
        super().__init__(unit, annotator, fields)
        self.extent = fields[2]
        self.offsets = ':'.join(fields[1].split()[1:])
        pairs = [pair for pair in self.offsets.split(';')]
        self.positions = [[int(pos) for pos in pair.split(':')] for pair in pairs]
        self.start = self.positions[0][0]
        self.end = self.positions[-1][-1]
        self.left_context = self.get_left_context()
        self.right_context = self.get_right_context()

    def __str__(self):
        return ("<%s %s %s %s '%s' %s>"
                % (self.tag, self.unit.name, self.identifier, self.offsets,
                   self.extent, self.annotator))

    def __eq__(self, other):
        # for sorting purposes, we do not care about the identifier or annotator
        return self.positions == other.positions

    def __lt__(self, other):
        return self.positions[0][0] < other.positions[0][0]

    def get_key_string(self):
        # This just returns the string from the text."""
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

    def positions_as_string(self):
        return ';'.join([self.fragment_as_string(f) for f in self.positions])

    def fragment_as_string(self, fragment):
        return "%s:%s" % (fragment[0], fragment[1])


class RelationAnnotation(Annotation):

    def __init__(self, unit, annotator, fields):
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


class ExtentKwic(object):

    """Implements a context for an ExtentAnnotation. The main task of this class
    is to generate a line in the KWIC for the extent annotation.

    annotation     -  Annotations that the KWIC is being created for
    annotator      -  the annotator (str)
    identifier     -  identifier copied from the annotation 
    left_context   -  left context string copied from the Annotation
    right_context  -  right context string copied from the Annotation
    keyphrase      -  list of pairs for the keyphrase

    The keyphrase usually is something like

        [('KEY', 'emotional abuse')]

    But for discontinuous tags it is more complex since it includes a key for
    each fragment as weel as the internal contexts spans between the fragments:
    
        [('KEY', 'car'), ('CONTEXT', 'behind'), ('KEY', 'hit her')]

    """
    
    def __init__(self, annotation):
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
                substring = annotation.text[p2:annotation.positions[i+1][0]]
                self.keyphrase.append(('CONTEXT', annotation.text[p2:annotation.positions[i+1][0]]))

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
            location = "%s:%s" % (self.annotation.unit.name, self.annotation.start)
            dtrs.insert(0, TD(dtrs=[Text(location)]))
        return TR(dtrs=dtrs)

    def keyphrase_as_html(self):
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
    content   -  an other instance of html.Tag, used for the main content
    fname     -  the file that the HTML object will be written to

    """

    def __init__(self, title):
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
            ".count { text-align: right; width: 25px; color: darkgreen; }\n"
            ".pmi { text-align: right; width: 45px; color: darkgreen; }\n"
            "table.fancy { margin-top: 10px; padding: 5px; border: thin solid lightgray; }\n")
    
    def write(self):
        """Write the HTML file, using the fname variable as the destination."""
        print('Writing %s' % self.fname)
        with open(self.fname, 'w') as fh:
            fh.write(str(self.html_obj))



class IndexFile(HtmlFile):

    """Creates an HTML file for the index."""
    
    def __init__(self, corpus, units):
        title = "Annalysis of annotations in %s" % corpus.name
        super().__init__(title)
        self.units = units
        self.fname = 'html/%s/index.html' % corpus.name
        self._add_left_context_links()
        self._add_annotation_units_table()

    def _add_left_context_links(self):
        header = "Global observations on left context:\n"
        hrefs = Tag('ul')
        self.content.add(P(header), hrefs, NL)
        # the contexts for each tag
        for tag in sorted(EXTENT_TAGS):
            hrefs.add(Tag('li', dtrs=[Href('global-%s.html' % tag, tag)]))

    def _add_annotation_units_table(self):
        header = "Annotations per file:\n"
        table = Tag('table', attrs={'cellspacing': 0, 'cellpadding': 5, 'border': 1})
        self.content.add(P(header), NL, Tag('blockquote', dtrs=[table, NL]), NL)
        tds = [SPACE, Tag('td', dtrs=[Text('')]), NL,
               SPACE, Tag('td', dtrs=[Text('text size')]), NL,
               SPACE, Tag('td', dtrs=[Text('extents')]), NL,
               SPACE, Tag('td', dtrs=[Text('relations')]), NL,
               SPACE, Tag('td', dtrs=[Text('annotators')]), NL]
        tr = Tag('tr', dtrs=tds)
        table.add(tr, NL)
        for unit in self.units:
            name = unit.name
            annotators = ', '.join(unit.annotators)
            size = Text(unit.text_size())
            extents = Text(unit.get_number_of_extents())
            relations = Text(unit.get_number_of_relations())
            tds = [SPACE, Tag('td', dtrs=[Href('kwic-%s.html' % name, name)]), NL,
                   SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[size]), NL,
                   SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[extents]), NL,
                   SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[relations]), NL,
                   SPACE, Tag('td', dtrs=[Text(annotators)]), NL]
            table.add(TR(dtrs=tds), NL)


class UnitFile(HtmlFile):

    """Creates an HTML file for an AnnotationUnit, which is basically a
   concordance of all annotations."""
    
    def __init__(self, corpus, unit):
        super().__init__(unit.name)
        self.unit = unit
        self.fname = 'html/%s/kwic-%s.html' % (corpus.name, self.unit.name)
        self.grouped_extents = self.unit.get_grouped_extents()
        self._add_navigation()
        self._add_extents()
        self._add_relations()

    def _add_navigation(self):
        # only add those tags for which there are examples
        annotations_filtered = {k:v for k, v in self.grouped_extents.items() if v}
        relations = self.unit.get_relations()
        div = Tag('div', attrs={'class': 'navigation'})
        for i, tag in enumerate(sorted(annotations_filtered)):
            prefix = '\n| ' if i > 0 else '[ '
            div.add(Text(prefix), Href('#' + tag, tag + 's'))
        if relations:
            tag = 'Relation'
            div.add(Text('\n| '))
            div.add(Href('#' + tag, tag + 's'))
            div.add(NL)
        div.add(Text(']\n'))
        self.content.add(Text('\n\n'), div, NL)

    def _add_extents(self):
        for tag in sorted(self.grouped_extents):
            if self.grouped_extents[tag]:
                table = Tag('table', nl=True, attrs={'class': 'fancy', 'cellpadding': 5})
                self.content.add(NL, Anchor(tag), H3('%ss' % tag), NL, table)
                for a in remove_duplicates(self.grouped_extents[tag]):
                    table.add(ExtentKwic(a).as_html())

    def _add_relations(self):
        if self.unit.get_relations():
            table = Tag('table', nl=True, attrs={'class': 'fancy', 'cellpadding': 5})
            relations = self.unit.get_relations()
            for rel in sorted(relations, key=lambda x: x.locations()):
                tag = rel.tag
                arg1_id = rel.arg1.identifier
                arg2_id = rel.arg2.identifier
                args = rel.sorted_args()
                signature = "%s(%s, %s)" % (rel.identifier, arg1_id, arg2_id)
                signature = "%s&nbsp;&rightarrow;&nbsp;%s" % (arg1_id, arg2_id)
                td_ann = TD(attrs={'class': 'source', 'valign': 'top'},
                            dtrs=[Text(rel.annotator)])
                td_signature =  TD(attrs={'class': 'signature', 'valign': 'top'},
                                   dtrs=[Text(signature)])
                context = ([Text(args[0].left_context)]
                           + [green('{')]
                           + ExtentKwic(args[0]).keyphrase_as_html().dtrs
                           + [green('}'), red_subscript(args[0].identifier)]
                           + [Text(rel.text_between_args())]
                           + [green('{')]
                           + ExtentKwic(args[1]).keyphrase_as_html().dtrs
                           + [green('}'), red_subscript(args[1].identifier)]
                           + [Text(args[1].right_context)]
                           )
                td_context = TD(dtrs=context)
                table.add(TR(dtrs=[td_ann, td_signature, td_context]))
            self.content.add(NL, Anchor('Relation'), H3('%ss' % 'Relation'), NL, table)
        

class ContextFile(HtmlFile):

    def __init__(self, corpus, tag, contexts_for_tag, extents):
        super().__init__("%s left context observations" % tag)
        self.fname = 'html/%s/global-%s.html' % (corpus.name, tag)
        self.contexts_for_tag = contexts_for_tag
        self.extents = extents
        for extent in sorted(self.contexts_for_tag):
            self.contexts_for_extent = self.contexts_for_tag[extent]
            if self.contexts_for_extent['data']:
                total_count = str(self.contexts_for_extent['count'])
                for pair, pair_data in self.contexts_for_extent['data'].items():
                    self._add_summary(total_count, pair, pair_data, extent)
                    self._add_contexts(pair_data)

    def _add_summary(self, total_count, pair, pair_data, extent):
        table =  Tag('table', nl=True, attrs={'class': 'fancy', 'cellpadding': 8})
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
                           TD(attrs=self.attrs(40), dtrs=[Text(pmi)]),
                           TD(attrs=self.attrs(20), dtrs=[Text(wider_extent_count)]),
                           TD(dtrs=[text_left, SPACE, text_key])]))

    def _add_contexts(self, pair_data):
        table = Tag('table', nl=True, attrs={'cellpadding': 8})
        self.content.add(Tag('blockquote', dtrs=[table]))
        lextents = pair_data['data']
        for lextent, t1, t2, pmi in lextents:
            table.add(ExtentKwic(lextent).as_html(add_location=True))

    @staticmethod
    def attrs(width):
        return {'width': width, 'align': 'right'}

        

def remove_duplicates(annotations):
    """Remove the duplicates from a list of annotations. This is done for the
    pre-annotations which appear to all have multiple instances with the same
    annotator, tag and offset, but different identifiers."""
    seen = set()
    filtered_annotations = []
    for anno in annotations:
        tup = (anno.annotator, anno.tag, anno.offsets)
        if not tup in seen:
            seen.add(tup)
            filtered_annotations.append(anno)
    return filtered_annotations

def protect(text):
    return text.replace('\n', ' ').replace('&', '&amp;').replace('<', '&lt;')

def green(text):
    return Tag('span', attrs={'class': 'green'}, dtrs=[Text(text)])

def green_subscript(text):
    return Tag('sub', attrs={'class': 'green'}, dtrs=[Text(text)])

def red_subscript(text):
    return Tag('sub', attrs={'class': 'darkred'}, dtrs=[Text(text)])


## Main function

def run():
    date = sys.argv[1]
    corpus = Corpus(os.path.join(BRAT_BACKUP, date))
    corpus.write_reports()


if __name__ == '__main__':
    
    run()
