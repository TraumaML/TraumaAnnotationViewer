"""show.py

Show all the annotations in a series of concordances for each file and show some
global analyses on the left context.

$ python run.py <BACKUP>

This runs on a directory with all Brat backups as defined in the BRAT_BACKUP
global variable, the argument identifies the backup we want. As of the beginning
of April, three backups can be processed: 3_7_22, 3_31_22 and 4_4_22. The first
is somewhat different from the others and code is sensitive to the differences
by using an if-then-else statement that tests what functionality to use. If you
add a new kind of backup you may need to makes some changes whetever '3_7_22' is
used in a test.

Results can be viewed at html/<BACKUP>/index.html.


TODO:

- ✓ include grounded_to relations

- ✓ include pre-annotations
    (only printing them for the files without any manual annotation)

- ✓ round up all left contexts from all documents
    (maybe use a flag for whether to use manual annotations or pre-annotations)

- ✓ get mutual information of left token with token prior to extent

- ✓ when we have a promising left extent (for example "physical [abuse]") find
    cases where that left contents may have been included in the annotation

- ✓ some things like selecting the files and determining what are annotations
    files may change with the next backup, would like to make those things
    dependent on what backup we are dealing with.

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


class Corpus(object):

    """Basically a directory with files from Brat annotations where related sets
    of annotation and text files are turned into an annotation unit."""

    def __init__(self, directory):
        self.directory = directory
        self.name = os.path.basename(self.directory)
        self.fnames = self._get_filenames()
        self.units = self._get_units()
        self.vocabulary = Vocabulary(self)
        self.annotations = None      # list of extent annotations for each tag
        self.annotations_idx = None  # annotations indexed on tag and extent
        self.contexts = None         # contexts indexed on tag and extent
        self._collect_annotations()
        self._collect_contexts()
        print('>>>', self)

    def __str__(self):
        return "<Corpus '%s' with %d units and %d files>" % (self.name, len(self.units), len(self.fnames))

    def _get_filenames(self):
        if self.name == '3_7_22':
            return glob.glob("%s/[0-9]*_[0-9]*" % self.directory)
        elif self.name == '5_12_22':
            fnames = glob.glob("%s/EHR/*/[0-9]*_[0-9]*" % self.directory)
            annotators = set(['Ann', 'Mei', 'Phil'])
            filtered_fnames = []
            for f in fnames:
                annotator = os.path.normpath(f).split(os.sep)[-2]
                if annotator in annotators:
                    filtered_fnames.append(f)
            return filtered_fnames
        else:
            return glob.glob("%s/EHR/*/[0-9]*_[0-9]*" % self.directory)

    def _get_units(self):
        units = {}
        for unit_name, fnames in self.split_files().items():
            if False:
                print('>>>', unit_name)
                for fname in fnames:
                    print('        ', os.sep.join(fname.split(os.sep)[-2:]))
            units[unit_name] = AnnotationUnit(self, unit_name, fnames)
        return units

    def _collect_annotations(self):
        # get the list of annotations for all tags
        self.annotations = {tag: [] for tag in EXTENT_TAGS}
        for unit in self.units.values():
            for annotation in unit.new_annotations.extents:
                self.annotations[annotation.tag].append(annotation)
        # group the annotations on tag and extent
        self.annotations_idx = {tag: {} for tag in EXTENT_TAGS}
        for tag in self.annotations:
            for annotation in self.annotations[tag]:
                # ignoring the difficult cases with discontinuous annotations
                if len(annotation.positions) > 1:
                    continue
                text = annotation.get_key_string()
                self.annotations_idx[tag].setdefault(text, []).append(annotation)

    def _collect_contexts(self):
        """Collect triples of the last word of the left context and the first
        word of the key phrase and the pmi of those pairs."""
        self.contexts = {tag: {} for tag in EXTENT_TAGS}
        for tag in self.annotations_idx:
            for keyphrase in self.annotations_idx[tag]:
                self.contexts[tag][keyphrase] = []
                for extent in self.annotations_idx[tag][keyphrase]:
                    #print('>>>', extent.fname, extent)
                    t1, t2 = extent.previous_token(), extent.first_token()
                    pmi = self.vocabulary.get_pmi((t1, t2))
                    if None not in (t1, t2, pmi):
                        self.contexts[tag][keyphrase].append((t1, t2, pmi))
        #self.get_significant_contexts()

    def get_significant_contexts(self):
        """Returns a subset of the contexts, only including those where the
        frequency and PMI are above a certain threshold."""
        significant_contexts =  {tag: {} for tag in EXTENT_TAGS}
        for tag in self.contexts:
            for keyphrase in sorted(self.contexts[tag]):
                triples = sorted(self.contexts[tag][keyphrase])
                triples_bag = Counter(triples)
                significant = [(t, c) for (t, c) in triples_bag.most_common()
                               if c >= MINIMUM_FREQUENCY and t[2] >= MINIMUM_PMI]
                if significant:
                    significant_contexts[tag][keyphrase] = significant
        return significant_contexts

    def split_files(self):
        """Split the files into groups with the same file identifier."""
        files = {}
        for fname in self.fnames:
            base = os.path.basename(fname)[:11]
            files.setdefault(base, []).append(fname)
        return files

    def _get_unit_names(self):
        if self.name == '3_7_22':
            # manually defined the files for the pilot annotation
            return ('104812792_0', '105663850_0', '105968572_0',
                    '106041479_0', '107176272_0')
        else:
            # TODO: using this for now, it includes the cases with the
            # non-expected line, remove for final run
            #return sorted(('104350400_0', '106041479_0', '106355919_0',
            #              '103505751_0', '105663850_0', '106357503_0', '101602488_0'))
            return sorted(k for k in list(self.units.keys()))

    def write_reports(self):
        """Write all annotations from all units to files as well as all context
        observations. Also write an index that ties all files together."""
        html_dir = os.path.join('html', self.name)
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)
        units_written = []
        for unit_name in self._get_unit_names():
            print('Writing html/%s/%s.html' % (self.name, unit_name))
            unit = self.units[unit_name]
            unit.write_html()
            units_written.append(unit_name)
            break # >>>
        contexts = self.get_significant_contexts()
        for tag in self.contexts:
            fname = 'html/%s/global-%s.html' % (self.name, tag)
            print('Writing %s' % fname)
            contexts_for_tag = self.get_significant_contexts()
            #print('>>>', contexts[tag])
            write_contexts(fname, tag, contexts[tag], self.annotations_idx[tag])
        write_index('html/%s/index.html' % self.name, self, units_written)

    def print_units(self):
        for unit_name in sorted(self.units):
            unit = self.units[unit_name]
            print(unit)

    def print_annotations_summary(self):
        """Print count for each tag."""
        for tag in self.annotations:
            print("%5d  %s" % (len(self.annotations[tag]), tag))

    def print_annotations(self):
        for tag in self.annotations_idx:
            for key in self.annotations_idx[tag]:
                annotations = self.annotations_idx[tag][key]
                print(tag, len(annotations), key)


class Vocabulary(object):

    """Keeps track of counts of tokens and bigrams as well as the PMI between all
    token pairs."""

    def __init__(self, corpus):
        self.corpus = corpus
        self.vocabulary = {}
        self.tokens = []
        for unit in corpus.units.values():
            if unit.text is not None:
                tokens = [t.lower() for t in unit.text.split()]
                self.tokens.extend(tokens)
        self.vocabulary = Counter(self.tokens)
        self.number_of_tokens = len(self.tokens)
        self.number_of_types = len(self.vocabulary)
        self.bigrams = self.get_bigrams()
        self.number_of_bigrams = len(self.bigrams)
        self.bigrams_counter = Counter(self.bigrams)
        self.pmi = None
        self._calculate_pmi_scores()

    def _calculate_pmi_scores(self):
        self.pmi = {}
        for pair, count in self.bigrams_counter.items():
            p_x_y = count / self.number_of_bigrams
            p_x = self.vocabulary[pair[0]] / self.number_of_tokens
            p_y = self.vocabulary[pair[1]] / self.number_of_tokens
            pmi = log(p_x_y / (p_x * p_y))
            self.pmi[pair] = (p_x_y, p_x, p_y, pmi)

    def __len__(self):
        return len(self.vocabulary)

    def __str__(self):
        return ("<Vocabulary on %s tokens with %d elements>"
                % (self.number_of_tokens, len(self)))
    
    def get_bigrams(self):
        bigrams = []
        for i in range(self.number_of_tokens - 1):
            bigrams.append((self.tokens[i], self.tokens[i+1]))
        return bigrams

    def get_pmi(self, pair):
        return self.pmi.get(pair, (None, None, None, None))[3]

    def print_pmi(self, threshold, limit=sys.maxsize):
        printed = 0
        for pair, (p_x_y, p_x, p_y, pmi) in self.pmi.items():
            if printed <= limit:
                c_x_y = self.bigrams_counter[pair]
                c_x = self.vocabulary[pair[0]]
                c_y = self.vocabulary[pair[1]]
                if pmi > threshold and c_x_y > 1:
                    print("%2d %f  %3d %f  %3d %f  %7.4f  %s"
                          % (c_x_y, p_x_y, c_x, p_x, c_y, p_y, pmi, pair))
                    printed += 1


class AnnotationUnit(object):

    """All information associated with an annotated file, including all the
    annotations and the text of the file.

    The name of the unit is the part of the filename that all files in the unit
    have in common, for example, if the two files in the unit are
    "107176272_0.ann" and "107176272_0.txt" then the name is "107176272_0"."""

    def __init__(self, corpus, name, fnames):
        self.corpus = corpus
        self.name = name
        self.fnames = fnames
        self.text = None
        self._initialize_text()
        self.annotators = set(self.get_annotator(fname) for fname in self.fnames)
        self.new_annotations = Annotations(unit=self)
        self.new_annotations.add_unit_annotations()

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
        return self.new_annotations.get_extent(annotator, identifier)

    def get_extents(self):
        return self.new_annotations.get_extents()

    def get_annotator(self, fname):
        if self.corpus.name == '3_7_22':
            annotator = fname.split('_')[-1].split('.')[0]
            if annotator == '0':
                annotator = 'pre'
        else:
            # we expect paths like "brat_backup/3_31_22/EHR/Phil/105541561_2.ann"
            annotator = os.path.split(os.path.split(fname)[0])[1]
        return annotator

    def get_annotation_files(self):
        """Get all files with annotations, allowing for one special case."""
        if self.corpus.name == '3_7_22':
            # special case for one file in 3_7_22
            ext = '.txt.ann' if self.name == '105663850_0' else '.ann'
        else:
            ext = 'ann'
        return [f for f in sorted(self.fnames) if f.endswith(ext)]

    def print_files(self):
        for fname in self.fnames:
            print(os.path.basename(fname))

    def print_annotations(annotations):
        for annotator, tag, offsets, text in self.annotations:
            print("%s\t%s\t%s\t%s" % (annotator, tag, offsets, text))

    def print_annotations_index(self):
        for annotator in self.annotations_idx.keys():
            print(annotator, end=': ')
            for identifier in self.annotations_idx[annotator].keys():
                print(identifier, end=' ')
            print()

    def write_html(self):
        write_unit('html/%s/kwic-%s.html' % (self.corpus.name, self.name), self)


def adjust_offset(unit_name, p, annotator):
    """This deals with the first line which is added for some documents, which
    introduced an inconsistency in offsets for identical annotation offsets
    depending on the length of the name of the annotator."""
    # TODO: annotation offset adjustment should be done at the unit level, now
    # it is called for each offset
    if unit_name == '105663850_0':
        return p - 21 if annotator == 'Phil' else p - 20
    else:
        return p


class Annotations(object):

    """Maintains the Annotations for an AnnotationUnit or a Corpus.

    self.extents      -  List of ExtentAnnotations
    self.extents_idx  -  { annotator => { identifier => ExtentAnnotation }}
    self.relations    -  List of RelationAnnotations

    """

    def __init__(self, unit=None, corpus=None):
        self.unit = unit
        self.corpus = corpus
        self.extents = []
        self.extents_idx = {}
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
        pass

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


class Annotation(object):

    def __init__(self, unit, annotator, fields):
        self.fields = fields
        self.unit = unit
        self.text = unit.text
        self.fname = unit.name
        self.annotator = annotator
        self.identifier = fields[0]
        self.tag = fields[1].split()[0]


class ExtentAnnotation(Annotation):

    def __init__(self, unit, annotator, fields):
        super().__init__(unit, annotator, fields)
        self.extent = fields[2]
        self.offsets = ':'.join(fields[1].split()[1:])
        pairs = [pair for pair in self.offsets.split(';')]
        self.positions = [[adjust_offset(self.fname, int(pos), self.annotator)
                           for pos in pair.split(':')] for pair in pairs]
        self.start = self.positions[0][0]
        self.end = self.positions[-1][-1]
        self.keyphrase = self.get_keyphrase()
        self.left_context = self.get_left_context()
        self.right_context = self.get_right_context()

    def __str__(self):
        return ("<%s %s %s '%s' %s>"
                % (self.tag, self.identifier, self.offsets,
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

    def get_keyphrase(self):
        # The key phrase is a bit different from the contexts in that this is a
        # list of html objects (tags and texts) rather than a string
        parts = []
        for i, (p1, p2) in enumerate(self.positions):
            parts.append(Span(attrs={'class': 'keyword'},
                              dtrs=[Text('[%s]' % self.text[p1:p2])]))
            if i+1 < len(self.positions):
                substring = self.text[p2:self.positions[i+1][0]]
                parts.append(Span(dtrs=[Text(substring)]))
        return parts

    def get_left_context(self):
        return protect(self.text[self.start-CONTEXT:self.start])

    def get_right_context(self):
        return protect(self.text[self.end:self.end+CONTEXT])

    def positions_as_string(self):
        return ';'.join([self.fragment_as_string(f) for f in self.positions])

    def fragment_as_string(self, fragment):
        return "%s:%s" % (fragment[0], fragment[1])

    def kwic(self):
        return "%s [%s] %s" % (self.left_context, self.extent, self.right_context)


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


## Writing HTML

def write_index(html_file, corpus, unit_names):
    title = "Annalysis of annotations in %s" % corpus.name
    header1 = "Global observations on left context:\n"
    header2 = "Annotations per file:\n"
    html_obj, content = get_toplevel_html(title)
    hrefs = Tag('ul')
    content.add(P(header1), hrefs, NL)
    # the contexts for each tag
    for tag in sorted(EXTENT_TAGS):
        hrefs.add(Tag('li', dtrs=[Href('global-%s.html' % tag, tag)]))
    # all units with their annotations
    table = Tag('table', attrs={'cellspacing': 0, 'cellpadding': 5, 'border': 1})
    content.add(P(header2), NL, Tag('blockquote', dtrs=[table, NL]), NL)
    tds = [SPACE, Tag('td', dtrs=[Text('')]), NL,
           SPACE, Tag('td', dtrs=[Text('text size')]), NL,
           SPACE, Tag('td', dtrs=[Text('extents')]), NL,
           SPACE, Tag('td', dtrs=[Text('relations')]), NL,
           SPACE, Tag('td', dtrs=[Text('annotators')]), NL]
    tr = Tag('tr', dtrs=tds)
    table.add(tr, NL)
    for name in unit_names:
        unit = corpus.units[name]
        annotators = ', '.join(unit.annotators)
        size = Text(unit.text_size())
        extents = Text(unit.new_annotations.number_of_extents)
        relations = Text(unit.new_annotations.number_of_relations)
        tds = [SPACE, Tag('td', dtrs=[Href('kwic-%s.html' % name, name)]), NL,
               SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[size]), NL,
               SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[extents]), NL,
               SPACE, Tag('td', attrs={'align': 'right'}, dtrs=[relations]), NL,
               SPACE, Tag('td', dtrs=[Text(annotators)]), NL]
        table.add(TR(dtrs=tds), NL)
    with open(html_file, 'w') as fh:
        fh.write(str(html_obj))

def write_unit(html_file, unit, show_position=False):
    extents_sorted = { tag: [] for tag in EXTENT_TAGS }
    for extent in sorted(unit.new_annotations.extents):
        extents_sorted[extent.tag].append(extent)
    html_obj, content = get_toplevel_html(html_file)
    _add_navigation(content, extents_sorted, unit.new_annotations.relations)
    for tag in sorted(extents_sorted):
        if extents_sorted[tag]:
            table = Tag('table', nl=True, attrs={'class': 'fancy', 'cellpadding': 5})
            content.add(NL, Anchor(tag), H3('%ss' % tag), NL, table)
            _add_anotations(extents_sorted[tag], table)
    if unit.new_annotations.relations:
        table = Tag('table', nl=True, attrs={'class': 'fancy', 'cellpadding': 5})
        content.add(NL, Anchor('Relation'), H3('%ss' % 'Relation'), NL, table)
        _add_relations(unit.new_annotations.relations, table)
    with open(html_file, 'w') as fh:
        fh.write(str(html_obj))

def write_contexts(fname, tag, contexts_for_tag, annotations):
    html_obj, content = get_toplevel_html("%s annotations" % tag)
    # count cases indicating that the key phrase may need to be extended
    total_count = 0
    for extent in sorted(contexts_for_tag):
        contexts_for_extent = contexts_for_tag[extent]
        total_count += sum(c for t,c in contexts_for_extent)
        table = Tag('table', nl=True, attrs={'class': 'fancy', 'cellpadding': 5})
        content.add(table)
        for (triple, count) in contexts_for_extent:
            wider_extent = '%s %s' % (triple[0], extent)
            wider_extent_annos = annotations.get(wider_extent, [])
            wider_extent_count = len(wider_extent_annos)
            _add_context_row(table, count, extent, triple, wider_extent_count)
    _add_count_paragraph(content, total_count)
    with open(fname, 'w') as fh:
        fh.write(str(html_obj))

def styles():
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

def get_toplevel_html(title):
    style = Tag('style', dtrs=[Text(styles())])
    title = H1(title)
    content = Tag('div')
    html_obj = Tag('html', dtrs=[Tag('head', dtrs=[style, NL]), NL,
                                 Tag('body', dtrs=[title, content, NL]), NL])
    return html_obj, content

def _add_navigation(tables, annotations_sorted, relations):
    # only add those that are non-empty
    annotations_filtered = {k:v for k, v in annotations_sorted.items() if v}
    div = Tag('div', attrs={'class': 'navigation'})
    tables.add(Text('\n\n'), div)
    for i, tag in enumerate(sorted(annotations_filtered)):
        if i > 0:
            div.add(Text(' | '))
        div.add(Href('#' + tag, tag + 's'))
    if relations:
        tag = 'Relation'
        div.add(Text(' | '))
        div.add(Href('#' + tag, tag + 's'))

def _add_anotations(annotations, table):
    for a in remove_duplicates(annotations):
        td_source = TD(attrs={'class': 'source'}, dtrs=[Text(a.annotator)])
        span_id = Span(attrs={'class': 'tagid'}, dtrs=[Text('(%s)' % a.identifier)])
        td_left = TD(attrs={'align': 'right'}, dtrs=[Text(a.left_context)])
        td_key = TD(dtrs=a.keyphrase.copy())
        td_key.add(Text('&nbsp;&nbsp;' + a.right_context))
        td_key.add(SPACE, span_id)
        table.add(TR(dtrs=[td_source, td_left, td_key]))

def _add_relations(relations, table):
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
                   + args[0].keyphrase
                   + [green('}'), red_subscript(args[0].identifier)]
                   + [Text(rel.text_between_args())]
                   + [green('{')]
                   + args[1].keyphrase
                   + [green('}'), red_subscript(args[1].identifier)]
                   + [Text(args[1].right_context)]
                   )
        td_context = TD(dtrs=context)
        table.add(TR(dtrs=[td_ann, td_signature, td_context]))

def _add_context_row(table, count, extent, triple, wider_extent_count):
    tr = TR(dtrs=[TD(class_='count', dtrs=[Text("%d" % count)]),
                  TD(class_='pmi', dtrs=[Text("%.2f" % triple[2])]),
                  TD(class_='count', dtrs=[Text("%d" % wider_extent_count)]),
                  TD(dtrs=[Text('&nbsp;')]),
                  TD(attrs={'align': 'right'}, dtrs=[Text(protect(triple[0]))]),
                  TD(class_='keyword', dtrs=[Text(protect('[%s]' % extent))])])
    table.add(tr)

def _add_count_paragraph(content, count):
    content.insert(0, Tag('p', dtrs=[Text("Total cases with potentially useful"
                                          " left context: %d" % count)]))

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
