
import os
import typing

TAGS = ('Event', 'Perpetrator', 'Symptom', 'Substance', 'Temporal_Frame')

BLANK = ''


class BratAnnotation:

    def __init__(self, filename: str, fields: list):
        """Extract the identifier from the fields, the identifier is the only
        field that the three kinds of Brat annotations have in common."""
        self.filename = filename
        self.identifier = fields[0]

    def is_extent(self):
        return False

    def is_attribute(self):
        return False

    def is_relation(self):
        return False


class BratExtent(BratAnnotation):

    """A Brat extent annotation.

    Instance variables (not including the ones inherited from BratAnnotation):

    tag            -  normalized tag name, for example Event or Symptom
    raw_tag        -  tag name, but including suffixes like _red or _yellow
    extent         -  annotation extent (string, for example 'made threats')
    offsets        -  start and end offsets ("8284:8295", "8284:8295;8304:8309")
    positions      -  offsets as a list ([[828,8295],[8304,8309]])
    start          -  start of first fragment (integer)
    end            -  end of last fragment (integer)

    The more complicated case for offsets and positions is for discontinuous
    annotations (that is, annotations with more than one fragment). There is
    currently no middle context, but there may be a future need for it.

    """

    def __init__(self, filename, fields):
        super().__init__(filename, fields)
        self.tag = fields[1].split()[0]
        self.raw_tag = self.tag
        self._normalize_tag()
        self.extent = fields[2]
        self.offsets = ':'.join(fields[1].split()[1:])
        pairs = [pair for pair in self.offsets.split(';')]
        self.positions = [[int(pos) for pos in pair.split(':')] for pair in pairs]
        self.start = self.positions[0][0]
        self.end = self.positions[-1][-1]
        self.attributes = []

    def _normalize_tag(self):
        if self.tag.endswith('_yellow'):
            self.tag = self.tag[:-7]
        if self.tag.endswith('_red'):
            self.tag = self.tag[:-4]

    def __str__(self):
        return ("<BratExtent %s %s %s %s '%s'>"
                % (self.filename, self.tag, self.identifier, self.offsets, self.extent))

    def __eq__(self, other):
        # for sorting purposes, we do not care about the identifier or annotator
        return self.positions[0][0] == other.positions[0][0]

    def __lt__(self, other):
        return self.positions[0][0] < other.positions[0][0]

    def is_extent(self):
        return True

    def short_string(self, max_size=32) -> str:
        extent = self.extent
        if len(extent) > max_size:
            extent = extent[:max_size-3] + '...'
        return f'{self.identifier:4} {self.start:-5}  {extent}'

    def attributes_as_short_string(self) -> str:
        return f'{{ {", ".join([a.as_short_string() for a in self.attributes]) } }}'


class BratAttribute(BratAnnotation):

    def __init__(self, filename: str, fields: list):
        super().__init__(filename, fields)
        attribute, extent, *rest = fields[1].split()
        self.attribute = fields[1].split()[0]
        self.attribute = attribute
        self.extent = extent
        self.value = rest[0] if rest else None

    def __str__(self):
        return ("<BratAttribute %s %s %s %s>"
                % (self.attribute, self.identifier, self.extent, self.value))

    def __eq__(self, other):
        return self.attribute == other.attribute and self.value == other.value

    def as_short_string(self):
        return f'{self.attribute}={self.value}'

    def is_attribute(self):
        return True


class BratRelation(BratAnnotation):

    def __init__(self, filename: str, fields: list):
        super().__init__(filename, fields)
        self.tag = fields[1].split()[0]
        self.args = {}
        for arg in fields[1].split()[1:]:
            if arg.startswith('Arg') and ':' in arg:
                name, value = arg.split(':')
                self.args[name] = value
        self.arg1 = self.args.get('Arg1')
        self.arg2 = self.args.get('Arg2')
        # self.arg1 = unit.get_extent(annotator, self.args.get('Arg1'))
        # self.arg2 = unit.get_extent(annotator, self.args.get('Arg2'))

    def __str__(self):
        return ("<BratRelation %s %s %s %s>"
                % (self.tag, self.identifier, self.arg1, self.arg2))

    def locations(self):
        """Return the starting positions of the two arguments with the arguments
        ordered on text position."""
        return (self.sorted_args()[0].positions[0][0],
                self.sorted_args()[1].positions[0][0])

    def sorted_args(self):
        """Returns the arguments sorted on text position."""
        return sorted([self.arg1, self.arg2], key=lambda x: x.positions[0][0])

    def is_relation(self):
        return True


class BratAnnotations:

    @classmethod
    def is_brat_extent(cls, fields: list) -> bool:
        return fields[0].startswith('T')

    @classmethod
    def is_brat_attribute(cls, fields: list) -> bool:
        return fields[0].startswith('A')

    @classmethod
    def is_brat_relation(cls, fields: list) -> bool:
        return fields[0].startswith('R')

    def __init__(self, filename: str):
        self.filename = filename
        self.filename_base = os.path.basename(filename)
        self.extents = {}
        self.attributes = {}
        self.relations = {}
        self._read_annotations()
        self._add_attributes_to_extents()

    def _read_annotations(self):
        with open(self.filename) as fh:
            for line in fh:
                fields = line.strip().split('\t')
                if self.__class__.is_brat_extent(fields):
                    annotation = BratExtent(self.filename_base, fields)
                elif self.__class__.is_brat_relation(fields):
                    annotation = BratRelation(self.filename_base, fields)
                elif self.__class__.is_brat_attribute(fields):
                    annotation = BratAttribute(self.filename_base, fields)
                else:
                    print("WARNING: unexpected annotation: '%s'" % fields)
                    continue
                self.add(annotation)

    def _add_attributes_to_extents(self):
        for attribute in self.attributes.values():
            extent = self.extents.get(attribute.extent)
            extent.attributes.append(attribute)

    def __str__(self):
        return ("<BratAnnotations %s extents=%d attributes=%d relations=%d>"
                % (self.filename_base,
                   len(self.extents), len(self.attributes), len(self.relations)))

    def add(self, annotation: BratAnnotation):
        if annotation.is_extent():
            self.extents[annotation.identifier] = annotation
        elif annotation.is_attribute():
            self.attributes[annotation.identifier] = annotation
        elif annotation.is_relation():
            self.relations[annotation.identifier] = annotation

    def sorted_extents(self) -> dict:
        """Return a dictionary of all extends, with for each tag a list of extents
        sorted on their position in the text."""
        result = {tag: [] for tag in TAGS}
        for extent in self.extents.values():
            result[extent.tag].append(extent)
        for key, value in result.items():
            result[key] = list(sorted(value))
        return result

    def compare(self, fh: typing.TextIO, other):
        self.compare_extents(fh, other)
        self.compare_relations(fh, other)

    def compare_extents(self, fh, other):
        sorted_extents_self = self.sorted_extents()
        sorted_extents_other = other.sorted_extents()
        # fh.write(f'\n==== {self.filename}\n')
        for tag in TAGS:
            extents_self = sorted_extents_self[tag]
            extents_other = sorted_extents_other[tag]
            # print("%s %s (%d:%s)" % (self.filename, tag.upper(), len(extents_self), len(extents_other)))
            # fh.write("\n%s (%d:%s)\n\n" % (tag.upper(), len(extents_self), len(extents_other)))
            self.compare_extents_for_tag(tag, fh, extents_self, extents_other)

    def compare_extents_for_tag(self, tag: str, fh: typing.TextIO,
                                tag_extents_self: list, tag_extents_other: list):
        queue_self = tag_extents_self.copy()
        queue_other = tag_extents_other.copy()
        while queue_self or queue_other:
            # print(len(queue_self), len(queue_other), queue_self[:1], queue_other[:1])
            if queue_self and queue_other:
                e1: BratExtent = queue_self[0]
                e2: BratExtent = queue_other[0]
                # print(len(queue_self), len(queue_other), queue_self[0], queue_other[0])
                if e1 == e2:
                    marker = self._comparison_marker(e1, e2)
                    if marker != '====':
                        fh.write("%s  %s  %s  %s  %s\n" %
                                 (f'{e1.filename:15}',
                                  f'{tag:12}',
                                  f'{e1.short_string():44}',
                                  marker,
                                  f'{e2.short_string():44}'))
                    if '*' in marker:
                        fh.write("%s  %s        %s\n" %
                                 (f'{"":29}',
                                  f'{e1.attributes_as_short_string():44}',
                                  f'{e2.attributes_as_short_string():44}'))
                    queue_self.pop(0)
                    queue_other.pop(0)
                elif e1 < e2:
                    fh.write(f'{e1.filename}  {tag:12}  {e1.short_string():44}  <<<<  \n')
                    queue_self.pop(0)
                elif e2 < e1:

                    fh.write(f'{e1.filename}  {tag:12}  {BLANK:44}  >>>>  {e2.short_string()}\n')
                    queue_other.pop(0)
                else:
                    print('\nhuh?\ne1', e1)
                    print('e2', e2)
                    exit()
            elif queue_self:
                e1: BratExtent = queue_self[0]
                fh.write(f'{e1.short_string():44}  <<<<  \n')
                queue_self.pop(0)
            elif queue_other:
                e2: BratExtent = queue_other[0]
                fh.write(f'{BLANK:44}  >>>>  {e2.short_string():44}\n')
                queue_other.pop(0)

    @staticmethod
    def _comparison_marker(e1, e2):
        positions_match = True if e1.positions == e2.positions else False
        attributes_match = True if e1.attributes == e2.attributes else False
        if positions_match and attributes_match:
            return '===='
        elif positions_match:
            return '*==*'
        elif attributes_match:
            return '<==>'
        else:
            return '<**>'

    def compare_relations(self, fh: typing.TextIO, other):
        pass

    def print_extents(self):
        for extent in self.extents.values():
            print(extent)
            for attribute in extent.attributes:
                print('   ', attribute)
        extents = self.sorted_extents()
        for tag in sorted(extents):
            dashes = '-' * 25
            print(f'\n{dashes} {tag} {dashes}\n')
            for extent in extents[tag]:
                print(extent)
                for attribute in extent.attributes:
                    print('   ', attribute)
                print()
