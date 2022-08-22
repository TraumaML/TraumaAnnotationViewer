import re


def split(text):
    """Split a text in a list of triples, where each triple has a line number,
    a line, and an end-of-line or the empty string."""
    parts = re.split('(\n)', text)
    lines = parts[::2]
    eols = parts[1::2]
    eols_indexes = list(range(len(eols)))
    result = []
    for i, line in enumerate(lines):
        eol = eols[i] if i in eols_indexes else ''
        result.append((i, line, eol))
    return result


def index_characters(lines):
    """Take the result of splitting a text and build an index from character
    offsets to sentence numbers, with sentence numbers starting at 0 because
    line numbers from split() start at 0."""
    p = 0
    idx = {}
    for i, line, eol in lines:
        length = len(line) + len(eol)
        for x in range(length):
            idx[p] = i
            p += 1
    return idx


if __name__ == '__main__':

    test_text = open('backup.sh').read()
    test_lines = split(test_text)
    idx = index_characters(test_lines)

    for i, line, eol in test_lines:
        print("%3d [%s, %s]" % (i, repr(line), repr(eol)))

    current_sentence_no = -1
    for p in sorted(idx.keys()):
        sentence_no = idx[p]
        if current_sentence_no != sentence_no:
            print("\n%d -" % sentence_no, end='')
            current_sentence_no = sentence_no
        print(" %d" % p, end='')
    print()
