"""HTML Utitilies.

Adapted from ~/Dropbox/code/python/html.py.

Utilities to help create HTML code. Contains classes that help build and print
HTML tags.

The two main classes are Tag and Text. To create a Text simply hand it the
string:

>>> text = Text("hello")

To wrap it in a paragraph and print the paragraph:

>>> p = Tag('p')
>>> p.add(text)
>>> print(p)
<p>hello</p>

You can add more text

>>> p.add(Text(" and goodbye"))
>>> print(p)
<p>hello and goodbye</p>

Tags have daughters:

>>> div = Tag('div', dtrs=Tag('p', dtrs=[Text("hello"), Text(" world")]))
>>> print(div)
<div>
<p>hello world</p></div>

Adding attributes:

>>> div = Tag('div',
...           attrs={'class': 'example'},
...           dtrs=Tag('p', dtrs=[Text("hello"), Text(" world")]))
>>> print(div)
<div class='example'>
<p>hello world</p></div>

"""

import io


BLOCK_TAGS = {'html', 'head', 'style', 'body',
              'table', 'tr', 'div', 'blockquote', 'ol', 'ul'}


class HtmlObject(object):

    """Abstract class."""

    def __str__(self):
        buffer = io.StringIO()
        self.write(buffer)
        return buffer.getvalue()


class Text(HtmlObject):

    """Class that holds some text. The text can include any html."""

    def __init__(self, text):
        self.text = str(text)

    def write(self, buffer):
        buffer.write(self.text)


class Tag(HtmlObject):

    def __init__(self, tag, nl=False, class_=None, attrs=None, dtrs=None):
        self.tag = tag
        self.nl = nl
        self.attrs = {}
        self.dtrs = []
        if attrs is not None:
            self.attrs = attrs
        if class_ is not None:
            self.attrs['class'] = class_
        if dtrs is not None:
            if isinstance(dtrs, list):
                self.dtrs = dtrs
            else:
                self.dtrs = [dtrs]

    def is_block(self):
        return True if self.tag.lower() in BLOCK_TAGS else False

    def add(self, *dtrs):
        """Add daughters, each of which is either a Tag or Text instance."""
        for dtr in dtrs:
            self.dtrs.append(dtr)

    def insert(self, index, dtr):
        """Insert a daughter at the given index."""
        self.dtrs.insert(index, dtr)

    def write(self, buffer):
        attrs = self._attribute_string()
        if not self.dtrs:
            buffer.write("<%s%s/>" % (self.tag, attrs))
        else:
            buffer.write("<%s%s>" % (self.tag, attrs))
            if self.is_block():
                buffer.write("\n")
            for dtr in self.dtrs:
                dtr.write(buffer)
            #if self.is_block():
            #    buffer.write("\n")
            buffer.write("</%s>" % self.tag)
        if self.nl:
            buffer.write("\n")

    def _attribute_string(self):
        attrs = ''
        if self.attrs:
            pairs = self.attrs.items()
            attrs = ' ' + ' '.join(["%s='%s'" % (k, v) for k, v in pairs])
        return attrs


## Convenience classes to make tag creation less verbose

class Span(Tag):
    def __init__(self, **attrs):
        super().__init__('span', **attrs)

class P(Tag):
    """Simple paragraph with just some text."""
    def __init__(self, text, **attrs):
        super().__init__('p', **attrs)
        self.add(Text(text))

class Href(Tag):
    def __init__(self, href, content):
        super().__init__('a', nl=False,
                         attrs={'href': href}, dtrs=[Text(content)])

class Anchor(Tag):
    def __init__(self, name):
        super().__init__('a', nl=False, attrs={'name': name})
        
class TR(Tag):
    def __init__(self, **attrs):
        super().__init__('tr', **attrs)

class TD(Tag):
    def __init__(self, **attrs):
        super().__init__('td', **attrs)

class H(Tag):
    def __init__(self, tag, content):
        super().__init__(tag)
        self.nl = True
        self.dtrs = [Text(content)]

class H1(H):
    def __init__(self, content):
        super().__init__('h1', content)

class H2(H):
    def __init__(self, content):
        super().__init__('h2', content)

class H3(H):
    def __init__(self, content):
        super().__init__('h3', content)

class H4(H):
    def __init__(self, content):
        super().__init__('h4', content)


# Some constants for fixed html text

NL = Text('\n')
HR = Text('\n<hr/>\n')
BR = Text('\n<br/>\n')
SPACE = Text(' ')



if __name__ == '__main__':

    import doctest
    doctest.testmod()
