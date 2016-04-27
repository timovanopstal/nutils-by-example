# Copyright (c) 2014 Joost van Zwieten
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import re
import unicodedata


class Block:

    def __init__( self, lines ):

        lines = tuple( map( self._get_characters, lines ) )
        n = 0 if len( lines ) == 0 else max( map( len, lines ) )
        self._lines = tuple( self._ljust( line, n ) for line in lines )
        self.shape = ( len( self._lines ), n )

    @staticmethod
    def _get_characters( line ):

        if not isinstance( line, str ):
            return line

        # Convert unicode strings to a list of characters.
        # Combining characters are merged with the preceeding character.

        characters = []
        for ch in line:
            if unicodedata.combining( ch ) and len( characters ) > 0:
                characters[-1] = characters[-1] + ch
            else:
                characters.append( ch )
        return tuple( characters )

    @staticmethod
    def _ljust( line, n ):

        if len( line ) < n:
            return line + ( ' ', ) * ( n - len( line ) )
        else:
            return line

    def __getitem__( self, item ):

        r, c = item
        if not isinstance( r, slice ):
            if not isinstance( c, slice ):
                return self._lines[r][c]
            else:
                r = slice( r, r + 1 )
        elif not isinstance( c, slice ):
            c = slice( c, c + 1 )

        return Block( line[c] for line in self._lines[r] )

    def is_empty( self ):

        return all( ch == ' ' for line in self._lines for ch in line )

    def __str__( self ):

        return '\n'.join( ''.join( line ) for line in self._lines )

    def __repr__( self ):

        return '\n' + '\n'.join( '>' + ''.join( line ) for line in self._lines )

    def __iter__( self ):

        for line in self._lines:
            for ch in line:
                yield ch

    def count( self, ch ):

        return sum( 1 for a in self if a == ch )

    def as_string( self ):

        if self.shape[0] == 0:
            raise ValueError( 'Cannot convert block with zero rows.' )
        elif self.shape[0] > 1:
            raise ValueError( 'Cannot convert block with more than one row to a string.' )

        return ''.join( self._lines[0] )


def _determine_fraction_length( line ):

    for c in range( line.shape[1] ):
        if line[0,c] not in ( '-', '―' ):
            return c
    return line.shape[1]


def parse_sub( block, *, outer_line = False, equation_label_prefix = '', display_style_in_array ):

    if not isinstance( block, Block ):
        block = Block( block )

    latex = ''
    baseline = None
    fraction_symbols = ( '-', '―' )
    is_multiline = False
    has_custom_tag = False

    # skip white columns

    while block.shape[1] and block[:,0].is_empty():
        block = block[:,1:]

    while block.shape[1]:

        # test for blocks surrounded by vertical lines

        proceed = True
        check_limits = True

        if set( block[:,0] ) <= set( ( ' ', '_' ) ) and block[:,0].count( '_' ) > 1:

            # substack

            rows = []
            start = 0
            for stop in range( block.shape[0] ):
                if block[stop,0] == '_':
                    rows.append( slice( start, stop + 1 ) )
                    start = stop + 1
            if start < block.shape[0] and not block[start:,:].is_empty():
                raise ValueError( 'space below last substack entry should be empty' )

            latex = '\\substack{{{}}}'.format( '\\\\'.join( parse_sub( block[ row, 1: ], display_style_in_array = display_style_in_array ) for row in rows ) )
            proceed = False
            check_limits = False
            block = block[block.shape[0]:,:]

        if proceed:

            vert_top = 0
            while block[vert_top,0] == ' ':
                vert_top += 1
            vert_bottom = block.shape[0]
            while block[vert_bottom-1,0] == ' ':
                vert_bottom -= 1

            if vert_top + 1 < vert_bottom and len( set( block[vert_top:vert_bottom,0] ) ) == 1:
                proceed = False

                ch_left = block[vert_top,0]
                for j in range( 1, block.shape[1] ):
                    ch_right = block[vert_top,j]
                    if ch_right != ' ' and all( block[i,j] == ( ch_right if i in range( vert_top, vert_bottom ) else ' ' ) for i in range( block.shape[0] ) ):
                        break
                else:
                    raise ValueError( 'expected closing vertical line' )

                latex_matrix = parse_matrix( block[vert_top:vert_bottom,1:j], display_style_in_array = display_style_in_array )
                if ch_left in ( '{', '}' ):
                    ch_left = '\\' + ch_left
                if ch_right in ( '{', '}' ):
                    ch_right = '\\' + ch_right
                if ch_left != '.' or ch_right != '.':
                    latex += '\\left{}{}\\right{}'.format( ch_left, latex_matrix, ch_right )
                else:
                    latex += latex_matrix

                block = block[:,j+1:]
                limit_top = vert_top
                limit_bottom = vert_bottom

        if proceed:

            if baseline is None:

                # find baseline

                symbols = []
                while block.shape[1] > 0:
                    for r in range( block.shape[0] ):
                        if block[r,0] in fraction_symbols:
                            symbols.append( ( r, _determine_fraction_length( block[r,:] ) ) )
                        elif block[r,0] != ' ':
                            symbols.append( ( r, 1 ) )
                    if len( symbols ) > 0:
                        break
                    block = block[:,1:]
                else:
                    raise ValueError( 'empty block' )

                baseline, length = max( symbols, key = lambda item: ( item[1], block[item[0],0] in fraction_symbols ) )
                for r, l in symbols:
                    if r == baseline:
                        continue
                    if block[baseline,0] in fraction_symbols and block[r,0] not in fraction_symbols:
                        continue
                    if l == length:
                        raise ValueError( 'multiple possible baselines in {!r}'.format( block ) )

            limit_top = baseline
            limit_bottom = baseline + 1

            if block[baseline,0] == ' ':
                raise ValueError( 'unexpected symbol, not on baseline' )

        if proceed and outer_line == True:

            if block[baseline,0:2].as_string() == '\\\\' and set( block[:,0:2] ) <= set( ' \\' ):
                is_multiline = True
                proceed = False
                check_limits = False
                latex += '\\\\'
                block = block[:,2:]

        if proceed:

            while True:

                if not block[:baseline,:].is_empty() or not block[baseline+1:,:].is_empty():
                    break

                match = re.match( r'\(([^()@]*)@([^()@]*)\)[ ]*', block[baseline,:].as_string() )
                if not match:
                    break

                tag = None
                label = None
                if match.group(1):
                    tag = match.group(1)
                    label = match.group(1)
                if match.group(2):
                    label = match.group(2)

                if tag:
                    has_custom_tag = True
                    latex += '\\tag{{{}}}'.format( tag )
                if label:
                    latex += '\\label{{{}{}}}'.format( equation_label_prefix, label )
                block = block[:,block.shape[1]:]
                proceed = False
                break

#               if block.shape[1] < 2 or block[baseline,0] != '(' or block[baseline,1] != '@':
#                   break
#
#               for i in range( 2, block.shape[1] ):
#                   if block[baseline,i] == ')':
#                       break
#               else:
#                   break
#
#               if not block[baseline,i+1:].is_empty():
#                   break
#
#               if not block[:baseline,:].is_empty() or not block[baseline+1:,:].is_empty():
#                   break
#
#               latex += '\\label{{{}}}'.format( ''.join( block[baseline,j] for j in range( 2, i ) ) )
#               block = block[:,block.shape[1]:]
#               proceed = False
#               break

        if proceed:

            if block[baseline,0] in fraction_symbols and (
                        _determine_fraction_length( block[baseline,:] ) > 1
                    or
                        not block[:baseline,0].is_empty() and not block[baseline+1:,0].is_empty()
                    ):
                # fraction
                fraction_length = _determine_fraction_length( block[baseline,:] )
                numerator = parse_sub( block[:baseline,:fraction_length], display_style_in_array = display_style_in_array )
                denominator = parse_sub( block[baseline+1:,:fraction_length], display_style_in_array = display_style_in_array )
                latex += '\\frac{{{}}}{{{}}}'.format( numerator, denominator )
                block = block[:,fraction_length:]
                check_limits = False
                proceed = False

        if proceed:

            if block[baseline,0] == '"':
                for i in range( 1, block.shape[1] ):
                    if block[baseline,i] == '"':
                        break
                else:
                    raise ValueError( 'could not find matching "' )

                if not block[:baseline,0:i+1].is_empty():
                    raise ValueError( 'unexpected symbols above text' )
                if not block[baseline+1:,0:i+1].is_empty():
                    raise ValueError( 'unexpected symbols below text' )

                latex += '\\text{{{}}}'.format( block[baseline,1:i].as_string() )
                block = block[:,i+1:]
                proceed = False

        if proceed:
            if not block[:baseline,0].is_empty():
                raise ValueError( 'unexpected symbol above symbol on baseline' )
            if not block[baseline+1:,0].is_empty():
                raise ValueError( 'unexpected symbol below symbol on baseline' )
            if block[baseline,0] in ( '(', ')' ):
                if block[baseline,0] == '(':
                    latex += '\\left('
                    check_limits = False
                else:
                    latex += '\\right)'
            else:
                latex += block[baseline,0]
            block = block[:,1:]

        if check_limits:
            i = 0
            while i < block.shape[1]:
                if not block[limit_top:limit_bottom,i].is_empty():
                    break
                i += 1
            while i > 0 and block[:,i-1].is_empty():
                i -= 1
            if i > 0:
                superscript = block[:limit_top,:i]
                if not superscript.is_empty():
                    latex += '^{{{}}}'.format( parse_sub( superscript, display_style_in_array = display_style_in_array ) )
                subscript = block[limit_bottom:,:i]
                if not subscript.is_empty():
                    latex += '_{{{}}}'.format( parse_sub( subscript, display_style_in_array = display_style_in_array ) )
                block = block[:,i:]

        # remove white columns

        while block.shape[1] > 0 and block[:,0].is_empty():
            block = block[:,1:]
            latex += ' '

    if outer_line:
        return latex.strip(), is_multiline, has_custom_tag
    else:
        return latex.strip()


def parse_matrix( block, *, display_style_in_array ):

    column_white = 3
    row_white = 1

    rows = []

    for start in range( block.shape[0] ):
        if not block[start,:].is_empty():
            break
    else:
        raise ValueError( 'empty block' )

    stop = start + 1
    while stop < block.shape[0]:
        if not block[stop:stop+row_white,:].is_empty():
            stop += 1
            continue
        rows.append( slice( start, stop ) )
        start = stop + row_white
        while start < block.shape[0] and block[start,:].is_empty():
            start += 1
        stop = start + 1
    if not block[start:stop,:].is_empty():
        rows.append( slice( start, stop ) )

    columns = []

    for start in range( block.shape[1] ):
        if not block[:,start].is_empty():
            break
    else:
        raise ValueError( 'empty block' )

    stop = start + 1
    while stop < block.shape[1]:
        if not block[:,stop:stop+column_white].is_empty():
            stop += 1
            continue
        columns.append( slice( start, stop ) )
        start = stop + column_white
        while start < block.shape[1] and block[:,start].is_empty():
            start += 1
        stop = start + 1
    if not block[:,start:stop].is_empty():
        columns.append( slice( start, stop ) )

    if len( rows ) == 1 and len( columns ) == 1:
        return parse_sub( block, display_style_in_array = display_style_in_array )

    # If each cell in a column has an `&`, vertically aligned, (on the
    # baseline, but we don't check this yet), use right alignment left of `&`
    # and left alignment right of `&`.  Otherwise use centre alignment.

    alignment = []
    columns, _columns = [], columns
    for column in _columns:
        for i in range( column.start, column.stop ):
            if '&' in block[:,i] and set( block[:,i] ) <= { '&', ' ' }:
                alignment.append( 'r@{}l' )
                columns.extend([ slice( column.start, i ), slice( i+1, column.stop ) ])
                break
        else:
            alignment.append( 'c' )
            columns.append( column )

    cell_prefix = ''
    if display_style_in_array:
        cell_prefix='\\displaystyle '
    latex_elements = '\\\\'.join( '&'.join( cell_prefix + parse_sub( block[row,column], display_style_in_array = display_style_in_array ) for column in columns ) for row in rows )
    return '\\begin{{array}}{{{}}}{}\\end{{array}}'.format( ''.join( alignment ), latex_elements )


def parse( block, *, break_line_on_dots = True, equation_label_prefix = '', starred_custom_tag = False, display_style_in_array = False ):

    if not isinstance( block, Block ):
        block = Block( block )

    row_groups = []

    # find row groups separated by two white rows

    for start in range( block.shape[0] ):
        if not block[start,:].is_empty():
            break
    else:
        raise ValueError( 'empty block' )

    stop = start + 1
    while stop < block.shape[0]:
        if not block[stop:stop+1,:].is_empty():
            stop += 1
            continue
        row_groups.append( block[start:stop,:] )
        start = stop
        while start < block.shape[0] and block[start,:].is_empty():
            start += 1
        stop = start + 1
    if not block[start:stop,:].is_empty():
        row_groups.append( block[start:stop,:] )

    assert len( row_groups ) > 0

    is_multiline = False
    has_custom_tag = False
    latex_groups = []
    for i, row_group in enumerate( row_groups ):
        latex, row_is_multiline, row_has_custom_tag = parse_sub( row_group, outer_line = True, equation_label_prefix = equation_label_prefix, display_style_in_array = display_style_in_array )
        latex = latex.strip()
        if row_is_multiline:
            is_multiline = True
        if row_has_custom_tag:
            has_custom_tag = True
        if i + 1 < len( row_groups ):
            assert latex.endswith( '...' )
            latex = latex[:-3]
            if break_line_on_dots:
                latex += '\\\\'
                is_multiline = True
        if i > 0:
            assert latex.startswith( '...' )
            latex = latex[3:]
        latex_groups.append( latex )

    if not is_multiline:
        return '\\begin{{equation{0}}}{1}\\end{{equation{0}}}'.format( '*' if has_custom_tag and starred_custom_tag else '', '\n'.join( latex_groups ) )
    else:
        return '\\begin{{multline{0}}}{1}\\end{{multline{0}}}'.format( '*' if has_custom_tag and starred_custom_tag else '', '\n'.join( latex_groups ) )


try:
    import markdown
except ImportError:
    pass
else:

    class InlineMath(markdown.inlinepatterns.Pattern):

        def __init__(self, *args, **kwargs):

            super().__init__(r'(\$)([^\$]+)\2', *args, **kwargs)

        def handleMatch(self, m):

            el = markdown.util.etree.Element('script', dict(type='math/tex'))
            el.text = markdown.util.AtomicString(m.group(3))
            return el


    class DisplayMath(markdown.inlinepatterns.Pattern):

        def __init__(self, parse_kwargs, *args, **kwargs):

            self.parse_kwargs = parse_kwargs
            super().__init__(r'(\s*\${2})\n([^\$]+)\n\2', *args, **kwargs)

        def handleMatch(self, m):

            el = markdown.util.etree.Element(
                'script', dict(type='math/tex; mode=display'))
            el.text = markdown.util.AtomicString(
                parse(m.group(3).split('\n'), **self.parse_kwargs))
            return el


    class MarkdownExtension(markdown.Extension):

        def __init__(self, **kwargs):
            self.config = dict(
                equation_label_prefix=['', 'prefix for equation labels'],
                starred_custom_tag=
                    [False, 'use equation* if a custom tag is given'],
                display_style_in_array=
                    [False, 'use display style in arrays'],
            )
            super().__init__(**kwargs)

        def extendMarkdown(self, md, md_globals):

            md.inlinePatterns.add(
                'displaymath', DisplayMath(self.getConfigs(), md), '<escape')
            md.inlinePatterns.add('inlinemath', InlineMath(md), '<escape')


    def makeExtension(configs=[]):

        return MarkdownExtension(configs)


def test():

    def _parse( lines ):
        lines = lines.split('\n')
        for line in lines:
            print( '>>>', line )
        print()
        print( '<<<', parse( lines, break_line_on_dots = False, starred_custom_tag = True ) )
        print()

    _parse( r'''
    1    π
    - + e
    2           ∂f(x)
    ------ + ∫  ----- dλ(x)    (@label)
       2      Ω  ∂x

    ''' )

    _parse( r'''
                    ∂F   (q(x))
                      kij
    ∑  ∫  ∑   v (x) ----------- dλ(x) + ...
     E  E  ik  i       ∂x
                         k

                                               ∂φ (τ)
                                 1               j
    ... ∑      ∫     ∑    v (x) ∫  F    (q(x)) ------ n (x) dτ dλ(x) = ...
         edges  edge  ijk  i     0  kij          ∂τ    k


                                  ... ∑  ∫  ∑  v (x) s (q) dλ(x).   (H1@)
                                       E  E  i  i     i

    ''' )

    _parse( r'''
    ∑        f(i,j)  (H2@test)
     _ i ∈ I
     _ j ∈ J
    ''' )

    _parse( r'''
    ∑         f(i,j)
     _ i ∈ I
            2
     _ j ∈ J
    ''' )

    _parse( r'''
    [ 1   2   3 ]
    [           ]
    [ 4   5   6 ]
    [           ]
    [ 7   8   9 ]
    ''' )

    _parse( r'''
                                    2
                |         "low"    |
             ∫  | q (x)- q     (x) |  dλ(x)
              E |  i      i        |
    S(q,E) = ―――――――――――――――――――――――――――――― .
                               2
                      |       |
                   ∫  | q (x) |  dλ(x)
                    E |  i    |
    ''' )

    _parse( r'''
          z       z
         y      y
        x   ≠ {x }
    ''' )

    _parse( r'''
    { a &= "foo"   c &= "spam" .
    {                          .
    { b &= "bar"   d &= "eggs" .
    ''' )

    _parse( r'''
    y(x) = ∑  ŷ  L (x)
            j  j  j
    ''' )


if __name__ == '__main__':

    test()
