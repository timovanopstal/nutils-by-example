mdeqn
=====

A library for converting equations written in an easy-to-read format, e.g.

            ∂u (x)
              i
    ∫  v(x) ------ dx
     Ω       ∂x

into something parsable by latex,

    \begin{equation}
        ∫_{Ω} v(x) \frac{∂u_i(x)}{∂x} dx
    \end{equation}

Syntax
======

Fractions
---------

A fraction can be written as

    numerator
    -----------
    denominator

Make sure that the horizontal line is at least as long as the numerator and
denominator.  When writing nested fractions, the outer fraction should have the
longest horizontal line.  Example:

     a
     -
     b
    ---
     c
     -
     d

The horizontal line should be on the baseline.  Example:

        3
    2 + -
        4

Super- and subscripts
---------------------

Super- and subscripts should be written to the right and above and below,
respectively, the block of characters.  Example:

     c+d
    a
     b

It is not allowed to place characters in the space between super- and
subscripts:

     c+d
    a + 1  (incorrect!)
     b

This should be written as

     c+d
    a    + 1
     b

Super- and subscripts can be nested.  Example:

    a
     b
      c

becomes in latex

    a_{b_{c}}

Substacks
---------

Subscripts can be stacked.  An element of the stack is identified by an
underscore (`_`) in the first column at the lowest line of the element.  Lines
below the underscore belong to the next subscript, if any.  For readability,
blank lines may be added between subscripts.  Example:

    ∑
     _ i

       j
     _  i
     _ k

becomes in latex

    \begin{equation}
        ∑_\substack{i\\j_i\\k}
    \end{equation}

This requires the latex package `amsmath`.

Text
----

Characters between double quotation marks (`"`) are interpreted as text.
Example:

    ∂Ω
      "Neumann"

becomes in latex

    \begin{equation}
        ∂Ω_{\text{Neumann}}
    \end{equation}

Autobrackets
------------

Vertical lines of brackets (`(){}[]|.`) are interpreted as latex autobrackets,
which scale with the enclosed block.  Example:

    ( 1     )
    ( - + x )
    ( 2     )

The vertical lines should always be written in pairs.  A line of dots (`.`) can
be used as an invisble bracket.  Example:

    { 1     .
    { - + x .
    { 2     .

A super- and subscript should occur above and below, respectively, and to the
right of the vertical line.  Example:

             2
    ( 1     )
    ( - + x )
    ( 2     )

Arrays
------

When three blank columns or one blank row is detected, the block is interpreted
as an array.

    [ 1   2   3 ]
    [           ]
    [ 4   5   6 ]
    [           ]
    [ 7   8   9 ]

Elements of the array are centered by default.  Alternatively, an ampersand
(`&`) can be used to align on.  The ampersands should be aligned vertically.
Per column, every array element should have either no or one ampersand.
Example:

    {   a &= b        &"if" b < 1  .
    {                              .
    {             2                .
    { 2 a &= 2 + e    &"otherwise" .

Without enclosing brackets, arrays are not yet detected:

      a &= b
                    (not yet supported)
                2
    2 c &= d + e

Tags and labels
---------------

An equation label is added by appending `(@LABEL)` to the equation.  Example:

    a = b   (@example label)

becomes in latex

    \begin{equation}
        \label{example label}
        a = b
    \end{equation}

A custom tag should be specified as `(TAG@LABEL)`.  Example:

    a = b   (test@example label)

    \begin{equation}
        \tag{test}
        \label{example label}
        a = b
    \end{equation}

If the label part is omitted, the tag part is used as a label.

In this version all equations are numbered, but this might change in the future.

Multiline equations
-------------------

Line continuation is denoted by three dots, `...`.  Example:

    A very ...

    ... very very ...

    ... very ...

    ... long equation.  (@some label)

becomes in latex

    \begin{multline}
        \label{some label}
        A very \\ very very \\ very \\ long equation.
    \end{multline}

Wrappers
========

lualatex-wrapper
----------------

The `lualatex-wrapper` python script preprocesses the latex document and feeds
the result to `lualatex`.  Everything between `\begin{mdeqn}` and `\end{mdeqn}`
is converted to latex equations.  Example:

    \begin{mdeqn}
      a = b   (@foo)
    \end{mdeqn}

becomes

    \begin{equation}
        \label{foo}
        a = b
    \end{equation}

When using unicode math symbols, load the latex package `mathspec` or
`unicode-math`.  The `amsmath` package may be required for some features.  The
script only parses the document specified on the command line.  Source files
included via `input` or `include` are not supported.

pandoc-filter
-------------

The `pandoc-filter` python script can be used as a `pandoc` filter to parse
display equations marked as follows:

    $$mdeqn
                ∂u (x)
                  i
        ∫  v(x) ------ dx
         Ω       ∂x
    $$

To use this filter, add the option `--filter=/path/to/mdeqn/pandoc-filter` to
the pandoc command line.

python-markdown
---------------

Mdeqn includes a [python-markdown] extension for converting inline math
(latex-style) between single dollars and mdeqn-style equations between double
dollars, e.g.

    Weak form of Laplace's equation with basis $φ_i$:
    $$
        ∫  φ    φ    u  dΩ = 0
         Ω  i,k  j,k  j
    $$

to latex-style equations in a `script` tag, suitable for rendering by [mathjax].

The following example converts above markdown snippet to html using
[python-markdown]:

    import mdeqn
    import markdown
    md = markdown.Markdown(
        output_format='xhtml5', extensions=[mdeqn.MarkdownExtension()])

    html_body = md.convert('''
    Weak form of Laplace's equation with basis $φ_i$:
    $$
        ∫  φ    φ    u  dΩ = 0
         Ω  i,k  j,k  j
    $$
    ''')

    html_head = '''\
    <!doctype html>
    <html>
      <head>
        <meta charset='utf-8'/>
        <script src='https://cdn.mathjax.org/mathjax/2.6-latest/MathJax.js'>
        </script>
        <script>
          MathJax.Hub.Config({
            jax: ['input/TeX','output/CommonHTML'],
            TeX: {
              extensions: [
                'AMSmath.js','AMSsymbols.js','noErrors.js','noUndefined.js'],
              equationNumbers: {autoNumber: "AMS"}
            }
          });
        </script>
        <title>example</title>
      </head>
      <body>
    '''

    html_tail = '''\
      </body>
    </html>
    '''

    print(html_head + html_body + html_tail)

[python-markdown]: https://pythonhosted.org/Markdown/
[mathjax]: https://www.mathjax.org/
