#!/usr/bin/env python
# Copyright 2026 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Implementation of orbl --alignmentSets.
Create an html file containing all visible alnsets from CodAlignView plus any hidden ones
    for which there orblq is defined. Add a column "ORBLq" that specifies for
    each alnset whether orblq is defined. Then open it in a browser.
"""
from __future__ import division, print_function
import re, os, sys, tempfile, webbrowser
from .DownloadLocalAlignment import get_from_url
from .orblq import get_untrans_ORBLv_dict_name
try:
    # Py3
    from html import escape as _html_escape
    from html import unescape as _html_unescape
except ImportError:
    # Py2
    import cgi
    import HTMLParser
    _html_escape = lambda s: cgi.escape(s, quote=True)
    _html_unescape = HTMLParser.HTMLParser().unescape

# For unit testing, don't invoke browser; just print temporary html file name to stderr.
PRINT_HTML_FILE_NAME = False

AlnsetsURL = 'https://data.broadinstitute.org/compbio1/cav.php?Alnsets'

def show_alignment_sets() :
    alnsetsWithoutHiddenHtml = get_from_url(AlnsetsURL)
    visibleAlnsetSet = set(AlnsetsHtmlEditor(alnsetsWithoutHiddenHtml).get_row_names())
    alnsetsWithHiddenHtml = get_from_url(AlnsetsURL + '&sha')
    alnsetsWithHiddenHtmlEditor = AlnsetsHtmlEditor(alnsetsWithHiddenHtml)
    allAlnsets = alnsetsWithHiddenHtmlEditor.get_row_names()
    orblqAlnsetSet = set(alnset for alnset in allAlnsets
                         if os.path.exists(get_untrans_ORBLv_dict_name(alnset)))
    keepAlnsets = visibleAlnsetSet | orblqAlnsetSet
    for alnset in allAlnsets :
        if not alnset in keepAlnsets :
            alnsetsWithHiddenHtmlEditor.remove_row(alnset)
    alnsetsWithHiddenHtmlEditor.add_column(
        'ORBLq', 1, ['yes' if alnset in orblqAlnsetSet else ''
                     for alnset in alnsetsWithHiddenHtmlEditor.get_row_names()])
    fhandle, fname = tempfile.mkstemp(prefix = 'temp_', suffix = '.html', dir = None,
                                      text = True)
    os.write(fhandle, alnsetsWithHiddenHtmlEditor.get_html().encode('utf-8'))
    os.close(fhandle)

    if PRINT_HTML_FILE_NAME :
        # For unit testing, don't invoke browser; just print temporary html file name.
        print(fname, file = sys.stderr)
        return

    if os.name == 'nt':
        # Normalize backslashes to forward slashes for a URL.
        url = 'file:///' + fname.replace('\\', '/')
    else:
        url = 'file://' + fname
    ok = webbrowser.open(url, new=2)  # new=2 requests a new tab if possible
    # os.remove(fname) # Don't remove because we get here before browser opens it
    if not ok:
        print('Warning: webbrowser could not confirm it opened %s.' % fname,
              file=sys.stderr)
        raise SystemExit(1)

class AlnsetsHtmlEditor(object):
    # Class for manipulating the html text containing a CodAlignView Alignment Sets table.
    # Mostly written by ChatGPT 5.2
    def __init__(self, htmlText) :
        """
        Read html text in format produced by CodAlignView's "Alignment Sets" command.
        """
        self.html = htmlText
        # Parse and cache a structured representation of the table.
        self._parse()

    def _parse(self):
        """
        Internal: parse the single CodAlignView table into:
          - self._table_open, self._table_close
          - self._thead_full (original thead), and replacement points
          - self._headers: list of dicts {attrs, inner_html}
          - self._rows: list of dicts {tr_attrs, cells:[{tag, attrs, inner_html}]}
        """
        # Extract the first <table> ... </table> block
        m = re.search(r'(?is)(<table\b[^>]*>)(.*?)(</table>)', self.html)
        if not m:
            raise ValueError("No <table>...</table> found in HTML.")
        self._table_open, table_inner, self._table_close = (
            m.group(1), m.group(2),  m.group(3))

        # Extract <thead>...</thead>
        mt = re.search(r'(?is)(<thead\b[^>]*>)(.*?)(</thead>)', table_inner)
        if not mt:
            raise ValueError("No <thead>...</thead> found inside table.")
        thead_open, thead_inner, thead_close = mt.group(1), mt.group(2), mt.group(3)
        self._thead_open = thead_open
        self._thead_close = thead_close

        # Find the first header <tr> within thead
        mtr = re.search(r'(?is)(<tr\b[^>]*>)(.*?)(</tr>)', thead_inner)
        if not mtr:
            raise ValueError("No header <tr> found inside <thead>.")
        self._thead_prefix = thead_inner[:mtr.start(1)]
        self._thead_tr_open = mtr.group(1)
        header_tr_inner = mtr.group(2)
        self._thead_tr_close = mtr.group(3)
        self._thead_suffix = thead_inner[mtr.end(3):]

        # Parse header cells (<th ...>...</th>)
        headers = []
        for tag, attrs, inner in re.findall(r'(?is)<(th)\b([^>]*)>(.*?)</\1>',
                                            header_tr_inner):
            headers.append({'attrs': attrs, 'inner_html': inner})
        if not headers:
            raise ValueError("No <th> header cells found in header row.")
        self._headers = headers

        # Body is everything after </thead> within table_inner
        body = table_inner[mt.end(3):]

        # Parse each body row <tr ...>...</tr>
        rows = []
        for tr_attrs, tr_inner in re.findall(r'(?is)<tr\b([^>]*)>(.*?)</tr>', body):
            cells = []
            for tag, attrs, inner in re.findall(r'(?is)<(th|td)\b([^>]*)>(.*?)</\1>',
                                                tr_inner):
                cells.append({'tag': tag.lower(), 'attrs': attrs, 'inner_html': inner})
            # Skip rows with no cells (defensive)
            if cells:
                rows.append({'tr_attrs': tr_attrs, 'cells': cells})
        if not rows:
            raise ValueError("No body <tr> rows found after <thead>.")

        self._rows = rows

    def _serialize(self):
        """
        Internal: build updated HTML text from cached structures and replace the original
        table.
        """
        # Rebuild header <tr>...</tr>
        header_cells_html = []
        for h in self._headers:
            attrs = h['attrs'] or ''
            header_cells_html.append('<th{0}>{1}</th>'.format(attrs, h['inner_html']))
        header_tr_html = (self._thead_tr_open + '\n'.join(header_cells_html) +
                          self._thead_tr_close)

        thead_html = (
            self._thead_open +
            self._thead_prefix +
            header_tr_html +
            self._thead_suffix +
            self._thead_close
        )

        # Rebuild body rows
        body_rows_html = []
        for r in self._rows:
            tr_attrs = r['tr_attrs'] or ''
            cells_html = []
            for c in r['cells']:
                attrs = c['attrs'] or ''
                tag = c['tag']
                cells_html.append('<{0}{1}>{2}</{0}>'.format(tag, attrs, c['inner_html']))
            body_rows_html.append('<tr{0}>\n{1}\n</tr>'.format(tr_attrs,
                                                               '\n'.join(cells_html)))

        table_html = (self._table_open + thead_html + '\n'.join(body_rows_html) + '\n' +
                      self._table_close)

        m = re.search(r'(?is)(<table\b[^>]*>)(.*?)(</table>)', self.html)
        start, end = m.start(1), m.end(3)
        self.html = self.html[:start] + table_html + self.html[end:]

    def get_column_names(self):
        return [cell['inner_html'] for cell in self._headers]

    def get_row_names(self):
        """
        Return a list containing the row names (i.e., the first column of the table).
        """
        names = []
        for r in self._rows:
            first = r['cells'][0]
            # In your sample, first cell is <th align="left">ROWNAME</th>
            # Treat inner_html as text and unescape entities; strip whitespace.
            text = re.sub(r'(?is)<.*?>', '', first['inner_html'])
            text = _html_unescape(text).strip()
            names.append(text)
        return names

    def get_values_from_row_name(self, rowName):
        for r in self._rows :
            if r['cells'][0]['inner_html'] == rowName :
                return [cell['inner_html'] for cell in r['cells']]
        else :
            raise ValueError('No row with name %s.' % rowName)

    def add_column(self, columnName, beforeIndex, values):
        """
        Add a column to the table with specified columnName in the header, and the
            specified list of values. Check that the number of values equals the number
            of rows in the table (not counting the header).
        """
        if values is None:
            raise ValueError("values must be a list (got None).")
        if len(values) != len(self._rows):
            raise ValueError(
                "values length ({0}) does not match number of rows ({1})."
                .format(len(values), len(self._rows))
            )

        ncols = len(self._headers)
        if beforeIndex < 0 or beforeIndex > ncols:
            raise IndexError(
                "beforeIndex out of range: {0} (valid: 0..{1})"
                .format(beforeIndex, ncols)
            )

        # Insert header cell
        safe_header = _html_escape(str(columnName))
        self._headers.insert(beforeIndex, {'attrs': '', 'inner_html': safe_header})

        # Insert each row cell at the same index. Use <td> for inserted column.
        for i, r in enumerate(self._rows):
            safe_val = _html_escape(str(values[i]))
            new_cell = {'tag': 'td', 'attrs': '', 'inner_html': safe_val}

            # Defensive: if a row has fewer cells than headers (unexpected), pad with
            # empty <td>.
            while len(r['cells']) < ncols:
                r['cells'].append({'tag': 'td', 'attrs': '', 'inner_html': ''})

            r['cells'].insert(beforeIndex, new_cell)

        # Rebuild HTML
        self._serialize()

    def remove_row(self, rowName):
        """
        Remove the row with specified name (name is the value in the first column)
        """
        # Match against decoded/stripped first-cell text (same method as get_row_names).
        target = _html_unescape(str(rowName)).strip()

        new_rows = []
        removed = 0
        for r in self._rows:
            first = r['cells'][0]
            text = re.sub(r'(?is)<.*?>', '', first['inner_html'])
            text = _html_unescape(text).strip()
            if text == target:
                removed += 1
                continue
            new_rows.append(r)

        if removed == 0:
            raise KeyError("Row not found: {0}".format(rowName))
        if removed > 1:
            # Unusual for this table, but explicitly flag it rather than silently removing
            # multiples.
            raise ValueError(
                "Multiple rows matched rowName {0}.".format(rowName))

        self._rows = new_rows
        self._serialize()

    def get_html(self) :
        return self.html
