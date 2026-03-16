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
    for which orblq is defined. Add a column "ORBLq" that specifies for each alnset
    whether orblq is defined. Then open it in a browser.
"""
from __future__ import division, print_function
import re, os, sys, tempfile, webbrowser
from .DownloadLocalAlignment import get_from_url
from .orblq import get_orblq_alnsets
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

# The following flag is to be set in unit testing to force printing of the name of the
# html file and a text version of the Alignment Set list rather than opening a browser
# window.
SHOW_ALNSETS_AS_TEXT = False

AlnsetsURL = 'https://data.broadinstitute.org/compbio1/cav.php?Alnsets'
TreeURLprefix = 'https://data.broadinstitute.org/compbio1/CodAlignViewFiles'

def show_alignment_sets() :
    """
    Create a temporary html file similar to CodAlignView's Alignment Set page, except:
        - Include those hidden Alignment Sets for which ORBLq is implemented.
        - Add a column indicating for each Alignment Set whether ORBLq is implemented.
    Normally, open the html file in a browser, but if SHOW_ALNSETS_AS_TEXT then print the
        name of the html file and a text version of the Alignment Set list to stderr
        rather than opening a browser window. Also do that if the browser fails to open.
    """
    alnsetsWithoutHiddenHtml = get_from_url(AlnsetsURL)
    visibleAlnsetSet = set(AlnsetsHtmlEditor(alnsetsWithoutHiddenHtml).get_row_names())
    alnsetsWithHiddenHtml = get_from_url(AlnsetsURL + '&sha')
    alnsetsWithHiddenHtmlEditor = AlnsetsHtmlEditor(alnsetsWithHiddenHtml)
    allAlnsets = alnsetsWithHiddenHtmlEditor.get_row_names()
    orblqAlnsetSet = set(get_orblq_alnsets())
    for alnset in allAlnsets :
        # Show alnset only if it is normally visible or orblq is implemented
        if alnset not in visibleAlnsetSet and alnset not in orblqAlnsetSet :
            alnsetsWithHiddenHtmlEditor.remove_row(alnset)
    alnsetsWithHiddenHtmlEditor.add_column(
        'ORBLq', 1, ['yes' if alnset in orblqAlnsetSet else ''
                     for alnset in alnsetsWithHiddenHtmlEditor.get_row_names()])
    fhandle, fname = tempfile.mkstemp(prefix = 'temp_', suffix = '.html', dir = None,
                                      text = True)
    os.write(fhandle, alnsetsWithHiddenHtmlEditor.get_html().encode('utf-8'))
    os.close(fhandle)

    if os.name == 'nt': # pragma: no cover
        # Normalize backslashes to forward slashes for a URL.
        url = 'file:///' + fname.replace('\\', '/')
    else:
        url = 'file://' + fname

    if SHOW_ALNSETS_AS_TEXT or not webbrowser.open(url, new = 2) :
        print('Alignment Set list was written to %s.' % fname, file = sys.stderr)
        print('However, could not open web browser. Printing a text version instead.',
              file = sys.stderr)
        print(get_alnsets_as_text(alnsetsWithHiddenHtmlEditor),
              file = sys.stderr) # Prints a blank line after table
        print('Tree files in Newick and pdf format can be found in:\n' +
              '    ' + TreeURLprefix + '/TreeNHs/{Alignment Set}.nh\n' +
              '    ' + TreeURLprefix + '/TreePdfs/{Alignment Set}.pdf',
              file = sys.stderr)

def get_alignment_set_assembly(alnset) :
    """ Return the reference assembly for this alignment set. """
    alnsetsWithHiddenHtml = get_from_url(AlnsetsURL + '&sha')
    alnsetsWithHiddenHtmlEditor = AlnsetsHtmlEditor(alnsetsWithHiddenHtml)
    assert alnsetsWithHiddenHtmlEditor.get_column_names()[1] == 'Reference<br>Assembly', (
        'Assembly is not column 1.')
    return alnsetsWithHiddenHtmlEditor.get_values_from_row_name(alnset)[1]

def get_alnsets_as_text(alnsetsHtmlEditor) :
    """
    Return a text representation of an AlnsetsHtmlEditor as follows:
    - For each column determine the maximum number of characters of any value in that
      column.
    - Create a text string in which the row headers are centered in fields of width one
      more than the maximum width for that column.
    - For each non-header row, add a text string in which the value in each column is left
      justified in a field of width one more than the maximum width for that column.
    - Exclude the last two columns, because they would be too wide.
    """
    columnNames = alnsetsHtmlEditor.get_column_names()[:-2]
    rowNames = alnsetsHtmlEditor.get_row_names()

    rows = [columnNames]
    for rowName in rowNames:
        rows.append(alnsetsHtmlEditor.get_values_from_row_name(rowName)[:-2])

    colWidths = []
    for colNum in range(len(columnNames)):
        headerLines = columnNames[colNum].split('<br>')
        headerWidth = max(len(line) for line in headerLines)
        valueWidth = max(len(row[colNum]) for row in rows[1:]) if len(rows) > 1 else 0
        colWidths.append(max(headerWidth, valueWidth))

    headerLine1 = []
    headerLine2 = []
    for colNum in range(len(columnNames)):
        headerParts = columnNames[colNum].split('<br>', 1)
        if len(headerParts) == 2:
            part1, part2 = headerParts
        else:
            part1, part2 = headerParts[0], ''
        headerLine1.append(part1.center(colWidths[colNum] + 1))
        headerLine2.append(part2.center(colWidths[colNum] + 1))

    lines = []
    lines.append(''.join(headerLine1))
    lines.append(''.join(headerLine2))
    for row in rows[1:]:
        lines.append(''.join(row[colNum].ljust(colWidths[colNum] + 1)
                             for colNum in range(len(columnNames))))

    return '\n'.join(lines) + '\n'

class AlnsetsHtmlEditor(object):
    # Class for manipulating the html text containing a CodAlignView Alignment Sets table.
    # Mostly written by ChatGPT 5.2
    def __init__(self, htmlText) :
        """
        Read html text in format produced by CodAlignView's "Alignment Sets" command.
        """
        self.html = htmlText
        self._parse()

    def _parse(self):
        beforeTable, rest = self.html.split('<table>\n', 1)
        tableText, afterTable = rest.split('</table>', 1)
        self._beforeTable = beforeTable
        self._afterTable = afterTable
        headerText, bodyText = tableText.split('</thead>', 1)
        headerRowText = re.search(r'<tr>\s*(.*?)\s*</tr>', headerText, re.S).group(1)
        self._columnNames = re.findall(r'<th>(.*?)</th>', headerRowText, re.S)
        rowTexts = re.findall(r'<tr>\s*(.*?)\s*</tr>', bodyText, re.S)
        self._rows = []
        for rowText in rowTexts :
            cells = re.findall(r'<t[hd](?: align="left")?>(.*?)</t[hd]>', rowText, re.S)
            self._rows.append(cells)

    def _serialize(self):
        lines = [self._beforeTable + '<table>',
                 '                <thead>',
                 '                    <tr>']
        for columnName in self._columnNames :
            lines.append('                        <th>%s</th>' % columnName)
        lines.extend(['                    </tr>',
                      '                </thead>'])

        for row in self._rows :
            lines.append('                <tr>')
            lines.append('                    <th align="left">%s</th>' % row[0])
            for value in row[1:] :
                lines.append('                    <td>%s</td>' % value)
            lines.append('                </tr>')

        lines.append('            </table>')
        self.html = '\n'.join(lines) + self._afterTable

    def get_html(self):
        return self.html

    def get_column_names(self):
        return self._columnNames[:]

    def get_row_names(self):
        return [row[0] for row in self._rows]

    def get_values_from_row_name(self, rowName):
        for row in self._rows :
            if row[0] == rowName :
                return row[:]
        raise ValueError('No row with name %s.' % rowName)  # pragma: no cover

    def add_column(self, columnName, beforeIndex, values):
        """
        Add a column to the table with specified columnName in the header, and the
            specified list of values.
        """
        assert len(values) == len(self._rows)
        self._columnNames.insert(beforeIndex, _html_escape(str(columnName)))
        for row, value in zip(self._rows, values) :
            row.insert(beforeIndex, _html_escape(str(value)))
        self._serialize()

    def remove_row(self, rowName):
        """
        Remove the row with specified name (name is the value in the first column)
        """
        for i, row in enumerate(self._rows) :
            if row[0] == rowName :
                del self._rows[i]
                self._serialize()
                return
        raise ValueError('No row with name %s.' % rowName)  # pragma: no cover
