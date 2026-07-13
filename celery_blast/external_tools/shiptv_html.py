import argparse
import re
from pathlib import Path


TITLE = "shiptv: Standalone HTML Interactive Phylogenetic Tree Visualization"
BRANCH_LENGTH_PATTERN = re.compile(
    r":([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
)
ZERO_BRANCH_LENGTH_BEFORE_TREE_TOKEN_PATTERN = re.compile(
    r":([+-]?(?:0+(?:\.0*)?|\.0+)(?:[eE][+-]?\d+)?)(?=\s*[,);])"
)
CHROMA_CDN_PATTERN = re.compile(
    r"\s*<script[^>]+https://cdnjs\.cloudflare\.com/ajax/libs/chroma-js/[^>]+></script>\s*",
    re.IGNORECASE,
)


def make_tree_drawable(newick):
    branch_lengths = [
        float(match.group(1))
        for match in BRANCH_LENGTH_PATTERN.finditer(newick)
    ]
    all_zero = branch_lengths and all(length == 0 for length in branch_lengths)

    if not all_zero:
        return newick

    return ZERO_BRANCH_LENGTH_BEFORE_TREE_TOKEN_PATTERN.sub(":1", newick)


def patch_shiptv_html_file(input_path, output_path=None):
    input_path = Path(input_path)
    output_path = Path(output_path) if output_path is not None else input_path
    output_path.write_text(patch_shiptv_html(input_path.read_text(encoding="utf-8")), encoding="utf-8")


def patch_shiptv_html(html):
    if "shiptv: Standalone HTML Interactive Phylogenetic Tree Visualization" not in html:
        return html

    html = CHROMA_CDN_PATTERN.sub("\n", html)
    html = _patch_document_header(html)
    html = html.replace('style="overflow: none;"', 'style="overflow: hidden;"')
    html = html.replace(
        '<div class="ag-theme-balham" id="metadata-grid" style="height: 600px; width: 100%;"></div>',
        '<div id="metadata-grid" style="height: 600px; width: 100%;"></div>',
    )
    html = _patch_ag_grid_code(html)
    html = _patch_tree_loading_code(html)
    _assert_no_obsolete_shiptv_code(html)
    return html


def _patch_document_header(html):
    malformed_header = (
        '<!DOCTYPE html>\n'
        '<meta charset="utf-8"\n'
        f'      title="{TITLE}">\n'
        '<html>\n'
        '<head>'
    )
    valid_header = (
        '<!DOCTYPE html>\n'
        '<html>\n'
        '<head>\n'
        '  <meta charset="utf-8">\n'
        f'  <title>{TITLE}</title>'
    )
    if malformed_header in html:
        return html.replace(malformed_header, valid_header, 1)

    duplicate_head = f'  <title>{TITLE}</title>\n<head>'
    if duplicate_head in html:
        return html.replace(duplicate_head, f'  <title>{TITLE}</title>', 1)

    if f'<title>{TITLE}</title>' in html and '<meta charset="utf-8">' in html:
        return html

    return html


def _patch_ag_grid_code(html):
    html = html.replace(
        "  var grid;\n",
        "  var gridApi;\n  var metadataColumns = [];\n",
        1,
    )
    html = html.replace(
        "      return cols.reduce(function (acc, col) {",
        "      return metadataColumns.reduce(function (acc, col) {",
    )
    html = html.replace(
        "    const cols = Object.keys(Object.values(genomeMetadata)[0]);\n"
        "    cols.unshift('genome');\n"
        "    const columnDefs = cols.map(function (x) {",
        "    metadataColumns = Object.keys(Object.values(genomeMetadata)[0]);\n"
        "    metadataColumns.unshift('genome');\n"
        "    const columnDefs = metadataColumns.map(function (x) {",
    )
    html = html.replace(
        "      return cols.reduce(function (acc, col) {",
        "      return metadataColumns.reduce(function (acc, col) {",
    )
    html = html.replace(
        "      rowSelection: 'multiple'\n"
        "    };\n"
        "    grid = new agGrid.Grid(document.getElementById('metadata-grid'), gridOptions);",
        "      rowSelection: {\n"
        "        mode: 'multiRow',\n"
        "        checkboxes: false,\n"
        "        headerCheckbox: false,\n"
        "        enableClickSelection: true\n"
        "      },\n"
        "      theme: agGrid.themeBalham\n"
        "    };\n"
        "    gridApi = agGrid.createGrid(document.getElementById('metadata-grid'), gridOptions);",
    )
    html = html.replace(
        "    if (!_.isUndefined(grid)) {\n"
        "      grid.context.beans.gridApi.beanInstance.forEachNode(function (x) {\n"
        "        if (sel_leaves.hasOwnProperty(x.data.genome)) {\n"
        "          console.log(x);\n"
        "          x.selected = true;\n"
        "        } else {\n"
        "          x.selected = false;\n"
        "        }\n"
        "      });\n"
        "      grid.context.beans.gridApi.beanInstance.redrawRows();\n"
        "    }",
        "    if (gridApi) {\n"
        "      gridApi.forEachNode(function (node) {\n"
        "        node.setSelected(\n"
        "          Object.prototype.hasOwnProperty.call(sel_leaves, node.data.genome)\n"
        "        );\n"
        "      });\n"
        "      gridApi.redrawRows();\n"
        "    }",
    )
    html = html.replace(
        "    grid.context.beans.gridApi.beanInstance.setRowData(newRowData)",
        "    if (gridApi) {\n"
        "      gridApi.setGridOption('rowData', newRowData);\n"
        "    }",
    )
    return html


def _patch_tree_loading_code(html):
    if "function makeTreeDrawable(newick)" in html:
        return html.replace(
            "  tree.load(newick_string);",
            "  tree.load(makeTreeDrawable(newick_string));",
            1,
        )

    make_tree_drawable_js = (
        "  function makeTreeDrawable(newick) {\n"
        "    var branchLengths = Array.from(\n"
        "      newick.matchAll(\n"
        "        /:([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)/g\n"
        "      )\n"
        "    ).map(function (match) {\n"
        "      return Number(match[1]);\n"
        "    });\n"
        "\n"
        "    var allZero =\n"
        "      branchLengths.length > 0 &&\n"
        "      branchLengths.every(function (length) {\n"
        "        return Number.isFinite(length) && length === 0;\n"
        "      });\n"
        "\n"
        "    if (!allZero) {\n"
        "      return newick;\n"
        "    }\n"
        "\n"
        "    return newick.replace(\n"
        "      /:([+-]?(?:0+(?:\\.0*)?|\\.0+)(?:[eE][+-]?\\d+)?)(?=\\s*[,);])/g,\n"
        "      ':1'\n"
        "    );\n"
        "  }\n"
        "\n"
        "  tree.load(makeTreeDrawable(newick_string));"
    )
    return html.replace("  tree.load(newick_string);", make_tree_drawable_js, 1)


def _assert_no_obsolete_shiptv_code(html):
    app_code = _extract_shiptv_app_code(html)
    obsolete_fragments = [
        "new agGrid.Grid",
        ".context.beans",
        ".beanInstance",
        "rowSelection: 'multiple'",
    ]
    for fragment in obsolete_fragments:
        if fragment in app_code:
            raise ValueError("shiptv HTML still contains obsolete code: {}".format(fragment))

    obsolete_document_fragments = [
        "overflow: none",
        "https://cdnjs.cloudflare.com/ajax/libs/chroma-js",
    ]
    for fragment in obsolete_document_fragments:
        if fragment in html:
            raise ValueError("shiptv HTML still contains obsolete code: {}".format(fragment))


def _extract_shiptv_app_code(html):
    start_marker = "  var DEFAULT_GRID_COL_OPTS = {"
    if start_marker not in html:
        return html
    start = html.index(start_marker)
    end = html.rfind("</script>")
    return html[start:end] if end > start else html[start:]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Patch generated shiptv standalone HTML.")
    parser.add_argument("input_html")
    parser.add_argument("output_html", nargs="?")
    args = parser.parse_args(argv)

    patch_shiptv_html_file(args.input_html, args.output_html)


if __name__ == "__main__":
    main()
