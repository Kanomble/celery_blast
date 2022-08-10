#script copied from https://github.com/msrosenberg/Phy2HTML
"""
Phy2HTML
"""

from math import trunc
from typing import Tuple
from .tree_utils import Node, read_newick_tree

class Branch:
    def __init__(self, min_col: int = 0, col_span: int = 0, row: int = 0, label: str = "&nbsp;"):
        self.min_col = min_col
        self.col_span = col_span
        self.row = row
        self.label = label


class VLine:
    def __init__(self, min_row: int = 0, row_span: int = 0, col: int = 0, label: str = "&nbsp;"):
        self.min_row = min_row
        self.row_span = row_span
        self.col = col
        self.label = label


class Taxon:
    def __init__(self, node: Node, row: int = 0, col: int = 0):
        self.node = node
        self.row = row
        self.col = col


def start_html(outlist: list) -> None:
    outlist.append("<html>\n")
    outlist.append("  <head>\n")


def end_head_section(outlist: list) -> None:
    outlist.append("  </head>\n")
    outlist.append("  <body>\n")


def end_html(outlist: list) -> None:
    outlist.append("  </body>\n")
    outlist.append("</html>\n")


def write_style_to_head(outlist: list, nrows: int, ncols: int, taxa: list, branches: list, vlines: list,
                        col_width: str, row_height: str, name_width: str, prefix: str, scale_branches: bool) -> None:
    outlist.append("    <style>\n")
    outlist.append("      #{}phylogeny {{\n".format(prefix))
    outlist.append("                   display: grid;\n")
    outlist.append("                   grid-template-rows: repeat({}, {});\n".format(nrows, row_height))
    outlist.append("                   grid-template-columns: repeat({}, {}) {};\n".format(ncols-1, col_width,
                                                                                           name_width))
    outlist.append("                 }\n")
    outlist.append("      .{}taxon-name {{ align-self: center; padding-left: 10px }}\n".format(prefix))
    outlist.append("      .{}genus-species-name {{ font-style: italic }}\n".format(prefix))
    outlist.append("      .{}branch-line {{ border-bottom: solid black 1px; text-align: center }}\n".format(prefix))
    outlist.append("      .{}vert-line {{ border-right: solid black 1px; text-align: right }}\n".format(prefix))
    outlist.append("\n")
    for i, t in enumerate(taxa):
        if scale_branches:
            tcol = t.col + 1
            if tcol > ncols:
                tcol = ncols
            cspan = ncols - tcol + 1
        else:
            tcol = ncols
            cspan = 1
        outlist.append("      #{}taxon{} {{ grid-area: {} / {} / span 2 / span {} }}\n".format(prefix, i+1, t.row,
                                                                                               tcol, cspan))
    outlist.append("\n")
    for i, b in enumerate(branches):
        outlist.append("      #{}branch{} {{ grid-area: {} / {} / span 1 / span {} }}\n".format(prefix, i+1, b.row,
                                                                                                b.min_col, b.col_span))
    outlist.append("\n")
    for i, v in enumerate(vlines):
        outlist.append("      #{}vline{} {{ grid-area: {} / {} / span {} / span 1 }}\n".format(prefix,  i+1, v.min_row,
                                                                                               v.col, v.row_span))
    outlist.append("\n")
    outlist.append("    </style>\n")


def write_tree_to_body(outlist: list, taxa: list, branches: list, vlines: list, prefix: str) -> None:
    outlist.append("    <div id=\"{}unique_phylogeny_container\" class=\"phylogeny_container\">\n".format(prefix))
    outlist.append("      <div id=\"{}phylogeny\" class=\"phylogeny_grid\">\n".format(prefix))
    outlist.append("\n")
    for i, t in enumerate(taxa):
        outlist.append("        <div id=\"{0}taxon{1}\" "
                       "class=\"{0}genus-species-name {0}taxon-name\">{2}</div>\n".format(prefix, i+1, t.node.name))
    outlist.append("\n")
    for b, branch in enumerate(branches):
        outlist.append("        <div id=\"{0}branch{1}\" class=\"{0}branch-line\">{2}</div>\n".format(prefix, b+1,
                                                                                                      branch.label))
    outlist.append("\n")
    for v, vline in enumerate(vlines):
        outlist.append("        <div id=\"{0}vline{1}\" class=\"{0}vert-line\">{2}</div>\n".format(prefix, v+1,
                                                                                                   vline.label))
    outlist.append("\n")
    outlist.append("      </div>\n")
    outlist.append("    </div>\n")


def total_rows_per_node(n: int, rows_per_tip: int) -> int:
    return (2*n - 1) * rows_per_tip


def tree_recursion(tree, min_col: int, max_col: int, min_row: int, max_row: int, taxa: list, branches: list,
                   vlines: list, rows_per_tip: int, label_branches: bool, scale_branches: bool, scale: float) -> int:
    """
    calculate positions of taxa, branches, and vertical connectors on subtrees
    """

    """
    determine the number of columns for the branch connecting a node to its ancestor
    """
    if tree.ancestor is not None:
        if scale_branches:
            col_span = trunc(tree.branch_length * scale)
        else:
            col_span = tree.node_depth - tree.ancestor.node_depth
    else:
        col_span = 0

    """
    if the node has descendants, first draw all of the descendants in the box which starts in the column to the 
    right of this node, and the rows defined for the entire node
    """
    if tree.n_descendants() > 0:  # this is an internal node
        horizontal_connections = []
        vert_top_row = 0
        vert_bottom_row = 0
        top_row = min_row
        for i, d in enumerate(tree.descendants):
            """
            calculate the total rows for each descendant based on the number of tips of the descendant
            """
            ndd = d.n_tips()
            d_rows = total_rows_per_node(ndd, rows_per_tip)
            bottom_row = top_row + d_rows - 1
            # draw the descendant in its own smaller bounded box
            row = tree_recursion(d, min_col + col_span, max_col, top_row, bottom_row, taxa, branches, vlines,
                                 rows_per_tip, label_branches, scale_branches, scale)
            """
            the rows of the first and last descendants represent the positions to draw the vertical line 
            connecting all of the descendants
            """
            if i == 0:
                vert_top_row = row + 1
            elif i == tree.n_descendants() - 1:
                vert_bottom_row = row
            top_row = bottom_row + rows_per_tip + 1
            horizontal_connections.append(row)

        """
        the vertical position of the node should be the midpoint of the vertical line connecting the descendants
        """
        row = ((vert_bottom_row - vert_top_row) // 2) + vert_top_row
        horizontal_connections.append(row)

        """
        add the vertical line connecting the descendants at the horizontal position of the node
        """
        horizontal_connections.sort()
        for i in range(1, len(horizontal_connections)):
            if horizontal_connections[i] != horizontal_connections[i-1]:  # skip for lines of zero height
                new_line = VLine(horizontal_connections[i-1]+1, horizontal_connections[i]-horizontal_connections[i-1],
                                 min_col+col_span)
                vlines.append(new_line)
                if label_branches:
                    new_line.label = "vline" + str(len(vlines))

    else:  # this is a tip node
        """
        if the node has no descendants, add it to the taxon list
        """
        row = min_row
        new_taxon = Taxon(tree, row, min_col + col_span)
        taxa.append(new_taxon)

    # add the branch connecting the node to its ancestor
    if col_span > 0:
        new_branch = Branch(min_col+1, col_span, row)
        branches.append(new_branch)
        if label_branches:
            new_branch.label = "branch" + str(len(branches))

    return row


def calculate_tree(tree: Node, nrows: int, ncols: int, rows_per_tip: int,
                   label_branches: bool, scale_branches: bool = False) -> Tuple[list, list, list]:
    taxa = []
    branches = []
    vlines = []
    if scale_branches:
        tree_depth = tree.max_node_tip_length()
        scale = (ncols - 1) / tree_depth
    else:
        scale = 1
    tree_recursion(tree, 1, ncols, 1, nrows, taxa, branches, vlines, rows_per_tip, label_branches, scale_branches,
                   scale)
    return taxa, branches, vlines


def add_node_depth(tree: Node, max_depth: int) -> None:
    """
    add the column depth of each node on the tree, where the root is column 1 and the tips are column x - 1
    where x is the last column which will contain the tip names
    """
    tree.node_depth = max_depth - tree.max_node_tip_count()
    for d in tree.descendants:
        add_node_depth(d, max_depth)


def create_html_tree(inname: str, outname: str, col_width: str = "40px", row_height: str = "10px",
                     name_width: str = "200px", prefix: str = "", label_branches: bool = False,
                     scale_branches: bool = False, tree_cols: int = 1, rows_per_tip: int = 2,
                     ) -> list:
    with open(inname, "r") as infile:
        newick_str = infile.read()
    newick_str = newick_str[:newick_str.find(";")+1]
    newick_str.replace("\n", "")

    tree = read_newick_tree(newick_str)

    ntips = tree.n_tips()
    nrows = total_rows_per_node(ntips, rows_per_tip)
    if scale_branches:
        ncols = tree_cols + 1  # ools for tree plus one for tip labels
    else:
        ncols = tree.max_node_tip_count() + 1
        add_node_depth(tree, ncols+1)
    taxa, branches, vlines = calculate_tree(tree, nrows, ncols, rows_per_tip, label_branches, scale_branches)
    outlist = []
    start_html(outlist)
    write_style_to_head(outlist, nrows, ncols, taxa, branches, vlines, col_width, row_height, name_width, prefix,
                        scale_branches)
    end_head_section(outlist)
    write_tree_to_body(outlist, taxa, branches, vlines, prefix)
    end_html(outlist)
    if outname != "":  # if output file name is provided, write to file
        with open(outname, "w") as outfile:
            outfile.writelines(outlist)

    return outlist


def query_user(prompt: str, default: str) -> str:
    x = input("{} [default={}]: ".format(prompt, default))
    if x == "":
        x = default
    return x
