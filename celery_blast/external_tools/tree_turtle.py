#script copied from https://github.com/msrosenberg/Phy2HTML

"""
Simple module to draw a phylogeny using Turtle Graphics

mostly used as an extra test of the node class
"""

import math
import turtle


def tree_turtle(tree, minx: int, maxx: int, miny: int, maxy: int, scale: float, draw_labels: bool,
                draw_branch_lengths: bool) -> int:
    """
    use turtle graphics to draw the tree

    minx, maxx, miny, and maxy represent the drawing bounds for the target node, with minx representing the
    horizontal position of the ancestor of the node, maxx representing the right hand edge of the entire tree,
    miny and maxy representing the vertical positioning of the node and all of its descendants.

    scale is a pre-calculated value that converts branch lengths to pixels

    branch length drawing is not currently enabled
    """

    yadj = 5  # this is a shift for font height
    name_padding = 10  # space names away from tree tips

    """
    the following calcuates the number of pixels necessary for the horizontal line connecting the node to it's 
    ancestor
    """
    x = math.trunc(tree.branch_length * scale)

    """
    if the node has descendants, first draw all of the descendants in the box which goes from the right edge 
    of the horizontal line for this node to the right hand edge of the drawing window, and the vertical box 
    defined for the entire node
    """
    if tree.n_descendants() > 0:  # this is an internal node
        top_vert_line = 0
        bottom_vert_line = 0
        nd = tree.n_tips()
        topy = miny
        for i, d in enumerate(tree.descendants):
            """
            divide the vertical plotting area for each descendant proportional
            to the number of tips contained within that descendant
            """
            ndd = d.n_tips()
            bottomy = topy + math.trunc((ndd / nd) * (maxy - miny))
            # draw the descendant in its own smaller bounded box
            y = tree_turtle(d, minx + x, maxx, topy, bottomy, scale, draw_labels, draw_branch_lengths)
            """
            the vertical position of the first and last descendants represent
            the positions to draw the vertical line connecting all of the
            descendants
            """
            if i == 0:
                bottom_vert_line = y
            elif i == tree.n_descendants() - 1:
                top_vert_line = y
            topy = bottomy
        """
        draw the vertical line connecting the descendants at the horizontal
        position of the node
        """
        turtle.penup()
        turtle.goto(minx + x, bottom_vert_line)
        turtle.pendown()
        turtle.goto(minx + x, top_vert_line)
        """
        the vertical position of the node should be the midpoint of the
        vertical line connecting the descendants
        """
        y = ((top_vert_line - bottom_vert_line) // 2) + bottom_vert_line
    else:  # this is a tip node
        """
        if the node has no descendants, figure out the vertical position as the midpoint of the vertical bounds
        """
        y = ((maxy - miny) // 2) + miny
        if draw_labels:
            # if desired, label the node
            turtle.penup()
            turtle.goto(minx + x + name_padding, y - yadj)
            turtle.pendown()
            turtle.write(tree.name)

    # draw the horizontal line connecting the node to its ancestor
    turtle.penup()
    turtle.goto(minx, y)  # ancestral node location
    turtle.pendown()
    turtle.goto(minx + x, y)  # this node location

    # add branch lengths
    if draw_branch_lengths:
        pass
        # not yet enabled

    return y


def draw_tree_turtle(root):
    xwidth = 1000
    yheight = 1000
    margin = 10
    maxbranch = root.max_node_tip_length()
    scale = (xwidth - 100) / maxbranch
    win = turtle.Screen()
    win.screensize(xwidth + 2 * margin, yheight + 2 * margin)
    win.setworldcoordinates(-margin, -margin, xwidth + margin, yheight + margin)
    turtle.hideturtle()
    turtle.color("black")
    tree_turtle(root, 0, xwidth, 0, yheight, scale, True, False)
