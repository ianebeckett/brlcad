/*			R B _ D E L E T E . C
 *
 *	Written by:	Paul Tanenbaum
 *
 *  $Header$
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "machine.h"
#include "vmath.h"
#include "raytrace.h"
#include "redblack.h"
#define		RB_DELETE	1
#include "rb_internals.h"

static void rb_fixup (rb_tree *tree, struct rb_node *node, int order)
{
    int			direction;
    struct rb_node	*parent;
    struct rb_node	*w;

    RB_CKMAG(tree, RB_TREE_MAGIC, "red-black tree");
    RB_CKMAG(node, RB_NODE_MAGIC, "red-black node");
    RB_CKORDER(tree, order);

    while ((node != rb_root(tree, order))
        && (rb_get_color(node, order) == RB_BLACK))
    {
	parent = rb_parent(node, order);
	if (node == rb_left_child(parent, order))
	    direction = RB_LEFT;
	else
	    direction = RB_RIGHT;

	w = rb_other_child(parent, order, direction);
	if (rb_get_color(w, order) == RB_RED)
	{
	    rb_set_color(w, order, RB_BLACK);
	    rb_set_color(parent, order, RB_RED);
	    rb_rotate(parent, order, direction);
	    w = rb_other_child(parent, order, direction);
	}
	if ((rb_get_color(rb_child(w, order, direction), order) == RB_BLACK)
	 && (rb_get_color(rb_other_child(w, order, direction), order)
		 == RB_BLACK))
	{
	    rb_set_color(w, order, RB_RED);
	    node = parent;
	}
	else
	{
	    if (rb_get_color(rb_other_child(w, order, direction), order)
		    == RB_BLACK)
	    {
		rb_set_color(rb_child(w, order, direction), order, RB_BLACK);
		rb_set_color(w, order, RB_RED);
		rb_other_rotate(w, order, direction);
		w = rb_other_child(parent, order, direction);
	    }
	    rb_set_color(w, order, rb_get_color(parent, order));
	    rb_set_color(parent, order, RB_BLACK);
	    rb_set_color(rb_other_child(w, order, direction), order, RB_BLACK);
	    rb_rotate(parent, order, direction);
	    node = rb_root(tree, order);
	}
    }
    rb_set_color(node, order, RB_BLACK);
}

static void rb_free_node (struct rb_node *node)
{
    rb_tree	*tree;

    RB_CKMAG(node, RB_NODE_MAGIC, "red-black node");

    fprintf(stderr, "rb_free_node(%x)...\n", node);

    tree = node -> rbn_tree;
    if (rb_current(tree) == node)
	rb_current(tree) = rb_null(tree);

    rt_free((char *) node -> rbn_parent, "red-black parents");
    rt_free((char *) node -> rbn_left, "red-black left children");
    rt_free((char *) node -> rbn_right, "red-black right children");
    rt_free((char *) node -> rbn_color, "red-black colors");
    rt_free((char *) node -> rbn_package, "red-black packages");
    rt_free((char *) node, "red-black node");
}

static void rb_free_package (struct rb_package *package)
{
    RB_CKMAG(package, RB_PKG_MAGIC, "red-black package");

    fprintf(stderr, "rb_free_package(%x)...\n", package);

    rt_free((char *) package -> rbp_node, "red-black package nodes");
    rt_free((char *) package, "red-black package");
}

/*		        _ R B _ D E L E T E ( )
 *
 *	        Delete a node from a red-black tree
 *
 */
static void _rb_delete (rb_tree *tree, struct rb_node *node, int order)
{
    struct rb_node	*y;		/* The node to splice out */
    struct rb_node	*parent;
    struct rb_node	*only_child;

    RB_CKMAG(tree, RB_TREE_MAGIC, "red-black tree");
    RB_CKMAG(node, RB_NODE_MAGIC, "red-black node");
    RB_CKORDER(tree, order);

    fprintf(stderr, "_rb_delete(<%x>, <%x>, %d)...\n", tree, node, order);
    fflush(stderr);

    if ((rb_left_child(node, order) == rb_null(tree))
     || (rb_right_child(node, order) == rb_null(tree)))
	y = node;
    else
	y = _rb_neighbor(node, order, SENSE_MAX);
    
    fprintf(stderr, "_rb_delete(%x, %d) splicing out <%x>\n",
	node, order, y);

    if (rb_left_child(y, order) == rb_null(tree))
	only_child = rb_right_child(y, order);
    else
	only_child = rb_left_child(y, order);
    
    parent = rb_parent(only_child, order) = rb_parent(y, order);
    if (parent == rb_null(tree))
	rb_root(tree, order) = only_child;
    else if (y == rb_left_child(parent, order))
	rb_left_child(parent, order) = only_child;
    else
	rb_right_child(parent, order) = only_child;
    
    /*
     *	Splice y out if it's not node
     */
    if (y != node)
    {
	(node -> rbn_package)[order] = (y -> rbn_package)[order];
	((node -> rbn_package)[order] -> rbp_node)[order] = node;
    }
    fprintf(stderr,
	"_rb_delete() decrementing refs to <%x>... now equals %d\n",
	y, y -> rbn_pkg_refs - 1);
    if (--(y -> rbn_pkg_refs) == 0)
	rb_free_node(y);
    
    if (rb_get_color(y, order) == RB_BLACK)
	rb_fixup(tree, only_child, order);
}

/*		        R B _ D E L E T E ( )
 *
 *	        Applications interface to _rb_delete()
 *
 *	This function has one parameter: a tree and the data block
 *	to be deleted.  rb_delete() removes the data block from each
 *	order in the tree.  rb_delete() returns the number of orders
 *	for which the discovered node contains a data pointer that
 *	is not equal to the specified data block pointer.
 */
void rb_delete (rb_tree *tree, int order)
{
    int			nm_orders;
    struct rb_node	**node;		/* Nodes containing data */
    struct rb_package	*package;

    RB_CKMAG(tree, RB_TREE_MAGIC, "red-black tree");
    RB_CKORDER(tree, order);

    fprintf(stderr, "rb_delete(<%x>, %d)...\n", tree, order); fflush(stderr);

    nm_orders = tree -> rbt_nm_orders;
    package = (rb_current(tree) -> rbn_package)[order];

    node = (struct rb_node **)
	    rt_malloc(nm_orders * sizeof(struct rb_node *), "node list");
	
    for (order = 0; order < nm_orders; ++order)
	node[order] = (package -> rbp_node)[order];

    /*
     *	Do the deletion from each order
     */
    for (order = 0; order < nm_orders; ++order)
	_rb_delete(tree, node[order], order);

    rb_free_package(package);
    rt_free((char *) node, "node list");
}
