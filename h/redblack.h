/*			R E D B L A C K . H
 *
 *	The data structures, constants, and applications interface
 *	to LIBREDBLACK(3), the BRL-CAD red-black tree library.
 *
 *	Many of the routines in LIBREDBLACK(3) are based on the algorithms
 *	in chapter 13 of T. H. Cormen, C. E. Leiserson, and R. L. Rivest,
 *	_Introduction to algorithms_, MIT Press, Cambridge, MA, 1990.
 *
 *	Author:	Paul Tanenbaum
 *
 *  $Header$
 */

#ifndef REDBLACK_H
#define REDBLACK_H seen

/*
 *  Macros for providing function prototypes, regardless of whether
 *  the compiler understands them or not.
 *  It is vital that the argument list given for "args" be enclosed
 *  in parens.
 */
#if __STDC__ || USE_PROTOTYPES
#	define	RB_EXTERN(type_and_name,args)	extern type_and_name args
#	define	RB_ARGS(args)			args
#else
#	define	RB_EXTERN(type_and_name,args)	extern type_and_name()
#	define	RB_ARGS(args)			()
#endif


/*
 *			R B _ T R E E
 *
 *	This is the only data structure used in LIBREDBLACK
 *	to which application software need make any explicit
 *	reference.
 *
 *	The members of this structure are grouped into three
 *	classes:
 *	    Class I:	Reading is appropriate, when necessary,
 *			but applications should not modify.
 *	    Class II:	Reading and modifying are both appropriate,
 *			when necessary.
 *	    Class III:	All access should be through routines
 *			provided in LIBREDBLACK.  Touch these
 *			at your own risk!
 */
typedef struct
{
    /* CLASS I - Applications may read directly... */
    long	 	rbt_magic;	  /* Magic no. for integrity check */
    int			rbt_nm_nodes;	  /* Number of nodes */
    /* CLASS II - Applications may read/write directly... */
    void		(*rbt_print)();	  /* Data pretty-print function */
    int			rbt_debug;	  /* Debug bits */
    char		*rbt_description; /* Comment for diagnostics */
    /* CLASS III - Applications should not manipulate directly... */
    int		 	rbt_nm_orders;	  /* Number of simultaneous orders */
    int			(**rbt_order)();  /* Comparison functions */
    struct rb_node	**rbt_root;	  /* The actual trees */
    char		*rbt_unique;	  /* Uniqueness flags */
    struct rb_node	*rbt_current;	  /* Current node */
    struct rb_node	*rbt_empty_node;  /* Sentinel representing nil */
}	rb_tree;
#define	RB_TREE_NULL	((rb_tree *) 0)
#define	RB_TREE_MAGIC	0x72627472

/*
 *	Debug bit flags for member rbt_debug
 */
#define RB_DEBUG_INSERT	0x00000001	/* Insertion process */
#define RB_DEBUG_UNIQ	0x00000002	/* Uniqueness of inserts */
#define RB_DEBUG_ROTATE	0x00000004	/* Rotation process */

/*
 *			R B _ P A C K A G E
 *
 *		    Wrapper for application data
 *
 *	This structure provides a level of indirection between
 *	the application software's data and the red-black nodes
 *	in which the data is stored.  It is necessary because of
 *	the algorithm for deletion, which generally shuffles data
 *	among nodes in the tree.  The package structure allows
 *	the application data to remember which node "contains" it
 *	for each order.
 */
struct rb_package
{
    long		rbp_magic;	/* Magic no. for integrity check */
    struct rb_node	**rbp_node;	/* Containing nodes */
    void		*rbp_data;	/* Application data */
};
#define	RB_PKG_NULL	((struct rb_package *) 0)

/*
 *			    R B _ N O D E
 *
 *	For the most part, there is a one-to-one correspondence
 *	between nodes and chunks of application data.  When a
 *	node is created, all of its package pointers (one per
 *	order of the tree) point to the same chunk of data.
 *	However, subsequent deletions usually muddy this tidy
 *	state of affairs.
 */
struct rb_node
{
    long		rbn_magic;	/* Magic no. for integrity check */
    rb_tree		*rbn_tree;	/* Tree containing this node */
    struct rb_node	**rbn_parent;	/* Parents */
    struct rb_node	**rbn_left;	/* Left subtrees */
    struct rb_node	**rbn_right;	/* Right subtrees */
    char		*rbn_color;	/* Colors of this node */
    struct rb_package	**rbn_package;	/* Contents of this node */
    int			rbn_pkg_refs;	/* How many orders are being used? */
};
#define	RB_NODE_NULL	((struct rb_node *) 0)

/*
 *	Applications interface to rb_extreme()
 */
#define	SENSE_MIN	0
#define	SENSE_MAX	1
#define	rb_min(t,o)	rb_extreme((t), (o), SENSE_MIN)
#define	rb_max(t,o)	rb_extreme((t), (o), SENSE_MAX)
#define rb_pred(t,o)	rb_neighbor((t), (o), SENSE_MIN)
#define rb_succ(t,o)	rb_neighbor((t), (o), SENSE_MAX)

/*
 *	Applications interface to rb_walk()
 */
#define	PREORDER	0
#define	INORDER		1
#define	POSTORDER	2

/*
 *	Applications interface to LIBREDBLACK
 */
RB_EXTERN(rb_tree *rb_create,	(char		*description,
				 int 		nm_orders,
				 int		(**order_funcs)()
				));
RB_EXTERN(rb_tree *rb_create1,	(char		*description,
				 int		(*order_func)()
				));
RB_EXTERN(void *rb_curr,	(rb_tree	*tree,
				 int		order
				));
#define		rb_curr1(t)	rb_curr((t), 0)
RB_EXTERN(void rb_delete,	(rb_tree	*tree,
				 int		order
				));
#define		rb_delete1(t)	rb_delete((t), 0)
RB_EXTERN(void rb_diagnose_tree,(rb_tree	*tree,
				 int		order,
				 int		trav_type
				));
RB_EXTERN(void *rb_extreme,	(rb_tree	*tree,
				 int		order,
				 int		sense
				));
RB_EXTERN(int rb_insert,	(rb_tree	*tree,
				 void		*data
				));
RB_EXTERN(int rb_is_uniq,	(rb_tree	*tree,
				 int		order
				));
#define		rb_is_uniq1(t)	rb_is_uniq((t), 0)
RB_EXTERN(void *rb_neighbor,	(rb_tree	*tree,
				 int		order,
				 int		sense
				));
RB_EXTERN(void *rb_search,	(rb_tree	*tree,
				 int		order,
				 void	*data
				));
#define		rb_search1(t,d)	rb_search((t), 0, (d))
RB_EXTERN(void rb_set_uniqv,	(rb_tree	*tree,
				 bitv_t		vec
				));
RB_EXTERN(void rb_summarize_tree,(rb_tree	*tree
				 ));
RB_EXTERN(void rb_uniq_all_off,	(rb_tree	*tree
				));
RB_EXTERN(void rb_uniq_all_on,	(rb_tree	*tree
				));
RB_EXTERN(int rb_uniq_off,	(rb_tree	*tree,
				 int		order
				));
#define		rb_uniq_off1(t)	rb_uniq_off((t), 0)
RB_EXTERN(int rb_uniq_on,	(rb_tree	*tree,
				 int		order
				));
#define		rb_uniq_on1(t)	rb_uniq_on((t), 0)
RB_EXTERN(void rb_walk,		(rb_tree	*tree,
				 int		order,
				 void		(*visit)(),
				 int		trav_type
				));
#define		rb_walk1(t,v,d)	rb_walk((t), 0, (v), (d))

#endif /* REDBLACK_H */
