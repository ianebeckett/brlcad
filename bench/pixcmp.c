/*
 *			P I X C M P . C
 *
 *  Compute the difference between two .pix files.
 *  To establish context, a monochrome image is produced when there
 *  are no differences;  otherwise the channels that differ are
 *  highlighted on differing pixels.
 *
 *  This routine operates on a pixel-by-pixel basis, and thus
 *  is independent of the resolution of the image.
 *
 *  Author -
 *	Charles M. Kennedy
 *	Michael John Muuss
 *
 *  Source -
 *	SECAD/VLD Computing Consortium, Bldg 394
 *	The U. S. Army Ballistic Research Laboratory
 *	Aberdeen Proving Ground, Maryland  21005
 *  
 *  Copyright Notice -
 *	This software is Copyright (C) 1988 by the United States Army.
 *	All rights reserved.
 */
#ifndef lint
static char RCSid[] = "@(#)$Header$ (BRL)";
#endif

#include <stdio.h>

long	matching;
long	off1;
long	offmany;

char usage[] = "Usage: pixcmp f1.pix f2.pix\n";

main(argc, argv)
char **argv;
{
	register FILE *f1, *f2;

	if( argc != 3 )  {
		fprintf(stderr, "%s", usage);
		exit(0);
	}

	if( strcmp( argv[1], "-" ) == 0 )
		f1 = stdin;
	else if( (f1 = fopen( argv[1], "r" ) ) == NULL )  {
		perror( argv[1] );
		exit(1);
	}
	if( strcmp( argv[2], "-" ) == 0 )
		f2 = stdin;
	else if( (f2 = fopen( argv[2], "r" ) ) == NULL )  {
		perror( argv[2] );
		exit(1);
	}
	while(1)  {
		register int r1, g1, b1;
		int r2, g2, b2;

		r1 = fgetc( f1 );
		g1 = fgetc( f1 );
		b1 = fgetc( f1 );
		r2 = fgetc( f2 );
		g2 = fgetc( f2 );
		b2 = fgetc( f2 );
		if( feof(f1) || feof(f2) )  break;

		if( r1 != r2 || g1 != g2 || b1 != b2 )  {
			register int i;

			/* Highlight differing channels */
			if( r1 != r2 )  {
				if( (i = r1 - r2) < 0 )  i = -i;
				if( i > 1 )  {
					offmany++;
				} else {
					off1++;
				}
			} else {
				matching++;
			}
			if( g1 != g2 )  {
				if( (i = g1 - g2) < 0 )  i = -i;
				if( i > 1 )  {
					offmany++;
				} else {
					off1++;
				}
			} else {
				matching++;
			}
			if( b1 != b2 )  {
				if( (i = b1 - b2) < 0 )  i = -i;
				if( i > 1 )  {
					offmany++;
				} else {
					off1++;
				}
			} else {
				matching++;
			}
		}  else  {
			/* Common case:  equal.  Give B&W NTSC average */
			matching += 3;
		}
	}
	fprintf(stderr,
		"pixcmp bytes: %7d matching, %7d off by 1, %7d off by many\n",
		matching, off1, offmany );
	if( offmany) {
		exit(1);
	}
	exit(0);
}
