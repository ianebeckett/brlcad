/*
 *			F B P O I N T . C
 *
 *  Tool to identify X,Y coordinates on the screen.
 *
 *  Authors -
 *	Michael Johns Muuss
 *	Bob Suckling
 *
 *  Source -
 *	SECAD/VLD Computing Consortium, Bldg 394
 *	The U. S. Army Ballistic Research Laboratory
 *	Aberdeen Proving Ground, Maryland  21005-5066
 *
 *  Copyright Notice -
 *	This software is Copyright (C) 1986 by the United States Army.
 *	All rights reserved.
 */
#ifndef lint
static char RCSid[] = "@(#)$Header$ (BRL)";
#endif

#include <stdio.h>
#include "../h/fb.h"

FBIO *fbp;

int	JumpSpeed;		/* # pixels skiped with fast commands. */

RGBpixel curPix; 		/* Current pixel value */
int	curX, curY;		/* current position */
int	oldX, oldY;		/* previous position */

int	Run = 1;		/* Tells when to stop the main loop */

int	xflag, yflag;
char	*xprefix = NULL;
char	*yprefix = NULL;
char	null_str = '\0';

char usage[] = "\
Usage: fbpoint [-h] [-x[prefix]] [-y[prefix]] [-X init_x] [-Y init_y]\n";

static char *help = "\
Char:   Command:                                                \r\n\
h B	Left (1)\r\n\
j N	Down (1)\r\n\
k P	Up (1)\r\n\
l F	Right (1)\r\n\
H ^B	Left (many)\r\n\
J ^N	Down (many)\r\n\
K ^P	Up (many)\r\n\
L ^F	Right (many)\r\n\
q Q	QUIT\r\n\
";

SimpleInput()	/* ==== get keyboard input.	*/
{
	char ch;

	if( read( 0, &ch, 1) <= 0 ) {
		Run = 0;
		return;
	}
	switch( ch ) {
	default:
		fprintf( stderr,
		"Unknown command(%c:%o). Type '?' for help!           \r\n",
			ch, ch );
		break;

	case '?':
		fprintf( stderr, "%s", help );
		return;

	case 'q':
	case 'Q':
		Run = 0;
		return;
#define ctl(x)	('x'&037)

	case 'B':
	case 'h':
		--curX;		/* Go left.	*/
		return;
	case 'N':
	case 'j':
		--curY;		/* Go down.	*/
		return;
	case 'P':
	case 'k':
		++curY;		/* Go up.	*/
		return;
	case 'F':
	case 'l':
		++curX;		/* Go right.	*/
		return;
	case ctl(b):
	case 'H':
		curX -= JumpSpeed;	/* Go LEFT.	*/
		return;
	case ctl(n):
	case 'J':
		curY -= JumpSpeed;	/* Go DOWN.	*/
		return;
	case ctl(p):
	case 'K':
		curY += JumpSpeed;	/* Go UP.	*/
		return;
	case ctl(f):
	case 'L':
		curX += JumpSpeed;	/* Go RIGHT.	*/
		return;
	}
}

main(argc, argv)
char **argv;
{
	int width, height;
	char *malloc();

	setbuf( stderr, malloc( BUFSIZ ) );
	width = height = 0;
	curX = curY = -1;

	while( argc > 1 ) {
		if( strcmp( argv[1], "-h" ) == 0 ) {
			width = height = 1024;
		} else if( strncmp( argv[1], "-x", 2 ) == 0 ) {
			if( xflag++ != 0 )
				break;
			xprefix = &argv[1][2];
		} else if( strncmp( argv[1], "-y", 2 ) == 0 ) {
			if( yflag++ != 0 )
				break;
			yprefix = &argv[1][2];
		} else if( strncmp( argv[1], "-X", 2 ) == 0 ) {
			if( strlen(argv[1]) > 2 )
				curX = atoi( &argv[1][2] );
			else if( argc > 1 ) {
				curX = atoi( argv[2] );
				argc--; argv++;
			} else
				;	/* no value given */
		} else if( strncmp( argv[1], "-Y", 2 ) == 0 ) {
			if( strlen(argv[1]) > 2 )
				curY = atoi( &argv[1][2] );
			else if( argc > 1 ) {
				curY = atoi( argv[2] );
				argc--; argv++;
			} else
				;	/* no value given */
		} else
			break;
		argc--;
		argv++;
	}
	if( argc > 1 ) {
		fprintf( stderr, usage );
		exit( 1 );
	}

	/* fix up pointers for printf */
	if( xprefix == NULL )
		xprefix = &null_str;
	if( yprefix == NULL )
		yprefix = &null_str;

	if( (fbp = fb_open( NULL, width, height )) == NULL )
		exit(12);

	JumpSpeed = fb_getwidth(fbp)/16;
	if( JumpSpeed < 2 )  JumpSpeed = 2;
	/* check for default starting positions */
	if( curX < 0 )
		curX = fb_getwidth(fbp)/2;
	if( curY < 0 )
		curY = fb_getheight(fbp)/2;
	oldX = oldY = -1;

	/* Set RAW mode */
	save_Tty( 0 );
	set_Raw( 0 );
	clr_Echo( 0 );

	while( Run )  {
		if( curX < 0 )
			curX = 0;
		if( curX >= fb_getwidth(fbp) )
			curX = fb_getwidth(fbp) -1;
		if( curY < 0 )
			curY = 0;
		if( curY >= fb_getheight(fbp) )
			curY = fb_getheight(fbp) -1;

		if( oldX != curX || oldY != curY ) {
			/* get pixel value, move cursor */
			fb_read( fbp, curX, curY, curPix, 1 );
			fb_cursor( fbp, 1, curX, curY );
			oldX = curX;
			oldY = curY;
		}
		fprintf( stderr, "xy=(%4d,%4d)  [%3d,%3d,%3d]      \r",
			curX, curY, curPix[RED], curPix[GRN], curPix[BLU] );
		fflush( stderr );

		SimpleInput();			/* read and do keybord	*/
	}

	fb_cursor( fbp, 0, curX, curY );	/* turn off */

	fprintf( stderr, "\n" );
	fflush( stderr );

	reset_Tty( 0 );

	/* write final location on stdout */
	if( xflag != 0 && yflag == 0 )
		printf( "%s%d\n", xprefix, curX );
	else if( yflag != 0 && xflag == 0 )
		printf( "%s%d\n", yprefix, curY );
	else
		printf( "%s%d %s%d\n", xprefix, curX, yprefix, curY );

	fb_close( fbp );
	exit( 0 );
}
