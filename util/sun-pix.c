/*
 *			S U N - P I X . C
 *
 *  Program to take Sun bitmap files created with Sun's ``screendump''
 *  command, and convert them to pix(5) format files.
 *
 *  Authors -
 *	Phillip Dykstra
 *	Michael John Muuss
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

extern int	getopt();
extern char	*optarg;
extern int	optind;

/*
 * Description of Sun header for files containing raster images
 */
struct rasterfile {
	int	ras_magic;		/* magic number */
	int	ras_width;		/* width (pixels) of image */
	int	ras_height;		/* height (pixels) of image */
	int	ras_depth;		/* depth (1, 8, or 24 bits) of pixel */
	int	ras_length;		/* length (bytes) of image */
	int	ras_type;		/* type of file; see RT_* below */
	int	ras_maptype;		/* type of colormap; see RMT_* below */
	int	ras_maplength;		/* length (bytes) of following map */
	/* color map follows for ras_maplength bytes, followed by image */
} header;

char	inbuf[sizeof(struct rasterfile)];

#define	RAS_MAGIC	0x59a66a95

	/* Sun supported ras_type's */
#define RT_OLD		0	/* Raw pixrect image in 68000 byte order */
#define RT_STANDARD	1	/* Raw pixrect image in 68000 byte order */
#define RT_BYTE_ENCODED	2	/* Run-length compression of bytes */
#define RT_EXPERIMENTAL 0xffff	/* Reserved for testing */

	/* Sun registered ras_maptype's */
#define RMT_RAW		2
	/* Sun supported ras_maptype's */
#define RMT_NONE	0	/* ras_maplength is expected to be 0 */
#define RMT_EQUAL_RGB	1	/* red[ras_maplength/3],green[],blue[] */

/*
 * NOTES:
 * 	Each line of the image is rounded out to a multiple of 16 bits.
 *   This corresponds to the rounding convention used by the memory pixrect
 *   package (/usr/include/pixrect/memvar.h) of the SunWindows system.
 *	The ras_encoding field (always set to 0 by Sun's supported software)
 *   was renamed to ras_length in release 2.0.  As a result, rasterfiles
 *   of type 0 generated by the old software claim to have 0 length; for
 *   compatibility, code reading rasterfiles must be prepared to compute the
 *   true length from the width, height, and depth fields.
 */

int	pixout = 1;		/* 0 = bw(5) output, 1 = pix(5) output */
int	hflag;
int	inverted;
int	pure;			/* No Sun header */
int	verbose;

static char	*file_name;
static FILE	*fp = stdin;

char	usage[] = "\
Usage: sun-pix [-b -h -i -P -v] [sun.bitmap]\n";


#define NET_LONG_LEN	4	/* # bytes to network long */

unsigned long
getlong(msgp)
	char *msgp;
{
	register unsigned char *p = (unsigned char *) msgp;
	register unsigned long u;

	u = *p++; u <<= 8;
	u |= *p++; u <<= 8;
	u |= *p++; u <<= 8;
	return (u | *p);
}

get_args( argc, argv )
register char **argv;
{
	register int c;

	while ( (c = getopt( argc, argv, "bhiPv" )) != EOF )  {
		switch( c )  {
		case 'b':
			pixout = 0;	/* bw(5) */
			break;
		case 'h':
			hflag = 1;	/* print header */
			break;
		case 'i':
			inverted = 1;
			break;
		case 'P':
			pure = 1;
			break;
		case 'v':
			verbose = 1;
			break;

		default:		/* '?' */
			return(0);
		}
	}

	if( optind >= argc )  {
		if( isatty(fileno(stdin)) )
			return(0);
		file_name = "-";
		fp = stdin;
	} else {
		file_name = argv[optind];
		if( (fp = fopen(file_name, "r")) == NULL )  {
			(void)fprintf( stderr,
				"sun-pix: cannot open \"%s\" for reading\n",
				file_name );
			return(0);
		}
	}

	if ( argc > ++optind )
		(void)fprintf( stderr, "sun-pix: excess argument(s) ignored\n" );

	return(1);		/* OK */
}

main( argc, argv )
int argc;
char **argv;
{
	register int	i, c;

	if ( !get_args( argc, argv ) || (isatty(fileno(stdout)) && (hflag == 0)) ) {
		(void)fputs(usage, stderr);
		exit( 1 );
	}

	if( !pure )  {
		fread( inbuf, sizeof(struct rasterfile), 1, fp );

		header.ras_magic = getlong( &inbuf[NET_LONG_LEN*0] );
		header.ras_width = getlong( &inbuf[NET_LONG_LEN*1] );
		header.ras_height = getlong( &inbuf[NET_LONG_LEN*2] );
		header.ras_depth = getlong( &inbuf[NET_LONG_LEN*3] );
		header.ras_length = getlong( &inbuf[NET_LONG_LEN*4] );
		header.ras_type = getlong( &inbuf[NET_LONG_LEN*5] );
		header.ras_maptype = getlong( &inbuf[NET_LONG_LEN*6] );
		header.ras_maplength = getlong( &inbuf[NET_LONG_LEN*7] );

		if( header.ras_magic != RAS_MAGIC )  {
			fprintf(stderr,
				"sun-pix: bad magic number, was x%x, s/b x%x\n",
				header.ras_magic, RAS_MAGIC );
			exit(1);
		}

		if(verbose)  {
			fprintf( stderr, 
				"ras_width = %d, ras_height = %d\nras_depth = %d, ras_length = %d\n",
				header.ras_width, header.ras_height,
				header.ras_depth, header.ras_length );
			fprintf( stderr,
				"ras_type = %d, ras_maptype = %d, ras_maplength = %d\n",
				header.ras_type,
				header.ras_maptype,
				header.ras_maplength );
		}
		if( hflag ) {
			printf( "-w%d -n%d\n", header.ras_width, header.ras_height );
			exit( 0 );
		}
	} else {
		/* "pure" bitmap */
		header.ras_type = RT_STANDARD;
		header.ras_depth = 1;
	}

	switch( header.ras_type )  {
	case RT_OLD:		/* ??? */
	case RT_STANDARD:
		break;
	default:
		fprintf(stderr,"sun-pix:  Unable to process type %d images\n",
			header.ras_type );
		exit(1);
	}

	/*  Gobble colormap -- ought to know what to do with it */
	for( c=0; c<header.ras_maplength; c++)  {
		(void)getchar();
	}

	switch( header.ras_depth )  {
	case 1:
		/* 1-bit image */
		while( !feof(fp) ) {
			c = getc(fp);
			if( inverted ) {
				for( i = 0x80; i > 0; i >>= 1 )
					if( c & i ) {
						putchar( 0 );
						if(pixout){putchar(0);putchar(0);}
					} else {
						putchar( 255 );
						if(pixout){putchar(255);putchar(255);}
					}
			} else {
				for( i = 0x80; i > 0; i >>= 1 )
					if( c & i ) {
						putchar( 255 );
						if(pixout){putchar(255);putchar(255);}
					} else {
						putchar( 0 );
						if(pixout){putchar(0);putchar(0);}
					}
			}
		}
		break;
	default:
		fprintf(stderr,"sun-pix:  unable to handle depth=%d\n",
			header.ras_depth );
		exit(1);
	}
	exit(0);
}
