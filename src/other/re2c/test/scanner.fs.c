/* Generated by re2c */
#line 1 "scanner.fs.re"
/* $Id$ */
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include "scanner.h"
#include "parser.h"
#include "y.tab.h"
#include "globals.h"
#include "dfa.h"

extern YYSTYPE yylval;

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#define	BSIZE	8192

#define	YYCTYPE		char
#define	YYCURSOR	cursor
#define	YYLIMIT		lim
#define	YYMARKER	ptr
#define	YYFILL(n)	{cursor = fill(cursor);}

#define	RETURN(i)	{cur = cursor; return i;}

namespace re2c
{

Scanner::Scanner(std::istream& i, std::ostream& o)
	: in(i)
	, out(o)
	, bot(NULL), tok(NULL), ptr(NULL), cur(NULL), pos(NULL), lim(NULL)
	, top(NULL), eof(NULL), tchar(0), tline(0), cline(1), iscfg(0)
{
    ;
}

char *Scanner::fill(char *cursor)
{
	if(!eof)
	{
		uint cnt = tok - bot;
		if(cnt)
		{
			memcpy(bot, tok, lim - tok);
			tok = bot;
			ptr -= cnt;
			cursor -= cnt;
			pos -= cnt;
			lim -= cnt;
		}
		if((top - lim) < BSIZE)
		{
			char *buf = new char[(lim - bot) + BSIZE];
			memcpy(buf, tok, lim - tok);
			tok = buf;
			ptr = &buf[ptr - bot];
			cursor = &buf[cursor - bot];
			pos = &buf[pos - bot];
			lim = &buf[lim - bot];
			top = &lim[BSIZE];
			delete [] bot;
			bot = buf;
		}
		in.read(lim, BSIZE);
		if ((cnt = in.gcount()) != BSIZE )
		{
			eof = &lim[cnt]; *eof++ = '\0';
		}
		lim += cnt;
	}
	return cursor;
}

#line 95 "scanner.fs.re"


int Scanner::echo()
{
    char *cursor = cur;
    bool ignore_eoc = false;

    if (eof && cursor == eof) // Catch EOF
	{
    	return 0;
	}

    tok = cursor;
echo:

#line 96 "<stdout>"

	switch (YYGETSTATE()) {
	default: goto yy0;
	case 0: goto yyFillLabel0;
	case 1: goto yyFillLabel1;
	case 2: goto yyFillLabel2;
	case 3: goto yyFillLabel3;
	case 4: goto yyFillLabel4;
	case 5: goto yyFillLabel5;
	case 6: goto yyFillLabel6;
	case 7: goto yyFillLabel7;
	case 8: goto yyFillLabel8;
	case 9: goto yyFillLabel9;
	case 10: goto yyFillLabel10;
	case 11: goto yyFillLabel11;
	case 12: goto yyFillLabel12;
	case 13: goto yyFillLabel13;
	case 14: goto yyFillLabel14;
	case 15: goto yyFillLabel15;
	case 16: goto yyFillLabel16;
	case 17: goto yyFillLabel17;
	case 18: goto yyFillLabel18;
	case 19: goto yyFillLabel19;
	case 20: goto yyFillLabel20;
	case 21: goto yyFillLabel21;
	case 22: goto yyFillLabel22;
	case 23: goto yyFillLabel23;
	case 24: goto yyFillLabel24;
	case 25: goto yyFillLabel25;
	case 26: goto yyFillLabel26;
	case 27: goto yyFillLabel27;
	case 28: goto yyFillLabel28;
	case 29: goto yyFillLabel29;
	case 30: goto yyFillLabel30;
	case 31: goto yyFillLabel31;
	case 32: goto yyFillLabel32;
	case 33: goto yyFillLabel33;
	case 34: goto yyFillLabel34;
	case 35: goto yyFillLabel35;
	}
yy0:
	YYSETSTATE(0);
	if ((YYLIMIT - YYCURSOR) < 11) YYFILL(11);
yyFillLabel0:
	yych = *YYCURSOR;
	if (yych <= ')') {
		if (yych <= 0x00) goto yy7;
		if (yych == '\n') goto yy5;
		goto yy9;
	} else {
		if (yych <= '*') goto yy4;
		if (yych != '/') goto yy9;
	}
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych == '*') goto yy12;
yy3:
#line 141 "scanner.fs.re"
	{
					goto echo;
				}
#line 157 "<stdout>"
yy4:
	yych = *++YYCURSOR;
	if (yych == '/') goto yy10;
	goto yy3;
yy5:
	++YYCURSOR;
#line 130 "scanner.fs.re"
	{
					out.write((const char*)(tok), (const char*)(cursor) - (const char*)(tok));
					tok = pos = cursor; cline++;
				  	goto echo;
				}
#line 170 "<stdout>"
yy7:
	++YYCURSOR;
#line 135 "scanner.fs.re"
	{
					out.write((const char*)(tok), (const char*)(cursor) - (const char*)(tok) - 1); // -1 so we don't write out the \0
					if(cursor == eof) {
						RETURN(0);
					}
				}
#line 180 "<stdout>"
yy9:
	yych = *++YYCURSOR;
	goto yy3;
yy10:
	++YYCURSOR;
#line 121 "scanner.fs.re"
	{
					if (ignore_eoc) {
						ignore_eoc = false;
					} else {
						out.write((const char*)(tok), (const char*)(cursor) - (const char*)(tok));
					}
					tok = pos = cursor;
					goto echo;
				}
#line 196 "<stdout>"
yy12:
	yych = *++YYCURSOR;
	if (yych == '!') goto yy14;
yy13:
	YYCURSOR = YYMARKER;
	goto yy3;
yy14:
	yych = *++YYCURSOR;
	if (yych == 'm') goto yy15;
	if (yych == 'r') goto yy16;
	goto yy13;
yy15:
	yych = *++YYCURSOR;
	if (yych == 'a') goto yy21;
	goto yy13;
yy16:
	yych = *++YYCURSOR;
	if (yych != 'e') goto yy13;
	yych = *++YYCURSOR;
	if (yych != '2') goto yy13;
	yych = *++YYCURSOR;
	if (yych != 'c') goto yy13;
	++YYCURSOR;
#line 110 "scanner.fs.re"
	{ 
					out.write((const char*)(tok), (const char*)(&cursor[-7]) - (const char*)(tok));
					tok = cursor;
					RETURN(1);
				}
#line 226 "<stdout>"
yy21:
	yych = *++YYCURSOR;
	if (yych != 'x') goto yy13;
	yych = *++YYCURSOR;
	if (yych != ':') goto yy13;
	yych = *++YYCURSOR;
	if (yych != 'r') goto yy13;
	yych = *++YYCURSOR;
	if (yych != 'e') goto yy13;
	yych = *++YYCURSOR;
	if (yych != '2') goto yy13;
	yych = *++YYCURSOR;
	if (yych != 'c') goto yy13;
	++YYCURSOR;
#line 115 "scanner.fs.re"
	{
					out << "#define YYMAXFILL " << maxFill << std::endl;
					tok = pos = cursor;
					ignore_eoc = true;
					goto echo;
				}
#line 248 "<stdout>"
#line 144 "scanner.fs.re"

}


int Scanner::scan()
{
    char *cursor = cur;
    uint depth;

scan:
    tchar = cursor - pos;
    tline = cline;
    tok = cursor;
	if (iscfg == 1)
	{
		goto config;
	}
	else if (iscfg == 2)
	{
   		goto value;
    }

#line 272 "<stdout>"

	YYSETSTATE(1);
	if ((YYLIMIT - YYCURSOR) < 5) YYFILL(5);
yyFillLabel1:
	yych = *YYCURSOR;
	if (yych <= '/') {
		if (yych <= '!') {
			if (yych <= '\f') {
				if (yych <= 0x08) goto yy56;
				if (yych <= '\t') goto yy50;
				if (yych <= '\n') goto yy52;
				goto yy56;
			} else {
				if (yych <= '\r') goto yy54;
				if (yych == ' ') goto yy50;
				goto yy56;
			}
		} else {
			if (yych <= ')') {
				if (yych <= '"') goto yy37;
				if (yych <= '&') goto yy56;
				if (yych <= '\'') goto yy39;
				goto yy43;
			} else {
				if (yych <= '+') {
					if (yych <= '*') goto yy35;
					goto yy44;
				} else {
					if (yych <= '-') goto yy56;
					if (yych <= '.') goto yy48;
					goto yy33;
				}
			}
		}
	} else {
		if (yych <= 'Z') {
			if (yych <= '=') {
				if (yych == ';') goto yy43;
				if (yych <= '<') goto yy56;
				goto yy43;
			} else {
				if (yych == '?') goto yy44;
				if (yych <= '@') goto yy56;
				goto yy47;
			}
		} else {
			if (yych <= 'q') {
				if (yych <= '[') goto yy41;
				if (yych <= '\\') goto yy43;
				if (yych <= '`') goto yy56;
				goto yy47;
			} else {
				if (yych <= 'z') {
					if (yych <= 'r') goto yy45;
					goto yy47;
				} else {
					if (yych <= '{') goto yy31;
					if (yych <= '|') goto yy43;
					goto yy56;
				}
			}
		}
	}
yy31:
	yyaccept = 0;
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych <= '/') {
		if (yych == ',') goto yy97;
	} else {
		if (yych <= '0') goto yy94;
		if (yych <= '9') goto yy95;
	}
yy32:
#line 166 "scanner.fs.re"
	{ depth = 1;
				  goto code;
				}
#line 350 "<stdout>"
yy33:
	++YYCURSOR;
	if ((yych = *YYCURSOR) == '*') goto yy92;
yy34:
#line 196 "scanner.fs.re"
	{ RETURN(*tok); }
#line 357 "<stdout>"
yy35:
	++YYCURSOR;
	if ((yych = *YYCURSOR) == '/') goto yy90;
yy36:
#line 198 "scanner.fs.re"
	{ yylval.op = *tok;
				  RETURN(CLOSE); }
#line 365 "<stdout>"
yy37:
	yyaccept = 1;
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych != '\n') goto yy86;
yy38:
#line 183 "scanner.fs.re"
	{ fatal("unterminated string constant (missing \")"); }
#line 373 "<stdout>"
yy39:
	yyaccept = 2;
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych != '\n') goto yy81;
yy40:
#line 184 "scanner.fs.re"
	{ fatal("unterminated string constant (missing ')"); }
#line 381 "<stdout>"
yy41:
	yyaccept = 3;
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych == '\n') goto yy42;
	if (yych == '^') goto yy72;
	goto yy71;
yy42:
#line 194 "scanner.fs.re"
	{ fatal("unterminated range (missing ])"); }
#line 391 "<stdout>"
yy43:
	yych = *++YYCURSOR;
	goto yy34;
yy44:
	yych = *++YYCURSOR;
	goto yy36;
yy45:
	++YYCURSOR;
	if ((yych = *YYCURSOR) == 'e') goto yy62;
	goto yy61;
yy46:
#line 225 "scanner.fs.re"
	{ cur = cursor;
				  yylval.symbol = Symbol::find(token());
				  return ID; }
#line 407 "<stdout>"
yy47:
	yych = *++YYCURSOR;
	goto yy61;
yy48:
	++YYCURSOR;
#line 229 "scanner.fs.re"
	{ cur = cursor;
				  yylval.regexp = mkDot();
				  return RANGE;
				}
#line 418 "<stdout>"
yy50:
	++YYCURSOR;
	yych = *YYCURSOR;
	goto yy59;
yy51:
#line 234 "scanner.fs.re"
	{ goto scan; }
#line 426 "<stdout>"
yy52:
	++YYCURSOR;
yy53:
#line 236 "scanner.fs.re"
	{ if(cursor == eof) RETURN(0);
				  pos = cursor; cline++;
				  goto scan;
	    			}
#line 435 "<stdout>"
yy54:
	++YYCURSOR;
	if ((yych = *YYCURSOR) == '\n') goto yy57;
yy55:
#line 241 "scanner.fs.re"
	{ std::ostringstream msg;
				  msg << "unexpected character: ";
				  prtChOrHex(msg, *tok);
				  fatal(msg.str().c_str());
				  goto scan;
				}
#line 447 "<stdout>"
yy56:
	yych = *++YYCURSOR;
	goto yy55;
yy57:
	yych = *++YYCURSOR;
	goto yy53;
yy58:
	++YYCURSOR;
	YYSETSTATE(2);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel2:
	yych = *YYCURSOR;
yy59:
	if (yych == '\t') goto yy58;
	if (yych == ' ') goto yy58;
	goto yy51;
yy60:
	++YYCURSOR;
	YYSETSTATE(3);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel3:
	yych = *YYCURSOR;
yy61:
	if (yych <= '@') {
		if (yych <= '/') goto yy46;
		if (yych <= '9') goto yy60;
		goto yy46;
	} else {
		if (yych <= 'Z') goto yy60;
		if (yych <= '`') goto yy46;
		if (yych <= 'z') goto yy60;
		goto yy46;
	}
yy62:
	yych = *++YYCURSOR;
	if (yych != '2') goto yy61;
	yych = *++YYCURSOR;
	if (yych != 'c') goto yy61;
	yyaccept = 4;
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych != ':') goto yy61;
yy65:
	++YYCURSOR;
	YYSETSTATE(4);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel4:
	yych = *YYCURSOR;
	if (yych <= '@') goto yy66;
	if (yych <= 'Z') goto yy67;
	if (yych <= '`') goto yy66;
	if (yych <= 'z') goto yy67;
yy66:
	YYCURSOR = YYMARKER;
	if (yyaccept <= 3) {
		if (yyaccept <= 1) {
			if (yyaccept <= 0) {
				goto yy32;
			} else {
				goto yy38;
			}
		} else {
			if (yyaccept <= 2) {
				goto yy40;
			} else {
				goto yy42;
			}
		}
	} else {
		if (yyaccept <= 5) {
			if (yyaccept <= 4) {
				goto yy46;
			} else {
				goto yy69;
			}
		} else {
			goto yy98;
		}
	}
yy67:
	yyaccept = 5;
	YYMARKER = ++YYCURSOR;
	YYSETSTATE(5);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel5:
	yych = *YYCURSOR;
	if (yych <= 'Z') {
		if (yych <= '9') {
			if (yych >= '0') goto yy67;
		} else {
			if (yych <= ':') goto yy65;
			if (yych >= 'A') goto yy67;
		}
	} else {
		if (yych <= '_') {
			if (yych >= '_') goto yy67;
		} else {
			if (yych <= '`') goto yy69;
			if (yych <= 'z') goto yy67;
		}
	}
yy69:
#line 218 "scanner.fs.re"
	{ cur = cursor;
				  tok+= 5; /* skip "re2c:" */
				  iscfg = 1;
				  yylval.str = new Str(token());
				  return CONFIG;
				}
#line 556 "<stdout>"
yy70:
	++YYCURSOR;
	YYSETSTATE(6);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel6:
	yych = *YYCURSOR;
yy71:
	if (yych <= '[') {
		if (yych == '\n') goto yy66;
		goto yy70;
	} else {
		if (yych <= '\\') goto yy74;
		if (yych <= ']') goto yy75;
		goto yy70;
	}
yy72:
	++YYCURSOR;
	YYSETSTATE(7);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel7:
	yych = *YYCURSOR;
	if (yych <= '[') {
		if (yych == '\n') goto yy66;
		goto yy72;
	} else {
		if (yych <= '\\') goto yy77;
		if (yych <= ']') goto yy78;
		goto yy72;
	}
yy74:
	++YYCURSOR;
	YYSETSTATE(8);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel8:
	yych = *YYCURSOR;
	if (yych == '\n') goto yy66;
	goto yy70;
yy75:
	++YYCURSOR;
#line 190 "scanner.fs.re"
	{ cur = cursor;
				  yylval.regexp = ranToRE(token());
				  return RANGE; }
#line 600 "<stdout>"
yy77:
	++YYCURSOR;
	YYSETSTATE(9);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel9:
	yych = *YYCURSOR;
	if (yych == '\n') goto yy66;
	goto yy72;
yy78:
	++YYCURSOR;
#line 186 "scanner.fs.re"
	{ cur = cursor;
				  yylval.regexp = invToRE(token());
				  return RANGE; }
#line 615 "<stdout>"
yy80:
	++YYCURSOR;
	YYSETSTATE(10);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel10:
	yych = *YYCURSOR;
yy81:
	if (yych <= '&') {
		if (yych == '\n') goto yy66;
		goto yy80;
	} else {
		if (yych <= '\'') goto yy83;
		if (yych != '\\') goto yy80;
	}
	++YYCURSOR;
	YYSETSTATE(11);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel11:
	yych = *YYCURSOR;
	if (yych == '\n') goto yy66;
	goto yy80;
yy83:
	++YYCURSOR;
#line 179 "scanner.fs.re"
	{ cur = cursor;
				  yylval.regexp = strToCaseInsensitiveRE(token());
				  return STRING; }
#line 643 "<stdout>"
yy85:
	++YYCURSOR;
	YYSETSTATE(12);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel12:
	yych = *YYCURSOR;
yy86:
	if (yych <= '!') {
		if (yych == '\n') goto yy66;
		goto yy85;
	} else {
		if (yych <= '"') goto yy88;
		if (yych != '\\') goto yy85;
	}
	++YYCURSOR;
	YYSETSTATE(13);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel13:
	yych = *YYCURSOR;
	if (yych == '\n') goto yy66;
	goto yy85;
yy88:
	++YYCURSOR;
#line 175 "scanner.fs.re"
	{ cur = cursor;
				  yylval.regexp = strToRE(token());
				  return STRING; }
#line 671 "<stdout>"
yy90:
	++YYCURSOR;
#line 172 "scanner.fs.re"
	{ tok = cursor;
				  RETURN(0); }
#line 677 "<stdout>"
yy92:
	++YYCURSOR;
#line 169 "scanner.fs.re"
	{ depth = 1;
				  goto comment; }
#line 683 "<stdout>"
yy94:
	yych = *++YYCURSOR;
	if (yych == ',') goto yy108;
	goto yy96;
yy95:
	++YYCURSOR;
	YYSETSTATE(14);
	if ((YYLIMIT - YYCURSOR) < 2) YYFILL(2);
yyFillLabel14:
	yych = *YYCURSOR;
yy96:
	if (yych <= '/') {
		if (yych == ',') goto yy101;
		goto yy66;
	} else {
		if (yych <= '9') goto yy95;
		if (yych == '}') goto yy99;
		goto yy66;
	}
yy97:
	++YYCURSOR;
yy98:
#line 216 "scanner.fs.re"
	{ fatal("illegal closure form, use '{n}', '{n,}', '{n,m}' where n and m are numbers"); }
#line 708 "<stdout>"
yy99:
	++YYCURSOR;
#line 204 "scanner.fs.re"
	{ yylval.extop.minsize = atoi((char *)tok+1);
				  yylval.extop.maxsize = atoi((char *)tok+1);
				  RETURN(CLOSESIZE); }
#line 715 "<stdout>"
yy101:
	yyaccept = 6;
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych <= '/') goto yy98;
	if (yych <= '9') goto yy104;
	if (yych != '}') goto yy98;
	++YYCURSOR;
#line 212 "scanner.fs.re"
	{ yylval.extop.minsize = atoi((char *)tok+1);
				  yylval.extop.maxsize = -1;
				  RETURN(CLOSESIZE); }
#line 727 "<stdout>"
yy104:
	++YYCURSOR;
	YYSETSTATE(15);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel15:
	yych = *YYCURSOR;
	if (yych <= '/') goto yy66;
	if (yych <= '9') goto yy104;
	if (yych != '}') goto yy66;
	++YYCURSOR;
#line 208 "scanner.fs.re"
	{ yylval.extop.minsize = atoi((char *)tok+1);
				  yylval.extop.maxsize = MAX(yylval.extop.minsize,atoi(strchr((char *)tok, ',')+1));
				  RETURN(CLOSESIZE); }
#line 742 "<stdout>"
yy108:
	yyaccept = 6;
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych <= '/') goto yy98;
	if (yych <= '9') goto yy104;
	if (yych != '}') goto yy98;
	++YYCURSOR;
#line 201 "scanner.fs.re"
	{ yylval.op = '*';
				  RETURN(CLOSE); }
#line 753 "<stdout>"
#line 247 "scanner.fs.re"


code:

#line 759 "<stdout>"

	YYSETSTATE(16);
	if ((YYLIMIT - YYCURSOR) < 2) YYFILL(2);
yyFillLabel16:
	yych = *YYCURSOR;
	if (yych <= '&') {
		if (yych <= '\n') {
			if (yych <= '\t') goto yy119;
			goto yy117;
		} else {
			if (yych == '"') goto yy121;
			goto yy119;
		}
	} else {
		if (yych <= '{') {
			if (yych <= '\'') goto yy122;
			if (yych <= 'z') goto yy119;
			goto yy115;
		} else {
			if (yych != '}') goto yy119;
		}
	}
	++YYCURSOR;
#line 251 "scanner.fs.re"
	{ if(--depth == 0){
					cur = cursor;
					yylval.token = new Token(token(), tline);
					return CODE;
				  }
				  goto code; }
#line 790 "<stdout>"
yy115:
	++YYCURSOR;
#line 257 "scanner.fs.re"
	{ ++depth;
				  goto code; }
#line 796 "<stdout>"
yy117:
	++YYCURSOR;
#line 259 "scanner.fs.re"
	{ if(cursor == eof) fatal("missing '}'");
				  pos = cursor; cline++;
				  goto code;
				}
#line 804 "<stdout>"
yy119:
	++YYCURSOR;
yy120:
#line 263 "scanner.fs.re"
	{ goto code; }
#line 810 "<stdout>"
yy121:
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych == '\n') goto yy120;
	goto yy128;
yy122:
	yych = *(YYMARKER = ++YYCURSOR);
	if (yych == '\n') goto yy120;
	goto yy124;
yy123:
	++YYCURSOR;
	YYSETSTATE(17);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel17:
	yych = *YYCURSOR;
yy124:
	if (yych <= '&') {
		if (yych != '\n') goto yy123;
	} else {
		if (yych <= '\'') goto yy119;
		if (yych == '\\') goto yy126;
		goto yy123;
	}
yy125:
	YYCURSOR = YYMARKER;
	goto yy120;
yy126:
	++YYCURSOR;
	YYSETSTATE(18);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel18:
	yych = *YYCURSOR;
	if (yych == '\n') goto yy125;
	goto yy123;
yy127:
	++YYCURSOR;
	YYSETSTATE(19);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel19:
	yych = *YYCURSOR;
yy128:
	if (yych <= '!') {
		if (yych == '\n') goto yy125;
		goto yy127;
	} else {
		if (yych <= '"') goto yy119;
		if (yych != '\\') goto yy127;
	}
	++YYCURSOR;
	YYSETSTATE(20);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel20:
	yych = *YYCURSOR;
	if (yych == '\n') goto yy125;
	goto yy127;
#line 264 "scanner.fs.re"


comment:

#line 870 "<stdout>"

	YYSETSTATE(21);
	if ((YYLIMIT - YYCURSOR) < 2) YYFILL(2);
yyFillLabel21:
	yych = *YYCURSOR;
	if (yych <= ')') {
		if (yych == '\n') goto yy135;
		goto yy137;
	} else {
		if (yych <= '*') goto yy132;
		if (yych == '/') goto yy134;
		goto yy137;
	}
yy132:
	++YYCURSOR;
	if ((yych = *YYCURSOR) == '/') goto yy140;
yy133:
#line 279 "scanner.fs.re"
	{ if(cursor == eof) RETURN(0);
				  goto comment; }
#line 891 "<stdout>"
yy134:
	yych = *++YYCURSOR;
	if (yych == '*') goto yy138;
	goto yy133;
yy135:
	++YYCURSOR;
#line 275 "scanner.fs.re"
	{ if(cursor == eof) RETURN(0);
				  tok = pos = cursor; cline++;
				  goto comment;
				}
#line 903 "<stdout>"
yy137:
	yych = *++YYCURSOR;
	goto yy133;
yy138:
	++YYCURSOR;
#line 272 "scanner.fs.re"
	{ ++depth;
				  fatal("ambiguous /* found");
				  goto comment; }
#line 913 "<stdout>"
yy140:
	++YYCURSOR;
#line 268 "scanner.fs.re"
	{ if(--depth == 0)
					goto scan;
				    else
					goto comment; }
#line 921 "<stdout>"
#line 281 "scanner.fs.re"


config:

#line 927 "<stdout>"

	YYSETSTATE(22);
	if ((YYLIMIT - YYCURSOR) < 2) YYFILL(2);
yyFillLabel22:
	yych = *YYCURSOR;
	if (yych <= 0x1F) {
		if (yych != '\t') goto yy148;
	} else {
		if (yych <= ' ') goto yy144;
		if (yych == '=') goto yy146;
		goto yy148;
	}
yy144:
	++YYCURSOR;
	yych = *YYCURSOR;
	goto yy153;
yy145:
#line 285 "scanner.fs.re"
	{ goto config; }
#line 947 "<stdout>"
yy146:
	++YYCURSOR;
	yych = *YYCURSOR;
	goto yy151;
yy147:
#line 286 "scanner.fs.re"
	{ iscfg = 2;
				  cur = cursor;
				  RETURN('='); 
				}
#line 958 "<stdout>"
yy148:
	++YYCURSOR;
#line 290 "scanner.fs.re"
	{ fatal("missing '='"); }
#line 963 "<stdout>"
yy150:
	++YYCURSOR;
	YYSETSTATE(23);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel23:
	yych = *YYCURSOR;
yy151:
	if (yych == '\t') goto yy150;
	if (yych == ' ') goto yy150;
	goto yy147;
yy152:
	++YYCURSOR;
	YYSETSTATE(24);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel24:
	yych = *YYCURSOR;
yy153:
	if (yych == '\t') goto yy152;
	if (yych == ' ') goto yy152;
	goto yy145;
#line 291 "scanner.fs.re"


value:

#line 989 "<stdout>"

	YYSETSTATE(25);
	if ((YYLIMIT - YYCURSOR) < 2) YYFILL(2);
yyFillLabel25:
	yych = *YYCURSOR;
	if (yych <= '&') {
		if (yych <= '\r') {
			if (yych <= 0x08) goto yy162;
			if (yych <= '\n') goto yy156;
			if (yych <= '\f') goto yy162;
		} else {
			if (yych <= ' ') {
				if (yych <= 0x1F) goto yy162;
			} else {
				if (yych == '"') goto yy164;
				goto yy162;
			}
		}
	} else {
		if (yych <= '/') {
			if (yych <= '\'') goto yy166;
			if (yych == '-') goto yy159;
			goto yy162;
		} else {
			if (yych <= '9') {
				if (yych <= '0') goto yy157;
				goto yy160;
			} else {
				if (yych != ';') goto yy162;
			}
		}
	}
yy156:
#line 300 "scanner.fs.re"
	{ cur = cursor;
				  yylval.str = new Str(token());
				  iscfg = 0;
				  return VALUE;
				}
#line 1029 "<stdout>"
yy157:
	++YYCURSOR;
	if ((yych = *YYCURSOR) <= '\r') {
		if (yych <= 0x08) goto yy162;
		if (yych <= '\n') goto yy158;
		if (yych <= '\f') goto yy162;
	} else {
		if (yych <= ' ') {
			if (yych <= 0x1F) goto yy162;
		} else {
			if (yych != ';') goto yy162;
		}
	}
yy158:
#line 295 "scanner.fs.re"
	{ cur = cursor;
				  yylval.number = atoi(token().to_string().c_str());
				  iscfg = 0;
				  return NUMBER;
				}
#line 1050 "<stdout>"
yy159:
	yych = *++YYCURSOR;
	if (yych <= '0') goto yy163;
	if (yych >= ':') goto yy163;
yy160:
	++YYCURSOR;
	YYSETSTATE(26);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel26:
	yych = *YYCURSOR;
	if (yych <= 0x1F) {
		if (yych <= '\n') {
			if (yych >= '\t') goto yy158;
		} else {
			if (yych == '\r') goto yy158;
		}
	} else {
		if (yych <= '9') {
			if (yych <= ' ') goto yy158;
			if (yych >= '0') goto yy160;
		} else {
			if (yych == ';') goto yy158;
		}
	}
yy162:
	++YYCURSOR;
	YYSETSTATE(27);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel27:
	yych = *YYCURSOR;
yy163:
	if (yych <= '\r') {
		if (yych <= 0x08) goto yy162;
		if (yych <= '\n') goto yy156;
		if (yych <= '\f') goto yy162;
		goto yy156;
	} else {
		if (yych <= ' ') {
			if (yych <= 0x1F) goto yy162;
			goto yy156;
		} else {
			if (yych == ';') goto yy156;
			goto yy162;
		}
	}
yy164:
	YYMARKER = ++YYCURSOR;
	YYSETSTATE(28);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel28:
	yych = *YYCURSOR;
	if (yych <= ' ') {
		if (yych <= '\n') {
			if (yych <= 0x08) goto yy164;
			if (yych <= '\t') goto yy174;
			goto yy156;
		} else {
			if (yych == '\r') goto yy174;
			if (yych <= 0x1F) goto yy164;
			goto yy174;
		}
	} else {
		if (yych <= ':') {
			if (yych == '"') goto yy162;
			goto yy164;
		} else {
			if (yych <= ';') goto yy174;
			if (yych == '\\') goto yy176;
			goto yy164;
		}
	}
yy166:
	YYMARKER = ++YYCURSOR;
	YYSETSTATE(29);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel29:
	yych = *YYCURSOR;
	if (yych <= ' ') {
		if (yych <= '\n') {
			if (yych <= 0x08) goto yy166;
			if (yych >= '\n') goto yy156;
		} else {
			if (yych == '\r') goto yy168;
			if (yych <= 0x1F) goto yy166;
		}
	} else {
		if (yych <= ':') {
			if (yych == '\'') goto yy162;
			goto yy166;
		} else {
			if (yych <= ';') goto yy168;
			if (yych == '\\') goto yy171;
			goto yy166;
		}
	}
yy168:
	++YYCURSOR;
	YYSETSTATE(30);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel30:
	yych = *YYCURSOR;
	if (yych <= '&') {
		if (yych != '\n') goto yy168;
	} else {
		if (yych <= '\'') goto yy172;
		if (yych == '\\') goto yy173;
		goto yy168;
	}
yy170:
	YYCURSOR = YYMARKER;
	goto yy156;
yy171:
	YYMARKER = ++YYCURSOR;
	YYSETSTATE(31);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel31:
	yych = *YYCURSOR;
	if (yych <= '\r') {
		if (yych <= '\t') {
			if (yych <= 0x08) goto yy166;
			goto yy168;
		} else {
			if (yych <= '\n') goto yy156;
			if (yych <= '\f') goto yy166;
			goto yy168;
		}
	} else {
		if (yych <= ' ') {
			if (yych <= 0x1F) goto yy166;
			goto yy168;
		} else {
			if (yych == ';') goto yy168;
			goto yy166;
		}
	}
yy172:
	yych = *++YYCURSOR;
	goto yy156;
yy173:
	++YYCURSOR;
	YYSETSTATE(32);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel32:
	yych = *YYCURSOR;
	if (yych == '\n') goto yy170;
	goto yy168;
yy174:
	++YYCURSOR;
	YYSETSTATE(33);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel33:
	yych = *YYCURSOR;
	if (yych <= '!') {
		if (yych == '\n') goto yy170;
		goto yy174;
	} else {
		if (yych <= '"') goto yy172;
		if (yych == '\\') goto yy177;
		goto yy174;
	}
yy176:
	YYMARKER = ++YYCURSOR;
	YYSETSTATE(34);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel34:
	yych = *YYCURSOR;
	if (yych <= '\r') {
		if (yych <= '\t') {
			if (yych <= 0x08) goto yy164;
			goto yy174;
		} else {
			if (yych <= '\n') goto yy156;
			if (yych <= '\f') goto yy164;
			goto yy174;
		}
	} else {
		if (yych <= ' ') {
			if (yych <= 0x1F) goto yy164;
			goto yy174;
		} else {
			if (yych == ';') goto yy174;
			goto yy164;
		}
	}
yy177:
	++YYCURSOR;
	YYSETSTATE(35);
	if (YYLIMIT <= YYCURSOR) YYFILL(1);
yyFillLabel35:
	yych = *YYCURSOR;
	if (yych == '\n') goto yy170;
	goto yy174;
#line 305 "scanner.fs.re"

}

void Scanner::fatal(uint ofs, const char *msg) const
{
	out.flush();
	std::cerr << "re2c: error: "
		<< "line " << tline << ", column " << (tchar + ofs + 1) << ": "
		<< msg << std::endl;
   	exit(1);
}

} // end namespace re2c

