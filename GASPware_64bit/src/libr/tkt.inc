	IMPLICIT NONE
	integer  STAT
	PARAMETER (STAT =5)
	integer HOMEC
	PARAMETER (HOMEC=12)
	integer STRTX
	PARAMETER (STRTX=26)
	integer PFIX
	PARAMETER (PFIX =27)
	integer GMODE
	PARAMETER (GMODE=29)
	integer AMODE
	PARAMETER (AMODE=31)
	integer XMIN
	PARAMETER (XMIN =0)
	integer XMAX
	PARAMETER (XMAX =1023)
	integer YMIN
	PARAMETER (YMIN =0)
	integer YMAX
	PARAMETER (YMAX =780)
	integer YHOME
	PARAMETER (YHOME=767)
	INTEGER XLOW  ,XHIG  ,YLOW  ,YHIG
	INTEGER XLOW_N,XHIG_N,YLOW_N,YHIG_N
	INTEGER XPOS,YPOS,IPUNT
	INTEGER ALLOCATED
	INTEGER TTCHAN
	integer BUFLEN
	PARAMETER (BUFLEN=128)
	BYTE BUF(BUFLEN)
	COMMON /TKTCOM/ XLOW  ,XHIG  ,YLOW  ,YHIG  ,
     $                  XLOW_N,XHIG_N,YLOW_N,YHIG_N,
     $                  XPOS,YPOS,ALLOCATED,TTCHAN,IPUNT,BUF
