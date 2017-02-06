/****************************************************************************
!C

!File: jdoy.c

!Description:
  Contains main routine for converting Julian day to day of month and vice versa. 

!Input Parameters: (none)

!Input/Output Parameters: (none)

!Output Parameters: (none)

!Revision History:
  Original January/February 1998. Version 1.0
  Modified March/April 1998. (Removed Fisher's Perl code completely)
                             (Wrote the meta.c to read metadata as global attribute)
                             (Incorporated new command line format)
                             (Uses routine from meta.c to read metadata)
  Modified June 1998 Version 1.1 (process multiple input files)

!Team-unique Header:
  This software was developed by:
    Land Data Operational Product Evaluation (LDOPE) Team for the
    National Aeronautics and Space Administration, Goddard Space Flight
    Center

!References and Credits:

    Sadashiva Devadiga
    LDOPE                             Science Systems and Applications Inc.
    devadiga@ltpmail.gsfc.nasa.gov    NASA/GSFC Code 922
    phone: 301-614-5449               Greenbelt, MD 20771

!Design Notes: (none)
  call the routine with *mm = 0 if the input is in jjj-yyyy format

!END
*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>

#include "qa_tool.h"

#define HELP \
"NAME \n" \
"	jdoy - convert julian day to day of month and vice versa \n" \
"\n" \
"SYNOPSIS \n" \
"	jody [-help] \n" \
"	jdoy [yyyy-mm-dd] [yyyy-jjj] \n" \
" \n" \
"	echo [yyyy-mm-dd] [yyyy-jjj] | jdoy \n" \
" \n" \
"DESCRIPTION \n" \
"	Convert the input date in Julian day number to date in month day and \n" \
" year format. If the date input is in month day and year format then the \n" \
" routine will return the date in julian day number. \n" \
" \n" \
"OPTIONS \n" \
"	-help                       Display this help message \n" \
"	yyyy-mm-dd                  Input date in day of month format \n" \
"	                            '-' can be omitted or replaced by a space\n" \
"	yyyy-jjj                    Input date in julian day number \n" \
"	                            '-' can be omitted or replaced by a space\n" \
"	                            Year can omitted for the current year \n" \
" \n" \
"Examples: \n" \
"	jdoy 1997-2-25 \n" \
"	jdoy 1997265 \n" \
" \n" \
"Last update: 05/25/1999 \n" \
"Please report problems to Sadashiva Devadiga (devadiga@ltpmail.gsfc.nasa.gov) \n"

#define USAGE \
"usage:		jdoy [yyyy-mm-dd] [yyyy-jjj] \n" \
" \n" \
"OPTIONS \n" \
"	-help                       Display this help message \n" \
"	yyyy-mm-dd                  Input date in day of month format \n" \
"	                            '-' can be omitted or replaced by a space\n" \
"	yyyy-jjj                    Input date in julian day number \n" \
"	                            '-' can be omitted or replaced by a space\n" \
"	                            Year can omitted for the current year \n" \
" \n"

int parse_cmd_jdoy(char *s, int *mm, int *dd, int *yy);
int conv_date(int *mm, int *dd, int yyyy);

int main(int argc, char *argv[])
{
  int st = 1;
  int mm, dd, yy;
  char date_str[50];

  if ((argc == 2) && ((strcmp(argv[1], "-h")==0) || (strcmp(argv[1], "-help")==0)))
  {
    fprintf(stderr, "%s\n", HELP);
    exit(0);
  }
  if (argc == 2) strcpy(date_str, argv[1]);
  else if (argc == 3) sprintf(date_str, "%s %s", argv[1], argv[2]);
  else if (argc == 4) sprintf(date_str, "%s %s %s", argv[1], argv[2], argv[3]);
  if ((argc == 1) || ((st = parse_cmd_jdoy(date_str, &mm, &dd, &yy)) == -1))
  {
    fprintf(stderr, "%s\n", USAGE);
    exit(1);
  }
  else
  {
    st = (mm == 0) ? 2 : 3;
    if (conv_date(&mm, &dd, yy) != -1)
    {
      if (st == 2) 
	fprintf(stdout, "\tJulian day: %s\tDate(y m d): %d %d %d\n", date_str, yy, mm, dd);
      else
        fprintf(stdout, "\tDate: %s\tJulian day: %d %d\n", date_str, yy, dd);
    }
  }
  return 0;
}

int parse_cmd_jdoy(char *s, int *mm, int *dd, int *yy)
{
  int st = 1;
  int cnt, len, len2;
  int p1, p2, k, i, p;
  time_t cur_time, tloc;
  char tmp[50], yy_str[10], *time_str;

  len = (int)strlen(s);
  if (len == 0) return -1;
  for (i=0, cnt=0; i<len; i++)
    if ((s[i] == '-') || (s[i] == ' ')) 
    {
      cnt++;
      if (cnt == 1) p1 = i;
      else if (cnt == 2) p2 = i;
    }
  if (cnt == 0)
  {
    if (len <= 3) 
    {
      *mm = 0;
      *dd = (int)atoi(s);
      tloc = 0;
      cur_time = time(&tloc);
      time_str = ctime(&cur_time);
      len2 = (int)strlen(time_str);
      p = len2 - 5;
      for (k=p; k<=len2; k++)
        yy_str[k-p] = time_str[k];
      *yy = (int)atoi(yy_str);
      free(time_str);
    }
    else
    {
      for (i=0; i<4; i++) tmp[i] = s[i]; tmp[4] = '\0'; 
      *yy = (int)atoi(tmp);
      if (len == 8)
      {
        tmp[0] = s[4]; tmp[1] = s[5]; tmp[2] = '\0'; 
        *mm = (int)atoi(tmp);
        tmp[0] = s[6]; tmp[1] = s[7]; tmp[2] = '\0'; 
        *dd = (int)atoi(tmp);
      }
      else if (len < 8)
      {
	*mm = 0;
	for (i=4; i<=len; i++) tmp[i-4] = s[i];
        *dd = (int)atoi(tmp);
      }
      else st = -1;
    }
  }
  else if ((cnt == 1) || (cnt == 2))
  {
    for(i=0; i<p1; i++) tmp[i] = s[i]; tmp[i] = '\0';
    *yy = (int)atoi(tmp);
    p1++;
    if (cnt == 1) *mm = 0;
    else 
    {
      for (i=p1; i<p2; i++) tmp[i-p1] = s[i]; tmp[i-p1] = '\0';
      *mm = (int)atoi(tmp);
      p1 = p2 + 1;
    }
    for (i=p1; i<=len; i++) tmp[i-p1] = s[i];
    *dd = (int)atoi(tmp);
  }
  else st = -1;
  return st;
}

int conv_date(int *mm, int *dd, int yyyy)
/*
   Note: *mm and *dd are input/output varialbes
   Input julian day number:
        input *mm = 0;
        input *dd = julian day
        output in mm-dd-yyyy format
   Input in mm-dd-yyyy format
        output *mm = 0;
        output *dd = julian day number;
*/
{
  int nm, im;
  int st = -1;
  int ndays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  if ((yyyy%400 == 0) || ((yyyy%4 == 0) && (yyyy%100 != 0))) ndays[1] = 29; 
  if (*mm == 0)
  {
    if (*dd <= 0)
      fprintf(stderr, "Error in input julian date: %d %d\n", *dd, yyyy);
    else
    {
      for (im=0; ((im<12)&&(*dd>0)); im++)
        *dd -= ndays[im];
      if ((im > 12) || ((im == 12) && (*dd > 0)))
        fprintf(stderr, "Error in input julian date: %d %d\n", *dd, yyyy);
      else
      {
        *mm = im;
        *dd += ndays[*mm - 1];
        st = 1;
      }
    }
  }
  else
  {
    if ((*mm <= 0) || (*dd <= 0))
      fprintf(stderr, "Error in input date: %d %d %d\n", *mm, *dd, yyyy);
    else
    {
      nm = *mm - 1;
      for (im=0; im<nm; im++)
        *dd += ndays[im];
      *mm = 0;
      st = 1;
    }
  }
  return st;
}

