/* convert MODIS LST from Celsius to Kelvin degrees */

#include <stdio.h>
#include <math.h>

#define FILLV -9999

int main(int argc, char *argv[])
{
  float c, k;
  FILE *in, *out;

  if(argc!=3) {
    printf("Usage: %s <in_C> <out_K>\n", argv[0]);
    return -1;
  }
  if((in=fopen(argv[1], "rb"))==NULL) {
    printf("open file %s error\n", argv[1]);
    return -1;
  }
  if((out=fopen(argv[2], "wb"))==NULL) {
    printf("open file %s error\n", argv[2]);
    return -1;
  }

  while(!feof(in)) {
    fread(&c, sizeof(float), 1, in);
    /* convert DN (C) to K */
    if(c == FILLV)
      k = FILLV;
    else
      k=c+273.15;
    fwrite(&k, sizeof(float), 1, out);
  }

  fclose(in);
  fclose(out);
}
