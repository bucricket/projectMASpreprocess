#include <stdio.h>
#include <math.h>
#define SCALE 0.01

int main(int argc, char *argv[])
{
  short int d;
  float f;
  FILE *in, *out;

  if(argc!=3) {
    printf("Usage: %s <in_int> <out_float>\n", argv[0]);
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
    fread(&d, sizeof(short int), 1, in);
    /* convert DN (C) to K */
    f=d*SCALE+273.15;
    fwrite(&f, sizeof(float), 1, out);
  }

  fclose(in);
  fclose(out);
}
