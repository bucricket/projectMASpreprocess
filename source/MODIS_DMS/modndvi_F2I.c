/* convert MODIS NDVI from float to short integer */

#include <stdio.h>
#include <math.h>
#define FILLV -9999

int main(int argc, char *argv[])
{
  float f;
  short int i;
  FILE *in, *out;

  if(argc!=3) {
    printf("Usage: %s <in_float> <out_int>\n", argv[0]);
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
    fread(&f, sizeof(float), 1, in);

    /* convert NDVI from float to integer type */
    if(f == FILLV) i = FILLV;
    else i=(int)(f*10000+0.5);

    fwrite(&i, sizeof(short int), 1, out);
  }

  fclose(in);
  fclose(out);
}
