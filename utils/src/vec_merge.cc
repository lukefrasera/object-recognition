#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "log.h"

#define FREE_CHECK_NULL(x) do {                                                \
  if ((x) != NULL {                                                            \
    free(x);                                                                   \
    x = NULL;                                                                  \
  }                                                                            \
} while(0)                                                                     \

void ReadVecHeader(FILE *file, int *count, int *size) {

}

void ReadVecSample(FILE *file, int size, unsigned char *data) {

}

void TransferSamples(FILE *file_input, FILE *file_output, int count, int size) {

}

void WriteCompbinedHeader(FILE *file_output,
    int count_1, int count_2, int size) {

}
int main(int argc, char *argv[]) {
  int error_flag = false, opt;
  char *filename_a = NULL, *filename_b = NULL, *output = NULL;
  while ((opt = getopt(argc, argv, "a:b:c:")) != -1) {
    switch (opt) {
      case 'a':
        filename_a = strdup(optarg);
        break;
      case 'b':
        filename_b = strdup(optarg);
        break;
      case 'c':
        output     = strdup(optarg);
        break;
      default:
        printf(
          "Usage: %s [-a <Input Vec file>]"
          " [-b <Input Vec file>]"
          " [-c <Output vec file>]\n", argv[0]);
        exit(EXIT_FAILURE);
    }
  }
  if (!filename_a) {
    printf("ERROR: Input a is NULL\n");
    error_flag = true;
  }
  if (!filename_b) {
    printf("ERROR: Input b is NULL\n");
    error_flag = true;
  }
  if (!output) {
    printf("ERROR: Output is NULL\n");
    error_flag = true;
  }
  if (error_flag) {
    FREE_CHECK_NULL(filename_a);
    FREE_CHECK_NULL(filename_b);
    FREE_CHECK_NULL(output);
    return -1;
  }

  printf("Input A: %s, Input B: %s, Output: %s\n",
    filename_a, filename_b, output);

  // Read Headers from the two input vec files
  FILE *file_a      = fopen(filename_a, "rb");
  FILE *file_b      = fopen(filename_b, "rb");
  FILE *output_file = fopen(output, "wb");
  int count_a, count_b, size_a, size_b;
  ReadVecHeader(file_a, &count_a, &size_a);
//  TODO: READ Header B
  LOG_INFO("Header B Doesn't exist yet, %s", filename_b);
  // Check if image size is the same
  if (size_a != size_b) {
    printf("ERROR: Input vec files have different image size\n");
    FREE_CHECK_NULL(filename_a);
    FREE_CHECK_NULL(filename_b);
    FREE_CHECK_NULL(output);
    return -1;
  }
  // Generate Output combined vec file header
  WriteCompbinedHeader(output_file, count_a, count_a, size_a);
  // Copy vec samples from first vec file
  TransferSamples(file_a, output_file, count_a, size_a);
  // copy vec samples from second vec file
  TransferSamples(file_b, output_file, count_b, size_b);

  fclose(file_a);
  fclose(file_b);
  fclose(output_file);

  FREE_CHECK_NULL(filename_a);
  FREE_CHECK_NULL(filename_b);
  FREE_CHECK_NULL(output);
  return 0;
}
