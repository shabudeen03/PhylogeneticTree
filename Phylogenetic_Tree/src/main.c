#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "debug.h"

int main(int argc, char **argv) {
    if(validargs(argc, argv)) {
        USAGE(*argv, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    if(global_options == HELP_OPTION) {
        USAGE(*argv, EXIT_SUCCESS);
    } else {
        int success = read_distance_data(stdin);
        if(success == 0) {
            FILE* fp;
            if(global_options == 0) {
                fp = stdout;
            }

            success = build_taxonomy(fp);
            if(success != 0) { //Some error in build_taxanomy
                return EXIT_FAILURE;
            }

            if(global_options != 0) {
                fp = stdout;
            }
            if(global_options == MATRIX_OPTION) {
                success = emit_distance_matrix(fp);
            } else if(global_options == NEWICK_OPTION) {
                success = emit_newick_format(fp);
            }

            if(success != 0) {
                return EXIT_FAILURE;
            }

            if(fp != NULL) {
                fclose(fp);
            }
        } else {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
