#include <stdlib.h>

#include "global.h"
#include "debug.h"

static int isOption(char* argument, char* option);

//This function basically checks if argument matches the defined options which are h, m, n
static int isOption(char* argument, char* option) {
    while(*argument != '\0' && *option != '\0') {
        if(*argument == *option) {
            argument++;
            option++;
        } else {
            return 0;
        }
    }

    return 1;
}

/**
 * @brief Validates command line arguments passed to the program.
 * @details This function will validate all the arguments passed to the
 * program, returning 0 if validation succeeds and -1 if validation fails.
 * Upon successful return, the various options that were specified will be
 * encoded in the global variable 'global_options', where it will be
 * accessible elsewhere in the program.  For details of the required
 * encoding, see the assignment handout.
 *
 * @param argc The number of arguments passed to the program from the CLI.
 * @param argv The argument strings passed to the program from the CLI.
 * @return 0 if validation succeeds and -1 if validation fails.
 * @modifies global variable "global_options" to contain an encoded representation
 * of the selected program options.
 */
int validargs(int argc, char **argv)
{
    // TO BE IMPLEMENTED.
    argv++;

    if(argc < 2) {
        return 0;
    }

    char* option = "-h\n";
    int isH = isOption(*argv, option);
    if(isH) {
        global_options = 0x1;
        return 0;
    }

    option = "-m\n";
    int isM = isOption(*argv, option);
    if(isM && argc == 2) {
        global_options = 0x4;
        return 0;
    }

    option = "-n\n";
    int isN = isOption(*argv, option);
    if(isN) {
        if(argc > 2) {
            argv++;
            option = "-o\n";

            int hasO = isOption(*argv, option);
            if(hasO) {
                if(argc == 4) {
                    outlier_name = *(++argv);
                    global_options = 0x2;
                    return 0;
                }
            }
        } else {
            if(argc == 2) {
                global_options = 0x2;
                return 0;
            }
        }
    }

    return -1;
}
