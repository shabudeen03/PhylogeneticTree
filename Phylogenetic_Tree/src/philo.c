#include <stdlib.h>

#include "global.h"
#include "debug.h"

//Helper function prototypes
static double decimalPortion(char* string);
static double stringToDouble(char* string);
static int namesMatch(char* name1, char* name2);
static int findParent(int root);
static int findChild(int root, int isLeft);
static void data(int a, int b, FILE* out);
static void newick(int root, int prevRoot, FILE* out);
static char* integerToString(char* str, int val);


/**
 * Converts decimal portion of string to double value
 */
static double decimalPortion(char* string) {
    double val = 0.0;
    double factor = 0.1;
    while(*string != '\0') {
        if(*string < '0' || *string > '9') {
            fprintf(stderr, "%s\n", "Invalid characters.");
            return -1; //Error, values in distance matrix not valid
        }

        val += ((double) (*string++ - '0')) * factor;
        factor /= 10.0;
    }

    return val;
}

/**
 * Converts a string to double
 * Returns double value
 */
static double stringToDouble(char* string) {
    double val = 0;
    char* ptr = string;

    if(*string == '.') { //Handling decimal starting with .
        string++;
        val = decimalPortion(string);
    } else if(*string == '0') { //Handling decimal starting with 0.
        string++;

        if(*string == '\0') {
            return val;
        } else if(*string == '.') {
            string++;
        } else if(*string == '0') {
            fprintf(stderr, "%s\n", "Only one leading zero allowed!");
            return -1; //Error, only one leading zero allowed
        } else if(*string > '0' && *string <= '9') {
            fprintf(stderr, "%s\n", "Leading 0 must be in front of decimal point only!");
            return -1; //Also Error, leading 0 only in front of decimal point
        }

        val = decimalPortion(string);
    } else if(*string > '0' && *string <= '9') {
        int factor = 1;
        string++;

        while(*string != '\0') {
            if(*string == '.') {
                string++;
                double dec = decimalPortion(string);
                if(dec < 0) {
                    return -1;
                }
                val += dec;
                break;
            }

            factor *= 10;
            string++;
        }

        while(*ptr != '\0' && *ptr != '.') {
            val += factor * (*ptr++ - '0');
            factor /= 10;
        }
    } else {
        fprintf(stderr, "%s\n", "Invalid characters in distance matrix!");
        return -1; //Error, invalid values in distance matrix
    }

    if(val < 0) {
        fprintf(stderr, "%s\n", "Error parsing decimal portion of distance matrix");
        return -1;
    }

    return val;
}

//See if name1 equals name2 contents
static int namesMatch(char* name1, char* name2) {
    char* ptr1 = name1;
    char* ptr2 = name2;
    while(*ptr1 != '\0' && *ptr2 != '\0') {
        if(*ptr1++ != *ptr2++) {
            return 0;
        }
    }

    return 1;
}

/**
 * @brief  Read genetic distance data and initialize data structures.
 * @details  This function reads genetic distance data from a specified
 * input stream, parses and validates it, and initializes internal data
 * structures.
 *
 * The input format is a simplified version of Comma Separated Values
 * (CSV).  Each line consists of text characters, terminated by a newline.
 * Lines that start with '#' are considered comments and are ignored.
 * Each non-comment line consists of a nonempty sequence of data fields;
 * each field iis terminated either by ',' or else newline for the last
 * field on a line.  The constant INPUT_MAX specifies the maximum number
 * of data characters that may be in an input field; fields with more than
 * that many characters are regarded as invalid input and cause an error
 * return.  The first field of the first data line is empty;
 * the subsequent fields on that line specify names of "taxa", which comprise
 * the leaf nodes of a phylogenetic tree.  The total number N of taxa is
 * equal to the number of fields on the first data line, minus one (for the
 * blank first field).  Following the first data line are N additional lines.
 * Each of these lines has N+1 fields.  The first field is a taxon name,
 * which must match the name in the corresponding column of the first line.
 * The subsequent fields are numeric fields that specify N "distances"
 * between this taxon and the others.  Any additional lines of input following
 * the last data line are ignored.  The distance data must form a symmetric
 * matrix (i.e. D[i][j] == D[j][i]) with zeroes on the main diagonal
 * (i.e. D[i][i] == 0).
 *
 * If 0 is returned, indicating data successfully read, then upon return
 * the following global variables and data structures have been set:
 *   num_taxa - set to the number N of taxa, determined from the first data line
 *   num_all_nodes - initialized to be equal to num_taxa
 *   num_active_nodes - initialized to be equal to num_taxa
 *   node_names - the first N entries contain the N taxa names, as C strings
 *   distances - initialized to an NxN matrix of distance values, where each
 *     row of the matrix contains the distance data from one of the data lines
 *   nodes - the "name" fields of the first N entries have been initialized
 *     with pointers to the corresponding taxa names stored in the node_names
 *     array.
 *   active_node_map - initialized to the identity mapping on [0..N);
 *     that is, active_node_map[i] == i for 0 <= i < N.
 *
 * @param in  The input stream from which to read the data.
 * @return 0 in case the data was successfully read, otherwise -1
 * if there was any error.  Premature termination of the input data,
 * failure of each line to have the same number of fields, and distance
 * fields that are not in numeric format should cause a one-line error
 * message to be printed to stderr and -1 to be returned.
 */

int read_distance_data(FILE *in) {
    char ch = fgetc(in);
    if(ch == EOF) {
        fprintf(stderr, "%s\n", "Error: File does not contain taxa");
        return -1;
    }

    while(ch == '#') { //Skip comments until not a comment
        while(ch != '\n') { //Skip to last character which is new line character
            ch = fgetc(in);
            if(ch == EOF) {
                fprintf(stderr, "%s\n", "Error: File does not contain taxa");
                return -1; //IO error or no taxon provided
            }
        }

        ch = fgetc(in); //Skip to next line
        if(ch == EOF) {
            fprintf(stderr, "%s\n", "Error: File does not contain taxa");
            return -1; //IO error or no taxon provided
        }
    }

    //Now, it must be the first data line
    ch = fgetc(in); //skip first ,
    while(ch != '\n') {
        if(ch == '#') {
            fprintf(stderr, "%s\n", "Error: Taxon names can not start with #");
            return -1;
        }

        if(ch == EOF) {
            fprintf(stderr, "%s\n", "Error: File ended without providing distances");
            return -1;
        }

        if(num_taxa > MAX_TAXA) {
            fprintf(stderr, "%s %d%s %d\n", "Max # of taxa allowed:", MAX_TAXA, "; # of Taxa provided:", num_taxa);
            return -1; //Error
        }

        char* ptr = input_buffer;
        int numChar = 0;
        while(ch != ',' && ch != '\n') {
            if(ch == '#' && numChar == 0) {
                fprintf(stderr, "%s\n", "Error: Taxon names can not start with #");
                return -1;
            }

            if(numChar > INPUT_MAX) {
                fprintf(stderr, "%s %d %s\n", "The length of the taxon name must be less than", INPUT_MAX, "characters long!");
                return -1; //error
            }

            if(ch == EOF) {
                fprintf(stderr, "%s\n", "Error: File does not contain distance data");
                return -1;
            }

            *ptr++ = ch;
            numChar++;
            ch = fgetc(in);
        }

        if(numChar == 0) {
            fprintf(stderr, "%s\n", "Error: Taxon name can not be empty");
            return -1;
        }

        *ptr = '\0';
        ptr = input_buffer;

        //Check for if this name repeats in data line
        for(int i=0; i<num_taxa; i++) {
            if(namesMatch(ptr, *(node_names + i))) {
                fprintf(stderr, "%s\n", "Taxon names must be unique.");
                return -1;
            }
        }

        char* ptr2 = *(node_names + num_taxa);
        while(*ptr != '\0') {
            *ptr2++ = *ptr++;
        }

        *ptr2 = '\0';
        if(ch == ',') ch = fgetc(in); //skip comma
        num_taxa++;
    }

    if(ch == EOF) {
        fprintf(stderr, "%s\n", "Error: Incomplete Distances data");
        return -1;
    }

    //Now, rest of lines are either data lines or comments, so just get the N data lines
    ch = fgetc(in); //Skip newline char of first data line

    if(ch == EOF) {
        fprintf(stderr, "%s\n", "Error: Incomplete Distances data");
        return -1;
    }

    int lineCtr = 0; //# of data lines read, <= N
    while(lineCtr < num_taxa) { //Either line is comment or data line
        if(ch == EOF) {
            fprintf(stderr, "%s\n", "Error: Incomplete Distances data");
            return -1;
        }

        if(ch == '#') {
            while(ch != '\n') { //skip to newline char
                if(ch == EOF) {
                    fprintf(stderr, "%s\n", "Error: Incomplete Distances data");
                    return -1;
                }

                ch = fgetc(in);
            }

            ch = fgetc(in); //skip newline char to next line beginning
        } else {
            (*(nodes + lineCtr)).name = *(node_names + lineCtr); //initialize nodes[i].name
            *(active_node_map + lineCtr) = lineCtr; //active_node_map[i] = i for 0 <= i < N
            int entryCtr = 0; //Column # of D matrix, Linectr --> row #
            char* ptr = *(node_names + lineCtr);

            while(ch != ',') {
                if(ch == EOF) {
                    fprintf(stderr, "%s\n", "Error: Incomplete Distances data");
                    return -1;
                }

                if(*ptr++ != ch) {
                    fprintf(stderr, "%s\n", "The Node names do NOT MATCH from the names provided in the first data line!");
                    return -1; //Node names do not match
                }

                ch = fgetc(in); //skip name
            }

            ch = fgetc(in); //skip comma, now pointing to first entry
            while(entryCtr < num_taxa) {
                if(ch == '\n') {
                    fprintf(stderr, "%s\n", "There are not enough data lines to fill up the Distance matrix!");
                    return -1; //Matrix not fully filled up
                }

                ptr = input_buffer;
                int numValues = 0;
                while(ch != ',' && ch != '\n' && ch != EOF) {
                    if(numValues > INPUT_MAX) {
                        fprintf(stderr, "%s %d %s\n", "The values of distance matrix must be less than", INPUT_MAX, "characters long!");
                        return -1; //Error
                    }

                    *ptr++ = ch;
                    numValues++;
                    ch = fgetc(in);
                }

                *ptr = '\0';

                //Need to convert input_buffer to decimal value
                double dist = stringToDouble(input_buffer);
                if(dist < 0) {
                    return -1; //Invalid distance matrix entries
                }

                //Check for 0s on main diagonal
                if(entryCtr == lineCtr && dist > 0) {
                    fprintf(stderr, "%s\n", "The Distances Matrix has non-zero values on main diagonal!");
                    return -1;
                }

                //Check for symmetric matrix
                if(entryCtr < lineCtr && dist != *(*(distances + entryCtr) + lineCtr)) {
                    fprintf(stderr, "%s\n", "The Distances Matrix is not symmetric!");
                    return -1;
                }

                *(*(distances + lineCtr) + entryCtr) = dist;

                entryCtr++;

                if(ch == ',') {
                    ch = fgetc(in); //Skip ,
                }
            }

            lineCtr++;
            ch = fgetc(in);
        }
    }

    if(num_taxa == 0) {
        fprintf(stderr, "%s\n", "No taxons provided in input");
        return -1;
    }

    num_all_nodes = num_taxa;
    num_active_nodes = num_taxa;
    return 0;
}

//Return index of the parent of the root
static int findParent(int root) {
    NODE* parent = *((*(nodes + root)).neighbors);
    char* name = parent->name;
    for(int i=0; i<num_all_nodes; i++) {
        if(namesMatch(*(node_names + i), name)) {
            return i;
        }
    }

    return -1;
}

//Return neighbors[1] if isLeft = 0 and neighbors[2] if isLeft = 1
static int findChild(int root, int isLeft) {
    NODE* child = *((*(nodes + root)).neighbors + isLeft + 1);
    if(child == NULL) {
        return -1;
    }

    char* name = child->name;
    for(int i=0; i<num_all_nodes; i++) {
        if(namesMatch(*(node_names + i), name)) {
            return i;
        }
    }

    return -1;
}

//Output "name:distance"
static void data(int a, int b, FILE* out) {
    fprintf(out, "%s:%.2f", *(node_names + a), *(*(distances + a) + b));
}

//Recursively done to emit newick form of the provided tree
static void newick(int root, int prevRoot, FILE* out) {
    int parent = findParent(root);
    int left = findChild(root, 0);
    int right = findChild(root, 1);

    int a, b;
    if(parent == prevRoot) {
        a = left;
        b = right;
    } else if(left == prevRoot) {
        a = parent;
        b = right;
    } else {
        a = parent;
        b = left;
    }


    if(a == -1 && b == -1) {
        data(root, prevRoot, out);
        return;
    } else { //if(a >= 0 && b >= 0)
        fprintf(out, "%c", '(');
        newick(a, root, out);
        fprintf(out, "%c", ',');
        newick(b, root, out);
        fprintf(out, "%c", ')');
        data(root, prevRoot, out);
    }
}

/**
 * @brief  Emit a representation of the phylogenetic tree in Newick
 * format to a specified output stream.
 * @details  This function emits a representation in Newick format
 * of a synthesized phylogenetic tree to a specified output stream.
 * See (https://en.wikipedia.org/wiki/Newick_format) for a description
 * of Newick format.  The tree that is output will include for each
 * node the name of that node and the edge distance from that node
 * its parent.  Note that Newick format basically is only applicable
 * to rooted trees, whereas the trees constructed by the neighbor
 * joining method are unrooted.  In order to turn an unrooted tree
 * into a rooted one, a root will be identified according by the
 * following method: one of the original leaf nodes will be designated
 * as the "outlier" and the unique node adjacent to the outlier
 * will serve as the root of the tree.  Then for any other two nodes
 * adjacent in the tree, the node closer to the root will be regarded
 * as the "parent" and the node farther from the root as a "child".
 * The outlier node itself will not be included as part of the rooted
 * tree that is output.  The node to be used as the outlier will be
 * determined as follows:  If the global variable "outlier_name" is
 * non-NULL, then the leaf node having that name will be used as
 * the outlier.  If the value of "outlier_name" is NULL, then the
 * leaf node having the greatest total distance to the other leaves
 * will be used as the outlier.
 *
 * @param out  Stream to which to output a rooted tree represented in
 * Newick format.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.  If the global variable "outlier_name" is
 * non-NULL, then it is an error if no leaf node with that name exists
 * in the tree.
 */
int emit_newick_format(FILE *out) {
    if(out == NULL) {
        fprintf(stderr, "%s\n", "Output stream null in emit_newick_format");
        return -1;
    }

    if(num_all_nodes < 2) { //print nothing if less than 2 nodes
        return 0;
    } else if(num_all_nodes == 2) { //print the only edge if only 2 nodes
        fprintf(out, "%s:%.2f;\n", *node_names, *(*distances + 1));
        return 0;
    }

    //Otherwise, print the tree if more than 2 nodes
    int root = -1;
    if(outlier_name == NULL) { //default strategry
        double maxD = 0.0;
        double d;
        for(int i=0; i<num_taxa; i++) {
            d = 0.0;
            for(int j=0; j<num_taxa; j++) {
                d += *(*(distances + i) + j);
            }

            if(d > maxD) {
                root = i;
                maxD = d;
            }
        }
    } else {
        for(int i=0; i<num_taxa; i++) {
            if(namesMatch(*(node_names + i), outlier_name) == 1) {
                root = i;
                break;
            }
        }

        if(root < 0) { //Error if root has not been set to an index
            fprintf(stderr, "%s\n", "Outlier name provided does not exist.");
            return -1;
        }
    }

    int prevRoot = root;
    root = findParent(prevRoot);
    newick(root, prevRoot, out); //Modify this
    fprintf(out, "%s", ";\n");
    return 0;
}

/**
 * @brief  Emit the synthesized distance matrix as CSV.
 * @details  This function emits to a specified output stream a representation
 * of the synthesized distance matrix resulting from the neighbor joining
 * algorithm.  The output is in the same CSV form as the program input.
 * The number of rows and columns of the matrix is equal to the value
 * of num_all_nodes at the end of execution of the algorithm.
 * The submatrix that consists of the first num_leaves rows and columns
 * is identical to the matrix given as input.  The remaining rows and columns
 * contain estimated distances to internal nodes that were synthesized during
 * the execution of the algorithm.
 *
 * @param out  Stream to which to output a CSV representation of the
 * synthesized distance matrix.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */
int emit_distance_matrix(FILE *out) {
    if(out == NULL) {
        fprintf(stderr, "%s\n", "Output stream is null for emit_distance_matrix.");
        return -1;
    }

    int i=0;
    fprintf(out, "%c", ',');
    for(i=0; i<num_all_nodes-1; i++) {
        fprintf(out, "%s,", *(node_names + i));
    }

    fprintf(out, "%s\n", *(node_names + i));
    for(i=0; i<num_all_nodes; i++) {
        fprintf(out, "%s,", *(node_names + i));
        int j;
        for(j=0; j<num_all_nodes-1; j++) {
            fprintf(out, "%.2f,", *(*(distances + i) + j));
        }

        fprintf(out, "%.2f\n", *(*(distances + i) + j));
    }

    return 0;
}

char* integerToString(char* str, int val) {
    char* ptr = str;
    *ptr++ = '#';
    int dummy = val;

    while(dummy > 9) {
        *ptr++ = '0' + (dummy % 10);
        dummy /= 10;
    }

    *ptr++ = '0' + dummy;
    *ptr = '\0';

    //If val >= 10
    if(val >= 10) {
        char* ptr2 = ptr - 1; //Make it point just to the last digit before null terminator
        char* ptr3 = str + 1; //Make it point after the first character which is '#'
        while(ptr3 < ptr2) {
            char ch = *ptr2;
            *ptr2 = *ptr3;
            *ptr3 = ch;
            ptr2--;
            ptr3++;
        }
    }

    return str;
}

/**
 * @brief  Build a phylogenetic tree using the distance data read by
 * a prior successful invocation of read_distance_data().
 * @details  This function assumes that global variables and data
 * structures have been initialized by a prior successful call to
 * read_distance_data(), in accordance with the specification for
 * that function.  The "neighbor joining" method is used to reconstruct
 * phylogenetic tree from the distance data.  The resulting tree is
 * an unrooted binary tree having the N taxa from the original input
 * as its leaf nodes, and if (N > 2) having in addition N-2 synthesized
 * internal nodes, each of which is adjacent to exactly three other
 * nodes (leaf or internal) in the tree.  As each internal node is
 * synthesized, information about the edges connecting it to other
 * nodes is output.  Each line of output describes one edge and
 * consists of three comma-separated fields.  The first two fields
 * give the names of the nodes that are connected by the edge.
 * The third field gives the distance that has been estimated for
 * this edge by the neighbor-joining method.  After N-2 internal
 * nodes have been synthesized and 2*(N-2) corresponding edges have
 * been output, one final edge is output that connects the two
 * internal nodes that still have only two neighbors at the end of
 * the algorithm.  In the degenerate case of N=1 leaf, the tree
 * consists of a single leaf node and no edges are output.  In the
 * case of N=2 leaves, then no internal nodes are synthesized and
 * just one edge is output that connects the two leaves.
 *
 * Besides emitting edge data (unless it has been suppressed),
 * as the tree is built a representation of it is constructed using
 * the NODE structures in the nodes array.  By the time this function
 * returns, the "neighbors" array for each node will have been
 * initialized with pointers to the NODE structure(s) for each of
 * its adjacent nodes.  Entries with indices less than N correspond
 * to leaf nodes and for these only the neighbors[0] entry will be
 * non-NULL.  Entries with indices greater than or equal to N
 * correspond to internal nodes and each of these will have non-NULL
 * pointers in all three entries of its neighbors array.
 * In addition, the "name" field each NODE structure will contain a
 * pointer to the name of that node (which is stored in the corresponding
 * entry of the node_names array).
 *
 * @param out  If non-NULL, an output stream to which to emit the edge data.
 * If NULL, then no edge data is output.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */
int build_taxonomy(FILE *out) {
    FILE* fp = out;

    //while there are more than 2 active nodes
    while(num_active_nodes > 2) {
        //figure out vector s(i) w/ active nodes only
        int i=0, j, r, c;
        double sum;
        for(i=0; i<num_active_nodes; i++) {
            sum = 0.0;
            r = *(active_node_map + i);
            for(j=0; j<num_active_nodes; j++) {
                c = *(active_node_map + j);
                sum += *(*(distances + r) + c);
            }

            *(row_sums + r) = sum;
        }


        //Q(i, j) where I have to determine f and g s.t. Q(f, g) is least
        int f, g;
        int fIdx, gIdx; //Keep track of where f and g are stored
        double val, leastVal;
        int firstRun = 0; //First time leastVal is just that Q(i, j), next time compare for lower value
        for(i=0; i<num_active_nodes; i++) {
            r = *(active_node_map + i);
            for(j=0; j<num_active_nodes; j++) {
                c = *(active_node_map + j);
                if(r != c) {
                    double dVal = *(*(distances + r) + c);
                    double si = *(row_sums + r);
                    double sj = *(row_sums + c);
                    val = ((num_active_nodes - 2) * dVal) - si - sj;

                    if(firstRun == 0 || val < leastVal) {
                        f = r;
                        g = c;
                        leastVal = val;
                        fIdx = i;
                        gIdx = j;

                        if(firstRun == 0) {
                            firstRun++;
                        }
                    }
                }
            }
        }

        //New node u, connect w/ f and g
        int u = num_all_nodes;
        char* str = *(node_names + u);
        str = integerToString(str, u);
        (*(nodes + u)).name = str;
        *((*(nodes + f)).neighbors) = (nodes + u);
        *((*(nodes + g)).neighbors) = (nodes + u);
        *((*(nodes + u)).neighbors + 1) = (nodes + f);
        *((*(nodes + u)).neighbors + 2) = (nodes + g);

        //Extend matrix D w/ non actives nodes being 0.0 distance from new node u
        for(i=0; i<u; i++) {
            *(*(distances + i) + u) = 0.0;
            *(*(distances + u) + i) = 0.0;
        }

        int k;
        double fg = *(*(distances + f) + g);
        for(i=0; i<num_active_nodes; i++) {
            k = *(active_node_map + i);
            if(k == f) {
                double sf = *(row_sums + f);
                double sg = *(row_sums + g);
                val = (fg + ((sf - sg) / (num_active_nodes - 2))) / 2.0;

                if(fp != NULL) { //Default operation - emit edge data
                    if(f < u) {
                        fprintf(fp, "%d,%d,%.2f\n", f, u, val);
                    } else {
                        fprintf(fp, "%d,%d,%.2f\n", u, f, val);
                    }
                }
            } else if(k == g) {
                double sf = *(row_sums + f);
                double sg = *(row_sums + g);
                val = (fg + ((sg - sf) / (num_active_nodes - 2))) / 2.0;

                if(fp != NULL) { //Default operation - emit edge data
                    if(g < u) {
                        fprintf(fp, "%d,%d,%.2f\n", g, u, val);
                    } else {
                        fprintf(fp, "%d,%d,%.2f\n", u, g, val);
                    }
                }
            } else {
                double fk = *(*(distances + f) + k);
                double gk = *(*(distances + g) + k);
                val = (fk + gk - fg) / 2.0;
            }

            *(*(distances + k) + u) = val;
            *(*(distances + u) + k) = val;
        }

        *(*(distances + u) + u) = 0.0; //the if k = u condition

        //Inactivate nodes f and g, update N
        *(active_node_map + fIdx) = num_all_nodes;
        *(active_node_map + gIdx) = *(active_node_map + num_active_nodes - 1);
        num_all_nodes++;
        num_active_nodes--;

        //Add edge b/w remaining 2 nodes
        if(num_active_nodes == 2) {
            f = *(active_node_map + fIdx);
            g = *(active_node_map + gIdx);
            *((*(nodes + f)).neighbors) = (nodes + g);
            *((*(nodes + g)).neighbors) = (nodes + f);
            if(fp != NULL) {
                if(g < f) {
                    fprintf(fp, "%d,%d,%.2f\n", g, f, *(*(distances + f) + g));
                } else {
                    fprintf(fp, "%d,%d,%.2f\n", f, g, *(*(distances + f) + g));
                }
            }
        }
    }

    return 0;
}