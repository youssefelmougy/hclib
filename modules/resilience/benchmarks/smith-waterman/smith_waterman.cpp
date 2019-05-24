#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include <iostream>
#include <cstring>

//#define DEBUG

#ifdef DEBUG
#include <sstream>
#endif

#define GAP_PENALTY -1
#define TRANSITION_PENALTY -2
#define TRANSVERSION_PENALTY -4
#define MATCH 2

//#define USE_REPLAY
//#define USE_REPLICATION

enum TASK_STATE {NON_LEAF, LEAF};

typedef int int_t;

typedef struct {
    // default values defined below; can be overridden on command line
    int_t tile_width = 32;
    int_t tile_height = 32;
    bool do_checksums = true;
    bool do_injection = false;
    double injection_rate = 0.0;
    std::string input1;
    std::string input2;
    std::string failfile;
} params_t;

enum Nucleotide {GAP=0, ADENINE, CYTOSINE, GUANINE, THYMINE};

signed char char_mapping ( char c ) {
    signed char to_be_returned = -1;
    switch(c) {
        case '_': to_be_returned = GAP; break;
        case 'A': to_be_returned = ADENINE; break;
        case 'C': to_be_returned = CYTOSINE; break;
        case 'G': to_be_returned = GUANINE; break;
        case 'T': to_be_returned = THYMINE; break;
    }
    return to_be_returned;
}

void print_matrix ( int** matrix, int n_rows, int n_columns ) {
    int i, j;
    for ( i = 0; i < n_rows; ++i ) {
        for ( j = 0; j < n_columns; ++j ) {
            fprintf(stdout, "%d ", matrix[i][j]);
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout,"--------------------------------\n");
}

static char alignment_score_matrix[5][5] =
{
    {GAP_PENALTY,GAP_PENALTY,GAP_PENALTY,GAP_PENALTY,GAP_PENALTY},
    {GAP_PENALTY,MATCH,TRANSVERSION_PENALTY,TRANSITION_PENALTY,TRANSVERSION_PENALTY},
    {GAP_PENALTY,TRANSVERSION_PENALTY, MATCH,TRANSVERSION_PENALTY,TRANSITION_PENALTY},
    {GAP_PENALTY,TRANSITION_PENALTY,TRANSVERSION_PENALTY, MATCH,TRANSVERSION_PENALTY},
    {GAP_PENALTY,TRANSVERSION_PENALTY,TRANSITION_PENALTY,TRANSVERSION_PENALTY, MATCH}
};

size_t clear_whitespaces_do_mapping ( signed char* buffer, long lsize ) {
    size_t non_ws_index = 0, traverse_index = 0;

    while ( traverse_index < lsize ) {
        char curr_char = buffer[traverse_index];
        switch ( curr_char ) {
            case 'A': case 'C': case 'G': case 'T':
                /*this used to be a copy not also does mapping*/
                buffer[non_ws_index++] = char_mapping(curr_char);
                break;
        }
        ++traverse_index;
    }
    return non_ws_index;
}

signed char* read_file( FILE* file, size_t* n_chars ) {
    fseek (file, 0L, SEEK_END);
    long file_size = ftell (file);
    fseek (file, 0L, SEEK_SET);

    signed char *file_buffer = (signed char *)malloc((1+file_size)*sizeof(signed char));

    size_t n_read_from_file = fread(file_buffer, sizeof(signed char), file_size, file);
    file_buffer[file_size] = '\n';

    /* shams' sample inputs have newlines in them */
    *n_chars = clear_whitespaces_do_mapping(file_buffer, file_size);
    return file_buffer;
}

class IntVec: public hclib::resilience::obj,
          public std::vector<int> {
public:
  virtual ~IntVec() {}

  bool equals(obj* obj2){
    IntVec *o2 = (IntVec*)obj2;
    return *this == *o2;
  }
#ifdef DEBUG
  void atomic_print(unsigned long long task_idx){
    std::stringstream ss;
    ss << task_idx << ": ";
    for (auto it=cbegin(); it!=cend(); ++it){
      ss << *it << " ";
    }
    ss << std::endl;
    std::cout << ss.str();
  }
#endif


};

typedef struct {
#if defined (USE_REPLAY)
    hclib::resilience::replay::promise_t<IntVec*>* bottom_row;
    hclib::resilience::replay::promise_t<IntVec*>* right_column;
    hclib::resilience::replay::promise_t<IntVec*>* bottom_right;
#elif defined (USE_REPLICATION)
    hclib::resilience::diamond::promise_t<IntVec*>* bottom_row;
    hclib::resilience::diamond::promise_t<IntVec*>* right_column;
    hclib::resilience::diamond::promise_t<IntVec*>* bottom_right;
#else
    hclib::ref_count::promise_t<IntVec*>* bottom_row;
    hclib::ref_count::promise_t<IntVec*>* right_column;
    hclib::ref_count::promise_t<IntVec*>* bottom_right;
#endif
    int id;
} Tile_t;

#if defined (USE_REPLAY)
typedef struct {
    int x;
    int y;
} Point;

int check(void *args)
{
  bool succeeded = *(bool*)args;
  return (succeeded ? 1 : 0);
}
#endif


void get_params(int argc, char *argv[], params_t &params)
{
    params_t defaults(params);
    for (int i=1; i < argc; ++i){
        if (strcmp(argv[i], "-tile_height") == 0){
            if (i+1 == argc){
                std::cerr << "Error: didn't specify value for ntiles\n";
                exit(EXIT_FAILURE);
            }
            params.tile_height = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-tile_width") == 0){
            if (i+1 == argc){
                std::cerr << "Error: didn't specify value for tilesize\n";
                exit(EXIT_FAILURE);
            }
            params.tile_width = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-checksum") == 0){
            params.do_checksums = true;
        } else if (strcmp(argv[i], "-nochecksum") == 0){
            params.do_checksums = false;
        } else if (strcmp(argv[i], "-inject") == 0){
            if (i+1 == argc){
                std::cerr << "Error: didn't specify failure file for injection\n";
                exit(EXIT_FAILURE);
            }
            params.do_injection = true;
            params.failfile = std::string(argv[++i]);
        } else if (strcmp(argv[i], "-inject_rate") == 0){
            if (i+1 == argc){
                std::cerr << "Error: didn't specify failure rate \n";
                exit(EXIT_FAILURE);
            }
            params.injection_rate = atof(argv[++i]);
        } else if (strcmp(argv[i], "-input1") == 0) {
            if (i+1 == argc){
                std::cerr << "Error: didn't specify input file 1\n";
                exit(EXIT_FAILURE);
            }
            params.input1 = std::string(argv[++i]);
        } else if (strcmp(argv[i], "-input2") == 0) {
            if (i+1 == argc){
                std::cerr << "Error: didn't specify input file 2\n";
                exit(EXIT_FAILURE);
            }
            params.input2 = std::string(argv[++i]);
        } else if ((strcmp(argv[i], "-help") == 0)
                   || (strcmp(argv[i], "-h") == 0)){
            std::cout << "Options:\n"
                      << "  -tile_width <int>: Width per tile"
                      << " (default " << defaults.tile_width << ")\n"
                      << "  -tile_height <int>: Height per tile"
                      << " (default " << defaults.tile_height << ")\n"
                      << "  -[no]checksum: Whether or not to use checksums"
                      << " (default "
                      << (defaults.do_checksums ? "checksum" : "nochecksum") << ")\n"
                      << "  -inject <filename>: Failure list for failure injection"
                      << " (default disabled)\n"
                      << "  -inject_rate: the probability of failures"
                      << "  -input1 <filename>: Input file 1\n"
                      << "  -input2 <filename>: Input file 2\n";
            exit(EXIT_SUCCESS);
        } else {
            std::cerr << "Unexpected option: " << argv[i] << "; use -h for syntax\n";
            exit(EXIT_FAILURE);
        }
    }
}

bool should_fail(double rate, int_t tileID, int_t nTiles) {
    //    fprintf(stderr, "rate: %lf, ID: %d, ntiles: %d\n", rate, tileID, nTiles);
    if (rate != 0.0) {
        return (tileID % (int)ceil(nTiles * 1.0 / (nTiles * rate))) == 0;
    } else {
        return false;
    }
}

int main ( int argc, char* argv[] ) {
    const char *deps[] = { "system" };
    hclib::launch(deps, 1, [&]() {
        int i, j;

        params_t params;
        get_params(argc, argv, params);

        if (params.tile_width < 1) {
            fprintf(stderr, "Fatal error: tile width must be positive\n");
            exit(EXIT_FAILURE);
        }

        if (params.tile_height < 1) {
            fprintf(stderr, "Fatal error: tile height must be positive\n");
            exit(EXIT_FAILURE);
        }

        if (params.input1.empty() || params.input2.empty()) {
            fprintf(stderr, "Fatal error: input file(s) is not specified\n");
            exit(EXIT_FAILURE);
        }

        if (params.injection_rate > 1) {
            fprintf(stderr, "Fatal error: inject rate should be less than 1");
            exit(EXIT_FAILURE);
        }

        signed char* string_1;
        signed char* string_2;

        FILE* file_1 = fopen(params.input1.c_str(), "r");
        if (!file_1) {
            fprintf(stderr, "could not open file %s\n", params.input1);
            exit(EXIT_FAILURE);
        }
        size_t n_char_in_file_1 = 0;
        string_1 = read_file(file_1, &n_char_in_file_1);
        fprintf(stdout, "Size of input string 1 is %lu\n", n_char_in_file_1 );

        FILE* file_2 = fopen(params.input2.c_str(), "r");
        if (!file_2) {
            fprintf(stderr, "could not open file %s\n", params.input2);
            exit(EXIT_FAILURE);
        }
        size_t n_char_in_file_2 = 0;
        string_2 = read_file(file_2, &n_char_in_file_2);
        fprintf(stdout, "Size of input string 2 is %lu\n", n_char_in_file_2 );

        int num_workers = hclib::get_num_workers();
        int n_tiles_width = n_char_in_file_1 / params.tile_width;
        int n_tiles_height = n_char_in_file_2 / params.tile_height;

        fprintf(stdout, "num_workers=%d\n", num_workers);
        fprintf(stdout, "n_tiles=%dx%d\n", n_tiles_width, n_tiles_height);
        fprintf(stdout, "totaltiles=%d\n", n_tiles_width*n_tiles_height);
        fprintf(stdout, "tilesize=%dx%d\n", params.tile_width, params.tile_height);
        fprintf(stdout, "checksums=%s\n", (params.do_checksums ? "yes" : "no"));
        fprintf(stdout, "injection=%s\n", (params.do_injection ? "yes" : "no"));
        fprintf(stdout, "injection_rate=%lf\n", params.injection_rate);
#if defined(USE_REPLAY)
        fprintf(stdout, "resilience=replay\n");
#elif defined(USE_REPLICATION)
        fprintf(stdout, "resilience=replication\n");
#else
        fprintf(stdout, "resilience=none\n");
#endif

        fprintf(stdout, "Allocating tile matrix\n");

        int_t tile_width = params.tile_width;
        int_t tile_height = params.tile_height;

        int_t id = 0;
        Tile_t** tile_matrix = (Tile_t **) malloc(sizeof(Tile_t*)*(n_tiles_height+1));
        for ( i = 0; i < n_tiles_height+1; ++i ) {
            tile_matrix[i] = (Tile_t *) malloc(sizeof(Tile_t)*(n_tiles_width+1));
            for ( j = 0; j < n_tiles_width+1; ++j ) {
#if defined (USE_REPLAY)
                tile_matrix[i][j].bottom_row = new hclib::resilience::replay::promise_t<IntVec*>(3);
                tile_matrix[i][j].right_column = new hclib::resilience::replay::promise_t<IntVec*>(3);
                tile_matrix[i][j].bottom_right = new hclib::resilience::replay::promise_t<IntVec*>(3);
#elif defined (USE_REPLICATION)
                tile_matrix[i][j].bottom_row = new hclib::resilience::diamond::promise_t<IntVec*>(3);
                tile_matrix[i][j].right_column = new hclib::resilience::diamond::promise_t<IntVec*>(3);
                tile_matrix[i][j].bottom_right = new hclib::resilience::diamond::promise_t<IntVec*>(3);
#else
                tile_matrix[i][j].bottom_row = new hclib::ref_count::promise_t<IntVec*>(3);
                tile_matrix[i][j].right_column = new hclib::ref_count::promise_t<IntVec*>(3);
                tile_matrix[i][j].bottom_right = new hclib::ref_count::promise_t<IntVec*>(3);
#endif
                tile_matrix[i][j].id = id++;
            }
        }

        fprintf(stdout, "Allocated tile matrix\n");

        IntVec* allocated = new IntVec;
        assert(allocated);
        allocated->resize(1);
        (*allocated)[0] = 0;
        tile_matrix[0][0].bottom_right->put(allocated);

        for ( j = 1; j < n_tiles_width + 1; ++j ) {
            allocated = new IntVec;
            assert(allocated);
            allocated->resize(tile_width);
            for( i = 0; i < tile_width ; ++i ) {
                (*allocated)[i] = GAP_PENALTY*((j-1)*tile_width+i+1);
            }
            tile_matrix[0][j].bottom_row->put(allocated);

            allocated = new IntVec;
            assert(allocated);
            allocated->resize(1);
            (*allocated)[0] = GAP_PENALTY*(j*tile_width); //sagnak: needed to handle tilesize 2
            tile_matrix[0][j].bottom_right->put(allocated);
        }

        for ( i = 1; i < n_tiles_height + 1; ++i ) {
            allocated = new IntVec;
            assert(allocated);
            allocated->resize(tile_height);
            for ( j = 0; j < tile_height ; ++j ) {
                (*allocated)[j] = GAP_PENALTY*((i-1)*tile_height+j+1);
            }
            tile_matrix[i][0].right_column->put(allocated);

            allocated = new IntVec;
            assert(allocated);
            allocated->resize(1);
            (*allocated)[0] = GAP_PENALTY*(i*tile_height); //sagnak: needed to handle tilesize 2
            tile_matrix[i][0].bottom_right->put(allocated);
        }

        struct timeval begin,end;
        gettimeofday(&begin,0);

        hclib::finish([=]() {
            for (int i = 1; i < n_tiles_height+1; ++i ) {
                for (int j = 1; j < n_tiles_width+1; ++j ) {
                    auto futures = new std::vector<hclib_future_t*>(3);
                    (*futures)[0] = tile_matrix[i][j-1].right_column->get_future();
                    (*futures)[1] = tile_matrix[i-1][j].bottom_row->get_future();
                    (*futures)[2] = tile_matrix[i-1][j-1].bottom_right->get_future();
#if defined (USE_REPLAY)
                    bool *succeeded = new bool;
                    hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();
                    hclib::resilience::replay::async_await_check<LEAF, 2>([=]{
#elif defined (USE_REPLICATION)
                    hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();
                    hclib::resilience::diamond::async_await_check([=]{
#else
                    hclib::ref_count::async_await([=] {
#endif
                        int index, ii, jj;
                        IntVec* above_tile_bottom_row = (IntVec *)tile_matrix[i-1][j  ].bottom_row->get_future()->get();
                        IntVec* left_tile_right_column = (IntVec *)tile_matrix[  i][j-1].right_column->get_future()->get();
                        IntVec* diagonal_tile_bottom_right = (IntVec *)tile_matrix[i-1][j-1].bottom_right->get_future()->get();

                        int  * curr_tile_tmp = (int*) malloc(sizeof(int)*(1+tile_width)*(1+tile_height));
                        int ** curr_tile = (int**) malloc(sizeof(int*)*(1+tile_height));
                        for (index = 0; index < tile_height+1; ++index) {
                            curr_tile[index] = &curr_tile_tmp[index*(1+tile_width)];
                        }

                        curr_tile[0][0] = (*diagonal_tile_bottom_right)[0];
                        for ( index = 1; index < tile_height+1; ++index ) {
                            curr_tile[index][0] = (*left_tile_right_column)[index-1];
                        }

                        for ( index = 1; index < tile_width+1; ++index ) {
                            curr_tile[0][index] = (*above_tile_bottom_row)[index-1];
                        }

#if defined (USE_REPLAY)
                        Point *p_tile = nullptr;
                        if (params.do_checksums) {
                            p_tile = (Point*)calloc((tile_height+1)*(tile_width+1), sizeof(Point));
                            for ( index = 1; index < tile_height+1; ++index ) {
                                p_tile[index * (tile_width+1) + 0].y = index - 1;
                            }
                            for ( index = 1; index < tile_width+1; ++index ) {
                                p_tile[0 * (tile_width+1) + index].x = index - 1;
                            }
                        }
#endif

                        for ( ii = 1; ii < tile_height+1; ++ii ) {
                            for ( jj = 1; jj < tile_width+1; ++jj ) {
                                signed char char_from_1 = string_1[(j-1)*tile_width+(jj-1)];
                                signed char char_from_2 = string_2[(i-1)*tile_height+(ii-1)];

                                int diag_score = curr_tile[ii-1][jj-1] + alignment_score_matrix[char_from_2][char_from_1];
                                int left_score = curr_tile[ii  ][jj-1] + alignment_score_matrix[char_from_1][GAP];
                                int  top_score = curr_tile[ii-1][jj  ] + alignment_score_matrix[GAP][char_from_2];

                                if (left_score > top_score) {
                                    if (left_score > diag_score) {
                                        curr_tile[ii][jj] = left_score;
#if defined (USE_REPLAY)
                                        if (params.do_checksums) {
                                            p_tile[ii * (tile_width+1) + jj].x = p_tile[ii *  (tile_width+1) + jj - 1].x;
                                            p_tile[ii * (tile_width+1) + jj].y = p_tile[ii *  (tile_width+1) + jj - 1].y;
                                        }
#endif
                                    } else {
                                        curr_tile[ii][jj] = diag_score;
#if defined (USE_REPLAY)
                                        if (params.do_checksums) {
                                            p_tile[ii * (tile_width+1) + jj].x = p_tile[(ii-1) * (tile_width+1) + jj - 1].x;
                                            p_tile[ii * (tile_width+1) + jj].y = p_tile[(ii-1) * (tile_width+1) + jj - 1].y;
                                        }
#endif
                                    }
                                } else {
                                    if (top_score > diag_score) {
                                        curr_tile[ii][jj] = top_score;
#if defined (USE_REPLAY)
                                        if (params.do_checksums) {
                                            p_tile[ii * (tile_width+1) + jj].x = p_tile[(ii-1) * (tile_width+1) + jj].x;
                                            p_tile[ii * (tile_width+1) + jj].y = p_tile[(ii-1) * (tile_width+1) + jj].y;
                                        }
#endif
                                    } else {
                                        curr_tile[ii][jj] = diag_score;
#if defined (USE_REPLAY)
                                        if (params.do_checksums) {
                                            p_tile[ii * (tile_width+1) + jj].x = p_tile[(ii-1) * (tile_width+1) + jj - 1].x;
                                            p_tile[ii * (tile_width+1) + jj].y = p_tile[(ii-1) * (tile_width+1) + jj - 1].y;
                                        }
#endif
                                    }
                                }
                            }
                        }

#if defined(USE_REPLAY)
                        *succeeded = true;
                        if (params.do_checksums) {
                            for ( ii = 0; ii < tile_height+1; ++ii ) {
                                for ( jj = 0; jj < tile_width+1; ++jj ) {
                                    if ( p_tile[ii * (tile_width+1) + jj].x != 0
                                         &&  p_tile[ii * (tile_width+1) + jj].y != 0) {
                                        *succeeded = false;
                                    }
                                }
                            }
                            free(p_tile);
                        }
                        if (should_fail(params.injection_rate, tile_matrix[i][j].id, n_tiles_width * n_tiles_height)
                            && hclib::resilience::replay::get_index() == 0) {
                            // fprintf(stderr, "tile %d failed\n", tile_matrix[i][j].id);
                            *succeeded = false;
                        }
#elif defined(USE_REPLICATION)
                        if (should_fail(params.injection_rate, tile_matrix[i][j].id, n_tiles_width * n_tiles_height)
                            && hclib::resilience::diamond::get_index() == 0) {
                            // fprintf(stderr, "tile %d failed\n", tile_matrix[i][j].id);
                            curr_tile[tile_height][tile_width] = -1000;
                        }
#endif

                        IntVec* curr_bottom_right = new IntVec;
                        assert(curr_bottom_right);
                        curr_bottom_right->resize(1);
                        (*curr_bottom_right)[0] = curr_tile[tile_height][tile_width];
                        tile_matrix[i][j].bottom_right->put(curr_bottom_right);

                        IntVec* curr_right_column = new IntVec;
                        assert(curr_right_column);
                        curr_right_column->resize(tile_height);
                        for ( index = 0; index < tile_height; ++index ) {
                            (*curr_right_column)[index] = curr_tile[index+1][tile_width];
                        }
                        tile_matrix[i][j].right_column->put(curr_right_column);

                        IntVec* curr_bottom_row= new IntVec;
                        assert(curr_bottom_row);
                        curr_bottom_row->resize(tile_width);
                        for ( index = 0; index < tile_width; ++index ) {
                            (*curr_bottom_row)[index] = curr_tile[tile_height][index+1];
                        }
                        tile_matrix[i][j].bottom_row->put(curr_bottom_row);

                        free(curr_tile);
                        free(curr_tile_tmp);
                        },
#if defined (USE_REPLAY)
                        prom_res, check, (void*)(succeeded),
#elif defined (USE_REPLICATION)
                        prom_res,
#endif
                        futures);

#if defined(USE_REPLAY) || defined(USE_REPLICATION)
                    hclib::async_await([=]{
                        int res = prom_res->get_future()->get();
                        if (res == 0){
                            fprintf(stderr,"Fatal error: resilience failed to save you\n");
                            exit(EXIT_FAILURE);
                        }
#if defined(USE_REPLAY)
                        delete succeeded;
#endif
                        delete futures;
                        }, prom_res->get_future());
#endif
                }
            }
        }); // finish

        gettimeofday(&end,0);
        fprintf(stdout, "The computation took %f seconds\n",((end.tv_sec - begin.tv_sec)*1000000+(end.tv_usec - begin.tv_usec))*1.0/1000000);
        int score = (*(IntVec *)(tile_matrix[n_tiles_height][n_tiles_width].bottom_row->get_future()->get()))[tile_width-1];
        fprintf(stdout, "score: %d\n", score);

    }); // launch

    return 0;
}
