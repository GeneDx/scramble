#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <unistd.h>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/cram.h"

/*
  ./soft_clipped_clusters file.bam
  ./soft_clipped_clusters -m 10 -s 5 -r chr9:1-20000 file.bam
*/

#define AUTHOR "Kevin Galens"
#define VERSION "0.1.1"
#define PROG "cluster_identifier"

// Holds the settings for this run
typedef struct sc_settings {
    int min_sc_bases;   // Minimum soft clipped bases to be put in a cluster [default: 10]
    int min_sc_reads;   // Minimum soft clipped reads to be a cluster [default: 5]
    char *region;       // Region [default: all]
    char *fn;           // File name
    bool debug;         // Will print debug information to stderr
} sc_settings_t;

// Enum for determine which side of the read
// is clipped.
typedef enum sc_side {
    left,
    right,
} sc_side_t;

// struct for representing a soft clipped read
typedef struct SCRead {
    char *read_name;
    uint32_t *cigar;
    int n_cigar;
    sc_side_t side; // enum of which side is softclipped
    int32_t mapped_pos; // Position the read mapped.
    int32_t sc_pos;  // Position where softclipped
    int32_t tid;     // chromosome
    int32_t l_qseq; // Length of the read
    int32_t sc_len; // Length of soft clipped portion
    uint8_t *seq;
} sc_read_t;

typedef struct SCCluster {
    char *clipped_consensus;
    char *anchored_consensus;
    int n_seqs;
    sc_side_t side;
    int tid;
    int pos;
} sc_seq_t;

// Function declarations.
static void process_region( htsFile *in_fp, hts_idx_t *idx, bam_hdr_t *header, char *region, sc_settings_t );
static bool is_soft_clipped( bam1_t* );
static void populate_cigar_info( bam1_t*, sc_read_t* );
static int sort_by_sc_pos( const void*, const void* );
static char* get_side_str( sc_side_t );
static void destroy_read( sc_read_t* );
static void handle_cluster( sc_read_t**, int, bam_hdr_t* );
static char seq_get_base( int );
static char *get_consensus( char*, int, char**, int, sc_side_t );
    
// These aren't use currently.
static void print_read( sc_read_t* ) __attribute__ ((unused));
static void print_left_align( char**, int ) __attribute__ ((unused));
static void print_right_align( int, char**, int ) __attribute__ ((unused));

static void usage(const char *argv0) {
    fprintf(stderr, "Usage: %s [-m min soft clipped bases] [-s min soft clipped reads] [-r region] cram/bam\n", argv0 );
    exit(EXIT_FAILURE);
}

int main( int argc, char** argv) {
    
    int c;
    int min_sc_bases = 10;
    int min_sc_reads = 5;
    bool debug = false;
    char *region = "all";

    while( (c = getopt(argc, argv, "m:s:r:v")) != -1 ) {
        switch(c) {
            case 'm':
                min_sc_bases = atoi(optarg);
                break;
            case 's':
                min_sc_reads = atoi(optarg);
                break;
            case 'r':
                region = optarg;
                break;
            case 'v':
                fprintf(stderr, "%s: version %s\n", PROG, VERSION);
                exit(0);
            case '?':
                fprintf(stderr, "Unknown option character -%c", optopt);
                usage(argv[0]);
                break;
            default:
                usage(argv[0]);
        }
    }

    // Make sure we have 1 non option argument;
    if( optind + 1 > argc ) {
        usage(argv[0]);
    }

    // Grab the file name.
    char *fn = argv[optind];

    // Keep the settings.
    sc_settings_t settings = {
        .min_sc_bases = min_sc_bases,
        .min_sc_reads = min_sc_reads,
        .region = region,
        .fn = fn,
        .debug = debug
    };

    // Open file and index and create header object.
    htsFile *in_fp = hts_open( fn, "r");
    hts_idx_t *idx = sam_index_load( in_fp, fn ); 

    if( idx == NULL ) {
        fprintf( stderr, "Could not find index for input alignment file %s\n", fn);
        exit(2);
    }

    bam_hdr_t *header = sam_hdr_read( in_fp );

    // If no region was provided, run for all regions found in header.
    if( strcmp(region,"all") == 0 ) {
        int i;
        for(i=0; i<header->n_targets; i++) {
            if( debug ) {
                fprintf(stderr, "Processing region %s\n", header->target_name[i]);
            }
            process_region( in_fp, idx, header, header->target_name[i], settings );
        }

    } else {
        if( debug ) {
            fprintf(stderr, "Processing region %s\n", region );
        }
        process_region( in_fp, idx, header, region, settings );
    }

    bam_hdr_destroy( header );
    hts_idx_destroy( idx );
    hts_close( in_fp);

    return 0;
}

static void process_region( htsFile *in_fp, hts_idx_t *idx, bam_hdr_t *header, char *region, sc_settings_t opts ) {
    hts_itr_t *iter = sam_itr_querys( idx, header, region );

    // This is used to store information from the iterator.
    bam1_t *b = bam_init1();

    int count = 0;
    int32_t prev_map_pos = 0;

    int sc_reads_n = 1000;
    sc_read_t **sc_reads = (sc_read_t**)calloc( sc_reads_n, sizeof(sc_read_t*));
    if( opts.debug ) {
        fprintf(stderr, "\tAllocated sc_reads array [%d]\n", sc_reads_n );
    }

    while( sam_itr_next( in_fp, iter, b) >= 0 ) {

        // Assert that the bam was sorted by position.
        assert( b->core.pos >= prev_map_pos && 
                "Input bam/cram should be sorted in position order" );
        prev_map_pos = b->core.pos;

        if( is_soft_clipped( b ) ) {

            // Filter 
            if( b->core.qual == 0 || b->core.flag & BAM_FDUP ) {
                continue;
            }

            // Allocate space for the read struct
            sc_read_t *read = (sc_read_t*)malloc(sizeof(sc_read_t));

            // Store some information about the read
            read->read_name = (char*)malloc(sizeof(char)*b->core.l_qname);
            strncpy( read->read_name, bam_get_qname(b), b->core.l_qname);
            read->l_qseq = b->core.l_qseq;

            // Add Cigar info.
            populate_cigar_info( b, read );

            // Filter any reads which have less than 10 soft clipped bases.
            if( read->sc_len < opts.min_sc_bases ) {
                destroy_read( read );
                continue;
            }

            read->mapped_pos = b->core.pos;
            read->tid = b->core.tid;
            read->seq = (uint8_t*)malloc(sizeof(uint8_t) * b->core.l_qseq);
            memcpy( read->seq, bam_get_seq(b), sizeof(uint8_t) * b->core.l_qseq );
            if( read->side == right ) {
                read->sc_pos = read->mapped_pos + (read->l_qseq - read->sc_len); 
            } else {
                read->sc_pos = read->mapped_pos;
            }

            count++;
            sc_reads[count-1] = read;
            if( count >= sc_reads_n ) {
                sc_reads_n *= 2;
                sc_reads = (sc_read_t**)realloc(sc_reads, sizeof(sc_read_t*) * sc_reads_n );
            }
        }

    }    
    if( opts.debug ) {
        fprintf(stderr, "\tFound %i soft clipped reads\n", count);
        fprintf(stderr, "\tFinal array size: %d\n", sc_reads_n );
    }

    int cluster_count = 0;

    if( count != 0 ) {

        // Sort the reads by soft clipped position.
        qsort(sc_reads, count, sizeof(sc_read_t*), sort_by_sc_pos);

        int j;
        int cur_pos = sc_reads[0]->sc_pos;
        sc_side_t cur_side = sc_reads[0]->side;
        int clust_n = 0;

        int cluster_arr_size = 10;
        sc_read_t **cluster = (sc_read_t**)calloc(cluster_arr_size, sizeof(sc_read_t*));

        for(j = 0; j < count; j++ ) {

            // If we have a new cluster
            if( (sc_reads[j]->sc_pos != cur_pos || sc_reads[j]->side != cur_side) && clust_n > 0 ) {
                if( clust_n >= opts.min_sc_reads ) {
                    cluster_count++;
                    handle_cluster( cluster, clust_n, header );
                }

                cur_pos = sc_reads[j]->sc_pos;
                cur_side = sc_reads[j]->side;
                clust_n = 1;


                // Clear the last cluster and create some space for the new.
                // I think we should also be able to free the reads here instead
                // of after the loop. TODO
                free(cluster);
                cluster_arr_size = 20;
                cluster = (sc_read_t**)calloc(cluster_arr_size, sizeof(sc_read_t*));
                cluster[0] = sc_reads[j];


                // We found a member of the cluster.
            } else {
                cluster[clust_n] = sc_reads[j];
                clust_n++;
                if( clust_n >= cluster_arr_size ) {
                    cluster_arr_size *= 2;
                    cluster = (sc_read_t**)realloc(cluster, sizeof(sc_read_t*) * cluster_arr_size);
                }
            }

        }

        if( clust_n >= opts.min_sc_reads ) {
            cluster_count++;
            handle_cluster( cluster, clust_n, header );
        }
        free(cluster);
    }

    // Free the memory allocated for the soft clipped reads
    int i;
    for( i = 0; i < count; i++ ) {
        free( sc_reads[i]->seq );
        destroy_read( sc_reads[i] );
    }

    free(sc_reads);
    hts_itr_destroy( iter );
    bam_destroy1( b );

    if( opts.debug ) {
        fprintf(stderr, "\tFound %d clusters\n\n", cluster_count );
    }
}

static void destroy_read( sc_read_t *r ) {
    free( r->cigar );
    free( r->read_name );
    free( r );
}

static void handle_cluster( sc_read_t **cluster, int cluster_n, bam_hdr_t *header ) {
    int max_sc_len = 0;
    int max_an_len = 0;
    char *clipped_seqs[cluster_n];
    char *anchor_seqs[cluster_n];

    int tid = -1;
    int sc_pos = -1;
    sc_side_t side = left;
    int i;
    for(i=0; i<cluster_n; i++) {
        sc_read_t *r = cluster[i];

        if( tid == -1 ) {
            tid = r->tid;
            side = r->side;
            sc_pos = r->sc_pos;
        }

        if( r->sc_len > max_sc_len ) {
            max_sc_len = r->sc_len;
        }
        if( (r->l_qseq - r->sc_len) > max_an_len ) {
            max_an_len = r->l_qseq - r->sc_len;
        }

        char *clipped_seq = (char*)malloc( sizeof(char) * r->sc_len + 1 );
        char *anchor_seq = (char*)malloc( sizeof(char) * (r->l_qseq - r->sc_len) + 1 );

        int j;

        for(j=0; j< r->l_qseq; j++) {
            int junct = r->sc_len;
            if( side == right ) {
                junct = (r->l_qseq - r->sc_len);
            }
            
            char base = seq_get_base( bam_seqi( r->seq, j ) );
            if( j < junct ) {
                if( r->side == left ) {
                    clipped_seq[j] = base;
                } else {
                    anchor_seq[j] = base;
                }
            } else {
                if( r->side == left ) {
                    anchor_seq[(j - junct)] = base;
                } else {
                    clipped_seq[(j - junct)] = base;
                }
            }
        }

        clipped_seq[r->sc_len] = '\0';
        anchor_seq[(r->l_qseq-r->sc_len)] = '\0';
        clipped_seqs[i] = clipped_seq;
        anchor_seqs[i] = anchor_seq;
    }


    // Get the consensus sequence
    char *clipped_consensus = (char*)malloc(sizeof(char)*max_sc_len + 1);
    char *anchored_consensus = (char*)malloc(sizeof(char)*max_an_len + 1);

    //printf("Clipped:\n");
    if( side == left ) { 
        get_consensus( clipped_consensus, max_sc_len, clipped_seqs, cluster_n, right );
        //print_right_align(max_sc_len, clipped_seqs, cluster_n );
        //printf("\n%s\n", clipped_consensus );
    } else {
        //print_left_align( clipped_seqs, cluster_n );
        get_consensus( clipped_consensus, max_sc_len, clipped_seqs, cluster_n, left );
        //printf("\n%s\n", clipped_consensus );
    }

    //printf("\nAnchor:\n");
    if( side == left ) {
        //print_left_align( anchor_seqs, cluster_n );
        get_consensus( anchored_consensus, max_an_len, anchor_seqs, cluster_n, left );
        //printf("\n%s\n", anchored_consensus );
    } else {
        //print_right_align( max_an_len, anchor_seqs, cluster_n );
        get_consensus( anchored_consensus, max_an_len, anchor_seqs, cluster_n, right );
        //printf("\n%s\n", anchored_consensus );
    }

    printf("%s:%i\t%s\t%d\t%s\t%s\n", header->target_name[tid], sc_pos, get_side_str(side), cluster_n, clipped_consensus, anchored_consensus);

    for(i=0; i<cluster_n; i++) {
        free(clipped_seqs[i]);
        free(anchor_seqs[i]);
    }

    free(clipped_consensus);
    free(anchored_consensus);
}

static char *get_consensus( char *retval, int max_len, char **seqs, int n_seqs, sc_side_t align ) {
    int i;
    for(i=0; i<max_len; i++) {
        int max=0, j, ac=0,cc=0,gc=0,tc=0,nc=0;
        char maxval = 'n';
        int char_ind = i;
        int retval_ind = i;
        if( align == right ) {
            retval_ind = (max_len - i) - 1;
        }
        for(j=0; j<n_seqs; j++) {
            int seqlen = strlen(seqs[j]);

            if( align == right ) {
                if( seqlen <= i ) {
                    continue;
                }
                char_ind = (seqlen - i) - 1;
            } else {
                if( i >= strlen(seqs[j]) ) {
                    continue;
                }
            }
            switch(seqs[j][char_ind]) {
                case 'a':
                    ac++;
                    if( ac > max ) {
                        max = ac;
                        maxval = 'a';
                    }
                    break;
                case 'c':
                    cc++;
                    if( cc > max ) {
                        max = cc;
                        maxval = 'c';
                    }
                    break;
                case 'g':
                    gc++;
                    if( gc > max ) {
                        max = gc;
                        maxval = 'g';
                    }
                    break;
                case 't':
                    tc++;
                    if( tc > max ) {
                        max = tc;
                        maxval = 't';
                    }
                    break;
                case 'n':
                    nc++;
                    if( nc > max ) {
                        max = nc;
                        maxval = 'n';
                    }
                    break;

            }
        }
        retval[retval_ind] = maxval;
    }
    retval[max_len] = '\0';

    return retval;
}

static void print_left_align( char **strs, int n_strs ) {
    int i;
    for(i=0;i<n_strs;i++) {
        printf("%s\n", strs[i]);
    }
}

static void print_right_align( int max, char **strs, int n_strs ) {
    int i;
    for(i=0;i<n_strs;i++) {
        printf("%*s\n", max, strs[i]);
    }
}

static char seq_get_base( int i ) {
    char base;
    switch( i ) {
        case 1:
            base = 'a';
            break;
        case 2:
            base = 'c';
            break;
        case 4:
            base = 'g';
            break;
        case 8:
            base = 't';
            break;
        default: 
            base = 'n';
            break;
    }
    return base;
}


static char* get_side_str( sc_side_t side ) {
    char* retval = "left";
    if( side == right ) {
        retval = "right";
    }
    return retval;
}

static int sort_by_sc_pos( const void *a, const void *b ) {
    sc_read_t *reada = *(sc_read_t**)a;
    sc_read_t *readb = *(sc_read_t**)b;
    int retval = (reada->sc_pos - readb->sc_pos);
    if( retval == 0 && reada->side != readb->side ) {
        if( reada->side == left ) {
            retval = 1;
        } else {
            retval = -1;
        }
    }

    return retval;
}

static void populate_cigar_info( bam1_t *b, sc_read_t *read ) {

    read->cigar = (uint32_t*)malloc(sizeof(uint32_t) * b->core.n_cigar);
    memcpy( read->cigar, bam_get_cigar(b), sizeof(uint32_t) * b->core.n_cigar );
    read->n_cigar = b->core.n_cigar;

    // We already know that this is soft clipped. If the first cigar
    // entry is soft clipped, that means left. Otherwise right.
    // ** Should have already filtered out reads with soft clipping
    // on both ends.
    if( bam_cigar_op(read->cigar[0]) == BAM_CSOFT_CLIP ) {
        read->side = left;
        read->sc_len = bam_cigar_oplen(read->cigar[0]);
    } else {
        read->side = right;
        read->sc_len = bam_cigar_oplen(read->cigar[b->core.n_cigar-1]);
    }
}
static void print_read( sc_read_t *read ) { 
    fprintf(stdout, "%s\t%d\t", read->read_name, read->l_qseq );
    int k;
    for( k = 0; k < read->n_cigar; k++ ) { 
        fprintf(stdout, "%u%c",         
                bam_cigar_oplen(read->cigar[k]), 
                bam_cigar_opchr(read->cigar[k]) 
               );
    }
    printf("\t%i\t", read->sc_len);   
    switch( read->side ) {
        case left :
            printf("left");
            break;

        case right :
            printf("right");
            break;
    }

    // Pretty sure this is now how you're supposed to get
    // the chromosome length. Maybe there is an array lookup
    // to the header somehow?
    printf("\tchr%i:%i\t", read->tid, read->mapped_pos);
    printf("chr%i:%i", read->tid, read->sc_pos);
    
    fprintf( stdout, "\n" );
}

static bool is_soft_clipped( bam1_t *b ) {
    uint32_t *cigar = bam_get_cigar(b);
    int k;
    for (k = 0; k < b->core.n_cigar; ++k) {
        if (bam_cigar_op(cigar[k]) == BAM_CSOFT_CLIP)  {
            return true;
        }
    }
    return false;
}
