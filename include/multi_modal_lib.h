#include "multi_modal.hpp"
#include <vector>

typedef multi_modal<std::vector<float>> mm_vector_float;

extern "C" {

typedef struct {
    mm_vector_float * ds;
    unsigned long dimensions;
} multi_modal_wrapper;

typedef struct {
    float * mean;
    unsigned long mean_size;
    float standard_deviation;
    unsigned long sample_count;
    unsigned long id;
} distribution_wrapper;

multi_modal_wrapper * mm_create(unsigned long dimensions, unsigned long maximum_nodes);
void mm_destroy(multi_modal_wrapper * wrapper);

void mm_insert(multi_modal_wrapper * wrapper, float * sample, unsigned long dimensions);
unsigned long mm_get_count(multi_modal_wrapper * wrapper);
void mm_extract_peaks(multi_modal_wrapper * wrapper, distribution_wrapper ** wrappers, unsigned long * wrapper_count);
void mm_destroy_peaks(multi_modal_wrapper * wrapper, distribution_wrapper * wrappers, unsigned long wrapper_count);


}