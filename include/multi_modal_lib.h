#include "multi_modal.hpp"
#include <vector>

typedef multi_modal<std::vector<float>> mm_vector_float;

extern "C" {

typedef struct {
    mm_vector_float * ds;
    size_t dimensions;
} multi_modal_wrapper;

typedef struct {
    float * mean;
    size_t mean_size;
    float standard_deviation;
    size_t sample_count;
} distribution_wrapper;

multi_modal_wrapper * mm_create(size_t dimensions, size_t maximum_nodes);
void mm_destroy(multi_modal_wrapper * wrapper);

void mm_insert(multi_modal_wrapper * wrapper, float * sample, size_t dimensions);
size_t mm_get_count(multi_modal_wrapper * wrapper);
void mm_extract_peaks(multi_modal_wrapper * wrapper, distribution_wrapper ** wrappers, size_t * wrapper_count);
void mm_destroy_peaks(multi_modal_wrapper * wrapper, distribution_wrapper * wrappers, size_t wrapper_count);


}