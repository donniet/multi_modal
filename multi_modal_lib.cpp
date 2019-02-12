#include "multi_modal_lib.h"

extern "C" {

    // typedef struct {
    //     mm_vector_float * ds;
    //     size_t dimensions;
    // } multi_modal_wrapper;

    // typedef struct {
    //     float * mean;
    //     size_t mean_size;
    //     float standard_deviation;
    //     size_t sample_count;
    // } distribution_wrapper;

    multi_modal_wrapper * mm_create(size_t dimensions, size_t maximum_nodes) {
        mm_vector_float * ds = new mm_vector_float(maximum_nodes);
        return new multi_modal_wrapper{ds, dimensions};
    }

    void mm_destroy(multi_modal_wrapper * wrapper) {
        delete wrapper->ds;
        delete wrapper;
    }

    void mm_insert(multi_modal_wrapper * wrapper, float * sample, size_t dimensions) {
        wrapper->ds->insert(std::vector<float>(sample, sample + dimensions));
    }

    size_t mm_get_count(multi_modal_wrapper * wrapper) {
        return wrapper->ds->get_count();
    }
    void mm_extract_peaks(multi_modal_wrapper * wrapper, distribution_wrapper ** wrappers, size_t * wrapper_count) {
        auto peaks = wrapper->ds->extract_peaks();

        *wrappers = new distribution_wrapper[peaks.size()];
        *wrapper_count = peaks.size();

        for(size_t i = 0; i < peaks.size(); i++) {
            (*wrappers)[i].mean = new float[wrapper->dimensions];
            std::copy(peaks[i].mean.begin(), peaks[i].mean.end(), (*wrappers)[i].mean);
            (*wrappers)[i].mean_size = wrapper->dimensions;
            (*wrappers)[i].standard_deviation = peaks[i].standard_deviation();
            (*wrappers)[i].sample_count = peaks[i].count;
        }
    }
    void mm_destroy_peaks(multi_modal_wrapper * wrapper, distribution_wrapper * wrappers, size_t wrapper_count) {
        for(size_t i = 0; i < wrapper_count; i++) {
            delete [] wrappers[i].mean;
        }
        delete [] wrappers;
    }
}