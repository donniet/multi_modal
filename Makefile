test: test.cpp include/multi_modal.hpp
	g++ test.cpp -g -o test -Iinclude -std=c++17

multi_modal_lib.o: multi_modal_lib.cpp include/multi_modal.hpp include/multi_modal_lib.h
	g++ -c multi_modal_lib.cpp -g -o multi_modal_lib.o -Iinclude -std=c++17
