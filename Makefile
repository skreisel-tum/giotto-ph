all: ripser_gredux_trivial ripser_gredux_tbb

FLAGS=-pthread -MMD -MP

-include *.d

ripser_gredux_trivial: main.cpp
	c++ -std=c++14 -Wall main.cpp -o $@ -Ofast ${FLAGS} -DUSE_TRIVIAL_CONCURRENT_HASHMAP

ripser_gredux_tbb: main.cpp
	c++ -std=c++14 -Wall main.cpp -o $@ -Ofast ${FLAGS} -DUSE_TBB_HASHMAP -ltbb


clean:
	rm -f ripser_gredux_trivial ripser_gredux_tbb

