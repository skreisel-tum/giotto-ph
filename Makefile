all: ripser_giotto_redux

ripser_giotto_redux: main.cpp
	c++ -std=c++14 -Wall -pthread main.cpp -o ripser_giotto_redux -g

clean: rm -f ripser_giotto_redux
