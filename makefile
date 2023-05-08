# compiler options
#--------------------------------------------
COMPILER = g++
FLAGS = -std=c++11 -fopenmp -O3 -g -Wno-unused-result -lstdc++ -lm -lz

# project files
#--------------------------------------------
PROGRAM = kestane
OBJECTS = build/encode.o build/fastcluster.o build/index.o build/io.o build/kestane.o build/lsh.o build/common.o

# rules
#--------------------------------------------
all: clean $(PROGRAM)

# generic rule for compiling *.cpp -> *.o
build/%.o: src/%.cpp
	@mkdir -p build
	$(COMPILER) $(FLAGS) -c src/$*.cpp -o build/$*.o

$(PROGRAM): $(OBJECTS)
	$(COMPILER) $(FLAGS) -o $@ $+

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@echo "Succesfully cleaned."
