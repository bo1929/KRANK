# compiler options
#--------------------------------------------
COMPILER = g++
LDLIBS = -lm -lz -lstdc++
CXXFLAGS = -fopenmp -std=c++11 -g -O3
WFLAGS = -Wno-unused-result

# project files
#--------------------------------------------
PROGRAM = kestane
OBJECTS = build/encode.o build/table.o build/io.o build/kestane.o build/lsh.o build/common.o

# rules
#--------------------------------------------
all: clean $(PROGRAM)

# generic rule for compiling *.cpp -> *.o
build/%.o: src/%.cpp
	@mkdir -p build
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) -c src/$*.cpp -o build/$*.o

$(PROGRAM): $(OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(LDLIBS) -o $@

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@if [ -d build ]; then rmdir build; fi
	@echo "Succesfully cleaned."
