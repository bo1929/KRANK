# compiler options
#--------------------------------------------
COMPILER = g++
LDLIBS = -lm -lz -lstdc++
CXXFLAGS = -g -std=c++11 -O3 -fopenmp
WFLAGS = -Wno-unused-result -Wno-unused-command-line-argument

# project files
#--------------------------------------------
PROGRAM = kestane
OBJECTS = build/common.o build/library.o build/taxonomy.o build/encode.o build/table.o build/assess.o build/io.o build/lsh.o build/kestane.o

# rules
#--------------------------------------------
all: clean $(PROGRAM)

# generic rule for compiling *.cpp -> *.o
build/%.o: src/%.cpp
	@mkdir -p build
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -c src/$*.cpp -o build/$*.o

$(PROGRAM): $(OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(LDLIBS) $(CPPFLAGS) $(LDFLAGS) -o $@

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@if [ -d build ]; then rmdir build; fi
	@echo "Succesfully cleaned."
