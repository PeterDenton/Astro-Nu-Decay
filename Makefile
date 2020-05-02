CXX=g++
OMPFlags=-fopenmp

CFlags=-c -Wall -O3 -pedantic -MMD -std=c++11 -Werror
CFlags+=-ffast-math
CFlags+=$(OMPFlags)
GSLFlags=`gsl-config --libs`

Sources=$(wildcard src/*.cpp)
IncludeDir=-I./include
AllObjects=$(addprefix obj/,$(notdir $(Sources:.cpp=.o)))
Executables=main
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))

all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	@mkdir -p data obj fig
	$(CXX) $(Objects) $(addprefix obj/,$@.o) $(GSLFlags) $(OMPFlags) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFlags) $(IncludeDir) $< -o $@

-include $(AllObjects:.o=.d)

test: $(Executables)
	$(foreach exe,$(Executables),./$(exe);)

clean:
	rm -rf obj/*.o obj/*.d $(Executables)

