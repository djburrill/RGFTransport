# Variables
HOME = ../
SOURCE = $(HOME)/src
EXEC = $(HOME)/exec
VPATH = . $(HOME) $(SOURCE) $(EXEC)

# C Flags
CXX=clang++
CFLAGS = -O3
CFLAGS += -std=c++11 -I$(SOURCE) -I$(EXEC)

# Libraries
LIBS:= -larmadillo

# Find source files
SRCFILES:= $(notdir $(wildcard $(SOURCE)/*.cpp $(EXEC)/*.cpp))

# Object & dependency files
odir = $(SOURCE)/o
ddir = $(SOURCE)/d
OBJS:=$(patsubst %.cpp,$(odir)/%.o,$(SRCFILES))
DEPS:=$(patsubst $(odir)/%.o,$(ddir)/%.d, $(OBJS))

# Make directories for object and dependency files
$(odir)/%.o:%.cpp Makefile
	mkdir -p $(odir);$(CXX) -c $(CFLAGS) $< -o $@
	mkdir -p $(ddir);$(CXX) -MM $(CFLAGS) $<  > $*.d;mv $*.d $(ddir)

# Make targets
all: Makefile $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) -o RGFTransport.exe

-include $(DEPS)

clean:
	rm -r *.exe $(odir) $(ddir)
