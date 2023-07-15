PROJ_NAME=mk_test_data

LINK=$(CXX)
#LINK=$(CC)

SRCS = mk_test_data.cpp
# SRCS += xxx.cpp

###################################################

CFLAGS_ADD :=

DEPSDIR=.deps
OBJDIR=.objs


ALLFLAGS := -g2 -O2 -grecord-gcc-switches -Wall -Wextra -Wfloat-equal
ALLFLAGS += -pthread
ALLFLAGS += -march=native

ALLFLAGS += $(CFLAGS_ADD)

CFLAGS   = $(ALLFLAGS)
CXXFLAGS = $(ALLFLAGS) -std=gnu++23 -fgnu-keywords -Weffc++ -Woverloaded-virtual # -std=c++17 -std=c++20  -std=c++23 -std=gnu++23 -fconcepts

# EXTLIBS += -lm

# ALLFLAGS += -I.

SRCPATHS =  .

vpath %.c   $(SRCPATHS)
vpath %.cpp $(SRCPATHS)
vpath %.s $(STMSRC)
vpath %.o $(OBJDIR)
vpath %.d $(DEPSDIR)


OBJS0a = $(SRCS:.cpp=.o)
OBJS0 = $(OBJS0a:.c=.o)
OBJS  = $(OBJS0:.s=.o)
OBJS1 = $(addprefix $(OBJDIR)/,$(OBJS))

###################################################

.PHONY: proj

all: proj dirs

dirs:
	mkdir -p $(DEPSDIR) $(OBJDIR)

proj:  dirs $(PROJ_NAME)

$(OBJDIR)/*.o:  Makefile


$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c -MMD -o $@ $<
	mv $(OBJDIR)/$*.d $(DEPSDIR)

$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -MMD -o $@ $<
	mv $(OBJDIR)/$*.d $(DEPSDIR)

$(OBJDIR)/%.o: %.s
	$(CC) $(CFLAGS) -c -o $@ $<

$(PROJ_NAME): $(OBJS1)
	$(LINK) $(CFLAGS) $(EXTLIBS) $^ -o $@



clean:
	rm -f *.o *.d $(OBJDIR)/*.o $(DEPSDIR)/*.d

include $(wildcard $(DEPSDIR)/*.d)
#

