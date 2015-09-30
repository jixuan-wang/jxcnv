#######################
#Source files
#######################

#basic dirs:
SRC_DIR =./sources
BLD_DIR =./build
INC_DIR =$(SRC_DIR)/include

LIB_INCLUDE_DIRS =$(INC_DIR)
LIB_INCLUDE_FLAGS =$(LIB_INCLUDE_DIRS:%=-I%)

#library files:
LIB_SUB_DIR =lib
LIB_SRC_DIR =$(SRC_DIR)/$(LIB_SUB_DIR)
LIB_SRC =$(wildcard $(LIB_SRC_DIR)/*.cpp $(LIB_SRC_DIR)/*.c)

BLD_LIB_DIR =$(BLD_DIR)/$(LIB_SUB_DIR)

#Library Object files:
LIB_OBJECTS =$(patsubst %.c,%.o,$(LIB_SRC:%.cpp=%.o))
LIB_OBJECTS_FULL_PATH =$(subst $(SRC_DIR),$(BLD_DIR),$(LIB_OBJECTS))

# All source files for which dependencies should be generated:
DEPENDENCIES_FILES =$(subst $(SRC_DIR),$(BLD_DIR),$(patsubst %.c,%.dep,$(SRC_DEPENDS:%.cpp=%.dep)))


# All directories in which the build will take place:
ALL_BLD_DIRS = $(BLD_LIB_DIR)

# All source files for which dependencies should be generated:
SRC_DEPENDS = $(LIB_SRC)

# Files to be removed by 'make clean':
REMOVES = $(BLD_DIR)

# Files to be included when making the 'zip' file:
ZIP_INCLUDES =


########################
#Static Libraries
########################

LIB_FULL_PATH =$(BLD_LIB_DIR)/lib$(MAKE_NAME).a


##########################################
#Compiler, Compilation flags:
##########################################
CXX =g++

#For debugging:
DBGCXXFLAGS :=
USE_DEBUG :=

ifneq ($(DEBUG_GDB),)
DBGCXXFLAGS :=-g
endif

ifneq ($(DEBUG),)
DBGCXXFLAGS :=-g
USE_DEBUG :=-DDEBUG=1
endif

ifneq ($(DEBUG_AND_PROF),)
DBGCXXFLAGS :=-g -pg
USE_DEBUG :=-DDEBUG=1
endif


CXXFLAGS =-Wall -Werror -O3 $(DBGCXXFLAGS)
ALL_CXXFLAGS =$(CXXFLAGS) $(LIB_INCLUDE_FLAGS) $(USE_DEBUG)


GCC =gcc
ALL_CFLAGS =-O3 -DSQLITE_OMIT_LOAD_EXTENSION $(LIB_INCLUDE_FLAGS)
