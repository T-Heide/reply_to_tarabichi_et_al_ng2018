# compiler:
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -O3 -march=native -Wall

# General variables:
SRC_DIR=src
CLASSES = CellType Cell Phylogeny Universe
CLASSES_HPP = $(addsuffix .hpp, $(CLASSES))
CLASSES_CPP = $(addsuffix .cpp, $(CLASSES))
CLASSES_O = $(addsuffix .o, $(CLASSES))
MAIN_O_FILES_PATH = $(addprefix $(SRC_DIR)/, $(CLASSES_O))
MAIN_HPP_FILES_PATH = $(addprefix $(SRC_DIR)/, $(CLASSES_HPP))

# Targets:
TARGET = gillespie_simulation

all: $(TARGET)

$(TARGET): $(SRC_DIR)/main.cpp $(MAIN_O_FILES_PATH) $(MAIN_HPP_FILES_PATH)
	$(CC) $(CFLAGS) $(MAIN_O_FILES_PATH) $< -o $@

%.o : %.cpp %.hpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET) $(OTHER_TARGETS) $(MAIN_O_FILES_PATH)
