CC=gcc
MPICC=mpicc
CFLAGS=-O3 -lm -Iheaders
SRC_DIR=src
BIN_DIR=bin
DATA_DIR=data
GRAPH_DIR=graph

# lists of sources and exevutable
all: directories $(BIN_DIR)/infect_seq $(BIN_DIR)/infect_central 

# create a folder
directories:
	@mkdir -p $(BIN_DIR) $(DATA_DIR) $(GRAPH_DIR)

# Compilation Seq
$(BIN_DIR)/infect_seq: $(SRC_DIR)/infect_seq.c
	$(CC) -o $@ $< $(CFLAGS)

# Compilation MPI Central (Algather)
$(BIN_DIR)/infect_central: $(SRC_DIR)/infect_central.c
	$(MPICC) -o $@ $< $(CFLAGS)



# to start tests
test: all
	@chmod +x run_test_infect.sh
	./run_test_infect.sh

# to create the graphs 
plot:
	python3 plots_infect.py

clean:
	rm -rf $(BIN_DIR) $(DATA_DIR) $(GRAPH_DIR)