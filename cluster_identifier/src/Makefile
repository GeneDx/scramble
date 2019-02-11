## Compilation options/libs
CC = gcc 
CFLAGS = -Wall -O2
INCLUDES = -I/usr/local/include/htslib
LIBS = /usr/local/lib/libhts.a -lz -lpthread -llzma -lbz2 -lcurl -lcrypto

## source and exe names
VERSION = 0.1.0
SRC = cluster_identifier.c
BASEMAIN = cluster_identifier
BUILD_DIR = build
MAIN = $(BUILD_DIR)/$(BASEMAIN)

.PHONY: clean

all: $(MAIN)
	@echo $(MAIN)

$(MAIN):
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(SRC) $(LIBS)

clean:
	$(RM) $(BUILD_DIR)/$(BASEMAIN)
