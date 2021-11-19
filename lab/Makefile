# Makefile for genomikon demo programs

LIB = -L../../genomikon -lgenomikon -lm
INC = -I../../genomikon
CC = gcc
CFLAGS = -O2 -Wall -Werror

APP1 = hello
OBJ1 = hello.o

APP2 = train
OBJ2 = train.o

APP3 = sw_v2
OBJ3 = sw_v2.o

APP4 = scrap
OBJ4 = scrap.o


default:
	make $(APP1)
	make $(APP2)
	make $(APP3)
	make $(APP4)
#	make $(APP5)

$(APP1): $(OBJ1)
	$(CC) -o $(APP1) $(OBJ1) $(LIB)

$(APP2): $(OBJ2)
	$(CC) -o $(APP2) $(OBJ2) $(LIB)

$(APP3): $(OBJ3)
	$(CC) -o $(APP3) $(OBJ3) $(LIB)

$(APP4): $(OBJ4)
	$(CC) -o $(APP4) $(OBJ4) $(LIB)

$(APP5): $(OBJ5)
	$(CC) -o $(APP5) $(OBJ5) $(LIB)

clean:
	rm -f *.o $(APP1) $(APP2) $(APP3) $(APP4) $(APP5)

test: $(APP3)
	./sw_v2 a.fa b.fa
%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
