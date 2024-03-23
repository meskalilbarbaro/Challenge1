# EXEC = "name of the executable"
EXEC=main
# OBJECTS = all .o object files compiled but not linked
OBJECTS= *.o
# SCR = list of all .cpp source files
SRC=main.cpp FunctionWrap.cpp MultFunctionWrap.cpp GradientDescent.cpp
# object files --> executable (linking)
$(EXEC): $(OBJECTS)
	g++ $(OBJECTS) -L/home/meskal_01/pacs-examples/Examples/lib -lmuparser -o $(EXEC)
# source files --> object files (compilation)
$(OBJECTS): $(SRC)
	g++ $(SRC) -I/home/meskal_01/pacs-examples/Examples/include -lmuParser -c
# cleaning routine (optional)
clean:
	@rm -f  *.o ./main
run:
	@LD_LIBRARY_PATH=/home/meskal_01/pacs-examples/Examples/lib  ./$(EXEC)

