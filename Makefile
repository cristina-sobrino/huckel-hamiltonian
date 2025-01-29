FC = gfortran
LIBS = -llapack

TARGET = huckel
SRC = polyenes.f90 

$(TARGET): $(SRC)
	$(FC) $(FLAGS) -o $(TARGET) $(SRC) $(LIBS)

clean:
	rm -f $(TARGET) *.o *.mod

