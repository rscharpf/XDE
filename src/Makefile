CXX = /usr/bin/g++ 
BIN = ./
LIBS = ./

flags = -I. -O3


all:	$(BIN)/xdeTest Rinterface_v2.o


Rinterface_v2.o:	Rinterface_v2.cpp Update_v2.h
			$(R_HOME)/bin/R CMD SHLIB -o Rinterface_v2.so Rinterface_v2.cpp

$(BIN)/xdeTest:		mainTest_v2.o Update_v2.o Potential_v2.o \
			Utility_v2.o Random.o Cholesky.o Matrix.o
			$(CXX) $(flags) -o $(BIN)xdeTest mainTest_v2.cpp \
			Update_v2.o Potential_v2.o Utility_v2.o \
			Random.o Cholesky.o Matrix.o -lm

mainTest_v2.o:		mainTest_v2.cpp Potential_v2.h Random.h Utility_v2.h \
			Update_v2.h
			$(CXX) $(flags) -c mainTest_v2.cpp

Update_v2.o:		Update_v2.cpp Random.h Cholesky.h Potential_v2.h \
			Utility_v2.h
			$(CXX) $(flags) -c Update_v2.cpp

Potential_v2.o:		Potential_v2.cpp Random.h Matrix.h Potential_v2.h \
			Utility_v2.h
			$(CXX) $(flags) -c Potential_v2.cpp

Utility_v2.o:		Random.h Utility_v2.cpp Matrix.h Utility_v2.h \
			Potential_v2.h
			$(CXX) $(flags) -c Utility_v2.cpp

Random.o:		Random.cpp Random.h Cholesky.h
			$(CXX) $(flags) -c Random.cpp

Cholesky.o:		Cholesky.cpp Cholesky.h
			$(CXX) $(flags) -c Cholesky.cpp

Matrix.o:		Matrix.cpp Matrix.h
			$(CXX) $(flags) -c Matrix.cpp

Cholesky.h:		Matrix.h
			touch Cholesky.h

Potential_v2.h:		Random.h Matrix.h
			touch Potential_v2.h

Update_v2.h:		Random.h Cholesky.h
			touch Update_v2.h

Utility_v2.h:		Matrix.h
			touch Utility_v2.h