CPPC=g++ -O

SimplRay: Main.cxx SimplRay.hxx 
	$(CPPC) Main.cxx -lm -o SimplRay

clean:
	rm  SimplRay
