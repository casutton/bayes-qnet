SOURCES= pwfun.pyx sampling.pyx qnet.pyx distributions.pyx queues.pyx arrivals.pyx

all: gqnet.so


gqnet.so: $(SOURCES)
		python setup.py build_ext --inplace

clean:
	rm $(SOURCES:.pyx=.c) $(SOURCES:.pyx=.so)
