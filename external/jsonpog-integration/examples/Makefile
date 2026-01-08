all: jercExample muonExample photonExample btvExample

%: %.C
	g++ $^ -o $@ $(shell correction config --cflags --ldflags --rpath)

clean:
	@rm -r jercExample btvExample

.PHONY: clean
