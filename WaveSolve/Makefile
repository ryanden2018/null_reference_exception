clean:
	rm -rf out

build-tests:
	mkdir -p out
	mkdir -p out/test
	mkdir -p out/test/wasm
	emcc tests/main.cpp -o out/test/wasm/main.wasm

test:
	wasmtime out/test/wasm/main.wasm

build:
	mkdir -p out
	mkdir -p out/build
	mkdir -p out/build/wasm

dist:
	mkdir -p out
	mkdir -p out/dist
	mkdir -p out/dist/wasm
	mkdir -p out/dist/asmjs
