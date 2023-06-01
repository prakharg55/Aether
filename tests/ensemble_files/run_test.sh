#!/bin/sh

rm -rf run.test

cp -R ../../share/run ./run.test
cd run.test
# cp ../aether1.json ./aether.json
# mpirun -np 3 ./aether

# cp ../aether2.json ./aether.json
# mpirun -np 2 ./aether

cp ../aether3.json ./aether.json
mpirun -np 4 ./aether