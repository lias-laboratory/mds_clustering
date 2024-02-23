#!/bin/bash

set -euxo pipefail

: Installing dependencies...
python3 -m pip install -r requirements.txt

cd src/algorithms

: Compiling EMOS...
pushd EMOS
gcc -O3 main-emos-mds.c -o emos-mds -DREP
popd

: Getting and compiling the MinimumDominatingSet heuristic...
if [ -d MinimumDominatingSet ]; then
  rm -rf MinimumDominatingSet
fi
git clone https://github.com/AlejandraCasado/MinimumDominatingSet.git
pushd MinimumDominatingSet
git checkout b83429c4e2cfd0aafa2964c0e9f46f63977d4a23
git apply ../MinimumDominatingSet.patch
cd Code/src
javac -cp . heuristic/Main.java
popd

: Installing protoclust...
Rscript -e 'install.packages("protoclust", ".")'

: Testing the installation...
cd ..  # mds_clustering/src
python3 bench_radius.py -c

: Installation completed successfully
