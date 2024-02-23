## Requirements

- Python 3 (tested with 3.11)
- Java Development Kit (tested with OpenJDK 17)
- R (tested with 4.2.2)
- GCC (tested with 12.2.0)
- Git

## Quick install

Once you have verified that you have the required software (see above), you can install the required dependencies by using the provided setup script:

```bash
./setup.sh
```

## Manual installation

### Required packages

Python requirements are listed in `requirements.txt`. To install via pip:

```bash
python3 -m pip install -r requirements.txt
```

Before any further installation step, you need to verify what the software requirements are already installed on your computer. You can do so by running the following commands in your terminal:

```bash
python3 src/bench_radius.py --check-requirements
```

Each failing test will indicate a missing software. To install the missing ones, please refer to the official documentation of the software you need to install.

If you prefer to install the tools and packages manually, you can follow the instructions below.

### Compiling EMOS

EMOS is coded in C, and compiled using GCC:

```bash
pushd src/algorithms/EMOS
gcc -O3 main-emos-mds.c -o emos-mds -DREP
popd
```

### Minimum dominating set heuristic

The MDS heuristic is coded in Java. The original source code is retrieved using git, and patched to fulfill the benchmark requirements:

```bash
pushd src/algorithms
git clone https://github.com/AlejandraCasado/MinimumDominatingSet.git
cd MinimumDominatingSet
git checkout b83429c4e2cfd0aafa2964c0e9f46f63977d4a23
git apply ../MinimumDominatingSet.patch
cd Code/src
javac -cp . heuristic/Main.java
popd
```

#### R interpreter and protoclust.R

Protoclust is available as an R library. The following command will attempt to install it in the default R library location, which might require admin privileges.

```bash
Rscript -e 'install.packages("protoclust")'
```

Alternatively, you can install protoclust in the `src/algorithms` subdirectory:

```bash
Rscript -e 'install.packages("protoclust", "src/algorithms")'
```

### Testing the installation

You can test the installation by running the following command. It will run a sample run on the iris dataset.

```bash
python3 src/bench_radius.py -c
```

If no error is raised, the installation is successful.
