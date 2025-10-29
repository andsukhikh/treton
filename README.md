## Treton

Treton is a software package for numerical calculations of the cores of fast reactors, based on the anisotropic porous solid model. The core code was developed by a team of faculty members at MEPhI.

### Requirements

  * Compiler: g++, Clang
  * C++20 or higher   
  * СMake 3.16 or higher

### Installation

Below are the instructions for installing the project.

1. Clone the repo
   ```
   git clone https://github.com/andsukhikh/scanner
   ```
2. Creating a build folder
   ```
   mkdir build
   ```
3. Сreating assembly instructions and preparing packages wihtout additional flags
   ```
   cmake ..
   ```
4. Assembling
   ```
	cmake --build .
   ```
Additional flags that can be set during the package preparation stage.

1. Flag for O3 optimization
   ```
    cmake -DCMAKE_BUILD_TYPE=Release
   ```
2. Assignment of analytical energy release in the active zone
   ```
    cmake -DANALYTICS=ON
   ```
3. Checking the correctness of the coordinates of the computational grid in the input file Thechyco.ini
   ```
    cmake -DTEST_CRD=ON 
   ```

## Usage

1. Prepare two input files in the `input` folder:
   - `T_IN.txt` — contains temperature data.
   - `Thechyco.ini` — sets the parameters of the active zone.

2. Start the executable file in the folder with the built project.

3. Follow the instructions displayed in the console.

4. Every 20 iterations, the program automatically saves intermediate data to a binary file named `data` in the `input` folder.



