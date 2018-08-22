
## ILDbench_ttbb4q

ttbar -> bb 4j analysis

See [https://confluence.desy.de/pages/viewpage.action?pageId=87344671](https://confluence.desy.de/pages/viewpage.action?pageId=87344671).

### Installation

Explain here:

- what are the package dependencies (iLCSoft, others ?)
- how to compile your package. Should normally be something like:

```shell
source /path/to/ilcsoft/init_ilcsoft.sh
mkdir build
cd build
cmake -C $ILCSOFT/ILCSoft.cmake ..
make install
```

### How to run the analysis

Explain here:

- where to find data needed for your analysis or how to produce them
- how to run you analysis: 
   - Marlin processors to run ?
   - ROOT macros to run ?
   - Shell scripts ?
   - Run the analysis on grid if you provide scripts for that

Example:

```shell
export MARLIN_DLL=./lib/libILDbench_ttbb4q.so
Marlin ./scripts/ExampleProcessor.xml
```

If you want to provide a lot of details on your analysis, use the doc/Readme.md and point to it from this Readme.md file:

More documentation available here in [doc/Readme.md](doc/Readme.md) !

### Issues and contact

Explain here how can people reach you:

- via the Github issue interface. For the skeleton package: https://github.com/ILDAnaSoft/ILDbench_ttbb4q/issues
- **not mandatory**:
    - email address
    - working institute


