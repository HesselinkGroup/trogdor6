# Trogdor 6 FDTD

Trogdor 6 FDTD
Paul C Hansen

Copyright 2018 Stanford University.  All rights reserved.

## Introduction

Forward and adjoint FDTD solver including dispersive materials.  Side-product of Paul Hansen's Ph.D. research at Stanford University.  Abandoned ca. January 2013 and then cleaned up for release in 2018.

This C++ code is simply an E&M solver.  It reads simulation descriptions in XML format and also reads and writes raw binary data, and writes the occasional .m file.  Its interface is written in Matlab.

### Dependencies

Boost and CMake.

This code contains the tinyxml library, copyright 2002-2003 Lee Thomason (www.grinninglizard.com).

### Installation

From the root-level directory:

```
mkdir build
cd build
cmake ..
make
```

The FDTD executable will be build/FDTD/trogdor6.

### Usage

```
trogdor6 parameters.xml
```

The XML format is undocumented.  This code should be used through its Matlab interface, which generates the XML.

## License

#### (The MIT License)

MIT License

Copyright (c) 2018 Stanford University

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.






