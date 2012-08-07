#!/bin/bash

# example of compiling with OpenGL:
# fsc -I ../third-party/taoframework-2.1.0/bin -r Tao.OpenGl.dll -r Tao.FreeGlut.dll "$@"

fsc src/Utilities.fs src/MathCore.fs src/LinearAlgebra.fs src/LATest.fs
