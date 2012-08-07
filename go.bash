#!/bin/bash

# exit on error and don't allow the use of unset variables
set -o errexit
set -o nounset
set -x

fsc --nologo --out:LATest.exe src/Utilities.fs src/MathCore.fs src/LinearAlgebra.fs src/LATest.fs
mono LATest.exe

