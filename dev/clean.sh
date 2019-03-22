#!/bin/bash

echo Removing all .c, .so and .html files...

find raysect_mayavi -type f -name '*.c' -exec rm {} +
find raysect_mayavi -type f -name '*.so' -exec rm {} +
find raysect_mayavi -type f -name '*.html' -exec rm {} +
rm build -rf
