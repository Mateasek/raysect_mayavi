#!/bin/bash

echo Removing all .c, .so and .html files...

find raycanvas -type f -name '*.c' -exec rm {} +
find raycanvas -type f -name '*.so' -exec rm {} +
find raycanvas -type f -name '*.html' -exec rm {} +
rm build -rf
