#!/bin/sh
git submodule init
git submodule update
(cd barvinok; git submodule init polylib; git submodule update polylib; cd ..)
(cd isl; git submodule init imath; git submodule update imath; cd ..)
