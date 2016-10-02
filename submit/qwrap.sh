#!/usr/bin/env bash

# echo $1
# shift
# echo SHIFTING
# echo $*
umask 0022

cd $1
shift
set -o noglob
exec $*


