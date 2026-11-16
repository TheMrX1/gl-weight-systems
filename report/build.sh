#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
tectonic -X compile report/main.tex
