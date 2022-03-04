#!/usr/bin/env bash

main() {


[ -z "$DART" ] && echo "ERROR: Must set DART environment variable" && exit 9
source "$DART"/build_templates/buildfunctions.sh

# build and run preprocess before making any other DART executables
buildpreprocess

}

main "$@"
