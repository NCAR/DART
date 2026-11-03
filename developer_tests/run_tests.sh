#!/usr/bin/env bash
# Build and run DART developer tests; emit aggregate TAP (one line per
# executable), optionally piped through tapview.  See tests.list.
#
# Run with -h for usage.
set -u
here=$(cd "$(dirname "$0")" && pwd)
DART=$(cd "$here/.." && pwd)
export DART
build=1 run=1 list=0 plain=0 out="$here/test_results" template=""
usage() {
  cat <<'USAGE'
Build and run DART developer tests and print TAP (one line per executable),
piped through tapview if available.  Tests are listed in tests.list.

Usage: run_tests.sh [options] [suite | suite/exe ...]

With no arguments, runs every test that is not manual and needs nothing extra.
A bare suite name runs that suite's non-manual tests; suite/exe runs one test,
including manual ones.

Options:
  -l         list tests and their status, then exit
  -b         build only, do not run
  -n         do not build; run existing executables
  -p         plain TAP, do not pipe through tapview
  -o DIR     results directory (default: developer_tests/test_results);
             build logs, test output and per-test copies of work/ go here
  -T FILE    mkmf template to copy to build_templates/mkmf.template
             if that file is missing
  -h, --help show this help

Requirements:
  Builds are serial because suites share generated sources.

Examples:
  ./run_tests.sh -l
  ./run_tests.sh perturb window
  ./run_tests.sh -n -p utilities/sort_test
  ./run_tests.sh -T ../build_templates/mkmf.template.gfortran

Exit status: 0 all tests passed, 1 a test or build failed, 2 usage/setup error.
USAGE
}
[ "${1:-}" = "--help" ] && { usage; exit 0; }
while getopts "blnpo:T:h" o; do
  case $o in
    b) run=0 ;; n) build=0 ;; l) list=1 ;; p) plain=1 ;;
    o) out=$OPTARG ;; T) template=$OPTARG ;;
    h) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
  esac
done
shift $((OPTIND-1))
sel=("$@")
# make the results dir absolute: the script cd's into suite directories
case $out in /*) ;; *) out="$PWD/$out" ;; esac

# ---- read manifest ----
suites=() exes=() stat=() opts=()
while read -r s e st rest; do
  case $s in ''|\#*) continue ;; esac
  suites+=("$s"); exes+=("$e"); stat+=("$st"); opts+=("$rest")
done < "$here/tests.list"

selected() {  # $1=index; manual/needs= tests run only when named as suite/exe
  local s=${suites[$1]} e=${exes[$1]} x skip=0
  [ "${stat[$1]}" = manual ] && skip=1
  case ${opts[$1]} in *needs=*) skip=1 ;; esac
  if [ ${#sel[@]} -eq 0 ]; then [ $skip -eq 0 ]; return; fi
  for x in "${sel[@]}"; do
    [ "$x" = "$s/$e" ] && return 0
    [ "$x" = "$s" ] && [ $skip -eq 0 ] && return 0
  done
  return 1
}

if [ $list -eq 1 ]; then
  printf '%-18s %-30s %-7s %s\n' SUITE EXECUTABLE STATUS OPTIONS
  for i in "${!suites[@]}"; do
    printf '%-18s %-30s %-7s %s\n' "${suites[$i]}" "${exes[$i]}" "${stat[$i]}" "${opts[$i]}"
  done
  exit 0
fi

# ---- unique selected suites ----
want=()
for i in "${!suites[@]}"; do
  selected $i || continue
  case " ${want[*]-} " in *" ${suites[$i]} "*) ;; *) want+=("${suites[$i]}") ;; esac
done
[ ${#want[@]} -gt 0 ] || { echo "no tests selected" >&2; exit 2; }

# ---- build ----
failed_builds=" "   # suites whose build failed; their tests are failed below
if [ $build -eq 1 ]; then
  if [ ! -f "$DART/build_templates/mkmf.template" ]; then
    if [ -n "$template" ] && [ -f "$template" ]; then
      echo "# copying $template to build_templates/mkmf.template" >&2
      cp "$template" "$DART/build_templates/mkmf.template"
    else
      echo "ERROR: $DART/build_templates/mkmf.template is missing." >&2
      echo "  Copy a compiler template there (e.g. mkmf.template.gfortran)," >&2
      echo "  or pass -T <template> to have it copied for you." >&2
      exit 2
    fi
  fi
  build_one() {
    mkdir -p "$out/$1"
    if ( cd "$here/$1/work" && ./quickbuild.sh nompi > "$out/$1/build.log" 2>&1 ); then
      echo "built $1" >&2
    else
      echo "BUILD FAILED $1 (see $out/$1/build.log)" >&2
      failed_builds="$failed_builds $1 "
    fi
  }
  for s in "${want[@]}"; do build_one "$s"; done
fi
# build-only: a failed build must fail the run
if [ $run -eq 0 ]; then
  [ -z "${failed_builds// /}" ] || exit 1
  exit 0
fi

# ---- run, emit TAP ----
mkdir -p "$out"
runtests() {
  local n=0 bad=0 i s e st o d out_txt rc sp sb lines
  echo "TAP version 13"
  for i in "${!suites[@]}"; do
    selected $i || continue
    s=${suites[$i]} e=${exes[$i]} st=${stat[$i]} o=${opts[$i]}
    n=$((n+1)); d="$here/$s/work"
    # a failed build fails every test in the suite, even if a stale
    # executable from an earlier build is still in work/
    case "$failed_builds" in *" $s "*)
      echo "not ok $n - $s/$e # build failed"; echo "# see $out/$s/build.log"
      bad=$((bad+1)); continue ;;
    esac
    if [ ! -x "$d/$e" ]; then
      echo "not ok $n - $s/$e # executable not built"; echo "# see $out/$s/build.log"
      bad=$((bad+1)); continue
    fi
    # run in a fresh copy of work/ so outputs never land in the source tree
    sb="$out/$s/$e"; rm -rf "$sb"; mkdir -p "$sb"; cp -R "$d/." "$sb/"
    sp=""; case $o in *stdin=*) sp=${o#*stdin=}; sp=${sp%% *} ;; esac
    if [ -n "$sp" ]; then out_txt=$(cd "$sb" && ./"$e" < "$sp" 2>&1); rc=$?
    else out_txt=$(cd "$sb" && ./"$e" < /dev/null 2>&1); rc=$?; fi
    printf '%s\n' "$out_txt" > "$out/$s/$e.out"
    ok=1; why=""
    if [ $rc -ne 0 ]; then ok=0; why="exit code $rc"; fi
    if [ "$st" = tap ]; then
      lines=$(printf '%s\n' "$out_txt" | grep -c '^ok ')
      if printf '%s\n' "$out_txt" | grep -q '^not ok '; then ok=0; why="${why:+$why, }failed assertions"; fi
      printf '%s\n' "$out_txt" | grep -q '^1\.\.[0-9]' || { ok=0; why="${why:+$why, }no TAP plan"; }
      why="${why:+$why }(${lines} passing assertions)"
    fi
    if [ $ok -eq 1 ]; then echo "ok $n - $s/$e"
    else
      bad=$((bad+1)); echo "not ok $n - $s/$e${why:+ # $why}"
      printf '%s\n' "$out_txt" | tail -n 40 | sed 's/^/# /'
    fi
  done
  echo "1..$n"
  [ $bad -eq 0 ]
}

if [ $plain -eq 0 ] && command -v tapview >/dev/null 2>&1; then
  runtests | tee "$out/results.tap" | tapview; exit ${PIPESTATUS[0]}
elif [ $plain -eq 0 ] && [ -x "$HOME/bin/tapview" ]; then
  runtests | tee "$out/results.tap" | "$HOME/bin/tapview"; exit ${PIPESTATUS[0]}
else
  [ $plain -eq 1 ] || echo "# tapview not found; printing plain TAP" >&2
  runtests | tee "$out/results.tap"; exit ${PIPESTATUS[0]}
fi
