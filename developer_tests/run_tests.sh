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
jobs=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
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
  -j N       run builds and tests N at a time (default: number of CPUs)
  -T FILE    mkmf template to copy to build_templates/mkmf.template
             if that file is missing
  -h, --help show this help

Build strategy:
  DART core is built once per location module (oned, threed_sphere) into
  lib/<loc>/libdart.a (see lib/build_lib.sh; the two libraries use a union of
  obs types and quantities, set in lib/<loc>/input.nml).  Each test program is
  then compiled and linked against the library, in parallel.  Suites that list
  a test with needs= still use their own work/quickbuild.sh, serially.
  A suite's location is "loc=oned" in tests.list; the default is threed_sphere.

Examples:
  ./run_tests.sh -l
  ./run_tests.sh perturb window
  ./run_tests.sh -n -p utilities/sort_test
  ./run_tests.sh -T ../build_templates/mkmf.template.gfortran

A test can also be a driver script that runs a program many times with
different inputs (e.g. qceff/work/runall.sh). List it in tests.list with
status tap; it must print TAP (see tap.sh) and exit non-zero on any failure.

Exit status: 0 all tests passed, 1 a test or build failed, 2 usage/setup error.
USAGE
}
[ "${1:-}" = "--help" ] && { usage; exit 0; }
while getopts "blnpj:o:T:h" o; do
  case $o in
    b) run=0 ;; n) build=0 ;; l) list=1 ;; p) plain=1 ;;
    j) jobs=$OPTARG ;; o) out=$OPTARG ;; T) template=$OPTARG ;;
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

locof() {  # location module of suite $1 (loc=... in tests.list; default threed_sphere)
  local i
  for i in "${!suites[@]}"; do
    if [ "${suites[$i]}" = "$1" ]; then
      case ${opts[$i]} in *loc=*) l=${opts[$i]#*loc=}; echo "${l%% *}"; return ;; esac
    fi
  done
  echo threed_sphere
}

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
  build_one() {   # legacy per-suite build (suites needing special sources)
    mkdir -p "$out/$1"
    if ( cd "$here/$1/work" && ./quickbuild.sh nompi > "$out/$1/build.log" 2>&1 ); then
      echo "built $1" >&2
    else
      echo "BUILD FAILED $1 (see $out/$1/build.log)" >&2
      failed_builds="$failed_builds $1 "
    fi
  }
  waitslot() { while [ "$(jobs -rp | wc -l)" -ge "$jobs" ]; do sleep 0.2; done; }

  # which suites use the shared libraries, and which locations are needed
  libsuites=() legacy=() locs=()
  for s in "${want[@]}"; do
    uselib=1
    for i in "${!suites[@]}"; do
      [ "${suites[$i]}" = "$s" ] && selected $i && case ${opts[$i]} in *needs=*) uselib=0 ;; esac
    done
    if [ $uselib -eq 1 ]; then
      libsuites+=("$s"); l=$(locof "$s")
      case " ${locs[*]-} " in *" $l "*) ;; *) locs+=("$l") ;; esac
    else legacy+=("$s"); fi
  done

  libfail=" "
  if [ ${#locs[@]} -gt 0 ]; then
    mkdir -p "$out"
    "$here/lib/build_lib.sh" pp > "$out/preprocess_build.log" 2>&1 || true
    for l in "${locs[@]}"; do
      ( "$here/lib/build_lib.sh" "$l" "$jobs" > "$out/lib_$l.log" 2>&1 ) &
    done
    wait
    for l in "${locs[@]}"; do
      if [ -f "$here/lib/$l/libdart.a" ] && grep -q '^built ' "$out/lib_$l.log"; then
        echo "built lib/$l/libdart.a" >&2
      else
        echo "BUILD FAILED lib/$l (see $out/lib_$l.log and lib/$l/*.log)" >&2
        libfail="$libfail$l "
      fi
    done
    # compile every program named in each suite's quickbuild.sh against its library
    tmpl="$DART/build_templates/mkmf.template"
    flags=$(printf 'include %s\nshow:\n\t@echo $(FC)\n\t@echo $(FFLAGS)\n\t@echo $(LIBS)\n' "$tmpl" | make -f - show)
    FC=$(echo "$flags" | sed -n 1p); FFLAGS=$(echo "$flags" | sed -n 2p); LIBS=$(echo "$flags" | sed -n 3p)
    compile_one() {  # suite prog loc
      local s=$1 p=$2 l=$3 d="$out/$1/build"
      if $FC $FFLAGS -I"$here/lib/$l" "$here/$s/$p.f90" "$here/lib/$l/libdart.a" $LIBS \
            -o "$d/$p" -J "$d/mod_$p" > "$d/$p.log" 2>&1; then rm -rf "$d/mod_$p" "$d/$p.log"
      else echo "BUILD FAILED $s/$p (see $d/$p.log)" >&2; rm -rf "$d/mod_$p"; fi
    }
    for s in "${libsuites[@]}"; do
      l=$(locof "$s")
      case "$libfail" in *" $l "*) failed_builds="$failed_builds $s "; continue ;; esac
      rm -rf "$out/$s/build"; mkdir -p "$out/$s/build"
      for p in $(sed -n '/programs=(/,/^)/p' "$here/$s/work/quickbuild.sh" | grep -v '[=()]' | grep -v '^#'); do
        [ -f "$here/$s/$p.f90" ] || continue
        waitslot; mkdir -p "$out/$s/build/mod_$p"; compile_one "$s" "$p" "$l" &
      done
    done
    wait
  fi
  for s in "${legacy[@]-}"; do [ -n "$s" ] && build_one "$s"; done
fi
# build-only: a failed build must fail the run
if [ $run -eq 0 ]; then
  [ -z "${failed_builds// /}" ] || exit 1
  exit 0
fi

# ---- run, emit TAP ----
mkdir -p "$out"
res="$out/.results"; rm -rf "$res"; mkdir -p "$res"

# run one test (index $1); writes $res/$1: line1 "ok"|"not ok", line2 title
# and detail, remaining lines are diagnostics.  Numbering is added when merging.
run_one() {
  local i=$1 s e st o d out_txt rc sp sb lines ok why bad
  s=${suites[$i]} e=${exes[$i]} st=${stat[$i]} o=${opts[$i]}
  d="$here/$s/work"
  local r="$res/$i" bin="$out/$s/build"
  case "$failed_builds" in *" $s "*)
    printf 'not ok\n%s\n# see %s\n' "$s/$e # build failed" "$out/$s/build.log" > "$r"; return ;;
  esac
  # executable: shared-library build (out/suite/build) or legacy build (work/)
  if [ ! -x "$bin/$e" ] && [ ! -x "$d/$e" ]; then
    printf 'not ok\n%s\n# see %s\n' "$s/$e # executable not built" "$bin/$e.log" > "$r"; return
  fi
  # run in a fresh copy of work/ so outputs never land in the source tree
  sb="$out/$s/$e"; rm -rf "$sb"; mkdir -p "$sb"
  # copy only inputs: skip build leftovers (Makefile, preprocess link, objects,
  # logs, quickbuild.sh) and stale binaries; fresh executables come from build/
  local f b
  for f in "$d"/* "$d"/.[!.]*; do
    [ -e "$f" ] || [ -L "$f" ] || continue
    b=${f##*/}
    case $b in Makefile|.cppdefs|preprocess|quickbuild.sh|*.o|*.mod|dart_log.*) continue ;; esac
    if [ -f "$f" ] && [ -x "$f" ]; then case $b in *.sh) ;; *) continue ;; esac; fi
    cp -R "$f" "$sb/"
  done
  [ -d "$bin" ] && find "$bin" -maxdepth 1 -type f -perm -u+x -exec cp {} "$sb/" \;
  # tests may read ../data
  [ -d "$here/$s/data" ] && [ ! -e "$out/$s/data" ] && cp -R "$here/$s/data" "$out/$s/data"
  sp=""; case $o in *stdin=*) sp=${o#*stdin=}; sp=${sp%% *} ;; esac
  if [ -n "$sp" ]; then out_txt=$(cd "$sb" && ./"$e" < "$sp" 2>&1); rc=$?
  else out_txt=$(cd "$sb" && ./"$e" < /dev/null 2>&1); rc=$?; fi
  printf '%s\n' "$out_txt" > "$out/$s/$e.out"
  ok=1; why=""; lines=""
  if [ $rc -ne 0 ]; then ok=0; why="exit code $rc"; fi
  if [ "$st" = passfail ]; then
    # tests in runtime modules that print PASS / FAIL lines but exit 0 on failure
    lines=$(printf '%s\n' "$out_txt" | grep -c 'PASS')
    if printf '%s\n' "$out_txt" | grep -q 'FAIL'; then ok=0; why="${why:+$why, }FAIL in output"; fi
    [ "$lines" -gt 0 ] || { ok=0; why="${why:+$why, }no PASS in output"; }
    [ $ok -eq 0 ] && why="${why:+$why }(${lines} PASS lines)"
  fi
  if [ "$st" = tap ]; then
    lines=$(printf '%s\n' "$out_txt" | grep -c '^ok ')
    if printf '%s\n' "$out_txt" | grep -q '^not ok '; then ok=0; why="${why:+$why, }failed assertions"; fi
    printf '%s\n' "$out_txt" | grep -q '^1\.\.[0-9]' || { ok=0; why="${why:+$why, }no TAP plan"; }
    [ $ok -eq 0 ] && why="${why:+$why }(${lines} passing assertions)"
  fi
  if [ $ok -eq 1 ]; then printf 'ok\n%s\n' "$s/$e${lines:+ # $lines $([ "$st" = passfail ] && echo PASS lines || echo assertions)}" > "$r"
  else
    { printf 'not ok\n%s\n' "$s/$e${why:+ # $why}"
      printf '%s\n' "$out_txt" | tail -n 40 | sed 's/^/# /'; } > "$r"
  fi
}

runtests() {
  local n=0 bad=0 i
  local idx=()
  for i in "${!suites[@]}"; do selected $i && idx+=($i); done
  for i in "${idx[@]}"; do
    while [ "$(jobs -rp | wc -l)" -ge "$jobs" ]; do sleep 0.1; done
    run_one $i &
  done
  wait
  echo "TAP version 13"
  for i in "${idx[@]}"; do
    n=$((n+1))
    { read -r status; read -r title; } < "$res/$i"
    echo "$status $n - $title"
    [ "$status" = "not ok" ] && bad=$((bad+1))
    tail -n +3 "$res/$i"
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
