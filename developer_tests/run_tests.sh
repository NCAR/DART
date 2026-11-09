#!/usr/bin/env bash
# Build and run DART developer tests; emit aggregate TAP (one line per
# executable), optionally piped through tapview.  See tests.list.
#
# Run with -h for usage.
set -u
here=$(cd "$(dirname "$0")" && pwd)
DART=$(cd "$here/.." && pwd)
export DART
build=1 run=1 verbose=0 list=0 plain=0 cov=0 nogcovr=0 out="$here/test_results" template=""
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
  -v         for tap tests, also count each assertion as its own test
  -G         with -c, use plain gcov text summary even if gcovr is installed
  -c         gfortran/gcov coverage: build with -O0 -g --coverage and print a
             per-file line coverage summary after the run (HTML too if gcovr
             is installed).  Reports go to DIR/coverage.  Rebuilds the library;
             suites using needs= (built by quickbuild.sh) are not instrumented.
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

  Note: the function compile_one builds a test program from a suite, the quickbuild.sh is
        not run.  Extra sources in the quickbuild.sh EXTRA list (other than the
        models/template model_mod, already in the library) are compiled with the
        test program, e.g. cam_common's cam-common-code modules.

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
while getopts "blnpvcGj:o:T:h" o; do
  case $o in
    b) run=0 ;; v) verbose=1 ;; n) build=0 ;; l) list=1 ;; p) plain=1 ;; c) cov=1 ;; G) nogcovr=1 ;;
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

suite_extra() {  # sources in a suite's quickbuild.sh EXTRA that the library lacks
  # (the library already has models/template/{oned,threed}_model_mod.f90)
  sed -n '/^EXTRA=/,/[^\\]$/p' "$here/$1/work/quickbuild.sh" | tr -d '\\"' | tr ' \t' '\n\n' \
    | sed -e 's/^EXTRA=//' -e "s|\$DART|$DART|" -e 's|\${DART}|'"$DART"'|' \
    | grep '\.f90$' | grep -v 'models/template/' | tr '\n' ' '
}

# ---- build ----
t_start=$SECONDS
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
    export DART_COVERAGE=$cov
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
    if [ $cov -eq 1 ]; then
      FFLAGS="${FFLAGS//-O[0-3s]/} -O0 -g --coverage"; LIBS="$LIBS --coverage"
    fi
    compile_one() {  # suite prog loc
      local s=$1 p=$2 l=$3 d="$out/$1/build"
      local extra; extra=$(suite_extra "$s")   # suite-specific sources, compiled with the program
      if $FC $FFLAGS -I"$here/lib/$l" $extra "$here/$s/$p.f90" "$here/lib/$l/libdart.a" $LIBS \
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
t_build=$((SECONDS - t_start))
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
  if [ $verbose -eq 1 ] && [ "$st" = tap ]; then
    printf '%s\n' "$out_txt" | grep -E '^(not )?ok ' > "$res/$i.sub"
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
    { read -r status; read -r title; } < "$res/$i"
    if [ -s "$res/$i.sub" ]; then   # -v: each assertion is a numbered test
      local line anybad=0 sname=${title%% *}
      while IFS= read -r line; do
        n=$((n+1)); case $line in "not ok "*) anybad=1; bad=$((bad+1)) ;; esac
        echo "${line%% [0-9]*} $n - $sname: ${line#*[0-9] - }" | sed 's/ *$//'
      done < "$res/$i.sub"
      if [ "$status" = "not ok" ] && [ $anybad -eq 0 ]; then
        n=$((n+1)); bad=$((bad+1)); echo "not ok $n - $title"
      fi
      tail -n +3 "$res/$i"
      continue
    fi
    n=$((n+1))
    echo "$status $n - $title"
    [ "$status" = "not ok" ] && bad=$((bad+1))
    tail -n +3 "$res/$i"
  done
  echo "1..$n"
  [ $bad -eq 0 ]
}

# ---- coverage report (gcov) ----
coverage_report() {
  local cdir="$out/coverage" gc d f
  gc=${GCOV:-gcov}
  command -v "$gc" >/dev/null 2>&1 || { echo "# coverage: gcov not found" >&2; return; }
  rm -rf "$cdir"; mkdir -p "$cdir"
  # source directories to report: DART core, the template model_mod built into
  # the library, plus the directory of every EXTRA source in the selected suites
  local dirs=("$DART/assimilation_code" "$DART/models/template") sdir x
  for x in "${want[@]}"; do
    for f in $(suite_extra "$x"); do
      sdir=$(dirname "$f")
      case " ${dirs[*]} " in *" $sdir "*) ;; *) dirs+=("$sdir") ;; esac
    done
  done
  local fargs=() re=""
  for d in "${dirs[@]}"; do fargs+=(--filter "$d"); re="${re:+$re|}${d#$DART/}"; done
  if [ $nogcovr -eq 0 ] && command -v gcovr >/dev/null 2>&1; then
    # --include-internal-functions: gfortran names module procedures __mod_MOD_x,
    # which gcovr would otherwise drop as compiler-generated (giving 0 lines)
    gcovr -r "$DART" --gcov-executable "$gc" --include-internal-functions "${fargs[@]}" \
      --html-details "$cdir/index.html" \
      --txt "$cdir/summary.txt" "$here/lib" "$out" > "$cdir/gcovr.log" 2>&1
    # gcovr tags 0% as coverage-none but colors it like coverage-low; make it
    # stand out (edit the colors here)
    cat >> "$cdir/index.css" <<'CSS'

/* 0% coverage: distinct from low */
table.legend span.coverage-none,
table.coverage td.coverage-none,
table.file-list td.coverage-none { background-color: #b30000 !important; color: #fff !important; }
meter.coverage-none::-webkit-meter-optimum-value { background: #b30000; }
CSS
    echo "# coverage: $cdir/index.html, $cdir/summary.txt" >&2
    return
  fi
  # fallback: plain gcov on the library objects; sources are DART core only
  for d in "$here"/lib/*/ "$out"/*/build/; do
    ls "$d"*.gcda >/dev/null 2>&1 || continue
    n=$(basename "$(dirname "$d")")_$(basename "$d"); mkdir -p "$cdir/$n"
    ( cd "$cdir/$n" && for f in "$d"*.gcda; do "$gc" -o "$d" "$f"; done ) > "$cdir/$n.gcov.log" 2>&1
  done
  # per-source "Lines executed" summary (skip system headers and test support)
  awk -v re="$re" '/^File /{gsub(/\x27/,"",$2); file=$2}
       /^Lines executed:/{ if (file ~ re && file !~ /testanything/) {
         split($0,a,/[:% ]+/); print a[3] "%  of " a[5] " lines  " file } }' \
    "$cdir"/*.gcov.log | sort -u > "$cdir/summary.txt"
  echo "# coverage (line %, lines, file) -- $cdir/summary.txt" >&2
  sed 's/^/# /' "$cdir/summary.txt" >&2
}

t_run0=$SECONDS
[ $cov -eq 1 ] && find "$here/lib" "$out" -name '*.gcda' -delete 2>/dev/null
if [ $plain -eq 0 ] && command -v tapview >/dev/null 2>&1; then
  runtests | tee "$out/results.tap" | tapview; rc=${PIPESTATUS[0]}
elif [ $plain -eq 0 ] && [ -x "$HOME/bin/tapview" ]; then
  runtests | tee "$out/results.tap" | "$HOME/bin/tapview"; rc=${PIPESTATUS[0]}
else
  [ $plain -eq 1 ] || echo "# tapview not found; printing plain TAP" >&2
  runtests | tee "$out/results.tap"; rc=${PIPESTATUS[0]}
fi
[ $cov -eq 1 ] && coverage_report
# timing goes last, after tapview's summary
echo "# run time: $((SECONDS - t_run0))s (${jobs} jobs)$([ $build -eq 1 ] && echo ", build time: ${t_build}s")" >&2
exit $rc
