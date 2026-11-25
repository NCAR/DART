# Helper for driver scripts that emit TAP.  Source it from a bash script:
#   source "${DART:-$(git rev-parse --show-toplevel)}/developer_tests/tap.sh"
# then call tap_ok / tap_not_ok per case and end with: tap_done; exit $?
# run_tests.sh treats such a script (status "tap" in tests.list) as passing
# only if it exits 0, prints no "not ok" line, and prints a "1..N" plan.

#count up tap calls
tap_n=0
tap_failed=0

tap_ok()     { tap_n=$((tap_n+1)); echo "ok $tap_n - $1"; }
tap_not_ok() { tap_n=$((tap_n+1)); tap_failed=$((tap_failed+1)); echo "not ok $tap_n - $1"; }

# print the plan; return non-zero if any case failed
tap_done() { echo "1..$tap_n"; [ "$tap_failed" -eq 0 ]; }
