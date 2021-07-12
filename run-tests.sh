#!/usr/bin/env bash

shopt -s nullglob

# This script is located on the root of the repository:
basedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check that macrel is in the PATH
if ! which macrel ; then
    echo "which macrel failed."
    exit 1
fi

echo ">>> Running tests with: $(macrel --version) <<<"

ok="yes"
failed_tests=""
for testdir in tests/*; do
    if test -d "$testdir"; then
        cur_ok=yes
        if test -f "${testdir}/TRAVIS_SKIP" -a "x$TRAVIS" = xtrue; then
            echo "Skipping $testdir on Travis"
            continue
        fi
        echo "Running $testdir"
        cd "$testdir"
        if [[ -d out ]]; then
            rm -rf out
        fi
        mkdir -p temp
        export TMPDIR="$PWD/temp"
        chmod 777 command.sh
        ./command.sh >stdout.txt 2>stderr.txt
        macrel_exit=$?
        if [[ $testdir == tests/error-* ]] ; then
            if test $macrel_exit -eq "0"; then
                echo "Macrel exited with exit code 0, even though an error was expected in test"
                cur_ok=no
            fi
        else
            if test $macrel_exit -ne "0"; then
                echo "Error non-zero exit in test: $testdir"
                echo "STDOUT:"
                cat stdout.txt
                echo
                echo "STDERR:"
                cat stderr.txt
                cur_ok=no
            fi
        fi
        for f in expected.*; do
            out=out/macrel.out${f#expected}
            diff -u "$f" "$out"
            if test $? -ne "0"; then
               echo "ERROR in test $testdir: $out did not match $f"
               cur_ok=no
            fi
        done

        if test $cur_ok = "no"; then
            ok=no
            failed_tests="${failed_tests} ${testdir}"
        fi

        rm stdout.txt stderr.txt
        rm -rf temp/
        rm -rf out/
        cd "$basedir"
    fi
done

if test $ok = "yes"; then
    echo "All done."
else
    echo "The following tests failed:"
    for f in $failed_tests; do
        echo " - $f"
    done
    exit 1
fi
