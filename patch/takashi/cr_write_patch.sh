#! /bin/bash

#######################################
# cr_write_patch.sh [patch dir[:patch dir...]]
# creates write_patch.f90 with a subroutine that writes the
# content of the patch dir(s) (*.f90) to disk at run time.
#
# Same output as utils/scripts/cr_write_patch.sh, but the text is
# emitted as one write statement per NGROUP lines (joined with
# new_line) instead of one write statement per line. 35k write
# statements in a single subroutine take ~15 min with ifx; this
# version compiles in seconds.
#######################################

NGROUP=50
export LC_ALL=C   # byte-wise sed/awk: sources may contain UTF-8 comments

if [ $# == 0 ]
    then
    echo "no patches" > .tmp_output.txt
    echo " " >> .tmp_output.txt
    else
    PATCHDIRS=$1
    PATCHDIRS=(${PATCHDIRS//:/ })
    for PATCHDIR in "${PATCHDIRS[@]}"; do
	for filename in ${PATCHDIR}/*.f90; do
            echo "$filename"
            cat "$filename"
	done
    done > .tmp_output.txt
fi

# one Fortran string literal per source line: "..." ($ -> space, " -> ')
sed 's/\$/ /g' .tmp_output.txt | sed "s/\"/'/g" | cat -e | sed 's/\$/\"/' | sed 's/^/"/' > .test_lines.txt

cat << EOF2 > .test_before.f90
subroutine output_patch(filename)
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=30)::format
  integer::ilun

  ilun=11

  fileloc=TRIM(filename)
  format="(A)"
  open(unit=ilun,file=fileloc,form='formatted')
EOF2

awk -v G=$NGROUP '
  { if (n%G==0) { if (n>0) print "    " last; printf "  write(ilun,format) &\n" }
    else        { print "    " last "//new_line(\"a\")// &" }
    last=$0; n++ }
  END { if (n>0) print "    " last }' .test_lines.txt > .test_middle.f90

cat << EOF2 > .test_after.f90
  close(ilun)
end subroutine output_patch
EOF2

cat .test_before.f90 .test_middle.f90 .test_after.f90 > write_patch.f90

rm .tmp_output.txt .test_lines.txt .test_before.f90 .test_middle.f90 .test_after.f90
