#!/bin/bash
#
# Getopts pilot
#!/bin/bash 

###############################################################################
#                         Get command line options
# Function to display help text
usage() {
    echo "Usage: $(basename $0) [-r N,...] [-U] "
    echo "Options:"
    echo "  -h          Display this help message"
    echo "  -r          select replicates to process by postion, comma delimited"
    echo "                  no spaces. e.g. -r 2,3 Default is all in config file"
    echo "  -U          Unpaired read (single-ended). Default is paired."
    echo "  -V          very verbose: print every command line"
}
#defaults
endedness="paired"

# Parse options using getopts
while getopts "hr:Uv" option; do
    case "${option}" in
        r)  IFS=,
            repl_indices=($OPTARG)
            ;;
        U)  endedness="single"
            ;;
        h)  # Help option
            usage
            exit 0
            ;;
        V)  set -x
            ;;
        \?) # Invalid option
            echo "Invalid option: -$OPTARG"
            usage
            exit 1
            ;;
    esac
done
shift $(($OPTIND - 1))

###############################################################################
echo " endedness: ${endedness}"
echo " replicate indices: ${repl_indices[@]}"
echo "Remaining arguments are: $*"
###############################################################################
###############################################################################
aa=3