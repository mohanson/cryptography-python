set -e

for eachfile in $(ls *.py)
do
    echo python $eachfile
    python $eachfile
done
