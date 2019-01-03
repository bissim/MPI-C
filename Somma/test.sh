#!/bin/bash
if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ] ||[ "$4" = "" ] ; then
    echo "Use of script: $0 <sourceFileName> <binaryFileName> <numberOfValues> <iterations> <outputFileName>"
    exit
fi

THIRD_STRATEGY_SOURCE=$1
THIRD_STRATEGY_BINARY=$2

echo "Testing $THIRD_STRATEGY_SOURCE..."
echo

echo "Compiling $THIRD_STRATEGY_SOURCE in $THIRD_STRATEGY_BINARY..."

mpicc $THIRD_STRATEGY_SOURCE -o $THIRD_STRATEGY_BINARY
echo "$THIRD_STRATEGY_BINARY compiled!"
echo

NUM_PROC=2
N=$3
ITERATIONS=$4

OUTPUT_FILE=$5
rm $OUTPUT_FILE
while [ $NUM_PROC -le 8 ]; do
    echo "Running $THIRD_STRATEGY_BINARY over $NUM_PROC processors..."
    echo "TEST WITH $NUM_PROC PROCESSORS" >> $OUTPUT_FILE
    while [ $ITERATIONS -gt 0 ]; do
       # echo "$NUM_PROC processors, $N numbers, $ITERATIONS MORE ITERATIONS"
        echo "$N numbers"
        mpirun -np $NUM_PROC ./$THIRD_STRATEGY_BINARY $N >> $OUTPUT_FILE
        let N=$N*2
        let ITERATIONS=$ITERATIONS-1
        echo "-----" >> $OUTPUT_FILE
        echo "" >> $OUTPUT_FILE
    done
    echo "END TEST FOR $NUM_PROC PROCESSORS" >> $OUTPUT_FILE
    echo "" >> $OUTPUT_FILE
    let NUM_PROC=$NUM_PROC*2
    let N=$3
    let ITERATIONS=$4
    echo
done

echo "Output saved into $OUTPUT_FILE"

echo
echo "End of tests for $THIRD_STRATEGY_BINARY"
