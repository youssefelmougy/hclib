if [ $# -ne 2 ]; then
    echo "USAGE: ./run.sh <WORKLOAD> <APP>"
    echo "WORKLOAD=tiny, medium, large, huge"
    exit
fi

#export HCLIB_LOCALITY_FILE=../config/sandia.2socket.json
#export HCLIB_WORKERS=$1
SIZE=$1
#export HCLIB_STATS=1

#if [ $HCLIB_WORKERS -eq 32 ]; then
#    HWLOC_CMD="hwloc-bind -p socket:0.pu:0-15 socket:1.pu:16-31"
#elif [ $HCLIB_WORKERS -le 16 ]; then
#    HWLOC_CMD="hwloc-bind -p socket:0.pu:0-`expr $HCLIB_WORKERS - 1`"
#else
#    echo "HCLIB_WORKERS should be either 32, 16, or less than 16"
#    exit
#fi

echo "hwloc-bind option: $HWLOC_CMD"

INPUT_FILE_1="./input/string1-$SIZE.txt"
INPUT_FILE_2="./input/string2-$SIZE.txt"

if [ "$SIZE" == "tiny" ]; then
    TILE_WIDTH=4
    TILE_HEIGHT=4
    EXPECTED_RESULT=12
    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
else if [ "$SIZE" == "medium" ]; then
    TILE_WIDTH=232
    TILE_HEIGHT=240
    EXPECTED_RESULT=3640
else if [ "$SIZE" == "large" ]; then
    # 8x8
    TILE_WIDTH=2320
    TILE_HEIGHT=2400
    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
    # 16x16
    TILE_WIDTH=1160
    TILE_HEIGHT=1200
    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
    # 32x32
    TILE_WIDTH=580
    TILE_HEIGHT=600
    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
    # 64x64
    TILE_WIDTH=290
    TILE_HEIGHT=300
    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
    EXPECTED_RESULT=36472
else if [ "$SIZE" == "larger" ]; then
    # 64x64
    TILE_WIDTH=1450
    TILE_HEIGHT=1500
    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT} -nochecksum
    EXPECTED_RESULT=182392
else if [ "$SIZE" == "huge" ]; then
    # 8x8
#    TILE_WIDTH=23200
#    TILE_HEIGHT=24000
#    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
#    # 16x16
#    TILE_WIDTH=11600
#    TILE_HEIGHT=12000
#    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
#    # 32x32
#    TILE_WIDTH=5800
#    TILE_HEIGHT=6000
#    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT} -nochecksum
    # 64x64
    TILE_WIDTH=2900
    TILE_HEIGHT=3000
    $HWLOC_CMD ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT} -nochecksum
    EXPECTED_RESULT=364792
fi
fi
fi
fi
fi
