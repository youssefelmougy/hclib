if [ $# -ne 2 ]; then
    echo "USAGE: ./run.sh <WORKLOAD> <APP>"
    echo "WORKLOAD=tiny, medium, large, larger, huge"
    exit
fi

SIZE=$1

INPUT_FILE_1="./input/string1-$SIZE.txt"
INPUT_FILE_2="./input/string2-$SIZE.txt"

if [ "$SIZE" = "tiny" ]; then
    TILE_WIDTH=4
    TILE_HEIGHT=4
    EXPECTED_RESULT=12
    ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
else if [ "$SIZE" = "medium" ]; then
    TILE_WIDTH=232
    TILE_HEIGHT=240
    EXPECTED_RESULT=3640
    ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
else if [ "$SIZE" = "large" ]; then
#    # 8x8
#    TILE_WIDTH=2320
#    TILE_HEIGHT=2400
#    ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
#    # 16x16
#    TILE_WIDTH=1160
#    TILE_HEIGHT=1200
#    ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
#    # 32x32
#    TILE_WIDTH=580
#    TILE_HEIGHT=600
#    ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
    # 64x64
    TILE_WIDTH=290
    TILE_HEIGHT=300
    ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT}
    EXPECTED_RESULT=36472
else if [ "$SIZE" = "larger" ]; then
    # 64x64
    TILE_WIDTH=1450
    TILE_HEIGHT=1500
    ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT} -nochecksum
    EXPECTED_RESULT=182392
else if [ "$SIZE" = "huge" ]; then
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
    ./$2 -input1 ${INPUT_FILE_1} -input2 ${INPUT_FILE_2} -tile_width ${TILE_WIDTH} -tile_height ${TILE_HEIGHT} -nochecksum
    EXPECTED_RESULT=364792
fi
fi
fi
fi
fi
