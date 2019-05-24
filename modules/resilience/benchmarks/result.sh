
NUM=5
FIG=1

while getopts f:i:n:h option
do
  case "${option}"
  in
    f) FIG=${OPTARG};;
    i) FILE=${OPTARG};;
    n) NUM=${OPTARG};;
    h) echo "Options
             -f <1|2|3|4> default is 1 i.e. generate figure 1
             -i <input file>
             -n <1..inf> default is 5 i.e. run each experiment 5 times
             -h show options"; exit 0;;
  esac
done

echo file $FILE
echo figure $FIG
echo num $NUM

echo
if [ $FIG -eq 1 ]; then
    echo "#########"
    echo "Figure 1"
    echo "#########"

    count=$(cat $FILE |grep took|wc -l)
    total=$((NUM*10))
    if [ $count -ne $total ]; then
        echo
        echo "WARNING: Requires " $total " readings but only " $count " found"
        echo
    fi
    echo
    echo "Resilience type|"  "Time(Sec)"
    i=1
    for name in stencil1D smith-waterman cholesky
        do
        echo
        echo $name
        for type in baseline replay replication
        do 
            #echo -n "  "$type" "
            printf "  %-12s " $type
            cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}'
            i=$((i+1))
        done
    done
    #echo -n "abft        "
    type=abft
    printf "  %-12s " $type
    cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}'

elif [ $FIG -eq 2 ]; then
    echo "#########"
    echo "Figure 2"
    echo "#########"

    count=$(cat $FILE |grep took|wc -l)
    total=$((NUM*6))
    if [ $count -ne $total ]; then
        echo
        echo "WARNING: Requires " $total " readings but only " $count " found"
        echo
    fi
    echo
    echo "Replic(%)|"  "Time(Sec)"
    i=1
    for name in stencil1D #stencil3D
        do
        echo
        echo $name
        for frac in 0 20 40 60 80 100
        do 
            #echo -n $frac " "
            printf "  %-5s " $frac
            cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}'
            i=$((i+1))
        done
    done

elif [ $FIG -eq 3 ]; then
    echo "#########"
    echo "Figure 3"
    echo "#########"

    count=$(cat $FILE |grep took|wc -l)
    total=$((NUM*21))
    if [ $count -ne $total ]; then
        echo
        echo "WARNING: Requires " $total " readings but only " $count " found"
        echo
    fi
    echo
    echo "Error(%)|"  "Time(Sec)|"  "Time increase(%)"
    i=1
    for name in stencil1D smith-waterman cholesky
        do
        echo
        echo $name
        for type in replay replication
        do 
            echo " " $type
            err_rate=0
            error_base=$(cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}')
            printf "  %4d %9.5f %9.5f\n" $err_rate $error_base 0
            i=$((i+1))
            for err_rate in 1 10
            do
                #echo -n "  " $err " "
                error=$(cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}')
                percent=$( echo "100*($error - $error_base)/$error_base" | bc -l)
                printf "  %4d %9.5f %9.5f\n" $err_rate $error $percent
                i=$((i+1))
            done
        done
    done
    type=abft
    echo " " $type
    err_rate=0
    error_base=$(cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}')
    printf "  %4d %9.5f %9.5f\n" $err_rate $error_base 0
    i=$((i+1))
    for err_rate in 1 10
    do
        #echo -n "  " $err " "
        error=$(cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}')
        percent=$( echo "100*($error - $error_base)/$error_base" | bc -l)
        printf "  %4d %9.5f %9.5f\n" $err_rate $error $percent
        i=$((i+1))
    done

elif [ $FIG -eq 4 ]; then
    echo "#########"
    echo "Figure 4"
    echo "#########"

    count=$(cat $FILE |grep took|wc -l)
    total=$((NUM*6))
    if [ $count -ne $total ]; then
        echo
        echo "WARNING: Requires " $total " readings but only " $count " found"
        echo
    fi
    echo
    echo "Resilience type|"  "Time(Sec)"
    i=1
    echo
    echo stencil1D MPI
    for nodes in 2-nodes 4-nodes
    do
        echo " " $nodes
        for type in baseline replay replication
        do 
            #echo -n $type
            printf "    %-12s " $type
            cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}'
            i=$((i+1))
        done
    done

else
    echo "Unknown figure"
    exit
fi

echo
