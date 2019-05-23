
NUM=5

while getopts i:n:h option
do
  case "${option}"
  in
    i) FILE=${OPTARG};;
    n) NUM=${OPTARG};;
    h) echo "Options
             -i <input file>
             -n <1..inf> default is 5 i.e. run each experiment 5 times
             -h show options"; exit 0;;
  esac
done

echo file $FILE
echo num $NUM

echo
    echo "#########"
    echo "Figure 1"
    echo "#########"
    echo
    echo "Resilience type|"  "Time(Sec)"
    i=1
    for name in CG
        do
        echo
        echo $name
        for type in baseline replay replication
        do 
            #echo -n "  "$type" "
            printf "  %-12s " $type
            cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}'
            let i+=1
        done
    done

    echo 
    echo "#########"
    echo "Figure 3"
    echo "#########"
    echo
    echo "Error(%)|"  "Time(Sec)|"  "Time increase(%)"
    for name in CG
        do
        echo
        echo $name
        for type in replay replication
        do 
            echo " " $type
            err_rate=0
            error_base=$(cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}')
            printf "  %4d %9.5f %9.5f\n" $err_rate $error_base 0
            let i+=1
            for err_rate in 1 10
            do
                #echo -n "  " $err " "
                error=$(cat $FILE |grep took|head -$(($NUM *$i))|tail -$NUM| awk '{sum+=$4} END{print sum/NR}')
                percent=$( echo "100*($error - $error_base)/$error_base" | bc -l)
                printf "  %4d %9.5f %9.5f\n" $err_rate $error $percent
                let i+=1
            done
        done
    done
echo
