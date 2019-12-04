#FILES=`cat SRR_Acc_List.txt`

for FILE in $FILES
do
        echo $FILE
        $fastqdump --origfmt -I --split-files --gzip $FILE
done

## fasterq-dump
