#!/bin/bash
for message in `cut -d '_' -f 1,4,21,24 ../new.txt`
do
        first=`echo $message | cut -d '_' -f1`
        second=`echo $message | cut -d '_' -f2`
        third=`echo $message | cut -d '_' -f 3`
        forth=`echo $message | cut -d '_' -f 4`
        grep $first ../new.txt | grep $second | grep $third| grep $forth | head -n1 >> tmp_select_before_after.txt
done
wc -l tmp_select_before_after.txt
sort tmp_select_before_after.txt | uniq > select_before_after.txt
cat select_before_after.txt | sed 's/#/ /g' | sed 's/_/\t/g' > select_before_after_modify.txt 
rm tmp_select_before_after.txt select_before_after.txt
