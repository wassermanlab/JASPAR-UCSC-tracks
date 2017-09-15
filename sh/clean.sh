for i in {1..13}
do
   echo "compute-1-$i"
   ssh compute-1-$i "rm /state/partition1/oriol/*"
done
