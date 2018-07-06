# run ./stat.sh targetpath will do statistics of dat in that dir, resulting file# stored in a subdir named statistics.
localpath=$(pwd)
echo localpath $localpath
para1=$1
cd $1
targetpath=$(pwd)
cd $targetpath
mkdir statistics
cd statistics
resultpath=$(pwd)

cd $targetpath
for i in *Nwalkers*
do
echo $i
resultfilename=`echo $i.stat`
echo $resultfilename
cp $i $localpath/rawdat
cd $localpath
./stat  rawdat
mv $localpath/rawdat.stat $resultpath/$resultfilename
rm rawdat
cd $targetpath
done
