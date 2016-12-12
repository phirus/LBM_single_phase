numOfLines=50
echo check in blocks of $numOfLines lines

while [ `ls *.vtk | head -$numOfLines | wc -l` -gt 0 ]
do
fileList=`ls *.vtk | head -$numOfLines`
pvpython conversion.py $fileList
rm $fileList
done

echo all vtk files converted