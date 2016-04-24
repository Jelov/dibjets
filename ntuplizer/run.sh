#creating folder for console output
mkdir -p out

echo "Building data ntuples"
sh buildtuples1.sh
echo "Done!"

echo "Building qcd and b-jet ntuples"
sh buildtuples2.sh
echo "Done!"

echo "Building and merging FCR ntuples"
sh buildtuples3.sh
echo "Finished!"