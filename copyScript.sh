#!bin/bash
echo "enter name of copied folder:"
read name
cp -r /home/dpd4k4/CA1_model/output /home/dpd4k4/CA1_Results
mv /home/dpd4k4/CA1_Results/output /home/dpd4k4/CA1_Results/${name}
